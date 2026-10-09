"""
Safeguards against unexpectedly long or memory-hungry computations.

* Time: batch-separable kernels (eigh, SVD, QR, matrix exponentials) are routed through
  `batched_call`. For large problems it first runs the kernel on a small slice of the batch,
  extrapolates the time for the whole batch and issues a CostWarning with the estimate when
  it exceeds `warn_seconds`. The result is identical to a single call.
* Memory: `check_memory` compares the estimated size of a large allocation with the memory
  available on the device; above `warn_memory_fraction` it warns, above `max_memory_fraction`
  it raises InsufficientMemoryError (instead of an out-of-memory crash halfway through).
* Integrators (dynamics.rk4 and everything built on it) warn with an ETA after the first step
  and can show a progress bar.

    import torch_vn_algebra as tv
    tv.cost.set_limits(warn_seconds=30)          # thresholds
    with tv.cost.disabled():                     # switch the checks off locally
        ...
    import warnings; warnings.simplefilter('ignore', tv.cost.CostWarning)   # silence warnings only
"""
import contextlib
import time
import warnings
from typing import Callable, Optional

import torch


class CostWarning(UserWarning):
    """A computation is expected to take long or to use a large share of the memory."""


class InsufficientMemoryError(MemoryError):
    """The estimated memory of an operation exceeds what is available on the device."""


_config = dict(
    enabled=True,
    warn_seconds=10.0,           # warn if a single call is expected to take longer
    warn_memory_fraction=0.5,    # warn above this share of the available memory
    max_memory_fraction=0.9,     # raise above this share
    probe_items=8,               # size of the timed slice
    probe_flops=2e8,             # problems below this many flops are not probed
)


def set_limits(**kwargs):
    """Change thresholds: warn_seconds, warn_memory_fraction, max_memory_fraction, probe_items,
    probe_flops, enabled."""
    for key, value in kwargs.items():
        if key not in _config:
            raise KeyError(f"unknown option {key!r}; options: {sorted(_config)}")
        _config[key] = value


def get_limits() -> dict:
    return dict(_config)


@contextlib.contextmanager
def disabled():
    """Temporarily switch all cost checks off."""
    old = _config['enabled']
    _config['enabled'] = False
    try:
        yield
    finally:
        _config['enabled'] = old


def _fmt_bytes(n: float) -> str:
    for unit in ('B', 'KB', 'MB', 'GB', 'TB'):
        if n < 1024 or unit == 'TB':
            return f"{n:.1f} {unit}"
        n /= 1024


def _fmt_seconds(s: float) -> str:
    if s < 120:
        return f"{s:.0f} s"
    if s < 7200:
        return f"{s / 60:.0f} min"
    return f"{s / 3600:.1f} h"


def _sync(device):
    if device.type == 'cuda':
        torch.cuda.synchronize(device)


def available_memory(device) -> Optional[int]:
    """Free memory on `device` in bytes (CUDA: mem_get_info; CPU: MemAvailable on Linux, psutil otherwise)."""
    device = torch.device(device)
    if device.type == 'cuda':
        free, _ = torch.cuda.mem_get_info(device)
        return int(free)
    try:
        with open('/proc/meminfo') as f:
            for line in f:
                if line.startswith('MemAvailable:'):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    try:
        import psutil
        return int(psutil.virtual_memory().available)
    except ImportError:
        return None


def check_memory(nbytes: float, device, what: str):
    """Warn or raise if an allocation of about `nbytes` is large compared with the free memory."""
    if not _config['enabled'] or nbytes < 64 * 2 ** 20:          # ignore below 64 MB
        return
    avail = available_memory(device)
    if avail is None:
        return
    hint = "reduce the batch size or process the batch in chunks"
    if nbytes > _config['max_memory_fraction'] * avail:
        raise InsufficientMemoryError(
            f"{what} needs about {_fmt_bytes(nbytes)}, but only {_fmt_bytes(avail)} are free on "
            f"{torch.device(device)}; {hint} (or relax tv.cost.set_limits(max_memory_fraction=...)).")
    if nbytes > _config['warn_memory_fraction'] * avail:
        warnings.warn(f"{what} will use about {_fmt_bytes(nbytes)} of {_fmt_bytes(avail)} free on "
                      f"{torch.device(device)}; {hint} if this fails.", CostWarning, stacklevel=3)


def tensor_bytes(shape, dtype) -> int:
    n = 1
    for s in shape:
        n *= int(s)
    return n * torch.empty((), dtype=dtype).element_size()


def kernel_flops(kind: str, k: int, complex_: bool) -> float:
    """Rough flop count of one k x k kernel (only used to skip probing of small problems)."""
    c = {'matmul': 2, 'qr': 4, 'eigh': 10, 'svd': 22, 'expm': 30}.get(kind, 10)
    return c * k ** 3 * (4 if complex_ else 1)


def _cat(parts):
    first = parts[0]
    if isinstance(first, torch.Tensor):
        return torch.cat(parts, dim=0)
    if isinstance(first, tuple):           # incl. torch.return_types (eigh, svd, ...)
        return tuple(torch.cat([p[i] for p in parts], dim=0) for i in range(len(first)))
    raise TypeError(f"cannot concatenate results of type {type(first)}")


def batched_call(fn: Callable, x: torch.Tensor, what: str, kind: str = 'eigh'):
    """
    fn(x) for a function acting independently on the items x[0], x[1], ... (leading dimension).
    Large problems are probed on a small slice first; a CostWarning reports the expected time
    if it exceeds warn_seconds. The returned value equals fn(x).
    """
    if not _config['enabled'] or x.dim() < 3 or x.shape[0] == 0:
        return fn(x)
    n_items = x.shape[0]
    probe = _config['probe_items']
    k = x.shape[-1]
    per_item = 1
    for d in x.shape[1:-2]:
        per_item *= int(d)
    flops = n_items * per_item * kernel_flops(kind, k, torch.is_complex(x))
    if n_items < 4 * probe or flops < _config['probe_flops']:
        return fn(x)
    warm = fn(x[:1])                       # first call may include one-off set-up costs
    _sync(x.device)
    t0 = time.perf_counter()
    head = fn(x[1:1 + probe])
    _sync(x.device)
    dt = time.perf_counter() - t0
    eta = dt * (n_items - 1 - probe) / probe
    if eta > _config['warn_seconds']:
        warnings.warn(f"{what}: about {_fmt_seconds(eta)} expected for {n_items * per_item} blocks of size {k} on "
                      f"{x.device} (measured on {probe} blocks). Consider a smaller batch, double -> single "
                      f"precision or a GPU. Silence with warnings.simplefilter('ignore', "
                      f"torch_vn_algebra.cost.CostWarning).", CostWarning, stacklevel=3)
    rest = fn(x[1 + probe:])
    return _cat([warm, head, rest])


class StepTimer:
    """ETA warning (and optional progress bar) for loops with a known number of equal steps."""

    def __init__(self, total_steps: int, what: str, progress: bool = False):
        self.total, self.what, self.done = total_steps, what, 0
        self.t0 = time.perf_counter()
        self.bar = None
        if progress:
            try:
                from tqdm.auto import tqdm
                self.bar = tqdm(total=total_steps, desc=what, leave=False)
            except ImportError:
                self.bar = None
                self._last_print = 0.0

        self.progress = progress

    def step(self, device=None):
        self.done += 1
        if self.bar is not None:
            self.bar.update(1)
        elif self.progress:
            now = time.perf_counter()
            if now - self._last_print > 5 or self.done == self.total:
                el = now - self.t0
                print(f"{self.what}: {self.done}/{self.total}, ETA {_fmt_seconds(el / self.done * (self.total - self.done))}")
                self._last_print = now
        if self.done == 1 and _config['enabled']:
            if device is not None:
                _sync(torch.device(device))
            self.t1 = time.perf_counter()        # the first step includes one-off set-up costs
            if self.total == 2 and (self.t1 - self.t0) > _config['warn_seconds']:
                warnings.warn(f"{self.what}: about {_fmt_seconds(self.t1 - self.t0)} expected for the "
                              f"remaining step.", CostWarning, stacklevel=4)
        if self.done == 2 and _config['enabled']:
            if device is not None:
                _sync(torch.device(device))
            eta = (time.perf_counter() - self.t1) * (self.total - 2)
            if eta > _config['warn_seconds']:
                warnings.warn(f"{self.what}: about {_fmt_seconds(eta)} expected for {self.total} steps "
                              f"(pass progress=True for a progress bar).", CostWarning, stacklevel=4)

    def close(self):
        if self.bar is not None:
            self.bar.close()
