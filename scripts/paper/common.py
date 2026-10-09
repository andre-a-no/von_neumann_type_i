"""Shared helpers for the scripts that produce the numbers and figures of paper/v2."""
import argparse
import atexit
import json
import math
import os
import platform
import sys
import time
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

PAPER = ROOT / 'paper' / 'v2'


def parse_args(description):
    p = argparse.ArgumentParser(description=description)
    p.add_argument('--device', default='cuda' if torch.cuda.is_available() else 'cpu')
    p.add_argument('--mode', choices=['quick', 'full', 'check'], default='quick',
                   help='quick: minutes on a CPU; full: the sizes reported in the paper (GPU); check: the sizes '
                        'of full with minimal repetitions, to test a full run quickly (output not for the paper)')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--out', default=str(PAPER / 'generated'), help='directory for generated .tex/.json')
    p.add_argument('--figdir', default=str(PAPER / 'figures'))
    p.add_argument('--tf32', action='store_true',
                   help='allow TF32 tensor-core matmuls for float32 on Ampere+ GPUs (default: full FP32)')
    args = p.parse_args()
    if args.mode == 'check':          # never overwrite the paper's tables with a check run
        if args.out == str(PAPER / 'generated'):
            args.out = str(ROOT / '.check_output' / 'generated')
        if args.figdir == str(PAPER / 'figures'):
            args.figdir = str(ROOT / '.check_output' / 'figures')
    torch.manual_seed(args.seed)
    # float32 matrix products: 'highest' = true FP32; 'high' lets cuBLAS use TF32 tensor cores
    # (10-bit mantissa) on A100/H100/B200, faster but only ~1e-3 relative accuracy
    torch.set_float32_matmul_precision('high' if args.tf32 else 'highest')
    torch.backends.cuda.matmul.allow_tf32 = bool(args.tf32)
    torch.backends.cudnn.allow_tf32 = bool(args.tf32)
    Path(args.out).mkdir(parents=True, exist_ok=True)
    Path(args.figdir).mkdir(parents=True, exist_ok=True)
    atexit.register(_report_resources, args.device, time.time())
    return args


def _report_resources(device, t0):
    """Printed at exit into every log: wall time and peak memory (host RSS, and GPU if used)."""
    msg = f"[resources] wall {time.time() - t0:.0f} s"
    try:
        import resource
        msg += f", peak host RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20:.2f} GiB"
    except ImportError:
        pass
    if str(device).startswith('cuda') and torch.cuda.is_available():
        msg += f", peak GPU memory {torch.cuda.max_memory_allocated() / 2 ** 30:.2f} GiB"
    print(msg, flush=True)


def environment(device):
    info = {'torch': torch.__version__, 'python': platform.python_version(), 'device': str(device),
            'date': time.strftime('%Y-%m-%d'), 'float32_matmul_precision': torch.get_float32_matmul_precision()}
    if torch.cuda.is_available():
        info['cuda'] = torch.version.cuda
        info['cudnn'] = torch.backends.cudnn.version()
    if str(device).startswith('cuda'):
        info['gpu'] = torch.cuda.get_device_name(torch.device(device))
        info['compute_capability'] = '.'.join(map(str, torch.cuda.get_device_capability(torch.device(device))))
    else:
        info['cpu'] = platform.processor() or platform.machine()
        info['threads'] = torch.get_num_threads()
    return info


def _partial_prefix(args):
    return f"{Path(sys.argv[0]).stem}_{args.mode}_"            # per script: parallel runs do not interfere


def _partial_path(args, key):
    return Path(args.out) / 'partial' / f'{_partial_prefix(args)}{key}.json'


def _fingerprint(args):
    """What a partial result depends on: seed, TF32, device type and the code (script and library)."""
    import hashlib
    h = hashlib.sha256()
    for p in [Path(sys.argv[0]).resolve()] + sorted((ROOT / 'torch_vn_algebra').glob('*.py')) + [Path(__file__)]:
        h.update(p.read_bytes())
    for p in sorted((ROOT / 'results' / 'experiments').glob('*/tables/*.csv.gz')):   # input data of bounds_search
        h.update(p.name.encode() + p.read_bytes())
    return dict(seed=args.seed, tf32=bool(args.tf32), device=torch.device(args.device).type, code=h.hexdigest()[:16])


def load_partial(args, key):
    """
    A partial result saved by save_partial in an earlier, interrupted run (crash, reboot, lost
    connection) with the same seed, TF32 setting, device type and code; None otherwise. Parts of a
    run with other settings or older code are ignored, and finish_partial() deletes this script's
    parts when it completes.
    """
    p = _partial_path(args, key)
    if not p.exists():
        return None
    with open(p) as f:
        saved = json.load(f)
    if saved.get('fingerprint') != _fingerprint(args):
        print(f"[resume] ignoring {p.name}: it comes from a run with other settings or older code", flush=True)
        return None
    print(f"[resume] using {p.name} from an interrupted earlier run", flush=True)
    return saved['result']


def save_partial(args, key, obj):
    p = _partial_path(args, key)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix('.tmp')
    with open(tmp, 'w') as f:
        json.dump(dict(fingerprint=_fingerprint(args), result=obj), f, default=float)
    tmp.replace(p)                                    # atomic: an interrupted write leaves no broken file


def finish_partial(args):
    """Delete this script's partial results (only its own: other scripts may run in parallel)."""
    d = Path(args.out) / 'partial'
    for p in d.glob(f'{_partial_prefix(args)}*.json'):
        p.unlink()
    try:
        d.rmdir()                                       # only if empty
    except OSError:
        pass


def save_json(args, name, payload):
    payload = dict(payload, environment=environment(args.device), mode=args.mode)
    with open(Path(args.out) / f'{name}.json', 'w') as f:
        json.dump(payload, f, indent=2, default=float)


def write_tex(args, name, text):
    with open(Path(args.out) / f'{name}.tex', 'w') as f:
        f.write('% generated by scripts/paper/*.py -- do not edit by hand\n' + text)


def env_macro(args, prefix):
    e = environment(args.device)
    hw = e.get('gpu') or f"CPU ({e.get('threads')} thread{'' if e.get('threads') == 1 else 's'})"
    return (f"\\newcommand{{\\{prefix}Device}}{{{hw}}}\n"
            f"\\newcommand{{\\{prefix}Mode}}{{{args.mode}}}\n"
            f"\\newcommand{{\\{prefix}Torch}}{{{e['torch']}}}\n"
            f"\\newcommand{{\\{prefix}Matmul}}{{{'TF32' if e['float32_matmul_precision'] != 'highest' else 'FP32'}}}\n")


def sync(device):
    if str(device).startswith('cuda'):
        torch.cuda.synchronize()


def fmt(x, digits=4):
    return f"{x:.{digits}g}"


def num3(x):
    """Three significant digits without exponent notation (1234 -> 1230, 0.01234 -> 0.0123)."""
    if x == 0 or not math.isfinite(x):
        return f"{x:g}"
    d = max(0, 2 - int(math.floor(math.log10(abs(x)))))
    return f"{round(x, d):.{d}f}" if abs(x) >= 1e-3 else f"{x:.2e}"


def sci(x, digits=2):
    if x == 0:
        return '0'
    m, e = f"{x:.{digits - 1}e}".split('e')
    return f"{m}\\times10^{{{int(e)}}}"
