#!/usr/bin/env python3
"""
Paper man/trace_envelopes: the violation of trace-characterising inequalities as a function of the Michelson contrasts
of the density W AND of the operator variables, and the question whether projections suffice.

For the inequalities of scripts/trace_characterizations.py whose variables are not projections (items 1-4, 6-8,
10-13, 16, 24, 25), the positive variables (and |X| of the arbitrary and self-adjoint ones) are restricted in three ways:

  free  every spectrum (as in trace_characterizations.py);
  two   two-level spectra, X = c (l + (1 - l) P) with a projection P (or alpha + beta P for self-adjoint variables):
        affine images of projections;
  proj  scaled projections, l = 0 (the smallest l allowed by the item where an inverse or a logarithm is needed).

Experiment A ("do projections suffice?"): V(delta) on a grid of Delta(W) for the three classes on M_2 and M_3.
Experiment B ("both contrasts"): V(delta, delta') with every operator variable of prescribed contrast
Delta = delta' (general spectra and two-level spectra), on M_2 and M_3.

The best final point of every grid value is recomputed in mpmath (trace_characterizations.certify), so every reported
positive value is a certified violation. Writes generated/trace_contrast.{json,tex} and figures/trace_contrast*.pdf
in man/trace_envelopes.

    python scripts/trace_contrast.py --mode full --device cpu
"""
import argparse
import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
for _p in (str(ROOT), str(ROOT / 'scripts')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import trace_characterizations as tc               # noqa: E402
from torch_vn_algebra import optimize as opt       # noqa: E402

OPERATOR_ITEMS = ['gardner', 'triangle', 'realpart', 'square', 'sgm', 'alt', 'young', 'monotone', 'det', 'subadd',
                  'powers', 'pb', 'lp', 'wmon']
MODES = {
    'quick': dict(dims=[2], grid=[0.3, 0.6, 0.9], grid_b=[0.5], vgrid=[0.0, 0.3, 0.6, 0.9], starts=8, rounds=3,
                  steps=40, digits=25),
    'full': dict(dims=[2, 3], grid=[0.1, 0.3, 0.5, 0.7, 0.9], grid_b=[0.3, 0.6, 0.9],
                 vgrid=[0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 0.999], starts=64, rounds=10, steps=80,
                 digits=40),
}


# ----------------------------------------------------------------------------------------------------------------
# restricted spectra
# ----------------------------------------------------------------------------------------------------------------
class TwoLevelSpectrum:
    """Eigenvalues where(D, 1, l): D a fixed rank pattern (1 <= rank <= k - 1 on the first block, so that both levels
    occur and the contrast is exactly (1 - l)/(1 + l)); l prescribed per start, or free in [lmin, 1] ('two'), or lmin
    ('proj')."""

    def __init__(self, alg, B, kind, lmin=0.0, delta=None):
        dev = alg.hilbert.device
        self.alg, self.kind, self.lmin = alg, kind, lmin
        self.mask = alg._const('rows')
        k = alg.k_factors[0]
        r = torch.randint(1, k, (B,), device=dev) if k > 1 else torch.ones(B, dtype=torch.long, device=dev)
        idx = torch.arange(alg.k_max, device=dev)
        hi = torch.zeros(B, alg.C, alg.k_max, dtype=torch.bool, device=dev)
        hi[:, 0, :] = idx[None, :] < r[:, None]
        self.hi = hi
        self.fixed_l = None if delta is None else ((1 - delta) / (1 + delta)).to(torch.float64).to(dev)
        self.eta = torch.randn(B, dtype=torch.float64, device=dev).requires_grad_(kind == 'two' and delta is None)

    def parameters(self):
        return [self.eta] if self.eta.requires_grad else []

    def level(self):
        if self.fixed_l is not None:
            return self.fixed_l
        if self.kind == 'proj':
            return torch.full_like(self.eta, self.lmin)
        return self.lmin + (1 - self.lmin) * 0.5 * (1 + torch.sin(self.eta))

    def delta(self):
        l = self.level()
        return (1 - l) / (1 + l)

    def eigenvalues(self):
        l = self.level()[:, None, None]
        lam = torch.where(self.hi, torch.ones_like(l).expand(self.hi.shape), l.expand(self.hi.shape))
        return lam * self.mask


class TwoLevelSA:
    """Self-adjoint two-level spectra alpha (on D) and beta (off D), alpha, beta free reals ('two'), or (alpha, 0)
    ('proj'); |X| and the signs follow from them. Stands in for SelfAdjointParam.spectrum, with `signs` dynamic."""

    def __init__(self, alg, B, kind):
        dev = alg.hilbert.device
        self.mask = alg._const('rows')
        k = alg.k_factors[0]
        r = torch.randint(1, k + 1, (B,), device=dev)
        idx = torch.arange(alg.k_max, device=dev)
        hi = torch.zeros(B, alg.C, alg.k_max, dtype=torch.bool, device=dev)
        hi[:, 0, :] = idx[None, :] < r[:, None]
        self.hi = hi
        self.ab = torch.randn(B, 2, dtype=torch.float64, device=dev)
        if kind == 'proj':
            self.ab[:, 1] = 0.0
            self.ab[:, 0] = torch.where(torch.rand(B, device=dev) < 0.5, -1.0, 1.0)
        else:
            self.ab.requires_grad_(True)

    def parameters(self):
        return [self.ab] if self.ab.requires_grad else []

    def values(self):
        a, b = self.ab[:, 0, None, None], self.ab[:, 1, None, None]
        v = torch.where(self.hi, a.expand(self.hi.shape), b.expand(self.hi.shape))
        return v * self.mask

    def eigenvalues(self):
        return self.values().abs()

    def delta(self):
        lam = self.eigenvalues()
        return torch.ones(lam.shape[0], dtype=lam.dtype, device=lam.device)


class _DynamicSigns(opt.SelfAdjointParam):
    @property
    def signs(self):
        v = self.spectrum.values()
        return torch.where(v < 0, -1.0, 1.0).to(v.dtype)

    def recenter(self):
        self.unitary.recenter()


SETTINGS = dict(kind='free', delta=None, positive_sa=False)        # set per worker before build()
POSITIVE_SA = {'pb'}    # Peierls-Bogoliubov: positive H, K (Bikchentaev's form) when their contrast is prescribed


def _lmin(dmax):
    return 0.0 if dmax >= 1.0 else (1 - dmax) / (1 + dmax)


class PosVarR(tc.PosVar):
    def __init__(self, alg, B, dmax=1.0):
        kind, delta = SETTINGS['kind'], SETTINGS['delta']
        if kind == 'free' and delta is None:
            super().__init__(alg, B, dmax)
            return
        self.alg = alg
        if kind == 'free':                                # general spectra with prescribed contrast
            d = delta.clamp(max=dmax) if dmax < 1 else delta
            self.p = opt.PositiveParam(alg, B, delta=d, init='pairs', init_fraction=tc.STRUCTURED,
                                       spectrum_init='pairs')
        else:
            self.p = opt.PositiveParam(alg, B, delta=torch.zeros(B, dtype=torch.float64), init='pairs',
                                       init_fraction=tc.STRUCTURED)
            d = None if delta is None else (delta.clamp(max=dmax) if dmax < 1 else delta)
            self.p.spectrum = TwoLevelSpectrum(alg, B, kind, _lmin(dmax), d)
        self.params = [self.p]


class SAVarR(tc.SAVar):
    def __init__(self, alg, B):
        kind, delta = SETTINGS['kind'], SETTINGS['delta']
        if kind == 'free':
            if delta is None:
                super().__init__(alg, B)
            elif SETTINGS['positive_sa']:
                self.p = opt.SelfAdjointParam(alg, B, delta=delta, init='pairs', init_fraction=tc.STRUCTURED,
                                              spectrum_init='pairs')
                self.p.signs = torch.ones_like(self.p.signs)
                self.params = [self.p]
            else:
                self.p = opt.SelfAdjointParam(alg, B, delta=delta, init='pairs', init_fraction=tc.STRUCTURED,
                                              spectrum_init='pairs')
                self.params = [self.p]
            return
        if delta is not None and SETTINGS['positive_sa']:   # positive variables with prescribed contrast
            self.p = opt.SelfAdjointParam(alg, B, delta=delta, init='pairs', init_fraction=tc.STRUCTURED,
                                          spectrum_init='pairs')
            self.p.signs = torch.ones_like(self.p.signs)
            if kind != 'free':
                self.p.spectrum = TwoLevelSpectrum(alg, B, 'two', 0.0, delta)
            self.params = [self.p]
            return
        if delta is not None:                             # prescribed contrast of |X|: two-level |X|, random signs
            self.p = opt.SelfAdjointParam(alg, B, delta=torch.zeros(B, dtype=torch.float64), init='pairs',
                                          init_fraction=tc.STRUCTURED)
            self.p.spectrum = TwoLevelSpectrum(alg, B, 'two', 0.0, delta)
        else:
            self.p = opt.SelfAdjointParam(alg, B, delta=torch.zeros(B, dtype=torch.float64), init='pairs',
                                          init_fraction=tc.STRUCTURED)
            sp = TwoLevelSA(alg, B, kind)
            self.p.__class__ = _DynamicSigns
            del self.p.__dict__['signs']
            self.p.spectrum = sp
        self.params = [self.p]


tc.PosVar = PosVarR
tc.SAVar = SAVarR


# ----------------------------------------------------------------------------------------------------------------
def run_one(job):
    """One envelope: item `name` on M_n, density contrasts `grid`, restriction `kind`, variable contrast `vdelta`
    (None: free). Returns best values per grid point and the certified values of the best final points."""
    name, n, grid, kind, vdelta, cfg, seed = job
    torch.set_num_threads(1)
    starts = cfg['starts']
    B = len(grid) * starts
    SETTINGS['kind'] = kind
    SETTINGS['positive_sa'] = vdelta is not None and name in POSITIVE_SA
    SETTINGS['delta'] = None if vdelta is None else torch.full((B,), float(vdelta), dtype=torch.float64)
    t0 = time.time()
    best, final, points, _ = tc.envelope(name, [n], grid, starts, cfg['rounds'], cfg['steps'], seed)
    cert = []
    for g, pt in enumerate(points):
        try:
            c, d, _ = tc.certify(name, pt, [n], cfg['digits'])
            cert.append(float(c) if abs(d - grid[g]) < 1e-9 else float('nan'))
        except Exception:                                   # noqa: BLE001
            cert.append(float('nan'))
    return dict(name=name, n=n, kind=kind, vdelta=vdelta, grid=list(grid), best=[float(x) for x in best],
                cert=cert, seconds=time.time() - t0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mode', default='quick', choices=list(MODES))
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--items', default=','.join(OPERATOR_ITEMS))
    ap.add_argument('--out', default=str(ROOT / 'man' / 'trace_envelopes'))
    args = ap.parse_args()
    cfg = MODES[args.mode]
    items = args.items.split(',')
    jobs = []
    for name in items:
        for n in cfg['dims']:
            for kind in ('free', 'two', 'proj'):
                jobs.append(('A', (name, n, cfg['grid'], kind, None, cfg, 1000 + n)))
            for vd in cfg['vgrid']:
                for kind in ('free', 'two'):
                    if vd == 0.0 and kind == 'two':
                        continue
                    jobs.append(('B', (name, n, cfg['grid_b'], kind, vd, cfg, 2000 + n)))
    t0 = time.time()
    with ProcessPoolExecutor(args.workers) as ex:
        results = list(ex.map(run_one, [j for _, j in jobs]))
    out = dict(mode=args.mode, config={k: v for k, v in cfg.items()}, seconds=time.time() - t0,
               A=[r for (e, _), r in zip(jobs, results) if e == 'A'],
               B=[r for (e, _), r in zip(jobs, results) if e == 'B'])
    gen = Path(args.out) / 'generated'
    gen.mkdir(parents=True, exist_ok=True)
    (gen / 'trace_contrast.json').write_text(json.dumps(out, indent=1))
    print(f"done in {out['seconds']:.0f} s; wrote {gen / 'trace_contrast.json'}")


if __name__ == '__main__':
    main()
