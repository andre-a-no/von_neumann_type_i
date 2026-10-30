#!/usr/bin/env python3
"""
Certificates for item 14 of the trace catalogue (Hölder's inequality for projections, alpha = 3) on the delta grid
of the trace preprint (man/article/generated/trace_characterizations.json, delta = 0 excluded: V(0) = 0).

    python scripts/certify_item14.py [--deltas 1/10 1/2 9/10] [--rel-tol 1e-6] [--jobs 2]

Writes certificates/holder_projections_alpha3_delta_<p>_<q>.json.gz and certificates/holder_projections_summary.json.
Check them with:  python torch_vn_algebra/certify/verify.py certificates/holder_projections_*.json.gz
"""
import argparse
from fractions import Fraction
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
GRID = ["1/100", "1/50", "1/20", "1/10", "1/5", "3/10", "2/5", "1/2", "3/5", "7/10", "4/5", "9/10", "19/20"]


def path_of(delta):
    return os.path.join(ROOT, "certificates", f"holder_projections_alpha3_delta_{delta.replace('/', '_')}.json.gz")


def one(args):
    delta, rel_tol = args
    from torch_vn_algebra import certify
    c = certify.certify("holder_projections", delta, rel_tol=rel_tol, alpha=3)
    c.save(path_of(delta))
    return dict(delta=delta, lower=c.data["claim"]["lower"], upper=c.data["claim"]["upper"],
                lower_float=float(c.lower), upper_float=float(c.upper), boxes=c.data["upper"]["leaves"],
                seconds=c.data["seconds"], file=os.path.relpath(path_of(delta), ROOT))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deltas", nargs="*", default=GRID)
    ap.add_argument("--rel-tol", type=float, default=1e-6)
    ap.add_argument("--jobs", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(os.path.join(ROOT, "certificates"), exist_ok=True)
    summ_path = os.path.join(ROOT, "certificates", "holder_projections_summary.json")
    summary = json.load(open(summ_path)) if os.path.exists(summ_path) else {}
    with ProcessPoolExecutor(a.jobs) as ex:
        for r in ex.map(one, [(d, a.rel_tol) for d in a.deltas]):
            print(f"delta = {r['delta']:>6}: [{r['lower_float']:.12g}, {r['upper_float']:.12g}]  "
                  f"{r['boxes']} boxes, {r['seconds']} s", flush=True)
            summary[r["delta"]] = r
            json.dump(dict(sorted(summary.items(), key=lambda kv: Fraction(kv[0]))), open(summ_path, "w"), indent=1)


if __name__ == "__main__":
    main()
