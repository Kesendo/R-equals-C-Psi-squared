"""handshake R_k: the "only computable, not derivable" gate (2026-06-19, Tom's framing).

Tom: a not-derivable result is itself valuable -- it joins the realm of BirthCanal / IsSteril, the
"computable but not closed-form-derivable" category. This script makes the membership FALSIFIABLE for
the bonding-mode defect response: paint the per-site f-profile on the cheap (1,1) block (reused from
_handshake_rk_block), then throw a BATTERY of natural closed forms at the per-bond response and the
per-site location. If they ALL miss (low R^2 / high residual), the response is "only computable" -- there
is no simple closed form, only the exact (computable) painter / Liouville-mixing recipe.

This is NOT a loosening: it is the positive statement. A closed form that FIT would refute the membership.
Anchors: (i) the rate is rate-PROTECTED (re_shift(0,1)~0, U(1)+Pi, proven elsewhere) so f is eigenvector
mixing not a rate functional; (ii) the painter computes f exactly + reproducibly (carrier-independent to
0.999, shown in _handshake_carrier_compare). What remains untested is only whether a SIMPLE closed form
exists -- this gate answers that.
"""
import importlib.util
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, "simulations")
sys.path.insert(0, "simulations/carbon")
from incompleteness_survivor import bonds
from value_vector_felt_time import re_shift


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


blk = _load("blk", "simulations/_handshake_rk_block.py")
psi = blk.psi


def r2_scaled(y, x):
    """Best R^2 of y ~ a*x (single global scale, no per-channel freedom). x,y over the same mask."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if np.allclose(x, 0) or len(y) < 3:
        return float("nan")
    a = (x @ y) / (x @ x)
    ss_res = np.sum((y - a * x) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1.0 - ss_res / max(ss_tot, 1e-30)


def main():
    Q, J, dJ = 1.5, 1.0, 0.02
    g = 1.0 / Q
    print("=== handshake R_k: the 'only computable' gate -- do ANY simple closed forms fit f(b)? ===\n",
          flush=True)
    # anchor: the (0,1) single-excitation rate is protected (f is NOT a rate functional)
    rp = max(abs(re_shift(N, 0, 1, J, g, bonds(N, "chain"), bonds(N, "chain")[len(bonds(N, "chain")) // 2])[0])
             for N in (4, 5, 6))
    print(f"ANCHOR  max |re_shift(0,1)| over N=4,5,6 = {rp:.1e}  (rate PROTECTED -> f is eigenvector mixing)\n",
          flush=True)
    print(f"{'N':>3} | best R^2 of STRENGTH(b) vs a simple closed form           | best LOCATION R^2", flush=True)
    for N in (4, 5, 6, 7, 8, 9):
        bnds = bonds(N, "chain")
        nb = len(bnds)
        f, rel, re = blk.paint_block(N, J, g, dJ, bnds)
        strength = np.array([f[bi][rel[bi]].mean() if rel[bi].any() else np.nan for bi in range(nb)])
        sb = ~np.isnan(strength)
        bb = [b for bi, b in enumerate(bnds) if sb[bi]]
        s = strength[sb]
        ctr = (N - 1) / 2.0
        # candidate closed forms for the per-bond strength profile s(b), b=(j,j+1)
        cands = {
            "const": np.ones(len(bb)),
            "eshift=2psi1 psi1'": np.array([2 * psi(1, b[0], N) * psi(1, b[1], N) for b in bb]),
            "SE-grad^2 (felt_time-D)": np.array([(psi(1, b[0], N) ** 2 - psi(1, b[1], N) ** 2) ** 2 for b in bb]),
            "sin^2(pi(j+1)/N)": np.array([np.sin(np.pi * (b[0] + 1) / N) ** 2 for b in bb]),
            "|j-center|": np.array([abs((b[0] + b[1]) / 2.0 - ctr) for b in bb]),
            "(j-center)^2": np.array([((b[0] + b[1]) / 2.0 - ctr) ** 2 for b in bb]),
        }
        best_name, best_r2 = "", -np.inf
        for name, x in cands.items():
            # center both (allow an additive const) then scale -- so "const" is the null model
            r = r2_scaled(s - s.mean(), x - x.mean()) if name != "const" else 0.0
            if r > best_r2:
                best_r2, best_name = r, name
        # location: per-bond-centered deviation vs the H_1 eigenvector footprint and a seesaw ramp
        m = rel.copy()

        def center(M):
            out = np.zeros_like(M, float)
            for bi in range(nb):
                if m[bi].any():
                    out[bi] = M[bi] - M[bi][m[bi]].mean()
            return out

        fd = center(f)
        seesaw = np.zeros((nb, N))
        for bi, b in enumerate(bnds):
            for a in range(N):
                seesaw[bi, a] = (a - ctr) * np.sign((b[0] + b[1]) / 2.0 - ctr if (b[0]+b[1])/2.0 != ctr else 1)
        loc_r2 = r2_scaled(fd[m], center(seesaw)[m]) if m.sum() > 3 else float("nan")
        print(f"{N:>3} | best: {best_name:>24}  R^2={best_r2:+.3f}              | seesaw R^2={loc_r2:+.3f}",
              flush=True)
    print("\nREAD: R^2 -> 1 for some simple form would REFUTE 'only computable'. Persistent low/erratic R^2"
          "\nacross N = the positive result: f(b) is computable (painter, exact, reproducible) but has no"
          "\nsimple closed form -- the BirthCanal / IsSteril category. The exact content is the Liouville"
          "\nmixing recipe dM_s = sum_s' <W_s'|V_L|M_s>/(lam_s-lam_s') M_s', evaluated, not reduced.", flush=True)


if __name__ == "__main__":
    main()
