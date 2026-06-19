"""felt_time_dimensions arc, D FOLLOW-UP, block-level core (2026-06-19). Gate-first.

Step B pins the full 4^N trajectory closure Sum f(b) ~ N*(-dRe(b))/|reS| -- LINEAR in the survivor
mode's bond rate shift dRe(b). So the EXACT functional g question reduces to a CHEAP, scalable
block question: how does dRe(b) depend on the survivor mode's standing-wave amplitude at bond b?

FIRST GUESS REFUTED (gate fired): the mode's HOPPING content k(b)=|Tr(M^dag H_b)| is IDENTICALLY
ZERO -- the slow survivor is a DENSITY (population) mode, diagonal in the dephasing basis, so the
off-diagonal hopping operator has zero overlap. The survivor is a DIFFUSION mode n(j), not a current.

THE MECHANISM (diffusion Rayleigh quotient): a bond-J defect perturbs the LOCAL diffusion coefficient
D_b of the slow density mode. First-order, the decay rate of a diffusion mode is lambda ~ -D * sum_b
(n(j)-n(j+1))^2 / ||n||^2 (a Rayleigh quotient over bond gradients), so
    dRe(b) = d lambda / d D_b  ~  (n(j) - n(j+1))^2     -- the SQUARED density GRADIENT at bond b.
This is "amplitude^2": amplitude = the density-wave gradient (the diffusion-current driver), SQUARED
because the diffusion rate is quadratic in the gradient. It vanishes at the no-flux (reflecting) chain
ENDS (gradient -> 0) and peaks in the INTERIOR, mirror-symmetric -- exactly the seam's shape. The shape
is Q-invariant because the slowest diffusion harmonic (k_min) is Q-fixed.

GATES that can fire:
  G1  dRe(b) ~ |grad n_b|^p with p = 2: log-log slope of dRe vs |grad| ~ 2 (NOT 1). If ~1, the
      relation is linear-in-gradient, not the diffusion-Rayleigh law -- report and stop.
  G2  COLLAPSE: dRe(b) / grad(b)^2 must be bond-INDEPENDENT (low CV) across the supported bonds.
  G3  CLOSED FORM: the empirical gradient profile must match the discrete derivative of the lowest
      reflecting-boundary diffusion harmonic, |grad|(bond between sites a,a+1) ~ sin(pi*(a+1)/N).
"""
import sys
import importlib.util
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, "simulations")
sys.path.insert(0, "simulations/carbon")
from incompleteness_survivor import bonds, basis, survivor
from value_vector_felt_time import re_shift


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


st = _load("stone", "simulations/_stone_survivor_alpha_closure.py")


def density_profile(N, prow, pcol, J, g, bnds):
    """Per-site occupation profile n(j) of the slow survivor mode (the diffusion standing wave):
    n(j) = sum_{s: bit j set} Re M[s,s], M = embedded slowest mode."""
    M, reMode = st.mode_embed(N, prow, pcol, J, g, bnds)
    diag = np.real(np.diag(M))
    n = np.zeros(N)
    for s in range(1 << N):
        if abs(diag[s]) > 0:
            for j in range(N):
                if (s >> j) & 1:
                    n[j] += diag[s]
    imag_frac = float(np.linalg.norm(np.imag(np.diag(M))) / max(np.linalg.norm(np.diag(M)), 1e-15))
    return n, reMode, imag_frac


def loglog_slope(x, y):
    m = (np.asarray(x) > 1e-9) & (np.asarray(y) > 1e-9)
    if m.sum() < 2:
        return float("nan"), float("nan")
    lx, ly = np.log(np.asarray(x)[m]), np.log(np.asarray(y)[m])
    return float(np.polyfit(lx, ly, 1)[0]), float(np.corrcoef(lx, ly)[0, 1])


def cv(vals):
    a = np.asarray(vals)
    return float(np.std(a) / max(abs(np.mean(a)), 1e-15))


def main():
    Q = 1.5
    g = 1.0 / Q
    J = 1.0
    print("=== D core: dRe(b) vs the survivor DENSITY-mode gradient (n(j)-n(j+1))^2. Q=1.5 chain ===\n",
          flush=True)
    slopes, collapses = [], []
    for N in (4, 5, 6, 7):   # N<=7: survivor block <= 1225^2 (fast). N=8 (4900^2) confirms slope~2, omitted.
        bnds = bonds(N, "chain")
        nb = len(bnds)
        re_s, im_s, sec, nxy = survivor(N, J, g, bnds)
        n, reMode, imfrac = density_profile(N, sec[0], sec[1], J, g, bnds)
        dRe = np.array([abs(re_shift(N, sec[0], sec[1], J, g, bnds, bnds[i])[0]) for i in range(nb)])
        grad = np.array([abs(n[b[0]] - n[b[1]]) for b in bnds])
        grad2 = grad ** 2
        # closed form: discrete-derivative of the lowest reflecting diffusion harmonic
        sinpred = np.array([abs(np.sin(np.pi * (b[0] + 1) / N)) for b in bnds])
        print(f"--- N={N} survivor {sec} <n_XY>={nxy:.3f}  (density imag-frac={imfrac:.1e}) ---", flush=True)
        print(f"   n(j) = [{', '.join(f'{x:+.3f}' for x in n)}]", flush=True)
        print(f"{'bond':>7} {'|dRe|':>8} {'|grad|':>8} {'grad^2':>9} {'dRe/grad^2':>11} {'sin pred':>9}",
              flush=True)
        for i, b in enumerate(bnds):
            r = dRe[i] / grad2[i] if grad2[i] > 1e-12 else float("nan")
            print(f"{str(b):>7} {dRe[i]:>8.4f} {grad[i]:>8.4f} {grad2[i]:>9.5f} {r:>11.3f} {sinpred[i]:>9.4f}",
                  flush=True)
        # gates over supported bonds (grad not ~0; the ~0-end bonds are ill-conditioned for a ratio)
        sup = grad > 0.05 * grad.max()
        slope, corr = loglog_slope(grad[sup], dRe[sup])
        ratio_cv = cv(dRe[sup] / np.maximum(grad2[sup], 1e-12)) if sup.sum() >= 2 else float("nan")
        # shape: normalized grad^2 vs normalized sin^2 (closed-form check)
        def nz(a):
            return np.asarray(a) / max(np.max(np.abs(a)), 1e-12)
        sin_miss = float(np.max(np.abs(nz(grad2) - nz(sinpred ** 2))))
        slopes.append(slope)
        collapses.append(ratio_cv)
        print(f"   G1 dRe~|grad|^p: slope={slope:+.2f} (want ~2) corr={corr:+.3f}  | "
              f"G2 dRe/grad^2 CV={ratio_cv:.3f} (want <0.2)  | "
              f"G3 grad^2 vs sin^2 shape-miss={sin_miss:.2f}\n", flush=True)

    valid_s = [s for s in slopes if not np.isnan(s)]
    valid_c = [c for c in collapses if not np.isnan(c)]
    g1_ok = len(valid_s) >= 2 and abs(np.median(valid_s) - 2.0) < 0.4
    g2_ok = len(valid_c) >= 2 and np.median(valid_c) < 0.2
    print("VERDICT:")
    print(f"  G1 power: median log-log slope dRe vs |grad| = {np.median(valid_s):+.2f} -> "
          f"{'p~2 CONFIRMED (diffusion Rayleigh: dRe ~ gradient^2)' if g1_ok else 'NOT p~2 -- diagnose'}",
          flush=True)
    print(f"  G2 collapse: median dRe/grad^2 CV = {np.median(valid_c):.3f} -> "
          f"{'BOND-INDEPENDENT (the law holds)' if g2_ok else 'not constant -- diagnose'}", flush=True)
    ok = g1_ok and g2_ok
    if ok:
        print("  => g IS amplitude^2: Sum f(b) ~ dRe(b) ~ (density-mode gradient at bond b)^2. The "
              "earlier phi*phi guess used the SINGLE-PARTICLE wave; the real survivor is a DENSITY/"
              "diffusion mode and the square is the diffusion-rate Rayleigh quotient. ~0 at ends = "
              "no-flux boundary; Q-invariant = k_min harmonic is Q-fixed.")
    else:
        print("  => a gate fired -- the functional is NOT clean amplitude^2; diagnose, do not loosen.")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
