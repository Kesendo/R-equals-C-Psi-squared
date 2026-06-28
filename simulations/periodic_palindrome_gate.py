"""GATE-FIRST verifier (HARDENED): is the periodic-table "palindrome" a real F1
signature beyond smoothness, or a monotonic-smoothness artifact?

Companion to periodic_palindrome.py. The honest verdict is written up in
docs/carbon/PERIODIC_PALINDROME_HARDENED.md.

WHY. periodic_palindrome.py reports the F1 palindrome shape
(v_k + v_{N+1-k} ~ const) across a period as "significant" (shuffle null, p<1e-4).
But pair-sum-constant is satisfied EXACTLY by ANY linear ramp, and the shuffle
null only rejects "unordered". So the hit may be smoothness, not F1. This lets the
elements try to say NO fairly.

DECOMPOSITION. Center each period (w = v - mean). Split about the centre:
  antisym a_k = (w_k - w_{N-1-k})/2   <- F1-RESPECTING (pair-sum const <=> s=0)
  symm    s_k = (w_k + w_{N-1-k})/2   <- F1-BREAKING
Split a into linear ramp a_lin (proj on l_k=k-(N-1)/2; TRIVIAL: any monotone prop)
and non-linear a_non = a - a_lin (F1's NON-TRIVIAL content). The discriminating
residual is r = w - a_lin = a_non + s.

HARDENING -- a real null. R = E_non / (E_non + E_sym) = fraction of the non-ramp
residual that respects the mirror. Magnitude-preserving sign-flip null: for each
mirror pair independently flip the partner's sign (preserves |r_k| exactly,
randomizes only mirror relationship). Flipping a pair SWAPS its antisym<->symm
energy, so the null mean of R is EXACTLY 0.5. Any departure is measurable.

THE GATE (committed before running):
  Per period, R_obs vs the sign-flip null (mean 0.5):
    R_obs significantly > 0.5  -> non-ramp structure RESPECTS the mirror  (F1-leaning)
    R_obs ~ 0.5                -> mirror-agnostic; F1 adds nothing beyond the ramp
    R_obs significantly < 0.5  -> non-ramp structure BREAKS the mirror     (anti-F1)
  Honest verdict needs BOTH direction and significance (n_pairs is small: 4 for
  periods 2-3, 9 for 4-5, 16 for 6; the discrete null floor is ~2^-n_pairs).
  Pool light (2,3) vs heavy (4,5,6) for power, since the effect is period-dependent.
"""
import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from periodic_palindrome import periods, allen_en_periods  # real NIST data, reused

RNG = np.random.RandomState(0)
B = 40000


def residual_pairs(v):
    """Return per-mirror-pair antisym/symm energies of the post-ramp residual."""
    v = np.asarray(v, float)
    n = len(v)
    assert n % 2 == 0
    w = v - v.mean()
    k = np.arange(n, dtype=float)
    ell = k - (n - 1) / 2.0
    a = (w - w[::-1]) / 2.0
    a_lin = (a @ ell) / (ell @ ell) * ell
    r = w - a_lin                      # residual: a_non + s  (ramp removed)
    m = n // 2
    rk, rp = r[:m], r[n - 1::-1][:m]   # partners (k, N-1-k)
    A = (rk - rp) ** 2 / 2.0           # antisym energy per pair (F1-respecting)
    S = (rk + rp) ** 2 / 2.0           # symm    energy per pair (F1-breaking)
    return A, S


def signflip_null(A, S):
    """R_obs and sign-flip null distribution of R = E_non/(E_non+E_sym)."""
    Eres = float((A + S).sum())
    if Eres < 1e-12:
        return None
    R_obs = float(A.sum() / Eres)
    bits = RNG.randint(0, 2, size=(B, len(A))).astype(bool)
    anti = np.where(bits, A, S).sum(axis=1)     # flip swaps A<->S per pair
    R_null = anti / Eres
    p_resp = float(np.mean(R_null >= R_obs))     # one-sided: F1-respecting
    p_anti = float(np.mean(R_null <= R_obs))     # one-sided: anti-F1
    return R_obs, p_resp, p_anti, Eres


def verdict(R, p_resp, p_anti):
    if R > 0.5:
        return f"F1-leaning (p_resp={p_resp:.3f})" + ("  *" if p_resp < 0.05 else "")
    if R < 0.5:
        return f"ANTI-F1   (p_anti={p_anti:.3f})" + ("  *" if p_anti < 0.05 else "")
    return "mirror-agnostic"


def run(title, table, key, pool_groups=None):
    print("=" * 98)
    print(title)
    print("=" * 98)
    store = {}
    for name, p in table.items():
        A, S = residual_pairs(p[key])
        store[name] = (A, S)
        res = signflip_null(A, S)
        if res is None:
            print(f"  {name:<24s} no residual structure (pure ramp)")
            continue
        R, pr, pa, Eres = res
        print(f"  {name:<24s} pairs={len(A):2d}  E_resid={Eres:7.3f}  R={R:.3f}  -> {verdict(R, pr, pa)}")
    if pool_groups:
        print(f"  {'-'*92}")
        for label, names in pool_groups.items():
            A = np.concatenate([store[n][0] for n in names])
            S = np.concatenate([store[n][1] for n in names])
            R, pr, pa, Eres = signflip_null(A, S)
            print(f"  POOL {label:<19s} pairs={len(A):2d}  E_resid={Eres:7.3f}  R={R:.3f}  -> {verdict(R, pr, pa)}")
    print()


def controls():
    print("=" * 98)
    print("CALIBRATION CONTROLS (synthetic, N=8) -- null must centre at R=0.5 and flag both directions")
    print("=" * 98)
    n = 8
    k = np.arange(n, dtype=float)
    bump_anti = np.sin(2 * np.pi * k / n)
    bump_sym = np.cos(2 * np.pi * k / n)
    for name, v in {
        "pure linear ramp": 10 + 0.5 * k,
        "ramp + antisym bump (F1)": 10 + 0.5 * k + 0.8 * bump_anti,
        "ramp + symm bump (anti-F1)": 10 + 0.5 * k + 0.8 * bump_sym,
    }.items():
        A, S = residual_pairs(v)
        res = signflip_null(A, S)
        if res is None:
            print(f"  {name:<28s} no residual structure (pure ramp)  [correct]")
        else:
            R, pr, pa, _ = res
            print(f"  {name:<28s} R={R:.3f}  -> {verdict(R, pr, pa)}")
    print()


if __name__ == "__main__":
    controls()
    run("HARDENED F1 GATE on real NIST first ionization energies (eV), periods 2-6",
        periods, 'IE',
        pool_groups={"light (2,3)": list(periods)[:2], "heavy (4,5,6)": list(periods)[2:]})
    run("Secondary: Allen 1989 electronegativity (incl. noble gases)",
        allen_en_periods, 'EN')
