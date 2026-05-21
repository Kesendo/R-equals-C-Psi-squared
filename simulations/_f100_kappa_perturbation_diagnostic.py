#!/usr/bin/env python3
"""_f100_kappa_perturbation_diagnostic.py

First-order Lindblad perturbation data for F100's κ -- a look, not a derivation.

κ_b = ∂c₁(b)/∂s (s the anti-palindromic sweep amplitude, J = J_sym + s·J_anti)
has no closed form (PROOF_F100 κ-section). That section conjectures the large
N=5 κ-spread -- 143% relative spread of κ across J_sym profiles, against
76% / 62% at N=3 / 4 -- reflects a NEAR-DEGENERACY in the slow-mode sector of
L(J_sym): a small denominator λ_s − λ_s' in the first-order eigenvector-mixing
series

    δM_s = Σ_{s'≠s} [ ⟨W_s'|V_L|M_s⟩ / (λ_s − λ_s') ] · M_s'.

A near-degeneracy makes that denominator tiny, δM_s blows up, and κ = ∂c₁/∂s
with it. For the near-degeneracy to *explain* the spread it must itself be
J_sym-dependent: a small gap for some profiles at N=5, not others.

This diagnostic computes the raw first-order data, no α-fit, no IFT. For each
F100 case (N, J_sym) it builds L(J_sym) in the exact (N+1)²-dim popcount sector
block (the c₁-witness engine, where the biorthonormal left eigenvectors come
free from inv(V_R), sidestepping the PTF full-space biorthogonalisation
blocker), eigendecomposes, forms V_L = ∂L/∂s along J_anti, and reads off:
  - the Hellmann-Feynman shifts  δλ_s = ⟨W_s|V_L|M_s⟩,
  - the mixing coefficients      ⟨W_s'|V_L|M_s⟩ / (λ_s − λ_s'),
  - the minimal non-exact eigenvalue gap.
It then asks whether the near-degeneracy proxy (the largest mixing coefficient)
tracks |κ|, loaded from the committed c₁-witness JSON.

Result (run 2026-05-21): the conjecture is refuted. δλ_s along J_anti vanishes
for every mode (F92; κ is a pure eigenvector-rotation effect), but the
eigenvector-mixing coefficients are flat across N (34 / 31 / 35 % J_sym-spread,
against κ's 76 / 62 / 143 %) and the N=5 near-degeneracy does not carry κ. The
first-order data does not capture κ's N=5 growth. Recorded in PROOF_F100's
κ-section ("The perturbation route, tested"). Re-run prints the verdict live.

Usage:
  python -u _f100_kappa_perturbation_diagnostic.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _f71_nonuniform_j_verification as f71w  # noqa: E402

GAMMA_0 = f71w.GAMMA_0
EXACT_DEGEN_TOL = 1e-9       # exact-degeneracy cutoff: smaller gap = kernel/symmetry pair,
                             # skipped in the near-degeneracy scan; also the |λ| kernel cutoff
KAPPA_DIR = f71w.RESULTS_DIR  # simulations/results/f71_nonuniform_j_verification/


# ---------------------------------------------------------------------------
def perturbation_data(N, j_sym):
    """First-order Lindblad data for L(J_sym) perturbed along J_anti.

    L is built in the exact (N+1)²-dim popcount sector block. V_L = ∂L/∂s along
    J_anti is the exact unit-step difference (L is affine in the J-profile).
    Returns the slow-first eigenvalue list, the minimal non-exact gap, the
    maximal first-order mixing coefficient and the slowness-rank pair carrying
    it, and the δλ spread."""
    S = f71w.kept_basis(N)
    anti = f71w.j_anti_direction(N)
    j_sym = np.asarray(j_sym, dtype=float)

    L0 = f71w.build_L_sub(list(j_sym), GAMMA_0, N, S)
    evals, V_R, V_Linv = f71w.eig_and_inv(L0)     # V_Linv @ V_R = I (biorthonormal)

    # V_L = dL/ds along J_anti. L is affine in the J-profile, so the unit-step
    # difference is exactly the derivative; reuse L0 = L(J_sym).
    V_L = f71w.build_L_sub(list(j_sym + anti), GAMMA_0, N, S) - L0

    G = V_Linv @ V_L @ V_R                        # G[s',s] = ⟨W_s'|V_L|M_s⟩

    order = np.argsort(-evals.real)               # slow-first (Re λ nearest 0)
    ev = evals[order]
    Gs = G[np.ix_(order, order)]

    dlam = np.diag(Gs)                            # δλ_s Hellmann-Feynman shifts
    m = len(ev)
    min_gap = np.inf
    max_mix = 0.0
    max_mix_pair = (-1, -1)
    for s in range(m):
        for sp in range(m):
            if s == sp:
                continue
            gap = abs(ev[s] - ev[sp])
            if gap < EXACT_DEGEN_TOL:
                continue                          # exact degeneracy, not "near"
            min_gap = min(min_gap, gap)
            mix = abs(Gs[sp, s] / (ev[s] - ev[sp]))
            if mix > max_mix:
                max_mix = mix
                max_mix_pair = (s, sp)
    return {
        "evals_slow_first": ev,
        "min_nonexact_gap": float(min_gap),
        "max_mixing": float(max_mix),
        "max_mixing_pair": max_mix_pair,
        "dlam_abs_max": float(np.max(np.abs(dlam))),
        "n_kernel": int(np.sum(np.abs(evals) < EXACT_DEGEN_TOL)),
    }


def load_kappa(N):
    """Per-J_sym max |κ| and the overall J_sym κ-spread, from the c₁-witness
    JSON results/f71_nonuniform_j_verification/f71_nonuniform_j_N{N}.json."""
    path = KAPPA_DIR / f"f71_nonuniform_j_N{N}.json"
    if not path.exists():
        return None, None
    data = json.loads(path.read_text(encoding="utf-8"))
    by_sym = {}
    for r in data["analysis"]["records"]:
        by_sym.setdefault(r["sym"], []).append(abs(r["kappa"]))
    kappa_max = {sym: max(ks) for sym, ks in by_sym.items()}
    spread = data["analysis"]["criterion_d_max_rel_kappa_spread"]
    return kappa_max, spread


def rel_spread(values):
    """(max − min) / max of a list of magnitudes; 0 if degenerate."""
    vals = [abs(v) for v in values if np.isfinite(v)]
    if len(vals) < 2 or max(vals) == 0.0:
        return 0.0
    return (max(vals) - min(vals)) / max(vals)


# ---------------------------------------------------------------------------
def main():
    print("=" * 78, flush=True)
    print("F100 κ -- first-order Lindblad perturbation diagnostic", flush=True)
    print("Conjecture (PROOF_F100 κ-section): the N=5 κ-spread is a near-",
          flush=True)
    print("degeneracy in the slow-mode sector of L(J_sym).", flush=True)
    print(f"gamma_0 = {GAMMA_0}", flush=True)
    print("=" * 78, flush=True)

    per_N = {}            # N -> list of (name, perturbation_data, max|κ|)
    for N in (3, 4, 5):
        kappa_max, spread = load_kappa(N)
        shapes = f71w.j_sym_profiles(N)
        head = f"\nN = {N}   sector dim {(N + 1) ** 2}"
        if spread is not None:
            head += f"   c₁-witness κ J_sym-spread = {spread * 100:.0f}%"
        print(head, flush=True)
        print("-" * 78, flush=True)
        print(f"  {'J_sym profile':18} {'min gap':>11} {'max mixing':>12} "
              f"{'|δλ|max':>10} {'max|κ| (c₁)':>12}", flush=True)
        rows = []
        for name, j_sym in shapes.items():
            d = perturbation_data(N, j_sym)
            km = kappa_max.get(name, float("nan")) if kappa_max else float("nan")
            rows.append((name, d, km))
            print(f"  {name:18} {d['min_nonexact_gap']:11.5f} "
                  f"{d['max_mixing']:12.3f} {d['dlam_abs_max']:10.4f} "
                  f"{km:12.4f}", flush=True)
        per_N[N] = rows
        ref = next((d for nm, d, _ in rows if nm == "uniform_1.0"), rows[0][1])
        slow = ref["evals_slow_first"][:8]
        print(f"  slowest 8 eigenvalues (uniform_1.0): "
              + "  ".join(f"{e.real:+.4f}{e.imag:+.4f}i" for e in slow),
              flush=True)
        print(f"  kernel modes (|λ| < {EXACT_DEGEN_TOL:.0e}): "
              f"{ref['n_kernel']}", flush=True)

    # --- verdict -----------------------------------------------------------
    print(f"\n{'=' * 78}", flush=True)
    print("VERDICT", flush=True)
    print("=" * 78, flush=True)
    print(f"  {'N':>3} {'κ J_sym-spread':>16} {'max-mixing spread':>20} "
          f"{'min gap (min)':>15}", flush=True)
    pooled_mix, pooled_kappa = [], []
    for N in (3, 4, 5):
        rows = per_N[N]
        _, sp = load_kappa(N)
        mixes = [d["max_mixing"] for _, d, _ in rows]
        gaps = [d["min_nonexact_gap"] for _, d, _ in rows]
        mix_spread = rel_spread(mixes)
        for _, d, km in rows:
            if np.isfinite(km):
                pooled_mix.append(d["max_mixing"])
                pooled_kappa.append(km)
        sp_str = f"{sp * 100:.0f}%" if sp is not None else "n/a"
        print(f"  {N:>3} {sp_str:>16} {mix_spread * 100:>19.0f}% "
              f"{min(gaps):>15.5f}", flush=True)

    if len(pooled_mix) >= 3:
        corr = float(np.corrcoef(pooled_mix, pooled_kappa)[0, 1])
        print(f"\n  pooled correlation  max-mixing vs max|κ|  "
              f"(over {len(pooled_mix)} cases): r = {corr:+.3f}", flush=True)
        print("  conjecture supported iff: max-mixing tracks |κ| (r near +1)",
              flush=True)
        print("  AND its J_sym-spread jumps at N=5 like the κ-spread "
              "(76 / 62 / 143).", flush=True)
    print("=" * 78, flush=True)


if __name__ == "__main__":
    main()
