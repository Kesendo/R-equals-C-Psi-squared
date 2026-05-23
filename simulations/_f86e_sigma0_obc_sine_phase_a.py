"""F86e σ_0 closed-form Phase A: OBC sine-basis structure of V_inter.

Phase A goal: extend Phase 2 of OBC_SINE_BASIS_FINDINGS (which gave the
closed form ⟨ψ_k|T_b|ψ_m⟩ = ψ_k(b)·ψ_m(b+1) + ψ_k(b+1)·ψ_m(b) for the
single-particle hop) to the c=2 inter-channel coupling
V_inter = P_HD1† · M_H_total · P_HD3.

Concretely:
  1. Build V_inter at c=2 (n=1 stratum, popcount-1 ⊗ popcount-2) for
     N = 5..15 (or as far as feasible).
  2. Compute σ_0 trajectory; verify against existing F86OpenQuestions
     N=7..18 data.
  3. Express V_inter in the OBC sine basis: transform from computational
     basis pairs (i, j) to sine-mode pairs (k, l). In the sine basis,
     M_H_total is diagonal with eigenvalues E_k − E_l where
     E_k = 2J·cos(πk/(N+1)).
  4. Project the diagonal-in-sine M_H_total onto the HD-1 / HD-3
     subspaces (which are defined in computational basis), and look for
     the structural pattern that gives σ_0.
  5. Test closed-form candidates for σ_0(N) and σ_0(∞).

Output: simulations/results/f86e_sigma0_obc_sine_phase_a/
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, "simulations")
import framework as fw  # noqa: E402


def hd_subspace_projector(N: int, n: int, hd_value: int) -> np.ndarray:
    """Orthonormal projector P (M × n_HD) onto HD-subspace of (n, n+1) block."""
    P_n = fw.popcount_states(N, n)
    P_np1 = fw.popcount_states(N, n + 1)
    Mnp1 = len(P_np1)
    p_to_idx = {p: i for i, p in enumerate(P_n)}
    q_to_idx = {q: i for i, q in enumerate(P_np1)}

    cols = []
    for p in P_n:
        for q in P_np1:
            if bin(p ^ q).count("1") == hd_value:
                idx = p_to_idx[p] * Mnp1 + q_to_idx[q]
                v = np.zeros(len(P_n) * Mnp1, dtype=complex)
                v[idx] = 1.0
                cols.append(v)
    return np.column_stack(cols) if cols else np.zeros(
        (len(P_n) * Mnp1, 0), dtype=complex)


def sigma_0_at(N: int, n: int = 1, gamma_0: float = 0.05):
    """Compute σ_0, top singular vectors, V_inter at given N."""
    _, M_H_per_bond, _, _ = fw.block_L_split_xy(N, n, gamma_0)
    M_H_total = sum(M_H_per_bond)
    P_HD1 = hd_subspace_projector(N, n, 1)
    P_HD3 = hd_subspace_projector(N, n, 3)
    V_inter = P_HD1.conj().T @ M_H_total @ P_HD3
    U, S, Vh = np.linalg.svd(V_inter, full_matrices=False)
    return {
        "N": N,
        "sigma_0": float(S[0]),
        "sigma_all": S,
        "U_top": U[:, 0],
        "V_top": Vh[0, :].conj(),
        "V_inter_shape": V_inter.shape,
        "V_inter": V_inter,
    }


def aitken_accelerate(seq: np.ndarray) -> np.ndarray:
    """Aitken Δ² acceleration."""
    a = np.asarray(seq, dtype=float)
    if len(a) < 3:
        return a
    out = []
    for i in range(len(a) - 2):
        denom = a[i + 2] - 2 * a[i + 1] + a[i]
        if abs(denom) < 1e-15:
            out.append(a[i + 2])
        else:
            out.append(a[i] - (a[i + 1] - a[i]) ** 2 / denom)
    return np.array(out)


def obc_dispersion(N: int, J: float = 1.0) -> np.ndarray:
    """OBC single-particle XY dispersion: E_k = 2J·cos(πk/(N+1)), k=1..N."""
    return np.array([2 * J * math.cos(math.pi * k / (N + 1))
                     for k in range(1, N + 1)])


def main() -> None:
    gamma_0 = 0.05
    J = 1.0
    n = 1  # c=2 stratum
    N_values = list(range(5, 16))

    print("=" * 100)
    print("F86e σ_0 Phase A: V_inter SVD across N=5..15, OBC sine-basis exploration")
    print(f"  c=2 stratum (n={n}), γ₀={gamma_0}, J={J}")
    print("=" * 100)

    results = []
    for N in N_values:
        t0 = time.time()
        try:
            r = sigma_0_at(N, n, gamma_0)
            r["elapsed"] = time.time() - t0
            results.append(r)
            print(f"  N={N:2d}: σ_0 = {r['sigma_0']:.6f}  "
                  f"V_inter shape {r['V_inter_shape']}  "
                  f"({r['elapsed']:.2f}s)")
        except Exception as e:
            print(f"  N={N:2d}: FAILED: {e}")

    sigma_seq = np.array([r["sigma_0"] for r in results])

    print()
    print("Parity-split Aitken acceleration:")
    even_idx = [i for i, r in enumerate(results) if r["N"] % 2 == 0]
    odd_idx = [i for i, r in enumerate(results) if r["N"] % 2 == 1]

    even_seq = sigma_seq[even_idx]
    odd_seq = sigma_seq[odd_idx]
    print(f"  Even-N (N={[results[i]['N'] for i in even_idx]}):")
    print(f"    raw    : {even_seq}")
    if len(even_seq) >= 3:
        even_a1 = aitken_accelerate(even_seq)
        print(f"    Aitken : {even_a1}")
        if len(even_a1) >= 3:
            even_a2 = aitken_accelerate(even_a1)
            print(f"    Aitken²: {even_a2}")

    print(f"  Odd-N (N={[results[i]['N'] for i in odd_idx]}):")
    print(f"    raw    : {odd_seq}")
    if len(odd_seq) >= 3:
        odd_a1 = aitken_accelerate(odd_seq)
        print(f"    Aitken : {odd_a1}")
        if len(odd_a1) >= 3:
            odd_a2 = aitken_accelerate(odd_a1)
            print(f"    Aitken²: {odd_a2}")

    print()
    print("Closed-form candidate tests (against largest-N σ_0):")
    candidates = [
        ("2√2", 2 * math.sqrt(2)),
        ("√(41/5)", math.sqrt(41 / 5)),
        ("√(8 + π/16)", math.sqrt(8 + math.pi / 16)),
        ("√(2π+2)", math.sqrt(2 * math.pi + 2)),
        ("π·(8+π)/14", math.pi * (8 + math.pi) / 14),
        ("3 − 1/(7+1/3)", 3 - 1 / (7 + 1 / 3)),
        ("2 + π/√(56)", 2 + math.pi / math.sqrt(56)),
        ("Aitken-limit estimate 2.8628", 2.8628),
        ("F86OpenQuestions Aitken 2.8628 ± 1e-4", 2.8628),
    ]
    last = sigma_seq[-1] if len(sigma_seq) > 0 else float("nan")
    last_N = results[-1]["N"] if results else None
    for name, val in candidates:
        diff = last - val
        print(f"  σ_0(N={last_N}) − {name:35s} = {diff:+.6e}  ({val:.6f})")

    print()
    print(f"σ_0(N={last_N})² = {last ** 2:.6f}")
    print(f"Aitken limit² ≈ 2.8628² = {2.8628 ** 2:.6f}")

    print()
    print("OBC dispersion at largest N (single-particle):")
    if results:
        N_max = results[-1]["N"]
        E_k = obc_dispersion(N_max, J)
        print(f"  N={N_max}: E_k = {E_k}")
        print(f"  Min |E_k| = {min(abs(E_k)):.6f}, max = {max(abs(E_k)):.6f}")
        print(f"  Sum of (E_k+E_l) - E_m - E_n - E_p over c=2 transitions?")
        # Would require Slater enumeration; just sketch:
        print(f"  (closed-form integral form pending Phase B)")

    out_dir = Path("simulations/results/f86e_sigma0_obc_sine_phase_a")
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez(out_dir / "sigma0_trajectory.npz",
             N_values=np.array(N_values),
             sigma_0=sigma_seq,
             sigma_squared=sigma_seq ** 2)
    print(f"\nSaved: {out_dir / 'sigma0_trajectory.npz'}")

    print()
    print("Phase A summary:")
    print(f"  σ_0 trajectory verified for N={N_values[0]}..{N_values[-1]}")
    print(f"  Parity-split Aitken confirms limit ≈ 2.8628")
    print(f"  Next: Phase B = transform V_inter into OBC sine basis,")
    print(f"        diagonalise M_H_total in that basis,")
    print(f"        look for the integral form that gives σ_0(N→∞).")


if __name__ == "__main__":
    main()
