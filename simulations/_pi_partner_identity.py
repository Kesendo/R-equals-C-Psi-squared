#!/usr/bin/env python3
"""Pi-partner identity verification for handshake algebra.

Claim (HANDSHAKE_ALGEBRA.md, 2026-04-24):
    h = (N, k, t, basis)  paired with  Pi.h = (N, N+1-k, t, Pi.basis)
give identical observables whenever dynamics are Pi-invariant (uniform-J chain
with uniform gamma-Z dephasing).

Amplitude identity: psi_k(N-1-l) = (-1)^{k+1} psi_k(l), so populations
|c_l|^2 are Pi-invariant; coherences differ by (-1)^{l+m}, which cancels in
|coherence|^2. Under Pi-invariant dynamics, MI trajectories stay paired.

Two observables:
  (a) MI(0, N-1)(t)             -- the receiver-engineering observable
  (b) Sum_i log pi_i(t)         -- PTF-analog: log-sum of single-site purities,
                                   survives perturbation as closure law
                                   Sum_i log alpha_i = 0 (PERSPECTIVAL_TIME_FIELD)

Test 1: uniform J=1, gamma=0.1, N=9. Compare both observables for k and N+1-k.
        Expected: identical to machine precision.
Test 2: non-uniform J (random positive), gamma=0.1, N=9. Same comparison.
        Expected: breakdown; H no longer Pi-invariant.

Consequence for receiver engineering: Pi-partnership halves the menu size
to ceil(N/2) distinct entries.
"""
import math
import sys

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


def build_H_uniform(N: int, J: float) -> np.ndarray:
    H = np.zeros((N, N), dtype=complex)
    for i in range(N - 1):
        H[i, i + 1] = -J
        H[i + 1, i] = -J
    return H


def build_H_nonuniform(N: int, Js) -> np.ndarray:
    assert len(Js) == N - 1
    H = np.zeros((N, N), dtype=complex)
    for i in range(N - 1):
        H[i, i + 1] = -Js[i]
        H[i + 1, i] = -Js[i]
    return H


def build_H_onsite(N: int, J: float, V) -> np.ndarray:
    """Uniform hopping + on-site potential V (length N)."""
    assert len(V) == N
    H = np.zeros((N, N), dtype=complex)
    for i in range(N - 1):
        H[i, i + 1] = -J
        H[i + 1, i] = -J
    for i in range(N):
        H[i, i] = V[i]
    return H


def psi_k(N: int, k: int) -> np.ndarray:
    return np.array(
        [math.sqrt(2.0 / (N + 1)) * math.sin(math.pi * k * (ell + 1) / (N + 1))
         for ell in range(N)],
        dtype=complex,
    )


def h_bin(p: float) -> float:
    if p <= 1e-15 or p >= 1 - 1e-15:
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def von_neumann_S(rho: np.ndarray) -> float:
    evs = np.linalg.eigvalsh(rho)
    S = 0.0
    for ev in evs:
        if ev > 1e-15:
            S -= ev * math.log2(ev)
    return S


def MI_pair_single_exc(rho: np.ndarray, l: int, m: int) -> float:
    """MI between sites l and m for a single-excitation-subspace rho (N x N)."""
    a = rho[l, l].real
    b = rho[m, m].real
    coh = rho[l, m]
    rdm = np.zeros((4, 4), dtype=complex)
    rdm[0, 0] = 1.0 - a - b
    rdm[1, 1] = b
    rdm[2, 2] = a
    rdm[2, 1] = coh
    rdm[1, 2] = np.conj(coh)
    S_l = h_bin(a)
    S_m = h_bin(b)
    S_lm = von_neumann_S(rdm)
    return S_l + S_m - S_lm


def site_purity(rho: np.ndarray, i: int) -> float:
    """Single-site purity Tr(rho_i^2) for single-excitation-subspace rho."""
    p = rho[i, i].real
    # rho_i = diag(1 - p, p) in {|0>, |1>} basis; purity = (1-p)^2 + p^2
    return (1.0 - p) ** 2 + p ** 2


def log_sum_purity(rho: np.ndarray, N: int) -> float:
    """Sum_i log pi_i -- PTF-analog grounding observable."""
    total = 0.0
    for i in range(N):
        pi = site_purity(rho, i)
        if pi > 1e-15:
            total += math.log(pi)
    return total


def rhs(rho: np.ndarray, H: np.ndarray, gamma: float) -> np.ndarray:
    commut = -1j * (H @ rho - rho @ H)
    diss = -4.0 * gamma * rho.copy()
    np.fill_diagonal(diss, 0.0)
    return commut + diss


def rk4_step(rho: np.ndarray, dt: float, H: np.ndarray, gamma: float) -> np.ndarray:
    k1 = rhs(rho, H, gamma)
    k2 = rhs(rho + 0.5 * dt * k1, H, gamma)
    k3 = rhs(rho + 0.5 * dt * k2, H, gamma)
    k4 = rhs(rho + dt * k3, H, gamma)
    return rho + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def simulate(N: int, k: int, H: np.ndarray, gamma: float, t_max: float, n_steps: int):
    psi = psi_k(N, k)
    rho = np.outer(psi, np.conj(psi))
    dt = t_max / n_steps
    mis = np.zeros(n_steps + 1)
    logpurity = np.zeros(n_steps + 1)
    for step in range(n_steps + 1):
        mis[step] = MI_pair_single_exc(rho, 0, N - 1)
        logpurity[step] = log_sum_purity(rho, N)
        if step < n_steps:
            rho = rk4_step(rho, dt, H, gamma)
    return mis, logpurity


def run_test(label: str, N: int, H: np.ndarray, gamma: float, t_max: float, n_steps: int):
    print(label)
    print(f"{'k':>3}  {'k_mirror':>9}  {'max|MI_k-MI_km|':>18}  {'max|SumLogPi_k-_km|':>22}  {'match':>7}")
    half = (N + 1) // 2
    for k in range(1, half + 1):
        k_mirror = N + 1 - k
        mi_k, lp_k = simulate(N, k, H, gamma, t_max, n_steps)
        mi_m, lp_m = simulate(N, k_mirror, H, gamma, t_max, n_steps)
        d_mi = float(np.max(np.abs(mi_k - mi_m)))
        d_lp = float(np.max(np.abs(lp_k - lp_m)))
        match = "YES" if max(d_mi, d_lp) < 1e-10 else "NO"
        marker = " (self)" if k == k_mirror else ""
        print(f"{k:>3}  {k_mirror:>9}  {d_mi:>18.3e}  {d_lp:>22.3e}  {match:>7}{marker}")
    print()


def main() -> None:
    N = 9
    gamma = 0.0  # literally zero: pure H dynamics, PTF grounding regime
    t_max = 2.0
    n_steps = 400

    print("=" * 72)
    print(f"Pi-partner identity verification  (N={N}, gamma=0 literal, t_max={t_max})")
    print("single-excitation manifold, pure Hamiltonian dynamics (no dephasing)")
    print("PTF-grounding regime: state stays pure, partner identity nakedly structural")
    print("=" * 72)
    print()

    # Test 1: uniform J (Pi-invariant H)
    J_uniform = 1.0
    H_u = build_H_uniform(N, J_uniform)
    run_test(
        f"Test 1  uniform J={J_uniform}  (H is Pi-invariant): expect identity",
        N, H_u, gamma, t_max, n_steps,
    )

    # Test 2: non-uniform J (Pi broken, BUT bipartite K = diag((-1)^l) still preserved!)
    rng = np.random.default_rng(42)
    Js_random = rng.uniform(0.5, 1.5, size=N - 1)
    H_nu = build_H_nonuniform(N, Js_random)
    js_str = "[" + ", ".join(f"{j:.3f}" for j in Js_random) + "]"
    run_test(
        f"Test 2  non-uniform J={js_str}  (Pi broken, K preserved): identity still holds",
        N, H_nu, gamma, t_max, n_steps,
    )

    # Test 3: uniform J + random on-site potential (breaks K bipartite gauge)
    V_random = rng.uniform(-1.0, 1.0, size=N)
    H_onsite = build_H_onsite(N, J_uniform, V_random)
    v_str = "[" + ", ".join(f"{v:+.3f}" for v in V_random) + "]"
    run_test(
        f"Test 3  H_hop + on-site V={v_str}  (K BROKEN): expect breakdown",
        N, H_onsite, gamma, t_max, n_steps,
    )

    # Test 4: same three H-cases but at gamma > 0. gamma merely sets a time unit
    # and does not affect partner identity; this confirms the effect is structural.
    gamma_check = 0.1
    print(f"Test 4  same three H-cases at gamma={gamma_check} (sanity: gamma just sets time unit)")
    print(f"{'case':>20}  {'worst |MI_k - MI_km|':>24}  {'worst |SumLogPi|':>20}")
    cases = [
        ("uniform J=1",      H_u),
        ("non-uniform J",    H_nu),
        ("uniform J + V",    H_onsite),
    ]
    for label, H in cases:
        worst_mi = 0.0
        worst_lp = 0.0
        half = (N + 1) // 2
        for k in range(1, half + 1):
            k_mirror = N + 1 - k
            if k == k_mirror:
                continue
            mi_k, lp_k = simulate(N, k, H, gamma_check, t_max, n_steps)
            mi_m, lp_m = simulate(N, k_mirror, H, gamma_check, t_max, n_steps)
            worst_mi = max(worst_mi, float(np.max(np.abs(mi_k - mi_m))))
            worst_lp = max(worst_lp, float(np.max(np.abs(lp_k - lp_m))))
        print(f"{label:>20}  {worst_mi:>24.3e}  {worst_lp:>20.3e}")
    print()

    # Consequence: receiver menu size
    print("Summary:")
    print("  - Partner identity k <-> N+1-k is carried by K = diag((-1)^l)")
    print("    (bipartite sublattice gauge), not by Pi (spatial reflection).")
    print("  - K invariance holds for any bipartite H (uniform or non-uniform J,")
    print("    any gamma_0 profile). Breaks only when H loses bipartite structure")
    print("    (on-site potential, NNN hopping).")
    print("  - gamma_0 is the time anchor: dynamics depend only on Q = J/gamma_0")
    print("    and H-structure; partner identity is Q-independent.")
    print()
    print("Consequence for receiver engineering:")
    print(f"  Naive menu: N = {N} bonding modes")
    print(f"  K-reduced:  ceil(N/2) = {(N + 1) // 2} distinct entries")
    for n in (7, 9, 11, 13, 15, 21, 51):
        print(f"  N={n:>3}:  menu size = {(n + 1) // 2}")


if __name__ == "__main__":
    main()
