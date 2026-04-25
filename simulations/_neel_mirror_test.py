#!/usr/bin/env python3
"""Néel-zu-Spiegel test: Z⊗N-partnership of antiferromagnetic xpatterns.

The xpattern |+−+−+⟩ is the antiferromagnetic Néel state of the Heisenberg-XXZ
chain (in X-basis). Its Z⊗N-mirror |−+−+−⟩ is what a global π/2 Z-magnet-puls
would produce. Z⊗N anti-commutes with X_l, Y_l (each Z X Z = −X) and commutes
with Z_l. So Z⊗N is preserved by any operator with an even number of X/Y per
term: XXZ, ZZ-detuning, Z-dephasing, and T1 (σ⁻σ⁺ pairs). It is broken only
by transverse fields (single X or Y).

The right Néel order is M_X = (1/N) Σ (-1)^l ⟨X_l⟩ (anti-invariant under Z⊗N:
M_X(b) = −M_X(a) when symmetry holds). M_AB = (1/N) Σ (-1)^l ⟨Z_l⟩ is Z⊗N-
invariant and identically zero on X-basis states; wrong observable.

Tests:
1. Pure Heisenberg + Z: Z⊗N exact, MI(a) = MI(b) and M_X(b) = −M_X(a).
2. + T1: Z⊗N still preserved (σ⁻σ⁺ pairs cancel the sign).
3. + Z-detuning (uniform and non-uniform): Z⊗N still preserved (Z commutes).
4. + transverse X-field: Z⊗N breaks. MI(a) ≠ MI(b), and M_X(a)+M_X(b) ≠ 0.
"""
import math
import sys

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SM = np.array([[0, 1], [0, 0]], dtype=complex)
SP = SM.conj().T

PLUS = np.array([1, 1], dtype=complex) / math.sqrt(2)
MINUS = np.array([1, -1], dtype=complex) / math.sqrt(2)


def kron_n(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(N, i, op):
    ops = [I2] * N
    ops[i] = op
    return kron_n(*ops)


def xpattern_state(signs):
    """Build product state ⊗_l |s_l⟩ where s_l is + or -.

    signs is a list of +1/-1 of length N.
    """
    state = (PLUS if signs[0] > 0 else MINUS).copy()
    for s in signs[1:]:
        next_qubit = PLUS if s > 0 else MINUS
        state = np.kron(state, next_qubit)
    return state


def build_H_xxz(N, J, delta):
    """Open-chain Heisenberg-XXZ: H = (J/2) Σ (X_iX_{i+1} + Y_iY_{i+1} + δ Z_iZ_{i+1})."""
    D = 2**N
    H = np.zeros((D, D), dtype=complex)
    for i in range(N - 1):
        Xi, Yi, Zi = site_op(N, i, X), site_op(N, i, Y), site_op(N, i, Z)
        Xj, Yj, Zj = site_op(N, i+1, X), site_op(N, i+1, Y), site_op(N, i+1, Z)
        H = H + 0.5 * J * (Xi @ Xj + Yi @ Yj + delta * Zi @ Zj)
    return H


def build_H_detuning(N, deltas):
    """Site-detuning: H = Σ_l δ_l Z_l. Mini-magnetic-field profile."""
    D = 2**N
    H = np.zeros((D, D), dtype=complex)
    for l in range(N):
        H = H + deltas[l] * site_op(N, l, Z)
    return H


def lindblad_rhs(rho, H, gamma_z, gamma_t1, Z_ops, sm_ops, sp_ops):
    out = -1j * (H @ rho - rho @ H)
    for Zl in Z_ops:
        out = out + gamma_z * (Zl @ rho @ Zl - rho)
    if gamma_t1 > 0:
        for sm, sp in zip(sm_ops, sp_ops):
            spsm = sp @ sm
            out = out + gamma_t1 * (sm @ rho @ sp - 0.5 * (spsm @ rho + rho @ spsm))
    return out


def rk4_step(rho, dt, H, gamma_z, gamma_t1, Z_ops, sm_ops, sp_ops):
    k1 = lindblad_rhs(rho, H, gamma_z, gamma_t1, Z_ops, sm_ops, sp_ops)
    k2 = lindblad_rhs(rho + 0.5 * dt * k1, H, gamma_z, gamma_t1, Z_ops, sm_ops, sp_ops)
    k3 = lindblad_rhs(rho + 0.5 * dt * k2, H, gamma_z, gamma_t1, Z_ops, sm_ops, sp_ops)
    k4 = lindblad_rhs(rho + dt * k3, H, gamma_z, gamma_t1, Z_ops, sm_ops, sp_ops)
    return rho + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)


def partial_trace_keep_sites(rho, sites, N):
    sites = tuple(sorted(sites))
    shape_2N = [2] * (2 * N)
    out = rho.reshape(shape_2N)
    ket_axes = list(range(N))
    bra_axes = list(range(N, 2 * N))
    for j in range(N - 1, -1, -1):
        if j in sites:
            continue
        a_k = ket_axes[j]
        a_b = bra_axes[j]
        out = np.trace(out, axis1=a_k, axis2=a_b)
        lo, hi = sorted((a_k, a_b))
        for k in range(N):
            if k == j:
                continue
            if ket_axes[k] > hi:
                ket_axes[k] -= 2
            elif ket_axes[k] > lo:
                ket_axes[k] -= 1
            if bra_axes[k] > hi:
                bra_axes[k] -= 2
            elif bra_axes[k] > lo:
                bra_axes[k] -= 1
    k = len(sites)
    d_sub = 2**k
    remaining_axes_ket = [ket_axes[i] for i in sites]
    remaining_axes_bra = [bra_axes[i] for i in sites]
    current_order = remaining_axes_ket + remaining_axes_bra
    perm = np.argsort(current_order)
    out = np.transpose(out, perm)
    return out.reshape(d_sub, d_sub)


def von_neumann(rho):
    evs = np.linalg.eigvalsh(rho).real
    return float(-sum(ev * math.log2(ev) for ev in evs if ev > 1e-15))


def MI_pair(rho, a, b, N):
    rho_ab = partial_trace_keep_sites(rho, (a, b), N)
    rho_a = partial_trace_keep_sites(rho, (a,), N)
    rho_b = partial_trace_keep_sites(rho, (b,), N)
    return von_neumann(rho_a) + von_neumann(rho_b) - von_neumann(rho_ab)


def neel_x(rho, N):
    """M_X = (1/N) Σ (-1)^l ⟨X_l⟩. X-basis Néel order parameter (Z⊗N-anti-invariant)."""
    M = 0.0
    for l in range(N):
        Xl = site_op(N, l, X)
        ev = float(np.trace(rho @ Xl).real)
        M = M + ((-1)**l) * ev
    return M / N


def run_pair(N, signs_a, signs_b, H, gamma_z, gamma_t1, t_max, n_steps):
    """Evolve both xpatterns and return MI(0, N-1) and M_X trajectories."""
    Z_ops = [site_op(N, i, Z) for i in range(N)]
    sm_ops = [site_op(N, i, SM) for i in range(N)]
    sp_ops = [site_op(N, i, SP) for i in range(N)]

    psi_a = xpattern_state(signs_a)
    psi_b = xpattern_state(signs_b)
    rho_a = np.outer(psi_a, psi_a.conj())
    rho_b = np.outer(psi_b, psi_b.conj())

    dt = t_max / n_steps
    times = np.empty(n_steps + 1)
    mi_a = np.empty(n_steps + 1)
    mi_b = np.empty(n_steps + 1)
    M_a = np.empty(n_steps + 1)
    M_b = np.empty(n_steps + 1)
    for step in range(n_steps + 1):
        times[step] = step * dt
        mi_a[step] = MI_pair(rho_a, 0, N - 1, N)
        mi_b[step] = MI_pair(rho_b, 0, N - 1, N)
        M_a[step] = neel_x(rho_a, N)
        M_b[step] = neel_x(rho_b, N)
        if step < n_steps:
            rho_a = rk4_step(rho_a, dt, H, gamma_z, gamma_t1, Z_ops, sm_ops, sp_ops)
            rho_b = rk4_step(rho_b, dt, H, gamma_z, gamma_t1, Z_ops, sm_ops, sp_ops)
    return times, mi_a, mi_b, M_a, M_b


def build_H_transverse_x(N, h_l):
    """Transverse X-field: H = Σ_l h_l X_l. Each term anti-commutes with Z⊗N."""
    D = 2**N
    H = np.zeros((D, D), dtype=complex)
    for l in range(N):
        H = H + h_l[l] * site_op(N, l, X)
    return H


def main():
    N = 5
    t_max = 3.0
    n_steps = 200
    gamma_z = 0.1
    delta_xxz = 1.0   # full Heisenberg (XX+YY+ZZ)

    signs_neel = [+1, -1, +1, -1, +1]   # |+−+−+⟩
    signs_mirror = [-s for s in signs_neel]   # |−+−+−⟩

    print("=" * 78)
    print("Néel-zu-Spiegel test: Z⊗N partnership of xpatterns")
    print(f"N={N}, Heisenberg Δ={delta_xxz}, t_max={t_max}, γ_z={gamma_z}")
    print("Initial states: |+−+−+⟩ vs |−+−+−⟩ (Z⊗N-partners)")
    print("Z⊗N-invariant probe: MI(0,N-1).  Z⊗N-anti-invariant probe: M_X.")
    print("Symmetry preserved ⇒ MI(a)=MI(b) AND M_X(a)+M_X(b)=0.")
    print("=" * 78)

    H_heis = build_H_xxz(N, J=1.0, delta=delta_xxz)

    def report(label, mi_a, mi_b, M_a, M_b):
        d_mi = float(np.max(np.abs(mi_a - mi_b)))
        d_mx_anti = float(np.max(np.abs(M_a + M_b)))   # nonzero iff Z⊗N broken
        print(f"  {label}")
        print(f"    max |MI(a) − MI(b)|        = {d_mi:.3e}")
        print(f"    max |M_X(a) + M_X(b)|      = {d_mx_anti:.3e}  (Z⊗N-break signature)")
        print(f"    M_X(t=0):   a={M_a[0]:+.4f}, b={M_b[0]:+.4f}  (expect ±1)")
        print(f"    M_X(t=end): a={M_a[-1]:+.4f}, b={M_b[-1]:+.4f}")

    # Test 1: pure Z-dephasing (no T1, no detuning). Z⊗N exact.
    print("\nTest 1: pure Heisenberg + uniform Z-dephasing (Z⊗N exact)")
    print("-" * 78)
    times, mi_a, mi_b, M_a, M_b = run_pair(N, signs_neel, signs_mirror, H_heis,
                                            gamma_z, gamma_t1=0, t_max=t_max,
                                            n_steps=n_steps)
    report("baseline", mi_a, mi_b, M_a, M_b)

    # Test 2: + T1 amplitude damping. T1 has σ⁻σ⁺ pairs which preserve Z⊗N.
    print("\nTest 2: + T1 amplitude damping (σ⁻σ⁺ pairs ⇒ even X+Y count, Z⊗N preserved)")
    print("-" * 78)
    for ratio in [0.5, 1.0, 2.0]:
        gamma_t1 = ratio * gamma_z
        times, mi_a, mi_b, M_a, M_b = run_pair(N, signs_neel, signs_mirror, H_heis,
                                                gamma_z, gamma_t1=gamma_t1,
                                                t_max=t_max, n_steps=n_steps)
        d_mi = float(np.max(np.abs(mi_a - mi_b)))
        d_anti = float(np.max(np.abs(M_a + M_b)))
        print(f"  γ_T1/γ_z = {ratio:.1f}: |ΔMI| = {d_mi:.3e}, |M_X(a)+M_X(b)| = {d_anti:.3e}")

    # Test 3: + non-uniform Z-detuning. Z commutes with Z⊗N, so this still preserves it.
    print("\nTest 3: + non-uniform Z-detuning δ_l Z_l (Z commutes with Z⊗N, preserved)")
    print("-" * 78)
    rng = np.random.default_rng(20260425)
    deltas = rng.uniform(-0.3, 0.3, size=N)
    H_z_detune = H_heis + build_H_detuning(N, deltas)
    print(f"  δ_l = {[f'{d:+.3f}' for d in deltas]}")
    times, mi_a, mi_b, M_a, M_b = run_pair(N, signs_neel, signs_mirror, H_z_detune,
                                            gamma_z, gamma_t1=0, t_max=t_max,
                                            n_steps=n_steps)
    d_mi = float(np.max(np.abs(mi_a - mi_b)))
    d_anti = float(np.max(np.abs(M_a + M_b)))
    print(f"  |ΔMI| = {d_mi:.3e}, |M_X(a)+M_X(b)| = {d_anti:.3e}  (still preserved)")

    # Test 4: + transverse X-field (uniform). X anti-commutes with Z⊗N: BREAKS Z⊗N.
    print("\nTest 4: + uniform transverse X-field h Σ X_l (X anti-commutes ⇒ Z⊗N BROKEN)")
    print("-" * 78)
    for h_x in [0.05, 0.1, 0.2]:
        H_tx = H_heis + build_H_transverse_x(N, [h_x] * N)
        times, mi_a, mi_b, M_a, M_b = run_pair(N, signs_neel, signs_mirror, H_tx,
                                                gamma_z, gamma_t1=0, t_max=t_max,
                                                n_steps=n_steps)
        d_mi = float(np.max(np.abs(mi_a - mi_b)))
        d_anti = float(np.max(np.abs(M_a + M_b)))
        print(f"  h_x = {h_x:.2f}: |ΔMI| = {d_mi:.3e}, |M_X(a)+M_X(b)| = {d_anti:.3e}  (break)")

    # Test 5: + non-uniform transverse X-field (random profile). Same family, asymmetric.
    print("\nTest 5: + non-uniform transverse X-field h_l X_l (random profile)")
    print("-" * 78)
    h_l = rng.uniform(-0.1, 0.1, size=N)
    H_tx_nu = H_heis + build_H_transverse_x(N, h_l)
    print(f"  h_l = {[f'{h:+.3f}' for h in h_l]}")
    times, mi_a, mi_b, M_a, M_b = run_pair(N, signs_neel, signs_mirror, H_tx_nu,
                                            gamma_z, gamma_t1=0, t_max=t_max,
                                            n_steps=n_steps)
    d_mi = float(np.max(np.abs(mi_a - mi_b)))
    d_anti = float(np.max(np.abs(M_a + M_b)))
    print(f"  |ΔMI| = {d_mi:.3e}, |M_X(a)+M_X(b)| = {d_anti:.3e}")

    print("\nSummary:")
    print("  Z⊗N is preserved by ANY operator with even X+Y count per term:")
    print("    - XXZ Heisenberg (XX, YY, ZZ all even)")
    print("    - Z-dephasing (Z), Z-detuning uniform OR non-uniform (Z commutes)")
    print("    - T1 amplitude damping (σ⁻σ⁺ pairs cancel signs)")
    print("  Z⊗N is BROKEN only by operators with odd X+Y count:")
    print("    - Transverse X-field h_l X_l (single X per site)")
    print("    - Single transverse Y-field h_l Y_l")
    print()
    print("  Physical reading: bipartite K (sublattice gauge) is broken by non-uniform Z;")
    print("  Z⊗N (X-basis Néel mirror) is broken only by transverse magnetic fields.")
    print("  Two distinct symmetries detect two distinct asymmetries in the same chain.")


if __name__ == "__main__":
    main()
