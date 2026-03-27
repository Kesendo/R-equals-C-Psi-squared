"""
Thermal Emergence Experiment
=============================
What emerges from heat?

Part A: Decay rate comparison (universal 2x law)
Part B: Time evolution (structure survives, heat dies)
Part C: Thermal bath only (heat cannot create oscillation)
Part D: Effect of heat on waves (dephasing vs thermal excitation)

Findings 2 and 3 in ENERGY_PARTITION.md.

Script: simulations/thermal_emergence.py
Output: stdout (run with PYTHONIOENCODING=utf-8 on Windows)
"""
import numpy as np
from scipy.linalg import expm

# === Pauli matrices ===
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
sp = np.array([[0, 1], [0, 0]], dtype=complex)  # sigma+
sm = np.array([[0, 0], [1, 0]], dtype=complex)  # sigma-


def site_op(op, site, N):
    ops = [I2] * N
    ops[site] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def heisenberg_H(N, J=1.0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [sx, sy, sz]:
            H += J * site_op(P, i, N) @ site_op(P, i + 1, N)
    return H


def build_liouvillian(H, L_ops):
    """General Liouvillian with arbitrary Lindblad operators."""
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for Lop in L_ops:
        LdL = Lop.conj().T @ Lop
        L += np.kron(Lop.conj(), Lop)
        L -= 0.5 * np.kron(Id, LdL)
        L -= 0.5 * np.kron(LdL.T, Id)
    return L


def dephasing_ops(gamma, N):
    """Z-dephasing on each qubit."""
    return [np.sqrt(gamma) * site_op(sz, k, N) for k in range(N)]


def thermal_ops(Gamma, n_bar, N):
    """Thermal bath: emission + absorption on each qubit."""
    ops = []
    for k in range(N):
        if (1 + n_bar) > 0:
            ops.append(np.sqrt(Gamma * (1 + n_bar)) * site_op(sm, k, N))
        if n_bar > 0:
            ops.append(np.sqrt(Gamma * n_bar) * site_op(sp, k, N))
    return ops


def classify_modes_Sg(evals, Sg, tol=1e-6):
    """Classify using known Sg (for pure dephasing)."""
    n = len(evals)
    paired = np.zeros(n, dtype=bool)

    for k in range(n):
        target = -(evals[k] + 2 * Sg)
        dists = np.abs(evals - target)
        best = np.argmin(dists)
        if dists[best] < tol * max(1, abs(target)):
            paired[k] = True

    return paired, ~paired


# =================================================================
# PART A: Decay rate comparison (paired vs unpaired)
# =================================================================
print("=" * 70)
print("PART A: Decay rate comparison")
print("Do unpaired modes die faster than paired modes?")
print("=" * 70)

J = 1.0
gamma = 0.1

for N in [2, 3, 4, 5]:
    d = 2 ** N
    Sg = N * gamma
    H = heisenberg_H(N, J)
    L_ops = dephasing_ops(gamma, N)
    L = build_liouvillian(H, L_ops)
    evals = np.linalg.eigvals(L)

    ev = evals[np.abs(evals) > 1e-8]
    paired, unpaired = classify_modes_Sg(ev, Sg)

    decay_paired = np.abs(ev[paired].real)
    decay_unpaired = np.abs(ev[unpaired].real)
    freq_paired = np.abs(ev[paired].imag)
    freq_unpaired = np.abs(ev[unpaired].imag)

    print(f"\nN={N}  (paired: {np.sum(paired)}, unpaired: {np.sum(unpaired)})")
    if np.sum(paired) > 0:
        print(f"  Paired modes:")
        print(f"    Decay rates: min={np.min(decay_paired):.4f}, "
              f"max={np.max(decay_paired):.4f}, mean={np.mean(decay_paired):.4f}")
        print(f"    Frequencies: min={np.min(freq_paired):.4f}, "
              f"max={np.max(freq_paired):.4f}, mean={np.mean(freq_paired):.4f}")
    if np.sum(unpaired) > 0:
        print(f"  Unpaired modes:")
        print(f"    Decay rates: min={np.min(decay_unpaired):.4f}, "
              f"max={np.max(decay_unpaired):.4f}, mean={np.mean(decay_unpaired):.4f}")
        print(f"    Frequencies: min={np.min(freq_unpaired):.4f}, "
              f"max={np.max(freq_unpaired):.4f}, mean={np.mean(freq_unpaired):.4f}")
    if np.sum(paired) > 0 and np.sum(unpaired) > 0:
        ratio = np.mean(decay_unpaired) / np.mean(decay_paired)
        print(f"  Ratio (unpaired/paired mean decay): {ratio:.4f}")
        if ratio > 1:
            print(f"  --> Unpaired modes die {ratio:.1f}x FASTER")
        else:
            print(f"  --> Paired modes die {1/ratio:.1f}x FASTER")


# =================================================================
# PART B: Time evolution via matrix exponential
# =================================================================
print("\n\n" + "=" * 70)
print("PART B: Time evolution")
print("Starting from |down,up,up> (single excitation, NOT eigenstate)")
print("As heat dissipates, does palindromic structure survive?")
print("=" * 70)

N = 3
d = 2 ** N
J = 1.0
gamma = 0.1
Sg = N * gamma

H = heisenberg_H(N, J)
L_ops = dephasing_ops(gamma, N)
L = build_liouvillian(H, L_ops)

# Steady state: maximally mixed
rho_ss = np.eye(d, dtype=complex) / d
rho_ss_vec = rho_ss.flatten()

# Initial state: |100> (single excitation, not an eigenstate of H_XXX)
# |100> = index 4 in computational basis (1*4 + 0*2 + 0*1)
# The Heisenberg Hamiltonian hops this excitation, creating coherences.
rho0 = np.zeros((d, d), dtype=complex)
rho0[4, 4] = 1.0
rho0_vec = rho0.flatten()

# Eigenvalues for classification
evals_L = np.linalg.eigvals(L)
ev_nz = evals_L[np.abs(evals_L) > 1e-8]
paired_nz, unpaired_nz = classify_modes_Sg(ev_nz, Sg)

# Decay rates of paired vs unpaired
decay_p = np.abs(ev_nz[paired_nz].real)
decay_u = np.abs(ev_nz[unpaired_nz].real)
print(f"\nPaired mean decay:   {np.mean(decay_p):.4f} (slowest: {np.min(decay_p):.4f})")
print(f"Unpaired mean decay: {np.mean(decay_u):.4f} (slowest: {np.min(decay_u):.4f})")
print(f"Unpaired halflife:   {np.log(2)/np.min(decay_u):.2f}")
print(f"Slowest paired halflife: {np.log(2)/np.min(decay_p):.2f}")

# Time evolution: track distance from steady state and spectral content
print(f"\n{'t':>6} {'|delta_rho|':>12} {'coherence':>12} {'population':>12} "
      f"{'coh_frac':>10} {'comment':>20}")
print("-" * 75)

for t in [0, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0]:
    rho_t_vec = expm(L * t) @ rho0_vec
    delta = rho_t_vec - rho_ss_vec

    # Total deviation from equilibrium
    norm_delta = np.sqrt(np.real(np.dot(delta.conj(), delta)))

    rho_t = rho_t_vec.reshape(d, d)

    # Off-diagonal elements = coherences (oscillation carriers)
    off_diag = rho_t.copy()
    np.fill_diagonal(off_diag, 0)
    coh_power = np.sum(np.abs(off_diag) ** 2)

    # Diagonal deviation from equilibrium = population imbalance
    diag_dev = np.abs(np.diag(rho_t) - 1.0 / d)
    pop_dev = np.sum(diag_dev ** 2)

    total_dev = coh_power + pop_dev
    coh_frac = coh_power / total_dev if total_dev > 1e-30 else 0

    if norm_delta < 1e-6:
        comment = "equilibrium"
    elif coh_frac > 0.5:
        comment = "coherence dominant"
    elif coh_frac < 0.01:
        comment = "population only"
    else:
        comment = f"mixed ({100*coh_frac:.0f}% coh)"

    print(f"{t:6.1f} {norm_delta:12.6f} {coh_power:12.8f} {pop_dev:12.8f} "
          f"{coh_frac:10.4f} {comment:>20}")


# =================================================================
# PART C: Thermal bath only (no Hamiltonian)
# =================================================================
print("\n\n" + "=" * 70)
print("PART C: Thermal bath only (J=0)")
print("Can heat alone create palindromic oscillation?")
print("=" * 70)

N = 3
d = 2 ** N
Gamma = 0.5  # bath coupling
H_zero = np.zeros((d, d), dtype=complex)

print(f"\n{'n_bar':>6} {'n_modes':>8} {'n_osc':>6} {'max_freq':>9} {'Efreq':>10}")
print("-" * 45)

for n_bar in [0.0, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    L_ops = thermal_ops(Gamma, n_bar, N)
    L = build_liouvillian(H_zero, L_ops)
    evals = np.linalg.eigvals(L)

    ev = evals[np.abs(evals) > 1e-8]
    n_modes = len(ev)
    n_osc = np.sum(np.abs(ev.imag) > 1e-8)
    max_freq = np.max(np.abs(ev.imag)) if n_osc > 0 else 0
    E_freq = np.sum(np.abs(ev.imag))

    print(f"{n_bar:6.2f} {n_modes:8d} {n_osc:6d} {max_freq:9.4f} {E_freq:10.4f}")


# =================================================================
# PART D: Effect of heat on waves
# =================================================================
print("\n\n" + "=" * 70)
print("PART D: Does heat affect waves?")
print("Compare: Z-dephasing (phase noise) vs thermal excitation (energy noise)")
print("=" * 70)

N = 3
J = 1.0

print("\n--- Pure Z-dephasing: increasing gamma ---")
header = f"{'gamma':>8} {'n_osc':>6} {'max_freq':>10} {'Efreq':>10} {'pair%':>7}"
print(f"\n{header}")
print("-" * len(header))

for gamma_val in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    Sg = N * gamma_val
    H = heisenberg_H(N, J)
    L_ops = dephasing_ops(gamma_val, N)
    L_mat = build_liouvillian(H, L_ops)
    evals = np.linalg.eigvals(L_mat)
    ev = evals[np.abs(evals) > 1e-8]

    paired, unpaired = classify_modes_Sg(ev, Sg)
    n_osc = np.sum(np.abs(ev.imag) > 1e-8)
    max_freq = np.max(np.abs(ev.imag))
    E_freq = np.sum(np.abs(ev.imag))
    pct = 100 * np.sum(paired) / len(ev)

    print(f"{gamma_val:8.2f} {n_osc:6d} {max_freq:10.4f} {E_freq:10.4f} {pct:7.1f}")


print("\n--- Thermal excitation only (no dephasing): increasing n_bar ---")
header = f"{'n_bar':>8} {'Gamma':>8} {'n_osc':>6} {'max_freq':>10} {'Efreq':>10}"
print(f"\n{header}")
print("-" * len(header))

Gamma_th = 0.1
for n_bar in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    H = heisenberg_H(N, J)
    L_ops = thermal_ops(Gamma_th, n_bar, N)
    L_mat = build_liouvillian(H, L_ops)
    evals = np.linalg.eigvals(L_mat)
    ev = evals[np.abs(evals) > 1e-8]

    n_osc = np.sum(np.abs(ev.imag) > 1e-8)
    max_freq = np.max(np.abs(ev.imag)) if n_osc > 0 else 0
    E_freq = np.sum(np.abs(ev.imag))

    print(f"{n_bar:8.2f} {Gamma_th:8.2f} {n_osc:6d} {max_freq:10.4f} {E_freq:10.4f}")


print("\n--- Dephasing + thermal excitation combined ---")
print("Baseline: gamma=0.1 (dephasing), adding thermal Gamma=0.1")
header = f"{'n_bar':>8} {'n_osc':>6} {'max_freq':>10} {'Efreq':>10} {'Edecay':>10} {'freq/decay':>10}"
print(f"\n{header}")
print("-" * len(header))

gamma_base = 0.1
Gamma_th = 0.1
for n_bar in [0.0, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    H = heisenberg_H(N, J)
    L_ops = dephasing_ops(gamma_base, N)
    L_ops += thermal_ops(Gamma_th, n_bar, N)
    L_mat = build_liouvillian(H, L_ops)
    evals = np.linalg.eigvals(L_mat)
    ev = evals[np.abs(evals) > 1e-8]

    n_osc = np.sum(np.abs(ev.imag) > 1e-8)
    max_freq = np.max(np.abs(ev.imag)) if n_osc > 0 else 0
    E_freq = np.sum(np.abs(ev.imag))
    E_decay = np.sum(np.abs(ev.real))
    ratio = E_freq / E_decay if E_decay > 1e-12 else 0

    print(f"{n_bar:8.2f} {n_osc:6d} {max_freq:10.4f} {E_freq:10.4f} {E_decay:10.4f} {ratio:10.4f}")


print("\n\nDone.")
