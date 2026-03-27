"""
Thermal Emergence Experiment
=============================
What emerges from heat?

Part A: Do unpaired modes die faster? (decay rate comparison)
Part B: Time evolution - does palindromic fraction of activity increase?
Part C: Thermal bath only (no Hamiltonian) - can heat create palindromic structure?
Part D: Thermal bath + coupling - does heat amplify palindromic modes?
Part E: Effect of heat on waves - dephasing vs thermal excitation

Not for publication. Exploratory.
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


def classify_modes(evals, tol_rel=1e-4):
    """Classify palindromic pairing by searching for best center S.

    For Z-dephasing Lindbladian, the mirror formula is:
      partner of lambda = -(lambda + 2*Sg)
    where Sg = N*gamma. But for general Lindblad operators, we
    search for the center empirically.
    """
    decays = -evals.real
    n = len(evals)

    # Try palindromic pairing with empirical center
    # Center = (max_decay + min_decay) / 2
    S = (np.max(decays) + np.min(decays)) / 2

    paired = np.zeros(n, dtype=bool)
    used = set()

    for k in range(n):
        if k in used:
            continue
        target_decay = 2 * S - decays[k]
        target_freq = -evals[k].imag  # mirror frequency

        best_j = -1
        best_dist = np.inf
        for j in range(n):
            if j in used or j == k:
                continue
            d_decay = abs(decays[j] - target_decay)
            d_freq = abs(evals[j].imag - target_freq)
            dist = d_decay + d_freq
            if dist < best_dist:
                best_dist = dist
                best_j = j

        tol = tol_rel * max(1, abs(2 * S))
        if best_j >= 0 and best_dist < tol:
            paired[k] = True
            paired[best_j] = True
            used.add(k)
            used.add(best_j)
        elif abs(decays[k] - S) < tol:
            # Self-paired at center
            paired[k] = True
            used.add(k)

    return paired, ~paired, S


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
print("Starting from |up,up,up>, how does the spectrum evolve?")
print("As heat dissipates, does palindromic structure emerge?")
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

# Initial state: |000> (all spin up) - far from equilibrium
rho0 = np.zeros((d, d), dtype=complex)
rho0[0, 0] = 1.0
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
print(f"\n{'t':>6} {'|delta_rho|':>12} {'osc_power':>12} {'decay_frac':>12} {'comment':>20}")
print("-" * 65)

for t in [0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0]:
    rho_t_vec = expm(L * t) @ rho0_vec
    delta = rho_t_vec - rho_ss_vec

    # Total deviation from equilibrium
    norm_delta = np.sqrt(np.real(np.dot(delta.conj(), delta)))

    # Spectral content: FFT-like analysis via eigenmode projection
    # Direct: how oscillatory is rho(t)?
    rho_t = rho_t_vec.reshape(d, d)

    # Off-diagonal elements carry oscillation, diagonal carry population
    off_diag = rho_t.copy()
    np.fill_diagonal(off_diag, 0)
    osc_power = np.sum(np.abs(off_diag) ** 2)

    diag_dev = np.abs(np.diag(rho_t) - 1.0 / d)
    pop_dev = np.sum(diag_dev ** 2)

    total_dev = osc_power + pop_dev
    decay_frac = pop_dev / total_dev if total_dev > 1e-30 else 0

    # What remains?
    if norm_delta < 1e-6:
        comment = "equilibrium"
    elif decay_frac > 0.9:
        comment = "mostly population"
    elif decay_frac < 0.1:
        comment = "mostly coherence"
    else:
        comment = f"mixed"

    print(f"{t:6.1f} {norm_delta:12.6f} {osc_power:12.6f} {decay_frac:12.4f} {comment:>20}")


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

print(f"\n{'n_bar':>6} {'n_modes':>8} {'n_osc':>6} {'n_paired':>8} {'pair%':>7} {'max_freq':>9}")
print("-" * 55)

for n_bar in [0.0, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    L_ops = thermal_ops(Gamma, n_bar, N)
    L = build_liouvillian(H_zero, L_ops)
    evals = np.linalg.eigvals(L)

    ev = evals[np.abs(evals) > 1e-8]
    n_modes = len(ev)
    n_osc = np.sum(np.abs(ev.imag) > 1e-8)
    max_freq = np.max(np.abs(ev.imag)) if n_osc > 0 else 0

    paired, unpaired, S = classify_modes(ev)
    n_paired = np.sum(paired)
    pct = 100 * n_paired / n_modes if n_modes > 0 else 0

    print(f"{n_bar:6.2f} {n_modes:8d} {n_osc:6d} {n_paired:8d} {pct:7.1f} {max_freq:9.4f}")


# =================================================================
# PART D: Thermal bath + coupling
# =================================================================
print("\n\n" + "=" * 70)
print("PART D: Thermal bath + Heisenberg coupling")
print("Does heat + coupling create MORE palindromic structure than coupling alone?")
print("=" * 70)

N = 3
J = 1.0
gamma_deph = 0.1  # dephasing (baseline)
Gamma_th = 0.1    # thermal bath strength

print(f"\nBaseline: J={J}, dephasing gamma={gamma_deph}")
print(f"Adding thermal bath Gamma={Gamma_th}")

header = f"{'n_bar':>6} {'n_modes':>8} {'pair%':>7} {'Efreq_p%':>9} {'Efreq_tot':>10} {'Edecay_tot':>11} {'ratio_f/d':>10}"
print(f"\n{header}")
print("-" * len(header))

for n_bar in [0.0, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0]:
    H = heisenberg_H(N, J)
    L_ops = dephasing_ops(gamma_deph, N)
    L_ops += thermal_ops(Gamma_th, n_bar, N)
    L = build_liouvillian(H, L_ops)
    evals = np.linalg.eigvals(L)

    ev = evals[np.abs(evals) > 1e-8]
    n_modes = len(ev)

    paired, unpaired, S = classify_modes(ev)
    n_paired = np.sum(paired)
    pct = 100 * n_paired / n_modes if n_modes > 0 else 0

    freq = np.abs(ev.imag)
    decay = np.abs(ev.real)
    E_freq_paired = np.sum(freq[paired])
    E_freq_total = np.sum(freq)
    E_decay_total = np.sum(decay)
    fp = 100 * E_freq_paired / E_freq_total if E_freq_total > 1e-12 else 0
    ratio = E_freq_total / E_decay_total if E_decay_total > 1e-12 else 0

    print(f"{n_bar:6.2f} {n_modes:8d} {pct:7.1f} {fp:9.1f} {E_freq_total:10.4f} "
          f"{E_decay_total:11.4f} {ratio:10.4f}")


# =================================================================
# PART E: Effect of heat on waves
# =================================================================
print("\n\n" + "=" * 70)
print("PART E: Does heat affect waves?")
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
