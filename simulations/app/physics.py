"""
Physics engine for the Five Regulator Simulator.
3-qubit open Heisenberg star (S-A-B) with Lindblad dynamics.
"""
import numpy as np
from scipy.fft import fft, fftfreq

# === Pauli matrices & identity ===
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def tensor(*args):
    """Tensor (Kronecker) product of multiple matrices or vectors."""
    result = args[0]
    for m in args[1:]:
        result = np.kron(result, m)
    return result


def op_on_qubit(op, q, n=3):
    """Place single-qubit operator on qubit q (0-indexed) in n-qubit system."""
    ops = [I2] * n
    ops[q] = op
    return tensor(*ops)


def two_body(op_a, op_b, qa, qb, n=3):
    """Two-body operator: op_a on qubit qa, op_b on qubit qb."""
    ops = [I2] * n
    ops[qa] = op_a
    ops[qb] = op_b
    return tensor(*ops)


def heisenberg_pair(J, qa, qb, xy_ratio=1.0, n=3):
    """J * (xy_ratio*(XX + YY) + ZZ) between qubits qa, qb."""
    H = J * xy_ratio * (two_body(sx, sx, qa, qb, n) +
                         two_body(sy, sy, qa, qb, n))
    H += J * two_body(sz, sz, qa, qb, n)
    return H


# === Hamiltonian ===

def star_hamiltonian(J_SA, J_SB, xy_ratio=1.0):
    """3-qubit star: S(0)–A(1) and S(0)–B(2) couplings.
    xy_ratio=1: isotropic Heisenberg, xy_ratio=0: pure ZZ (Ising).
    """
    return heisenberg_pair(J_SA, 0, 1, xy_ratio) + \
           heisenberg_pair(J_SB, 0, 2, xy_ratio)


# === Jump operators ===

def build_jump_ops(gamma, eta=0.0, phi=0.0):
    """Lindblad jump operators.
    gamma: base noise strength
    eta: 0=local only, 1=fully correlated (A-B bath)
    phi: bath basis angle (0=ZZ, pi/2=XX)
    """
    ops = []

    # S qubit: full local dephasing
    ops.append(np.sqrt(gamma) * op_on_qubit(sz, 0))

    # A, B: mix of local + correlated
    local_g = (1 - eta) * gamma
    if local_g > 1e-15:
        s = np.sqrt(local_g)
        ops.append(s * op_on_qubit(sz, 1))
        ops.append(s * op_on_qubit(sz, 2))

    if eta > 1e-15:
        s = np.sqrt(eta * gamma / 2)
        bath_op = np.cos(phi) * sz + np.sin(phi) * sx
        # Collective operator on A + B
        ops.append(s * (op_on_qubit(bath_op, 1) + op_on_qubit(bath_op, 2)))

    return ops


# === Lindblad integrator ===

def precompute_lindblad(L_ops):
    """Precompute (L, L†, L†L) triples for the simulation loop."""
    return [(L, L.conj().T, L.conj().T @ L) for L in L_ops]


def lindblad_rhs(rho, H, L_triples):
    """dρ/dt from Lindblad master equation. L_triples = [(L, L†, L†L), ...]."""
    drho = -1j * (H @ rho - rho @ H)
    for L, Ld, LdL in L_triples:
        drho += L @ rho @ Ld - 0.5 * (LdL @ rho + rho @ LdL)
    return drho


def rk4_step(rho, H, L_triples, dt):
    """One RK4 step with numerical hygiene."""
    k1 = lindblad_rhs(rho, H, L_triples) * dt
    k2 = lindblad_rhs(rho + k1 / 2, H, L_triples) * dt
    k3 = lindblad_rhs(rho + k2 / 2, H, L_triples) * dt
    k4 = lindblad_rhs(rho + k3, H, L_triples) * dt
    rho_new = rho + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    # Hermiticity + trace normalization
    rho_new = (rho_new + rho_new.conj().T) / 2
    rho_new /= np.trace(rho_new)
    return rho_new


# === Partial trace ===

def partial_trace_keep(rho, keep, n_qubits=3):
    """Partial trace keeping specified qubits."""
    dims = [2] * n_qubits
    rho_r = rho.reshape(dims + dims)
    trace_out = sorted(set(range(n_qubits)) - set(keep))
    for offset, q in enumerate(trace_out):
        ax_ket = q - offset
        ax_bra = ax_ket + (n_qubits - offset)
        rho_r = np.trace(rho_r, axis1=ax_ket, axis2=ax_bra)
    d_keep = 2 ** len(keep)
    return rho_r.reshape(d_keep, d_keep)


# === Initial states ===

def _density(psi):
    """Pure-state density matrix from state vector."""
    return np.outer(psi, psi.conj())


def bell_sa_plus_b():
    """Bell(SA) ⊗ |+>_B = (|00⟩+|11⟩)/√2 ⊗ (|0⟩+|1⟩)/√2"""
    psi = np.zeros(8, dtype=complex)
    psi[0b000] = 1  # |000⟩
    psi[0b001] = 1  # |001⟩
    psi[0b110] = 1  # |110⟩
    psi[0b111] = 1  # |111⟩
    return _density(psi / np.linalg.norm(psi))


def w_state():
    """W state: (|001⟩ + |010⟩ + |100⟩)/√3"""
    psi = np.zeros(8, dtype=complex)
    psi[0b001] = 1
    psi[0b010] = 1
    psi[0b100] = 1
    return _density(psi / np.linalg.norm(psi))


def product_plus():
    """|+⟩⊗|+⟩⊗|+⟩"""
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, plus), plus)
    return _density(psi)


def ghz_state():
    """GHZ: (|000⟩ + |111⟩)/√2"""
    psi = np.zeros(8, dtype=complex)
    psi[0] = 1
    psi[7] = 1
    return _density(psi / np.linalg.norm(psi))


INITIAL_STATES = {
    "Bell_SA ⊗ |+>_B": bell_sa_plus_b,
    "W-state": w_state,
    "Product |+++>": product_plus,
    "GHZ": ghz_state,
}


# === Observables ===

# Precompute 4x4 operators in AB subspace
_YZ = np.kron(sy, sz)
_ZY = np.kron(sz, sy)
_XX = np.kron(sx, sx)
_YY = np.kron(sy, sy)


def compute_cplusminus(rho_ab):
    """c+ = (⟨YZ⟩ + ⟨ZY⟩)/√2,  c- = (⟨YZ⟩ - ⟨ZY⟩)/√2"""
    yz = np.real(np.trace(rho_ab @ _YZ))
    zy = np.real(np.trace(rho_ab @ _ZY))
    return (yz + zy) / np.sqrt(2), (yz - zy) / np.sqrt(2)


def xx_commutator_norm(rho_ab):
    """||[ρ_AB, X⊗X]||_F"""
    comm = rho_ab @ _XX - _XX @ rho_ab
    return np.linalg.norm(comm, 'fro')


def concurrence_2q(rho_ab):
    """Wootters concurrence for a 2-qubit density matrix."""
    rho_tilde = _YY @ rho_ab.conj() @ _YY
    R = rho_ab @ rho_tilde
    ev = np.sort(np.sqrt(np.maximum(np.real(np.linalg.eigvals(R)), 0)))[::-1]
    return max(0.0, ev[0] - ev[1] - ev[2] - ev[3])


def l1_coherence(rho):
    """Off-diagonal L1 norm."""
    d = rho.shape[0]
    return np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))


def psi_norm(rho):
    """Normalized coherence Ψ = l1/(d-1)."""
    d = rho.shape[0]
    return l1_coherence(rho) / (d - 1)


# === FFT analysis ===

def fft_analysis(signal, dt):
    """Dominant frequency, amplitude, and spectrum from time series."""
    n = len(signal)
    centered = signal - np.mean(signal)
    freqs = fftfreq(n, d=dt)
    spectrum = np.abs(fft(centered))
    pos = freqs > 0
    pf = freqs[pos]
    ps = spectrum[pos] * 2 / n
    if len(ps) == 0:
        return 0.0, 0.0, pf, ps
    idx = np.argmax(ps)
    return pf[idx], ps[idx], pf, ps


# === Main simulation ===

def run_simulation(J_SA, J_SB, xy_ratio, gamma, eta, phi,
                   initial_state_name, dt=0.02, t_max=20.0):
    """Run full Lindblad simulation. Returns results dict."""

    H = star_hamiltonian(J_SA, J_SB, xy_ratio)
    L_triples = precompute_lindblad(build_jump_ops(gamma, eta, phi))

    if initial_state_name not in INITIAL_STATES:
        raise ValueError(f"Unknown state: {initial_state_name!r}. "
                         f"Choose from: {list(INITIAL_STATES)}")
    rho = INITIAL_STATES[initial_state_name]()

    n_steps = int(t_max / dt)

    times = np.zeros(n_steps + 1)
    c_plus = np.zeros(n_steps + 1)
    c_minus = np.zeros(n_steps + 1)
    xx_comms = np.zeros(n_steps + 1)
    concs = np.zeros(n_steps + 1)
    psis = np.zeros(n_steps + 1)

    # Incremental skeleton: accumulate rho_sum and purity_sum
    d_ab = 4
    rho_ab_sum = np.zeros((d_ab, d_ab), dtype=complex)
    purity_sum = 0.0

    for i in range(n_steps + 1):
        times[i] = i * dt
        rho_ab = partial_trace_keep(rho, [1, 2])
        c_plus[i], c_minus[i] = compute_cplusminus(rho_ab)
        xx_comms[i] = xx_commutator_norm(rho_ab)
        concs[i] = concurrence_2q(rho_ab)
        psis[i] = psi_norm(rho_ab)

        rho_ab_sum += rho_ab
        purity_sum += np.real(np.trace(rho_ab @ rho_ab))

        if i < n_steps:
            rho = rk4_step(rho, H, L_triples, dt)

    # FFT
    f_cp, a_cp, fft_f, fft_cp = fft_analysis(c_plus, dt)
    f_cm, a_cm, _, fft_cm = fft_analysis(c_minus, dt)

    # Amplitude ratio
    amp_ratio = a_cp / a_cm if a_cm > 1e-10 else float('inf')

    # XX symmetry (mean of last 10% for robustness)
    xx_final = np.mean(xx_comms[-max(1, len(xx_comms) // 10):])

    # Skeleton fraction (incremental)
    n_samples = n_steps + 1
    rho_avg = rho_ab_sum / n_samples
    p_avg = np.real(np.trace(rho_avg @ rho_avg))
    p_mean = purity_sum / n_samples
    skel = p_avg / p_mean if p_mean > 1e-15 else 0.0

    return {
        'times': times,
        'c_plus': c_plus,
        'c_minus': c_minus,
        'concurrence': concs,
        'psi': psis,
        'cpsi': concs * psis,
        'f_cp': f_cp, 'f_cm': f_cm,
        'a_cp': a_cp, 'a_cm': a_cm,
        'amp_ratio': amp_ratio,
        'fft_freqs': fft_f,
        'fft_cp': fft_cp,
        'fft_cm': fft_cm,
        'xx_comms': xx_comms,
        'xx_final': xx_final,
        'xx_sym': xx_final < 1e-6,
        'skeleton': skel,
    }
