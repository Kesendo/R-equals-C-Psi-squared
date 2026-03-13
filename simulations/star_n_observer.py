"""
N-OBSERVER STAR TOPOLOGY
=========================
Extend star_topology_v2.py from 3 to N qubits.
S (qubit 0) coupled to observers 1..N-1.
Track ALL pair CΨ maxima over full trajectory.

Key finding: AB crossing SURVIVES at N=4 and N=5.
The agent claim of "scaling failure" was wrong
(concurrence bug: eigvalsh vs eigvals).

See: experiments/STAR_TOPOLOGY_OBSERVERS.md
"""
import numpy as np
from itertools import combinations

# Pauli matrices
I2 = np.eye(2, dtype=complex)
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)

def kron_n(ops):
    """Tensor product of list of operators."""
    result = ops[0]
    for op in ops[1:]:
        result = np.kron(result, op)
    return result

def star_hamiltonian_n(N, J_list):
    """
    Star Hamiltonian: S (qubit 0) coupled to observers 1..N-1.
    J_list[i] = coupling strength S <-> observer i+1.
    len(J_list) = N-1.
    """
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for obs_idx in range(N-1):
        J = J_list[obs_idx]
        for sigma in [sx, sy, sz]:
            ops = [I2]*N
            ops[0] = sigma
            ops[obs_idx+1] = sigma
            H += J * kron_n(ops)
    return H

def dephasing_ops_n(N, gamma_list):
    """Local sigma_z dephasing. gamma_list[i] = gamma for qubit i."""
    L_ops = []
    for i in range(N):
        if gamma_list[i] > 0:
            ops = [I2]*N
            ops[i] = sz
            L_ops.append(np.sqrt(gamma_list[i]) * kron_n(ops))
    return L_ops

def lindblad_rhs(rho, H, L_ops):
    drho = -1j * (H @ rho - rho @ H)
    for L in L_ops:
        Ld = L.conj().T
        drho += L @ rho @ Ld - 0.5*(Ld @ L @ rho + rho @ Ld @ L)
    return drho

def rk4_step(rho, H, L_ops, dt):
    k1 = lindblad_rhs(rho, H, L_ops)
    k2 = lindblad_rhs(rho + 0.5*dt*k1, H, L_ops)
    k3 = lindblad_rhs(rho + 0.5*dt*k2, H, L_ops)
    k4 = lindblad_rhs(rho + dt*k3, H, L_ops)
    return rho + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

def ptrace_pair(rho, N, keep):
    """Partial trace: keep pair (i,j), trace out rest."""
    d_keep = 4
    d_trace = 2**N // d_keep
    rho_r = rho.reshape([2]*2*N)
    kept = list(keep)
    traced = [q for q in range(N) if q not in kept]
    bra_order = kept + traced
    ket_order = [x+N for x in kept] + [x+N for x in traced]
    new_order = bra_order + ket_order
    rho_r = rho_r.transpose(new_order)
    rho_r = rho_r.reshape(d_keep, d_trace, d_keep, d_trace)
    rho_pair = np.trace(rho_r, axis1=1, axis2=3)
    return rho_pair

def concurrence_2q(rho):
    """Wootters concurrence for 2-qubit density matrix.
    CRITICAL: uses eigvals, NOT eigvalsh. R is not Hermitian."""
    sy2 = np.kron(sy, sy)
    rho_tilde = sy2 @ rho.conj() @ sy2
    R = rho @ rho_tilde
    ev = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(R), 0))))[::-1]
    return max(0, ev[0] - ev[1] - ev[2] - ev[3])

def purity(rho):
    return np.real(np.trace(rho @ rho))

def l1_coherence(rho):
    return np.sum(np.abs(rho)) - np.real(np.trace(rho))

def psi_norm(rho):
    d = rho.shape[0]
    l1 = l1_coherence(rho)
    return l1 / (d - 1) if d > 1 else 0

def make_bell_SA_plus_rest(N):
    """Bell state between S(0) and A(1), |+> for all others."""
    bell = np.zeros(4, dtype=complex)
    bell[0] = 1/np.sqrt(2)
    bell[3] = 1/np.sqrt(2)
    plus = np.array([1,1], dtype=complex)/np.sqrt(2)
    psi = bell
    for _ in range(N-2):
        psi = np.kron(psi, plus)
    return np.outer(psi, psi.conj())

def run_n_observer(N, J_list, gamma_list, dt=0.005, t_max=5.0):
    """Run N-qubit star topology, track all pair CΨ maxima."""
    H = star_hamiltonian_n(N, J_list)
    L = dephasing_ops_n(N, gamma_list)
    rho = make_bell_SA_plus_rest(N)
    steps = int(t_max/dt)
    pairs = list(combinations(range(N), 2))
    results = {}
    for p in pairs:
        name = f"({p[0]},{p[1]})"
        results[name] = {'cpsi_max': 0, 'crosses': False,
                         'first_cross': None, 'window': 0,
                         'above': False, 't_start': 0}
    sample = max(1, int(0.02/dt))
    for step in range(steps+1):
        t = step * dt
        if step % sample == 0:
            for p in pairs:
                name = f"({p[0]},{p[1]})"
                rp = ptrace_pair(rho, N, p)
                c = concurrence_2q(rp)
                pn = psi_norm(rp)
                cpsi = c * pn
                if cpsi > results[name]['cpsi_max']:
                    results[name]['cpsi_max'] = cpsi
                if cpsi >= 0.25 and not results[name]['above']:
                    results[name]['above'] = True
                    results[name]['crosses'] = True
                    results[name]['t_start'] = t
                    if results[name]['first_cross'] is None:
                        results[name]['first_cross'] = t
                elif cpsi < 0.25 and results[name]['above']:
                    results[name]['above'] = False
                    results[name]['window'] += t - results[name]['t_start']
        if step < steps:
            rho = rk4_step(rho, H, L, dt)
    for p in pairs:
        name = f"({p[0]},{p[1]})"
        if results[name]['above']:
            results[name]['window'] += t_max - results[name]['t_start']
    return results

if __name__ == "__main__":
    for N, J_list, label in [
        (3, [1.0, 2.0], "S+A+B"),
        (4, [1.0, 2.0, 1.0], "S+A+B+C"),
        (5, [1.0, 2.0, 1.0, 1.0], "S+A+B+C+D"),
    ]:
        dt = 0.005 if N <= 4 else 0.01
        print(f"\n{'='*60}")
        print(f"  N={N}: {label}")
        print(f"{'='*60}")
        r = run_n_observer(N, J_list, [0.05]*N, dt=dt)
        for name, data in sorted(r.items()):
            cross = "YES" if data['crosses'] else "NO"
            fc = f"t={data['first_cross']:.2f}" if data['first_cross'] else "---"
            print(f"  {name}: CPsi max={data['cpsi_max']:.4f} {cross:>4} {fc:>8} w={data['window']:.3f}")
