"""
R=CPsi^2 meets Quantum State Transfer
======================================
Bridge between our palindromic spectral results and the QST literature.

Questions:
1. How does our system compare to standard QST protocols?
2. Does the palindrome give us design advantages?
3. What predicts transfer quality: topology, coupling ratio, or spectral structure?

References:
- Bose, PRL 91, 207901 (2003) - original QST via spin chains
- Christandl et al., PRL 92, 187902 (2004) - perfect state transfer
- Wojcik et al., PRA 75, 022330 (2007) - weak end-coupling
- Kay, Int J Quantum Inf 8, 641 (2010) - review of PST

March 14, 2026
"""

import numpy as np
from scipy.linalg import expm
import sys, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex) 
Z = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)
minus = (up - dn) / np.sqrt(2)

def kron_list(ops):
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r

def build_H_general(N, bonds, couplings=None):
    """Build Heisenberg H for N qubits with specified bonds and optional per-bond couplings."""
    dim = 2**N
    H = np.zeros((dim, dim), dtype=complex)
    if couplings is None:
        couplings = [1.0] * len(bonds)
    for (i, j), J in zip(bonds, couplings):
        for p in [X, Y, Z]:
            ops = [I2]*N; ops[i] = p; ops[j] = p
            H += J * kron_list(ops)
    return H

def build_L(H, gammas, N):
    dim = 2**N; dim2 = dim**2; Id = np.eye(dim)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for i in range(N):
        g = gammas[i] if isinstance(gammas, (list, np.ndarray)) else gammas
        ops = [I2]*N; ops[i] = Z; Zi = kron_list(ops)
        L += g * (np.kron(Zi, Zi.conj()) - np.eye(dim2))
    return L

def ptrace(rho, N, keep):
    dim = 2**N; nk = len(keep); dk = 2**nk
    rho_r = np.zeros((dk, dk), dtype=complex)
    for i in range(dim):
        for j in range(dim):
            bi = [(i >> (N-1-k)) & 1 for k in range(N)]
            bj = [(j >> (N-1-k)) & 1 for k in range(N)]
            if all(bi[k] == bj[k] for k in range(N) if k not in keep):
                ki = sum(bi[keep[m]] << (nk-1-m) for m in range(nk))
                kj = sum(bj[keep[m]] << (nk-1-m) for m in range(nk))
                rho_r[ki, kj] += rho[i, j]
    return rho_r

def state_transfer_fidelity(psi_in, rho_out):
    """F = <psi|rho|psi> - standard QST fidelity."""
    return np.real(psi_in.conj() @ rho_out @ psi_in)

def avg_fidelity_6states(eLt, N, sender_qubit, receiver_qubit, idle_state=up):
    """Average fidelity over 6 cardinal states (MUB). Standard QST metric."""
    states = [up, dn, plus, minus, (up+1j*dn)/np.sqrt(2), (up-1j*dn)/np.sqrt(2)]
    fids = []
    for psi_A in states:
        # Build initial state: idle on all qubits except sender
        components = [idle_state] * N
        components[sender_qubit] = psi_A
        psi_full = components[0]
        for c in components[1:]:
            psi_full = np.kron(psi_full, c)
        rho0 = np.outer(psi_full, psi_full.conj())
        rho_t_vec = eLt @ rho0.flatten()
        rho_t = rho_t_vec.reshape(2**N, 2**N)
        rho_t = (rho_t + rho_t.conj().T) / 2
        rho_t /= np.trace(rho_t).real
        rho_B = ptrace(rho_t, N, [receiver_qubit])
        fids.append(state_transfer_fidelity(psi_A, rho_B))
    return np.mean(fids), fids

def check_palindrome(L, N, gamma):
    """Verify palindromic spectrum and return decay rate structure."""
    evals = np.linalg.eigvals(L)
    rates = -evals.real
    center = N * gamma
    # Check palindrome
    sorted_rates = np.sort(rates)
    n_paired = 0
    for r in sorted_rates:
        partner = 2*center - r
        if np.min(np.abs(sorted_rates - partner)) < 1e-6:
            n_paired += 1
    # Unique rates (rounded)
    unique_rates = sorted(set(np.round(rates, 6)))
    return n_paired == len(rates), unique_rates, center

# ================================================================
print("=" * 80)
print("R=CPsi^2 MEETS QUANTUM STATE TRANSFER")
print("=" * 80)

# ================================================================
print(f"\n{'='*80}")
print("TEST 1: Topology Comparison as QST Channels")
print("Same N=3, same gamma=0.05, different graphs")
print("Sender=qubit 1, Receiver=qubit 2, Mediator=qubit 0")
print(f"{'='*80}")

gamma = 0.05
N = 3
times = np.linspace(0.05, 5.0, 200)

topologies = {
    'star (0=hub)': {'bonds': [(0,1), (0,2)], 'couplings': [1.0, 1.0],
                     'sender': 1, 'receiver': 2},
    'chain 0-1-2':  {'bonds': [(0,1), (1,2)], 'couplings': [1.0, 1.0],
                     'sender': 0, 'receiver': 2},
    'triangle':     {'bonds': [(0,1), (1,2), (0,2)], 'couplings': [1.0, 1.0, 1.0],
                     'sender': 0, 'receiver': 2},
}

print(f"\n{'Topology':>20} {'Best F_avg':>10} {'t_opt':>6} {'Palindrome':>10} {'Center':>8} {'# rates':>8}")
print(f"{'-'*65}")

for name, topo in topologies.items():
    H = build_H_general(N, topo['bonds'], topo['couplings'])
    L = build_L(H, gamma, N)
    is_pal, unique_rates, center = check_palindrome(L, N, gamma)
    
    best_f = 0; best_t = 0
    for t in times:
        eLt = expm(L * t)
        f_avg, _ = avg_fidelity_6states(eLt, N, topo['sender'], topo['receiver'])
        if f_avg > best_f:
            best_f = f_avg; best_t = t
    
    print(f"{name:>20} {best_f:>10.4f} {best_t:>6.2f} {'YES' if is_pal else 'NO':>10} "
          f"{center:>8.3f} {len(unique_rates):>8}")

# ================================================================
print(f"\n{'='*80}")
print("TEST 2: Coupling Ratio Optimization (Wojcik-style)")
print("Star topology, J_SA fixed at 1.0, sweep J_SB")
print("Does the optimal ratio relate to palindromic rates?")
print(f"{'='*80}")

j_ratios = [0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0]
print(f"\n{'J_SB':>6} {'J_SB/J_SA':>10} {'Best F':>8} {'t_opt':>6} {'F(|0>)':>8} {'F(|1>)':>8} {'F(|+>)':>8}")
print(f"{'-'*60}")

for jsb in j_ratios:
    H = build_H_general(3, [(0,1), (0,2)], [1.0, jsb])
    L = build_L(H, gamma, 3)
    
    best_f = 0; best_t = 0; best_fids = None
    for t in times:
        eLt = expm(L * t)
        f_avg, fids = avg_fidelity_6states(eLt, 3, 1, 2)
        if f_avg > best_f:
            best_f = f_avg; best_t = t; best_fids = fids
    
    print(f"{jsb:>6.1f} {jsb/1.0:>10.2f} {best_f:>8.4f} {best_t:>6.2f} "
          f"{best_fids[0]:>8.4f} {best_fids[1]:>8.4f} {best_fids[2]:>8.4f}")

# ================================================================
print(f"\n{'='*80}")
print("TEST 3: Chain Length Scaling (standard QST benchmark)")
print("Linear chain, sender=0, receiver=N-1, all J=1, gamma=0.05")
print(f"{'='*80}")

print(f"\n{'N':>4} {'Best F_avg':>10} {'t_opt':>6} {'Palindrome':>10} {'# unique rates':>15}")
print(f"{'-'*50}")

for N_chain in [2, 3, 4, 5]:
    bonds = [(i, i+1) for i in range(N_chain-1)]
    H = build_H_general(N_chain, bonds)
    L = build_L(H, gamma, N_chain)
    is_pal, unique_rates, center = check_palindrome(L, N_chain, gamma)
    
    best_f = 0; best_t = 0
    t_range = np.linspace(0.05, 8.0, 200)
    for t in t_range:
        eLt = expm(L * t)
        f_avg, _ = avg_fidelity_6states(eLt, N_chain, 0, N_chain-1)
        if f_avg > best_f:
            best_f = f_avg; best_t = t
    
    print(f"{N_chain:>4} {best_f:>10.4f} {best_t:>6.2f} {'YES' if is_pal else 'NO':>10} {len(unique_rates):>15}")

# ================================================================
print(f"\n{'='*80}")
print("TEST 4: Palindromic Rate Decomposition as QST Diagnostic")
print("Which component of the transferred state degrades first?")
print("Star: S=0, A=1, B=2, J_SA=1, J_SB=2, gamma=0.05")
print(f"{'='*80}")

def concurrence_2q(rho):
    YY = np.kron(Y, Y)
    R = rho @ (YY @ rho.conj() @ YY)
    ev = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(R), 0))))[::-1]
    return max(0, ev[0] - ev[1] - ev[2] - ev[3])

H = build_H_general(3, [(0,1), (0,2)], [1.0, 2.0])
L = build_L(H, gamma, 3)

# Track all components over time for |+> input (coherent, most interesting)
psi_A = plus
psi_full = np.kron(np.kron(up, psi_A), up)
rho0 = np.outer(psi_full, psi_full.conj())
rho0_vec = rho0.flatten()

print(f"\nTime evolution of QST components for |+> input:")
print(f"{'t':>6} {'F_B':>8} {'C_SB':>8} {'Purity_B':>8} {'l1_B':>8} {'Conc_SB':>8} "
      f"{'CPsi_SB':>8} {'theta_SB':>8}")
print(f"{'-'*70}")

t_fine = np.linspace(0.05, 3.0, 120)
components = []

for t in t_fine:
    eLt = expm(L * t)
    rho_t_vec = eLt @ rho0_vec
    rho_t = rho_t_vec.reshape(8, 8)
    rho_t = (rho_t + rho_t.conj().T) / 2
    rho_t /= np.trace(rho_t).real
    
    # B (qubit 2) fidelity
    rho_B = ptrace(rho_t, 3, [2])
    F_B = state_transfer_fidelity(psi_A, rho_B)
    
    # SB pair (qubits 0,2)
    rho_SB = ptrace(rho_t, 3, [0, 2])
    conc_SB = concurrence_2q(rho_SB)
    
    # B purity
    purity_B = np.real(np.trace(rho_B @ rho_B))
    
    # B coherence (l1 norm)
    l1_B = np.sum(np.abs(rho_B)) - np.sum(np.abs(np.diag(rho_B)))
    
    # CPsi_SB = C * Psi where C = purity of SB, Psi = l1/(d-1) of SB
    purity_SB = np.real(np.trace(rho_SB @ rho_SB))
    l1_SB = np.sum(np.abs(rho_SB)) - np.sum(np.abs(np.diag(rho_SB)))
    cpsi_SB = purity_SB * l1_SB / 3
    
    theta = np.degrees(np.arctan(np.sqrt(max(0, 4*cpsi_SB - 1)))) if cpsi_SB > 0.25 else 0.0
    
    components.append({
        't': t, 'F_B': F_B, 'purity_B': purity_B, 'l1_B': l1_B,
        'conc_SB': conc_SB, 'cpsi_SB': cpsi_SB, 'theta_SB': theta
    })
    
    if len(components) % 6 == 1:
        print(f"{t:>6.2f} {F_B:>8.4f} {purity_SB:>8.4f} {purity_B:>8.4f} {l1_B:>8.4f} "
              f"{conc_SB:>8.4f} {cpsi_SB:>8.4f} {theta:>8.1f}")

# ================================================================
print(f"\n{'='*80}")
print("TEST 5: Noiseless vs Noisy QST — What Does Dephasing Cost?")
print("Same star, J_SA=1, J_SB=2, compare gamma=0 vs gamma=0.05")
print(f"{'='*80}")

print(f"\n{'gamma':>8} {'Best F':>8} {'t_opt':>6} {'F at t=1.33':>12} {'F loss at t_opt':>15}")
print(f"{'-'*55}")

for g in [0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]:
    H = build_H_general(3, [(0,1), (0,2)], [1.0, 2.0])
    L = build_L(H, g, 3)
    
    best_f = 0; best_t = 0; f_at_133 = 0
    for t in np.linspace(0.05, 5.0, 200):
        eLt = expm(L * t)
        f_avg, _ = avg_fidelity_6states(eLt, 3, 1, 2)
        if f_avg > best_f:
            best_f = f_avg; best_t = t
        if abs(t - 1.33) < 0.02:
            f_at_133 = f_avg
    
    loss = 1.0 - best_f if g == 0 else None
    noiseless_best = best_f if g == 0 else noiseless_best
    f_loss = noiseless_best - best_f
    
    print(f"{g:>8.3f} {best_f:>8.4f} {best_t:>6.2f} {f_at_133:>12.4f} {f_loss:>15.4f}")

# ================================================================
print(f"\n{'='*80}")
print("TEST 6: Can Palindromic Rates Predict Optimal Transfer Time?")
print(f"{'='*80}")
print("\nHypothesis: t_opt is related to pi/(4J_eff) where J_eff depends on coupling.")
print("The palindromic rates set the WINDOW, not the peak.")

for jsb in [0.5, 1.0, 1.5, 2.0, 3.0]:
    H = build_H_general(3, [(0,1), (0,2)], [1.0, jsb])
    L = build_L(H, gamma, 3)
    
    # Find t_opt
    best_f = 0; best_t = 0
    for t in np.linspace(0.05, 5.0, 200):
        eLt = expm(L * t)
        f_avg, _ = avg_fidelity_6states(eLt, 3, 1, 2)
        if f_avg > best_f:
            best_f = f_avg; best_t = t
    
    # Hamiltonian eigenvalues
    eig_H = np.sort(np.real(np.linalg.eigvalsh(H)))
    bohr_freqs = sorted(set(np.round(np.abs(np.diff(eig_H)), 4)))
    
    # Palindromic rates
    evals_L = np.linalg.eigvals(L)
    decay_rates = sorted(set(np.round(-evals_L.real, 6)))
    decay_rates = [r for r in decay_rates if r > 1e-8]
    
    # Window: time until fidelity drops below 2/3
    window_end = 0
    for t in np.linspace(best_t, 5.0, 100):
        eLt = expm(L * t)
        f_avg, _ = avg_fidelity_6states(eLt, 3, 1, 2)
        if f_avg < 2/3:
            window_end = t
            break
    
    window = window_end - best_t if window_end > 0 else float('inf')
    
    # Slowest palindromic rate predicts window?
    slowest_rate = min(decay_rates) if decay_rates else 0
    predicted_window = 1.0 / slowest_rate if slowest_rate > 0 else float('inf')
    
    print(f"\n  J_SB={jsb:.1f}: t_opt={best_t:.2f}, F={best_f:.4f}, "
          f"window={window:.2f}, 1/slowest_rate={predicted_window:.2f}")
    print(f"    Bohr freqs: {bohr_freqs[:4]}")
    print(f"    Decay rates: {decay_rates[:5]}")

# ================================================================
print(f"\n{'='*80}")
print("TEST 7: Mirror-Symmetric Couplings (Kay/Christandl PST condition)")
print("Does mirror-symmetric coupling imply palindromic spectrum?")
print("Test: chain with couplings [J1, J2, ..., J2, J1] (mirror)")
print(f"{'='*80}")

# N=4 chain: couplings [J1, J2, J1] (mirror symmetric)
# N=5 chain: couplings [J1, J2, J2, J1]
# Also test non-mirror: [J1, J2, J3]

configs = [
    ("N=4 mirror [1,2,1]", 4, [(0,1),(1,2),(2,3)], [1.0, 2.0, 1.0]),
    ("N=4 non-mirror [1,2,3]", 4, [(0,1),(1,2),(2,3)], [1.0, 2.0, 3.0]),
    ("N=4 uniform [1,1,1]", 4, [(0,1),(1,2),(2,3)], [1.0, 1.0, 1.0]),
    ("N=4 weak-end [0.5,1,0.5]", 4, [(0,1),(1,2),(2,3)], [0.5, 1.0, 0.5]),
    ("N=5 mirror [1,2,2,1]", 5, [(0,1),(1,2),(2,3),(3,4)], [1.0, 2.0, 2.0, 1.0]),
    ("N=5 PST-like [1,s2,s2,1]", 5, [(0,1),(1,2),(2,3),(3,4)], 
     [1.0, np.sqrt(2), np.sqrt(2), 1.0]),
]

print(f"\n{'Config':>30} {'Palindrome':>10} {'Best F(0->N-1)':>15} {'t_opt':>6}")
print(f"{'-'*65}")

for name, N_cfg, bonds, coups in configs:
    H = build_H_general(N_cfg, bonds, coups)
    L = build_L(H, gamma, N_cfg)
    is_pal, _, _ = check_palindrome(L, N_cfg, gamma)
    
    best_f = 0; best_t = 0
    for t in np.linspace(0.05, 8.0, 300):
        eLt = expm(L * t)
        f_avg, _ = avg_fidelity_6states(eLt, N_cfg, 0, N_cfg-1)
        if f_avg > best_f:
            best_f = f_avg; best_t = t
    
    print(f"{name:>30} {'YES' if is_pal else 'NO':>10} {best_f:>15.4f} {best_t:>6.2f}")

# ================================================================
print(f"\n{'='*80}")
print("TEST 8: Star vs Chain vs Wojcik for Same N")
print("N=4: Which topology gives the best A->B transfer?")
print("gamma=0.05, sweep coupling ratio")
print(f"{'='*80}")

print(f"\n{'Setup':>35} {'Best F':>8} {'t_opt':>6}")
print(f"{'-'*55}")

# Star: hub=0, A=1, B=2, extra leaf=3
for jsb in [1.0, 2.0]:
    H = build_H_general(4, [(0,1),(0,2),(0,3)], [1.0, jsb, 0.1])
    L = build_L(H, gamma, 4)
    best_f = 0; best_t = 0
    for t in np.linspace(0.05, 6.0, 200):
        eLt = expm(L * t)
        f_avg, _ = avg_fidelity_6states(eLt, 4, 1, 2)
        if f_avg > best_f: best_f = f_avg; best_t = t
    print(f"{'Star (weak 3rd leaf) J_SB='+str(jsb):>35} {best_f:>8.4f} {best_t:>6.2f}")

# Chain: 0-1-2-3, sender=0, receiver=3
for setup in [([1,1,1], "uniform"), ([0.5,1,0.5], "weak-end"), ([1,2,1], "mirror")]:
    coups, label = setup
    H = build_H_general(4, [(0,1),(1,2),(2,3)], coups)
    L = build_L(H, gamma, 4)
    best_f = 0; best_t = 0
    for t in np.linspace(0.05, 8.0, 200):
        eLt = expm(L * t)
        f_avg, _ = avg_fidelity_6states(eLt, 4, 0, 3)
        if f_avg > best_f: best_f = f_avg; best_t = t
    print(f"{'Chain '+label:>35} {best_f:>8.4f} {best_t:>6.2f}")

print(f"\n{'='*80}")
print("DONE")
print(f"{'='*80}")
