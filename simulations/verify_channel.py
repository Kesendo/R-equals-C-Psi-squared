"""
Verify GPT's One-Shot Channel Analysis
=======================================
GPT claims for 3-qubit star, J_SA=1.0, J_SB=2.0, gamma=0.05:
  - Average state-transfer fidelity A->B: 0.884 at t=1.30
  - Holevo lower bound (binary {|0>,|1>}): 0.529 bits/use
  - Coherent information lower bound: 0.185 qubits/use
  - Classical benchmark F_avg = 2/3 beaten on t in [0.60, 1.65]

March 14, 2026
"""

import numpy as np
from scipy.linalg import expm, logm
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

def build_H_star_asym(J_SA, J_SB):
    """3-qubit star: S=0, A=1, B=2. S-A coupling J_SA, S-B coupling J_SB."""
    N = 3; dim = 8
    H = np.zeros((dim, dim), dtype=complex)
    for p in [X, Y, Z]:
        # S-A bond
        ops = [I2]*3; ops[0] = p; ops[1] = p
        H += J_SA * kron_list(ops)
        # S-B bond
        ops = [I2]*3; ops[0] = p; ops[2] = p
        H += J_SB * kron_list(ops)
    return H

def build_L(H, gamma, N=3):
    dim = 2**N; dim2 = dim**2; Id = np.eye(dim)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for i in range(N):
        ops = [I2]*N; ops[i] = Z; Zi = kron_list(ops)
        L += gamma * (np.kron(Zi, Zi.conj()) - np.eye(dim2))
    return L

def ptrace(rho, N, keep):
    """Partial trace of N-qubit density matrix, keeping specified qubits."""
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

def von_neumann_entropy(rho):
    """S(rho) = -Tr(rho log2 rho)"""
    evals = np.real(np.linalg.eigvalsh(rho))
    evals = evals[evals > 1e-15]
    return -np.sum(evals * np.log2(evals))

def channel_map(psi_A, H, L, gamma, t, N=3):
    """
    One-shot channel: |0>_S |psi>_A |0>_B -> trace_SA[rho(t)] = rho_B(t)
    Also returns full rho(t) for coherent information.
    """
    # Initial state: |0>_S x |psi>_A x |0>_B
    psi_S = up
    psi_B = up
    psi_full = np.kron(np.kron(psi_S, psi_A), psi_B)
    rho0 = np.outer(psi_full, psi_full.conj())
    
    # Evolve via Liouvillian
    rho0_vec = rho0.flatten()
    rho_t_vec = expm(L * t) @ rho0_vec
    rho_t = rho_t_vec.reshape(2**N, 2**N)
    rho_t = (rho_t + rho_t.conj().T) / 2
    rho_t /= np.trace(rho_t).real
    
    # Trace out S and A, keep B (qubit 2)
    rho_B = ptrace(rho_t, N, [2])
    
    return rho_B, rho_t

def fidelity_pure(psi, rho):
    """Fidelity between pure state |psi> and mixed state rho."""
    return np.real(psi.conj() @ rho @ psi)

# ================================================================
print("=" * 80)
print("VERIFYING GPT's ONE-SHOT CHANNEL ANALYSIS")
print("=" * 80)

# GPT's parameters
J_SA = 1.0
J_SB = 2.0
gamma = 0.05
N = 3

H = build_H_star_asym(J_SA, J_SB)
L = build_L(H, gamma, N)

print(f"\nParameters: J_SA={J_SA}, J_SB={J_SB}, gamma={gamma}")
print(f"Liouvillian: {L.shape[0]}x{L.shape[0]}")

# ================================================================
# Test 1: Average state-transfer fidelity
# F_avg = integral over Bloch sphere of F(|psi>, channel(|psi>))
# For a qubit channel, F_avg = (2*F_e + 1)/3 where F_e is entanglement fidelity
# Or we can sample uniformly from the Bloch sphere
# ================================================================
print(f"\n{'='*80}")
print("TEST 1: Average State-Transfer Fidelity")
print(f"{'='*80}")
print("GPT claims: F_avg = 0.884 at t = 1.30")

# Sample 1000 random pure states on the Bloch sphere
np.random.seed(42)
n_samples = 1000

# Pre-compute matrix exponentials for time sweep
times = np.linspace(0.1, 3.0, 60)
fidelities_avg = []

# Generate random pure states (Haar measure on single qubit)
random_states = []
for _ in range(n_samples):
    # Random state on Bloch sphere
    theta = np.arccos(1 - 2*np.random.random())
    phi = 2 * np.pi * np.random.random()
    psi = np.array([np.cos(theta/2), np.exp(1j*phi)*np.sin(theta/2)])
    random_states.append(psi)

# Also the 6 cardinal states for quick check
cardinal_states = {
    '|0>': up, '|1>': dn,
    '|+>': plus, '|->': minus,
    '|+i>': (up + 1j*dn)/np.sqrt(2),
    '|-i>': (up - 1j*dn)/np.sqrt(2),
}

print(f"\nTime sweep with {n_samples} random states + 6 cardinal states...")
print(f"{'t':>6} {'F_avg(rand)':>12} {'F_avg(card)':>12} {'F(|0>)':>8} {'F(|1>)':>8} {'F(|+>)':>8}")
print(f"{'-'*60}")

best_t = 0
best_f = 0

for t in times:
    eLt = expm(L * t)
    
    # Cardinal states
    f_cardinal = []
    f_details = {}
    for name, psi_A in cardinal_states.items():
        psi_full = np.kron(np.kron(up, psi_A), up)
        rho0 = np.outer(psi_full, psi_full.conj())
        rho_t_vec = eLt @ rho0.flatten()
        rho_t = rho_t_vec.reshape(8, 8)
        rho_t = (rho_t + rho_t.conj().T) / 2
        rho_t /= np.trace(rho_t).real
        rho_B = ptrace(rho_t, 3, [2])
        f = fidelity_pure(psi_A, rho_B)
        f_cardinal.append(f)
        f_details[name] = f
    f_avg_card = np.mean(f_cardinal)
    
    # Random states (subsample for speed at each t)
    f_random = []
    for psi_A in random_states:
        psi_full = np.kron(np.kron(up, psi_A), up)
        rho0 = np.outer(psi_full, psi_full.conj())
        rho_t_vec = eLt @ rho0.flatten()
        rho_t = rho_t_vec.reshape(8, 8)
        rho_t = (rho_t + rho_t.conj().T) / 2
        rho_t /= np.trace(rho_t).real
        rho_B = ptrace(rho_t, 3, [2])
        f = fidelity_pure(psi_A, rho_B)
        f_random.append(f)
    f_avg_rand = np.mean(f_random)
    
    if f_avg_rand > best_f:
        best_f = f_avg_rand
        best_t = t
    
    # Print every 5th point
    if len(fidelities_avg) % 3 == 0:
        print(f"{t:>6.2f} {f_avg_rand:>12.4f} {f_avg_card:>12.4f} "
              f"{f_details['|0>']:>8.4f} {f_details['|1>']:>8.4f} {f_details['|+>']:>8.4f}")
    
    fidelities_avg.append((t, f_avg_rand, f_avg_card))

print(f"\n  BEST: F_avg = {best_f:.4f} at t = {best_t:.2f}")
print(f"  GPT:  F_avg = 0.884 at t = 1.30")
print(f"  Classical limit: F = 2/3 = 0.6667")

# Find window where F > 2/3
above = [(t, f) for t, f, _ in fidelities_avg if f > 2/3]
if above:
    print(f"  F > 2/3 window: t in [{above[0][0]:.2f}, {above[-1][0]:.2f}]")
    print(f"  GPT claims:     t in [0.60, 1.65]")

# ================================================================
print(f"\n{'='*80}")
print("TEST 2: Holevo Lower Bound (binary alphabet {|0>, |1>})")
print(f"{'='*80}")
print("GPT claims: 0.529 bits/use at t = 1.30")

# Holevo bound chi = S(rho_avg) - (1/2)[S(rho_0) + S(rho_1)]
# where rho_avg = (1/2)(channel(|0>) + channel(|1>))

print(f"\n{'t':>6} {'S(rho_avg)':>12} {'S(rho_0)':>10} {'S(rho_1)':>10} {'chi':>10}")
print(f"{'-'*55}")

best_chi = 0
best_chi_t = 0

for t in np.linspace(0.1, 3.0, 60):
    eLt = expm(L * t)
    
    rho_B_list = []
    for psi_A in [up, dn]:
        psi_full = np.kron(np.kron(up, psi_A), up)
        rho0 = np.outer(psi_full, psi_full.conj())
        rho_t_vec = eLt @ rho0.flatten()
        rho_t = rho_t_vec.reshape(8, 8)
        rho_t = (rho_t + rho_t.conj().T) / 2
        rho_t /= np.trace(rho_t).real
        rho_B = ptrace(rho_t, 3, [2])
        rho_B_list.append(rho_B)
    
    rho_avg = 0.5 * (rho_B_list[0] + rho_B_list[1])
    S_avg = von_neumann_entropy(rho_avg)
    S_0 = von_neumann_entropy(rho_B_list[0])
    S_1 = von_neumann_entropy(rho_B_list[1])
    chi = S_avg - 0.5 * (S_0 + S_1)
    
    if chi > best_chi:
        best_chi = chi
        best_chi_t = t
    
    if abs(t - 1.30) < 0.03 or abs(t - best_chi_t) < 0.03:
        print(f"{t:>6.2f} {S_avg:>12.4f} {S_0:>10.4f} {S_1:>10.4f} {chi:>10.4f}")

print(f"\n  BEST: chi = {best_chi:.4f} bits at t = {best_chi_t:.2f}")
print(f"  GPT:  chi = 0.529 bits at t = 1.30")

# ================================================================
print(f"\n{'='*80}")
print("TEST 3: Coherent Information Lower Bound")
print(f"{'='*80}")
print("GPT claims: 0.185 qubits/use at t = 1.30")
print("I_coh = S(rho_B) - S(rho_SB) for maximally mixed input on A")

# Coherent information: I_coh(A>B) = S(B) - S(SB)
# For maximally mixed input: rho_A = I/2
# Initial state: |0>_S (I/2)_A |0>_B -> mixed

print(f"\n{'t':>6} {'S(B)':>10} {'S(SB)':>10} {'I_coh':>10}")
print(f"{'-'*40}")

best_icoh = -10
best_icoh_t = 0

for t in np.linspace(0.1, 3.0, 60):
    eLt = expm(L * t)
    
    # Maximally mixed input: average over |0> and |1> on A
    rho_t_total = np.zeros((8, 8), dtype=complex)
    for psi_A in [up, dn]:
        psi_full = np.kron(np.kron(up, psi_A), up)
        rho0 = np.outer(psi_full, psi_full.conj())
        rho_t_vec = eLt @ rho0.flatten()
        rho_t = rho_t_vec.reshape(8, 8)
        rho_t = (rho_t + rho_t.conj().T) / 2
        rho_t_total += 0.5 * rho_t
    
    rho_t_total /= np.trace(rho_t_total).real
    
    # S(B) - trace out S,A keep B
    rho_B = ptrace(rho_t_total, 3, [2])
    S_B = von_neumann_entropy(rho_B)
    
    # S(SB) - trace out A, keep S,B  
    rho_SB = ptrace(rho_t_total, 3, [0, 2])
    S_SB = von_neumann_entropy(rho_SB)
    
    I_coh = S_B - S_SB
    
    if I_coh > best_icoh:
        best_icoh = I_coh
        best_icoh_t = t
    
    if abs(t - 1.30) < 0.03 or abs(t - best_icoh_t) < 0.03:
        print(f"{t:>6.2f} {S_B:>10.4f} {S_SB:>10.4f} {I_coh:>10.4f}")

print(f"\n  BEST: I_coh = {best_icoh:.4f} qubits at t = {best_icoh_t:.2f}")
print(f"  GPT:  I_coh = 0.185 qubits at t = 1.30")

# ================================================================
print(f"\n{'='*80}")
print("TEST 4: Gamma Dependence")
print(f"{'='*80}")
print("GPT's table:")
print("  gamma  F_avg   Holevo   I_coh")
print("  0.001  0.933   0.625    0.516")
print("  0.010  0.923   0.605    0.438")
print("  0.050  0.884   0.529    0.185")
print("  0.100  0.841   0.458   -0.040")
print("  0.200  0.771   0.363   -0.335")

print(f"\nOur results:")
print(f"{'gamma':>8} {'F_avg':>8} {'t_opt':>6} {'Holevo':>8} {'t_opt':>6} {'I_coh':>8} {'t_opt':>6}")
print(f"{'-'*55}")

for g in [0.001, 0.010, 0.050, 0.100, 0.200]:
    H_test = build_H_star_asym(1.0, 2.0)
    L_test = build_L(H_test, g, 3)
    
    best_f = 0; best_f_t = 0
    best_chi = 0; best_chi_t = 0
    best_ic = -10; best_ic_t = 0
    
    for t in np.linspace(0.1, 3.0, 60):
        eLt = expm(L_test * t)
        
        # Fidelity (6 cardinal states)
        f_list = []
        rho_B_01 = []
        rho_t_mixed = np.zeros((8,8), dtype=complex)
        
        for psi_A in [up, dn, plus, minus, (up+1j*dn)/np.sqrt(2), (up-1j*dn)/np.sqrt(2)]:
            psi_full = np.kron(np.kron(up, psi_A), up)
            rho0 = np.outer(psi_full, psi_full.conj())
            rho_t_vec = eLt @ rho0.flatten()
            rho_t = rho_t_vec.reshape(8, 8)
            rho_t = (rho_t + rho_t.conj().T) / 2
            rho_t /= np.trace(rho_t).real
            rho_B = ptrace(rho_t, 3, [2])
            f_list.append(fidelity_pure(psi_A, rho_B))
        
        f_avg = np.mean(f_list)
        if f_avg > best_f:
            best_f = f_avg; best_f_t = t
        
        # Holevo (binary |0>, |1>)
        rho_B_01 = []
        for psi_A in [up, dn]:
            psi_full = np.kron(np.kron(up, psi_A), up)
            rho0 = np.outer(psi_full, psi_full.conj())
            rho_t_vec = eLt @ rho0.flatten()
            rho_t = rho_t_vec.reshape(8, 8)
            rho_t = (rho_t + rho_t.conj().T) / 2
            rho_t /= np.trace(rho_t).real
            rho_B_01.append(ptrace(rho_t, 3, [2]))
        
        rho_avg = 0.5*(rho_B_01[0] + rho_B_01[1])
        chi = von_neumann_entropy(rho_avg) - 0.5*(von_neumann_entropy(rho_B_01[0]) + von_neumann_entropy(rho_B_01[1]))
        if chi > best_chi:
            best_chi = chi; best_chi_t = t
        
        # Coherent info
        rho_t_mix = np.zeros((8,8), dtype=complex)
        for psi_A in [up, dn]:
            psi_full = np.kron(np.kron(up, psi_A), up)
            rho0 = np.outer(psi_full, psi_full.conj())
            rho_t_vec = eLt @ rho0.flatten()
            rho_t = rho_t_vec.reshape(8, 8)
            rho_t = (rho_t + rho_t.conj().T) / 2
            rho_t_mix += 0.5 * rho_t
        rho_t_mix /= np.trace(rho_t_mix).real
        rho_B_m = ptrace(rho_t_mix, 3, [2])
        rho_SB_m = ptrace(rho_t_mix, 3, [0, 2])
        ic = von_neumann_entropy(rho_B_m) - von_neumann_entropy(rho_SB_m)
        if ic > best_ic:
            best_ic = ic; best_ic_t = t
    
    print(f"{g:>8.3f} {best_f:>8.4f} {best_f_t:>6.2f} {best_chi:>8.4f} {best_chi_t:>6.2f} {best_ic:>8.4f} {best_ic_t:>6.2f}")

# ================================================================
print(f"\n{'='*80}")
print("TEST 5: J_SB Dependence (J_SA=1.0, gamma=0.05)")
print(f"{'='*80}")
print("GPT's table:")
print("  J_SB   F_avg   Holevo   I_coh   t_opt")
print("  1.0    0.850   0.471    0.088   1.15")
print("  1.5    0.861   0.499    0.109   1.60")
print("  2.0    0.884   0.529    0.185   1.30")
print("  3.0    0.830   0.394    0.015   0.95")

print(f"\nOur results:")
print(f"{'J_SB':>6} {'F_avg':>8} {'t_opt':>6} {'Holevo':>8} {'t_opt':>6} {'I_coh':>8} {'t_opt':>6}")
print(f"{'-'*55}")

for jsb in [1.0, 1.5, 2.0, 3.0]:
    H_test = build_H_star_asym(1.0, jsb)
    L_test = build_L(H_test, 0.05, 3)
    
    best_f = 0; best_f_t = 0
    best_chi = 0; best_chi_t = 0
    best_ic = -10; best_ic_t = 0
    
    for t in np.linspace(0.1, 3.0, 60):
        eLt = expm(L_test * t)
        
        # Fidelity
        f_list = []
        for psi_A in [up, dn, plus, minus, (up+1j*dn)/np.sqrt(2), (up-1j*dn)/np.sqrt(2)]:
            psi_full = np.kron(np.kron(up, psi_A), up)
            rho0 = np.outer(psi_full, psi_full.conj())
            rho_t_vec = eLt @ rho0.flatten()
            rho_t = rho_t_vec.reshape(8, 8)
            rho_t = (rho_t + rho_t.conj().T) / 2
            rho_t /= np.trace(rho_t).real
            rho_B = ptrace(rho_t, 3, [2])
            f_list.append(fidelity_pure(psi_A, rho_B))
        f_avg = np.mean(f_list)
        if f_avg > best_f:
            best_f = f_avg; best_f_t = t
        
        # Holevo
        rho_B_01 = []
        for psi_A in [up, dn]:
            psi_full = np.kron(np.kron(up, psi_A), up)
            rho0 = np.outer(psi_full, psi_full.conj())
            rho_t_vec = eLt @ rho0.flatten()
            rho_t = rho_t_vec.reshape(8, 8)
            rho_t = (rho_t + rho_t.conj().T) / 2
            rho_t /= np.trace(rho_t).real
            rho_B_01.append(ptrace(rho_t, 3, [2]))
        rho_avg = 0.5*(rho_B_01[0] + rho_B_01[1])
        chi = von_neumann_entropy(rho_avg) - 0.5*(von_neumann_entropy(rho_B_01[0]) + von_neumann_entropy(rho_B_01[1]))
        if chi > best_chi:
            best_chi = chi; best_chi_t = t
        
        # Coherent info
        rho_t_mix = np.zeros((8,8), dtype=complex)
        for psi_A in [up, dn]:
            psi_full = np.kron(np.kron(up, psi_A), up)
            rho0 = np.outer(psi_full, psi_full.conj())
            rho_t_vec = eLt @ rho0.flatten()
            rho_t = rho_t_vec.reshape(8, 8)
            rho_t = (rho_t + rho_t.conj().T) / 2
            rho_t_mix += 0.5 * rho_t
        rho_t_mix /= np.trace(rho_t_mix).real
        rho_B_m = ptrace(rho_t_mix, 3, [2])
        rho_SB_m = ptrace(rho_t_mix, 3, [0, 2])
        ic = von_neumann_entropy(rho_B_m) - von_neumann_entropy(rho_SB_m)
        if ic > best_ic:
            best_ic = ic; best_ic_t = t
    
    print(f"{jsb:>6.1f} {best_f:>8.4f} {best_f_t:>6.2f} {best_chi:>8.4f} {best_chi_t:>6.2f} {best_ic:>8.4f} {best_ic_t:>6.2f}")

print(f"\n{'='*80}")
print("DONE")
print(f"{'='*80}")
