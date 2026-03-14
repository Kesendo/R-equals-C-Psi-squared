"""
Why does |+-+-⟩ cross on a ring but |0+0+⟩ doesn't?
=====================================================
Decompose initial states into Liouvillian eigenmodes.
Track which palindrome pairs get activated.
Find the selection rule.

March 14, 2026
"""

import numpy as np
from itertools import product as iproduct
import sys, io, time as timer

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

def build_H(N, bonds, J=1.0):
    dim = 2**N
    H = np.zeros((dim, dim), dtype=complex)
    for (i, j) in bonds:
        for p in [X, Y, Z]:
            ops = [I2]*N; ops[i] = p; ops[j] = p
            H += J * kron_list(ops)
    return H

def build_L(N, H, gamma):
    dim = 2**N; dim2 = dim**2; Id = np.eye(dim)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for i in range(N):
        ops = [I2]*N; ops[i] = Z; Zi = kron_list(ops)
        L += gamma * (np.kron(Zi, Zi.conj()) - np.eye(dim2))
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

def concurrence_2q(rho):
    YY = np.kron(Y, Y)
    R = rho @ (YY @ rho.conj() @ YY)
    ev = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(R), 0))))[::-1]
    return max(0, ev[0] - ev[1] - ev[2] - ev[3])

def l1_coh(rho):
    return np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))

# XY-weight for Pauli index classification
PI_PERM = {0: 1, 1: 0, 2: 3, 3: 2}
PI_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}

def xy_weight_pauli(indices):
    return sum(1 for i in indices if i in (1, 2))

# ================================================================
N = 4
bonds_ring = [(0,1), (1,2), (2,3), (3,0)]
H = build_H(N, bonds_ring)
gamma = 0.05
L = build_L(N, H, gamma)

print("=" * 80)
print("WHY DOES |+-+-⟩ CROSS ON A RING BUT |0+0+⟩ DOESN'T?")
print("=" * 80)

# Diagonalize L
print("\nDiagonalizing 256x256 Liouvillian...", flush=True)
evals, R_vecs = np.linalg.eig(L)
L_vecs = np.linalg.inv(R_vecs)  # Left eigenvectors
decay_rates = -evals.real
freqs = evals.imag
print(f"Done. {len(evals)} eigenvalues.")

# Identify palindrome pairs
center = N * gamma  # = 0.2
pairs_found = []
used = set()
for i in range(len(evals)):
    if i in used:
        continue
    partner_rate = 2 * center - decay_rates[i]
    partner_freq = -freqs[i]
    for j in range(len(evals)):
        if j in used or j == i:
            continue
        if (abs(decay_rates[j] - partner_rate) < 1e-8 and
            abs(freqs[j] - partner_freq) < 1e-8):
            pairs_found.append((i, j))
            used.add(i)
            used.add(j)
            break

print(f"\nPalindrome pairs found: {len(pairs_found)}")
print(f"Unpaired modes: {len(evals) - 2*len(pairs_found)}")

# ================================================================
# Decompose initial states into Liouvillian eigenmodes
# ================================================================
initial_states = {
    '|0+0+⟩': np.kron(np.kron(up, plus), np.kron(up, plus)),
    '|+-+-⟩': np.kron(np.kron(plus, minus), np.kron(plus, minus)),
    '|0+0-⟩': np.kron(np.kron(up, plus), np.kron(up, minus)),
    '|++++⟩': kron_list([plus]*4),
    'Bell01+Bell23': np.kron(
        (np.kron(up,up)+np.kron(dn,dn))/np.sqrt(2),
        (np.kron(up,up)+np.kron(dn,dn))/np.sqrt(2)),
}

for name, psi0 in initial_states.items():
    rho0 = np.outer(psi0, psi0.conj())
    rho0_vec = rho0.flatten()
    
    # Project onto eigenmodes: c_k = L_vecs[k] @ rho0_vec
    coeffs = L_vecs @ rho0_vec
    
    print(f"\n{'=' * 80}")
    print(f"State: {name}")
    print(f"{'=' * 80}")
    
    # Amplitude distribution across decay rate sectors
    # Group by decay rate (rounded)
    rate_groups = {}
    for k in range(len(evals)):
        rate = round(decay_rates[k], 6)
        if rate not in rate_groups:
            rate_groups[rate] = []
        rate_groups[rate].append((k, coeffs[k]))
    
    # Sort by decay rate
    sorted_rates = sorted(rate_groups.keys())
    
    # Show amplitude distribution
    print(f"\n  Amplitude by decay rate sector:")
    print(f"  {'Rate':>10} {'# modes':>8} {'Total |c|':>12} {'Max |c|':>12} {'Freq range':>20}")
    print(f"  {'-'*65}")
    
    sig_rates = []
    for rate in sorted_rates:
        modes = rate_groups[rate]
        amps = [abs(c) for _, c in modes]
        total_amp = sum(amps)
        max_amp = max(amps)
        mode_freqs = [freqs[k] for k, _ in modes]
        freq_range = f"[{min(mode_freqs):.3f}, {max(mode_freqs):.3f}]"
        
        if total_amp > 1e-6:
            sig_rates.append((rate, total_amp, len(modes)))
            print(f"  {rate:>10.6f} {len(modes):>8} {total_amp:>12.6f} {max_amp:>12.6f} {freq_range:>20}")
    
    # Which palindrome pairs are excited?
    print(f"\n  Palindrome pair activation:")
    print(f"  {'Pair':>6} {'Rate_1':>10} {'Rate_2':>10} {'|c_1|':>10} {'|c_2|':>10} {'Product':>12} {'Active':>8}")
    print(f"  {'-'*70}")
    
    active_pairs = 0
    total_product = 0
    for i_mode, j_mode in pairs_found:
        c1 = abs(coeffs[i_mode])
        c2 = abs(coeffs[j_mode])
        product = c1 * c2
        active = product > 1e-8
        if active:
            active_pairs += 1
            total_product += product
            if product > 1e-4:  # only show significant ones
                print(f"  {f'({i_mode},{j_mode})':>6} {decay_rates[i_mode]:>10.4f} "
                      f"{decay_rates[j_mode]:>10.4f} {c1:>10.6f} {c2:>10.6f} "
                      f"{product:>12.6f} {'YES' if active else '':>8}")
    
    print(f"\n  Active palindrome pairs: {active_pairs}/{len(pairs_found)}")
    print(f"  Total pair product: {total_product:.6f}")

    # Time evolution: track CΨ for all pairs
    print(f"\n  CΨ_max by pair:")
    times = np.linspace(0, 15, 1500)
    pairs_phys = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
    labels = ['01 ring', '02 diag', '03 ring', '12 ring', '13 diag', '23 ring']
    
    Vinv_rho0 = L_vecs @ rho0_vec  # pre-compute
    
    for pair, label in zip(pairs_phys, labels):
        cpsi_max = 0
        t_max = 0
        c_at_max = 0
        
        for t in times:
            rho_vec_t = R_vecs @ (np.exp(evals * t) * Vinv_rho0)
            rho_t = rho_vec_t.reshape(2**N, 2**N)
            rho_t = (rho_t + rho_t.conj().T) / 2
            rho_t /= np.trace(rho_t)
            
            rho_p = ptrace(rho_t, N, list(pair))
            C = concurrence_2q(rho_p)
            l1 = l1_coh(rho_p)
            cpsi = C * l1/3
            if cpsi > cpsi_max:
                cpsi_max = cpsi
                t_max = t
                c_at_max = C
        
        crosses = "CROSSES" if cpsi_max > 0.25 else f"gap={0.25-cpsi_max:.4f}"
        print(f"    {label}: CΨ_max={cpsi_max:.6f} (C={c_at_max:.4f}) at t={t_max:.3f}  {crosses}")

# ================================================================
# Physical analysis: what's different about |+-+-⟩?
# ================================================================
print(f"\n\n{'=' * 80}")
print("PHYSICAL ANALYSIS: What makes |+-+-⟩ special?")
print("=" * 80)

# |+-+-⟩ in the computational basis
psi_alt = np.kron(np.kron(plus, minus), np.kron(plus, minus))
psi_same = np.kron(np.kron(up, plus), np.kron(up, plus))

# Hamiltonian expectation values
print(f"\n  Energy analysis:")
for name, psi in [('|+-+-⟩', psi_alt), ('|0+0+⟩', psi_same), ('|++++⟩', kron_list([plus]*4))]:
    E = np.real(psi.conj() @ H @ psi)
    E2 = np.real(psi.conj() @ (H @ H) @ psi)
    var_E = E2 - E**2
    print(f"  {name:>10}: <H> = {E:>8.4f}, var(H) = {var_E:>8.4f}, std(H) = {np.sqrt(max(0,var_E)):>8.4f}")

# Bond-by-bond energy
print(f"\n  Bond energies (XX + YY + ZZ per bond):")
print(f"  {'State':>10} {'(0,1)':>8} {'(1,2)':>8} {'(2,3)':>8} {'(3,0)':>8} {'Total':>8}")

for name, psi in [('|+-+-⟩', psi_alt), ('|0+0+⟩', psi_same)]:
    bond_E = []
    for (i, j) in bonds_ring:
        E_bond = 0
        for p in [X, Y, Z]:
            ops = [I2]*N; ops[i] = p; ops[j] = p
            term = kron_list(ops)
            E_bond += np.real(psi.conj() @ term @ psi)
        bond_E.append(E_bond)
    print(f"  {name:>10} {bond_E[0]:>8.4f} {bond_E[1]:>8.4f} {bond_E[2]:>8.4f} {bond_E[3]:>8.4f} {sum(bond_E):>8.4f}")

# What does the Hamiltonian DO to each state in the first instant?
print(f"\n  Initial dynamics: |drho/dt| at t=0 decomposition")
for name, psi in [('|+-+-⟩', psi_alt), ('|0+0+⟩', psi_same)]:
    rho = np.outer(psi, psi.conj())
    
    # Hamiltonian part: -i[H, rho]
    drho_H = -1j * (H @ rho - rho @ H)
    
    # Dephasing part: sum_i gamma(Zi rho Zi - rho)
    drho_D = np.zeros_like(rho)
    for i in range(N):
        ops = [I2]*N; ops[i] = Z; Zi = kron_list(ops)
        drho_D += gamma * (Zi @ rho @ Zi - rho)
    
    # Track what happens to each pair's entanglement
    print(f"\n  {name}:")
    print(f"    ||dρ_H/dt|| = {np.linalg.norm(drho_H):.4f} (Hamiltonian)")
    print(f"    ||dρ_D/dt|| = {np.linalg.norm(drho_D):.4f} (Dephasing)")
    print(f"    Ratio H/D = {np.linalg.norm(drho_H)/max(np.linalg.norm(drho_D),1e-15):.2f}")
    
    # Which pair subspaces get driven?
    for pair, label in zip(pairs_phys, labels):
        rho_p = ptrace(rho, N, list(pair))
        drho_H_p = ptrace(drho_H, N, list(pair))
        drho_D_p = ptrace(drho_D, N, list(pair))
        
        # Initial state of pair
        C0 = concurrence_2q(rho_p)
        l1_0 = l1_coh(rho_p)
        
        # How much does H drive this pair?
        drive = np.linalg.norm(drho_H_p)
        damp = np.linalg.norm(drho_D_p)
        
        if drive > 1e-6 or C0 > 1e-6:
            print(f"    {label}: C(0)={C0:.4f}, l1(0)={l1_0:.4f}, "
                  f"||drive||={drive:.4f}, ||damp||={damp:.4f}")

# ================================================================
# The antiferromagnet connection
# ================================================================
print(f"\n\n{'=' * 80}")
print("THE ANTIFERROMAGNET CONNECTION")
print("=" * 80)

# |+-+-⟩ is the Neel-like state in the X basis
# On a ring with Heisenberg coupling, alternating phases create
# maximum exchange energy between neighbors

# Check: what is |+-+-⟩ in the Z basis?
print(f"\n  |+-+-⟩ computational basis decomposition:")
psi = psi_alt
for i in range(2**N):
    amp = psi[i]
    if abs(amp) > 1e-10:
        bits = format(i, f'0{N}b')
        print(f"    |{bits}⟩: {amp:.4f}")

print(f"\n  Key property: |+-+-⟩ is a uniform superposition of ALL 2^N = {2**N} basis states")
print(f"  with ALTERNATING SIGNS determined by the +/- pattern.")
print(f"  This is the X-basis Neel state: maximum staggered magnetization in X.")

# Compare magnetization profiles
print(f"\n  Local observables at t=0:")
print(f"  {'State':>10} {'<X0>':>6} {'<X1>':>6} {'<X2>':>6} {'<X3>':>6} "
      f"{'<Z0>':>6} {'<Z1>':>6} {'<Z2>':>6} {'<Z3>':>6}")
for name, psi in [('|+-+-⟩', psi_alt), ('|0+0+⟩', psi_same)]:
    xvals = []
    zvals = []
    for i in range(N):
        ops_x = [I2]*N; ops_x[i] = X
        ops_z = [I2]*N; ops_z[i] = Z
        xvals.append(np.real(psi.conj() @ kron_list(ops_x) @ psi))
        zvals.append(np.real(psi.conj() @ kron_list(ops_z) @ psi))
    print(f"  {name:>10} {xvals[0]:>6.2f} {xvals[1]:>6.2f} {xvals[2]:>6.2f} {xvals[3]:>6.2f} "
          f"{zvals[0]:>6.2f} {zvals[1]:>6.2f} {zvals[2]:>6.2f} {zvals[3]:>6.2f}")

# Nearest-neighbor correlations at t=0
print(f"\n  Nearest-neighbor correlations <σ_i σ_j> at t=0:")
print(f"  {'State':>10} {'<XX>_01':>8} {'<YY>_01':>8} {'<ZZ>_01':>8} {'<XX>_02':>8}")
for name, psi in [('|+-+-⟩', psi_alt), ('|0+0+⟩', psi_same)]:
    for tag, (i,j) in [('_01', (0,1)), ('_02', (0,2))]:
        vals = []
        for p in [X, Y, Z]:
            ops = [I2]*N; ops[i] = p; ops[j] = p
            vals.append(np.real(psi.conj() @ kron_list(ops) @ psi))
        if tag == '_01':
            xx01, yy01, zz01 = vals
        else:
            xx02 = vals[0]
    print(f"  {name:>10} {xx01:>8.4f} {yy01:>8.4f} {zz01:>8.4f} {xx02:>8.4f}")

print(f"\n  KEY INSIGHT:")
print(f"  |+-+-⟩ has <XX>_neighbors = -1 (ANTI-correlated in X)")
print(f"  This is maximum exchange energy for the XX part of Heisenberg.")
print(f"  The YY and ZZ parts contribute zero (product state).")
print(f"  But the Hamiltonian immediately converts XX correlation into YY and ZZ,")
print(f"  building ENTANGLEMENT from the initial classical correlation.")

# ================================================================
# Selection rule hunt: what predicts crossing?
# ================================================================
print(f"\n\n{'=' * 80}")
print("SELECTION RULE: What predicts whether a state crosses 1/4?")
print("=" * 80)

# Hypothesis: crossing requires sufficient STAGGERED correlation
# between ring neighbors, which provides the Hamiltonian "fuel" to
# build entanglement.

# Test many product states
print(f"\n  Systematic scan: |a b c d⟩ product states on N=4 ring")
print(f"  Only testing single-qubit pure states: |0⟩, |1⟩, |+⟩, |-⟩")

qubit_states = {
    '0': up, '1': dn, '+': plus, '-': minus,
}

# Pre-compute evolution
Vinv = L_vecs  # already computed
V = R_vecs

results = []
times = np.linspace(0, 15, 800)

for s0 in qubit_states:
    for s1 in qubit_states:
        for s2 in qubit_states:
            for s3 in qubit_states:
                label = f"|{s0}{s1}{s2}{s3}⟩"
                psi = np.kron(np.kron(qubit_states[s0], qubit_states[s1]),
                              np.kron(qubit_states[s2], qubit_states[s3]))
                rho = np.outer(psi, psi.conj())
                rho_vec = rho.flatten()
                
                # Energy and variance
                E = np.real(psi.conj() @ H @ psi)
                E2 = np.real(psi.conj() @ (H@H) @ psi)
                var_E = E2 - E**2
                
                # Staggered X magnetization: M_stagg = <X0> - <X1> + <X2> - <X3>
                m_stagg = 0
                for i in range(N):
                    ops = [I2]*N; ops[i] = X
                    sign = 1 if i % 2 == 0 else -1
                    m_stagg += sign * np.real(psi.conj() @ kron_list(ops) @ psi)
                
                # NN XX correlation
                xx_nn = 0
                for (i,j) in bonds_ring:
                    ops = [I2]*N; ops[i] = X; ops[j] = X
                    xx_nn += np.real(psi.conj() @ kron_list(ops) @ psi)
                
                # Best CΨ across all pairs
                proj = Vinv @ rho_vec
                best_cpsi = 0
                best_pair = ""
                for pair, plabel in zip(pairs_phys, labels):
                    cpsi_max = 0
                    for t in times:
                        rho_vec_t = V @ (np.exp(evals * t) * proj)
                        rho_t = rho_vec_t.reshape(2**N, 2**N)
                        rho_t = (rho_t + rho_t.conj().T) / 2
                        rho_t /= np.trace(rho_t)
                        rho_p = ptrace(rho_t, N, list(pair))
                        C = concurrence_2q(rho_p)
                        l1 = l1_coh(rho_p)
                        cpsi_max = max(cpsi_max, C * l1/3)
                    if cpsi_max > best_cpsi:
                        best_cpsi = cpsi_max
                        best_pair = plabel
                
                crosses = best_cpsi > 0.25
                results.append({
                    'label': label, 'E': E, 'var_E': var_E,
                    'm_stagg': m_stagg, 'xx_nn': xx_nn,
                    'cpsi_max': best_cpsi, 'best_pair': best_pair,
                    'crosses': crosses
                })

# Sort by CΨ_max
results.sort(key=lambda x: -x['cpsi_max'])

print(f"\n  Top 20 states by CΨ_max:")
print(f"  {'State':>12} {'<H>':>7} {'var(H)':>8} {'M_stagg':>8} {'<XX>_nn':>8} "
      f"{'CΨ_max':>8} {'Best pair':>10} {'Crosses':>8}")
print(f"  {'-'*82}")
for r in results[:20]:
    tag = "YES" if r['crosses'] else ""
    print(f"  {r['label']:>12} {r['E']:>7.2f} {r['var_E']:>8.2f} "
          f"{r['m_stagg']:>8.2f} {r['xx_nn']:>8.2f} "
          f"{r['cpsi_max']:>8.4f} {r['best_pair']:>10} {tag:>8}")

# Count crossings
n_cross = sum(1 for r in results if r['crosses'])
print(f"\n  Total: {n_cross}/{len(results)} product states cross 1/4 on ring")

# Correlation analysis
print(f"\n  Correlation of CΨ_max with predictors:")
cpsi_vals = np.array([r['cpsi_max'] for r in results])
for pred_name, pred_key in [('Energy <H>', 'E'), ('Energy variance', 'var_E'),
                             ('Staggered M_x', 'm_stagg'), ('<XX> neighbors', 'xx_nn')]:
    pred_vals = np.array([r[pred_key] for r in results])
    if np.std(pred_vals) > 1e-10 and np.std(cpsi_vals) > 1e-10:
        corr = np.corrcoef(pred_vals, cpsi_vals)[0, 1]
        print(f"    {pred_name:>20}: r = {corr:>7.4f}")

print(f"\n{'=' * 80}")
print("DONE")
print(f"{'=' * 80}")
