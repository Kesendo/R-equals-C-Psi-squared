"""
Orphaned Results Deep Exploration - March 14, 2026
====================================================
Run on Tom's machine (Intel Core Ultra 9 285k, 128GB RAM).
No timeouts. Go deep.

Three explorations:
1. u = C(Ψ+R) under Π - connection between Mandelbrot and palindrome
2. Ring near-miss - systematic scan across initial states and parameters  
3. Echo effect - full characterization with Fourier analysis

Uses eigendecomposition for fast time evolution.
"""

import numpy as np
from itertools import product as iproduct
import time as timer
import sys
import io

# Fix Windows console encoding
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ================================================================
# Infrastructure
# ================================================================
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)
minus = (up - dn) / np.sqrt(2)


def kron_list(ops):
    result = ops[0]
    for o in ops[1:]:
        result = np.kron(result, o)
    return result


def build_H(N, bonds, J=1.0):
    dim = 2**N
    H = np.zeros((dim, dim), dtype=complex)
    for (i, j) in bonds:
        for p in [X, Y, Z]:
            ops = [I2]*N; ops[i] = p; ops[j] = p
            H += J * kron_list(ops)
    return H


def build_L(N, H, gamma):
    """Build Liouvillian superoperator. gamma can be float or list per site."""
    dim = 2**N; dim2 = dim**2; Id = np.eye(dim, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    if isinstance(gamma, (int, float)):
        gamma = [gamma] * N
    for i in range(N):
        ops = [I2]*N; ops[i] = Z; Zi = kron_list(ops)
        L += gamma[i] * (np.kron(Zi, Zi.conj()) - np.eye(dim2))
    return L


class FastEvolver:
    """Diagonalize L once, then evolve to any time via e^{λt}."""
    def __init__(self, L):
        self.dim = int(np.sqrt(L.shape[0]))
        print(f"    Diagonalizing {L.shape[0]}x{L.shape[0]} Liouvillian...", end=" ", flush=True)
        t0 = timer.time()
        evals, evecs = np.linalg.eig(L)
        print(f"done in {timer.time()-t0:.2f}s")
        self.evals = evals
        self.V = evecs
        self.Vinv = np.linalg.inv(evecs)
    
    def evolve(self, rho0, t):
        rho_vec = self.Vinv @ rho0.flatten()
        rho_vec = self.V @ (np.exp(self.evals * t) * rho_vec)
        rho = rho_vec.reshape(self.dim, self.dim)
        rho = (rho + rho.conj().T) / 2
        rho /= np.trace(rho)
        return rho
    
    def evolve_many(self, rho0, times):
        """Batch evolution - much faster for many timepoints."""
        rho_vec0 = self.Vinv @ rho0.flatten()
        states = []
        for t in times:
            rho_vec = self.V @ (np.exp(self.evals * t) * rho_vec0)
            rho = rho_vec.reshape(self.dim, self.dim)
            rho = (rho + rho.conj().T) / 2
            rho /= np.trace(rho)
            states.append(rho)
        return states


def ptrace(rho, N, keep):
    """Partial trace: keep specified qubits, trace out rest."""
    dim = 2**N; n_keep = len(keep); dim_k = 2**n_keep
    rho_r = np.zeros((dim_k, dim_k), dtype=complex)
    for i in range(dim):
        for j in range(dim):
            bi = [(i >> (N-1-k)) & 1 for k in range(N)]
            bj = [(j >> (N-1-k)) & 1 for k in range(N)]
            if all(bi[k] == bj[k] for k in range(N) if k not in keep):
                ki = sum(bi[keep[m]] << (n_keep-1-m) for m in range(n_keep))
                kj = sum(bj[keep[m]] << (n_keep-1-m) for m in range(n_keep))
                rho_r[ki, kj] += rho[i, j]
    return rho_r


def concurrence_2q(rho):
    YY = np.kron(Y, Y)
    R = rho @ (YY @ rho.conj() @ YY)
    ev = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(R), 0))))[::-1]
    return max(0, ev[0] - ev[1] - ev[2] - ev[3])


def l1_coh(rho):
    return np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))


def purity(rho):
    return np.real(np.trace(rho @ rho))


def make_bell_plus(q1, q2, N):
    """Bell+ state on qubits q1,q2, |0> on rest."""
    psi = np.zeros(2**N, dtype=complex)
    for bits in iproduct([0,1], repeat=N):
        if bits[q1] == bits[q2] and all(bits[k] == 0 for k in range(N) if k not in (q1, q2)):
            idx = sum(b << (N-1-i) for i, b in enumerate(bits))
            psi[idx] = 1.0 / np.sqrt(2)
    return psi


def get_bonds(topo, N):
    if topo == 'star': return [(0, i) for i in range(1, N)]
    if topo == 'chain': return [(i, i+1) for i in range(N-1)]
    if topo == 'ring': return [(i, (i+1) % N) for i in range(N)]
    if topo == 'complete': return [(i,j) for i in range(N) for j in range(i+1,N)]
    return []


# ================================================================
# EXPLORATION 1: Π and physical observables - deep dive
# ================================================================
def exploration_1():
    print("=" * 80)
    print("EXPLORATION 1: Π, CΨ, z*, u - seeking the Mandelbrot-Palindrome bridge")
    print("=" * 80)
    
    # Strategy: The Liouvillian eigenvalues are palindromic (proven).
    # CΨ(t) is a nonlinear function of the density matrix.
    # Question: does the palindromic spectral structure constrain CΨ(t)?
    # 
    # Approach: decompose ρ(t) in the Liouvillian eigenbasis.
    # Each eigenmode contributes with amplitude × e^{λt}.
    # Palindromic pairs have λ and -(λ+2Σγ). Their joint contribution
    # to any observable has a specific structure.
    #
    # We test across multiple initial states and topologies.
    
    configs = [
        # (N, topology, initial_state_label, initial_state_builder)
        (3, 'star',  'Bell_SA+0_B', lambda: np.kron((np.kron(up,up)+np.kron(dn,dn))/np.sqrt(2), up)),
        (3, 'star',  'Bell_AB+0_S', lambda: np.kron(up, (np.kron(up,up)+np.kron(dn,dn))/np.sqrt(2))),
        (3, 'chain', 'Bell_01+0_2', lambda: np.kron((np.kron(up,up)+np.kron(dn,dn))/np.sqrt(2), up)),
        (3, 'star',  '0+0',        lambda: np.kron(np.kron(up, plus), up)),
        (4, 'star',  'Bell_01+00',  lambda: np.kron((np.kron(up,up)+np.kron(dn,dn))/np.sqrt(2), np.kron(up, up))),
        (4, 'ring',  '0+0+',       lambda: np.kron(np.kron(up, plus), np.kron(up, plus))),
        (4, 'chain', '0+0+',       lambda: np.kron(np.kron(up, plus), np.kron(up, plus))),
    ]
    
    gamma = 0.05
    t_max = 15.0
    n_points = 2000
    times = np.linspace(0, t_max, n_points)
    
    results = []
    
    for N, topo, label, state_fn in configs:
        print(f"\n  --- N={N}, {topo}, {label} ---")
        bonds = get_bonds(topo, N)
        H = build_H(N, bonds)
        L = build_L(N, H, gamma)
        evolver = FastEvolver(L)
        
        psi0 = state_fn()
        rho0 = np.outer(psi0, psi0.conj())
        states = evolver.evolve_many(rho0, times)
        
        # Which pair to track? For star: AB=[1,2]. For others: try all.
        if N == 3:
            track_pairs = [(1,2)] if topo == 'star' else [(0,1), (1,2)]
        else:
            track_pairs = [(0,1), (0,2), (1,2)]
        
        for pair in track_pairs:
            C_vals, Psi_vals, CΨ_vals, zstar_vals = [], [], [], []
            
            for rho in states:
                rho_p = ptrace(rho, N, list(pair))
                C = concurrence_2q(rho_p)
                l1 = l1_coh(rho_p)
                psi = l1 / (2**len(pair) - 1)
                cpsi = C * psi
                zs = (1 - np.sqrt(max(0, 1-4*cpsi))) / 2 if 0 < cpsi <= 0.25 else np.nan
                
                C_vals.append(C)
                Psi_vals.append(psi)
                CΨ_vals.append(cpsi)
                zstar_vals.append(zs)
            
            CΨ_arr = np.array(CΨ_vals)
            zs_arr = np.array(zstar_vals)
            cpsi_max = np.max(CΨ_arr)
            t_max_cpsi = times[np.argmax(CΨ_arr)]
            
            crosses = cpsi_max > 0.25
            
            print(f"    Pair {pair}: CΨ_max={cpsi_max:.6f} at t={t_max_cpsi:.3f}"
                  f"  {'CROSSES 1/4' if crosses else ''}")
            
            # If crosses, find crossing times
            if crosses:
                for i in range(1, len(CΨ_arr)):
                    if CΨ_arr[i-1] > 0.25 and CΨ_arr[i] <= 0.25:
                        print(f"      Crosses DOWN at t={times[i]:.4f}: CΨ={CΨ_arr[i]:.6f}, z*={zs_arr[i]:.6f}")
                        break
                for i in range(1, len(CΨ_arr)):
                    if CΨ_arr[i-1] < 0.25 and CΨ_arr[i] >= 0.25:
                        print(f"      Crosses UP at t={times[i]:.4f}: CΨ={CΨ_arr[i]:.6f}")
            
            # z* statistics where defined
            valid = ~np.isnan(zs_arr) & (CΨ_arr > 0.01)
            if np.any(valid):
                print(f"      z* range: [{np.min(zs_arr[valid]):.4f}, {np.max(zs_arr[valid]):.4f}]")
                print(f"      z* at CΨ_max (if <1/4): {zs_arr[np.argmax(CΨ_arr)] if cpsi_max <= 0.25 else 'N/A (above 1/4)'}")
            
            results.append({
                'N': N, 'topo': topo, 'label': label, 'pair': pair,
                'cpsi_max': cpsi_max, 'crosses': crosses,
                'CΨ': CΨ_arr, 'zstar': zs_arr, 'times': times
            })
    
    # Cross-config analysis: is there a universal CΨ_max envelope?
    print(f"\n  --- Cross-configuration summary ---")
    print(f"  {'Config':>30} {'Pair':>6} {'CΨ_max':>8} {'Crosses':>8}")
    print(f"  {'-'*56}")
    for r in results:
        print(f"  {r['label']+'('+r['topo']+')':>30} {str(r['pair']):>6} "
              f"{r['cpsi_max']:>8.4f} {'YES' if r['crosses'] else 'no':>8}")
    
    return results


# ================================================================
# EXPLORATION 2: Ring near-miss - systematic scan
# ================================================================
def exploration_2():
    print(f"\n{'=' * 80}")
    print("EXPLORATION 2: Ring CΨ near-miss - systematic parameter scan")
    print("=" * 80)
    
    # The audit claims CΨ_max = 0.247 for ring neighbors.
    # We scan: initial states, gamma, J, N, topologies.
    # Question: is there a hard ceiling below 0.25 for ring neighbors?
    
    N = 4
    bonds_ring = get_bonds('ring', N)
    n_points = 1000
    
    # === Scan 1: Initial states ===
    print(f"\n  --- Scan 1: Different initial states (N=4 ring, γ=0.05) ---")
    
    initial_states = {
        '|0+0+⟩':      np.kron(np.kron(up, plus), np.kron(up, plus)),
        '|+0+0⟩':      np.kron(np.kron(plus, up), np.kron(plus, up)),
        '|++00⟩':      np.kron(np.kron(plus, plus), np.kron(up, up)),
        '|0000⟩+|1111⟩': (kron_list([up]*4) + kron_list([dn]*4)) / np.sqrt(2),  # GHZ
        'W':           (kron_list([dn,up,up,up]) + kron_list([up,dn,up,up]) + 
                        kron_list([up,up,dn,up]) + kron_list([up,up,up,dn])) / 2,
        '|+-+-⟩':      np.kron(np.kron(plus, minus), np.kron(plus, minus)),
        '|0+0-⟩':      np.kron(np.kron(up, plus), np.kron(up, minus)),
        'Bell01+Bell23': np.kron((np.kron(up,up)+np.kron(dn,dn))/np.sqrt(2),
                                 (np.kron(up,up)+np.kron(dn,dn))/np.sqrt(2)),
        'Bell01+00':    np.kron((np.kron(up,up)+np.kron(dn,dn))/np.sqrt(2),
                                np.kron(up, up)),
        '|++++⟩':      kron_list([plus]*4),
    }
    
    gamma = 0.05
    H = build_H(N, bonds_ring)
    L = build_L(N, H, gamma)
    evolver = FastEvolver(L)
    times = np.linspace(0, 20, n_points)
    
    print(f"  {'State':>20} {'01(ring)':>10} {'02(diag)':>10} {'12(ring)':>10} {'13(diag)':>10} {'Best':>10}")
    print(f"  {'-'*65}")
    
    global_max = 0
    global_max_config = ""
    
    for name, psi0 in initial_states.items():
        rho0 = np.outer(psi0, psi0.conj())
        states = evolver.evolve_many(rho0, times)
        
        pair_maxes = {}
        for pair, label in [((0,1),'01(ring)'), ((0,2),'02(diag)'), ((1,2),'12(ring)'), ((1,3),'13(diag)')]:
            cpsi_max = 0
            for rho in states:
                rho_p = ptrace(rho, N, list(pair))
                C = concurrence_2q(rho_p)
                l1 = l1_coh(rho_p)
                cpsi_max = max(cpsi_max, C * l1/3)
            pair_maxes[label] = cpsi_max
            if cpsi_max > global_max:
                global_max = cpsi_max
                global_max_config = f"{name}, pair {pair}"
        
        best = max(pair_maxes.values())
        best_pair = max(pair_maxes, key=pair_maxes.get)
        mark = " ← CROSSES" if best > 0.25 else ""
        print(f"  {name:>20} {pair_maxes['01(ring)']:>10.4f} {pair_maxes['02(diag)']:>10.4f} "
              f"{pair_maxes['12(ring)']:>10.4f} {pair_maxes['13(diag)']:>10.4f} "
              f"{best:>10.4f}{mark}")
    
    print(f"\n  Global max: {global_max:.6f} ({global_max_config})")
    print(f"  Gap to 1/4: {0.25 - global_max:.6f}")

    # === Scan 2: gamma sweep with best initial state ===
    print(f"\n  --- Scan 2: γ sweep (N=4 ring, |0+0+⟩) ---")
    print(f"  {'γ':>8} {'CΨ_max(01)':>12} {'CΨ_max(02)':>12} {'CΨ_max(13)':>12} {'Best':>10}")
    print(f"  {'-'*58}")
    
    for g in [0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.0]:
        L_g = build_L(N, H, g)
        ev_g = FastEvolver(L_g)
        psi0 = initial_states['|0+0+⟩']
        rho0 = np.outer(psi0, psi0.conj())
        states = ev_g.evolve_many(rho0, np.linspace(0, max(30, 5/g), n_points))
        
        maxes = {}
        for pair, label in [((0,1),'01'), ((0,2),'02'), ((1,3),'13')]:
            cpsi_max = 0
            for rho in states:
                rho_p = ptrace(rho, N, list(pair))
                cpsi_max = max(cpsi_max, concurrence_2q(rho_p) * l1_coh(rho_p)/3)
            maxes[label] = cpsi_max
        
        best = max(maxes.values())
        mark = " ← CROSSES" if best > 0.25 else ""
        print(f"  {g:>8.3f} {maxes['01']:>12.6f} {maxes['02']:>12.6f} "
              f"{maxes['13']:>12.6f} {best:>10.6f}{mark}")
    
    # === Scan 3: Compare ring vs chain vs star at N=4 ===
    print(f"\n  --- Scan 3: Topology comparison (N=4, γ=0.05, |0+0+⟩) ---")
    
    for topo in ['ring', 'chain', 'star', 'complete']:
        bonds = get_bonds(topo, N)
        H_t = build_H(N, bonds)
        L_t = build_L(N, H_t, 0.05)
        ev_t = FastEvolver(L_t)
        psi0 = initial_states['|0+0+⟩']
        rho0 = np.outer(psi0, psi0.conj())
        states = ev_t.evolve_many(rho0, np.linspace(0, 20, n_points))
        
        all_pairs = [(i,j) for i in range(N) for j in range(i+1, N)]
        print(f"\n    {topo}:")
        for pair in all_pairs:
            cpsi_max = 0; t_best = 0
            for k, rho in enumerate(states):
                rho_p = ptrace(rho, N, list(pair))
                cp = concurrence_2q(rho_p) * l1_coh(rho_p)/3
                if cp > cpsi_max:
                    cpsi_max = cp; t_best = times[k] if k < len(times) else 0
            gap = 0.25 - cpsi_max
            tag = "CROSSES" if cpsi_max > 0.25 else f"gap={gap:.4f}"
            print(f"      pair {pair}: CΨ_max={cpsi_max:.6f}  {tag}")
    
    # === Scan 4: N=5,6 ring ===
    print(f"\n  --- Scan 4: Larger rings (N=5,6) ---")
    for N_big in [5, 6]:
        bonds = get_bonds('ring', N_big)
        H_big = build_H(N_big, bonds)
        
        # |0+0+0...⟩ alternating
        qubit_states = []
        for i in range(N_big):
            qubit_states.append(up if i % 2 == 0 else plus)
        psi0 = qubit_states[0]
        for q in qubit_states[1:]:
            psi0 = np.kron(psi0, q)
        
        L_big = build_L(N_big, H_big, 0.05)
        ev_big = FastEvolver(L_big)
        rho0 = np.outer(psi0, psi0.conj())
        states = ev_big.evolve_many(rho0, np.linspace(0, 20, 500))
        
        print(f"\n    N={N_big} ring, |0+0+...⟩, γ=0.05:")
        # Check neighbor pairs and diagonal pairs
        for pair in [(0,1), (0,2), (1,2), (0, N_big//2)]:
            if pair[1] >= N_big:
                continue
            cpsi_max = 0
            for rho in states:
                rho_p = ptrace(rho, N_big, list(pair))
                cp = concurrence_2q(rho_p) * l1_coh(rho_p)/3
                cpsi_max = max(cpsi_max, cp)
            dist = min(abs(pair[1]-pair[0]), N_big - abs(pair[1]-pair[0]))
            gap = 0.25 - cpsi_max
            print(f"      pair {pair} (dist={dist}): CΨ_max={cpsi_max:.6f}, "
                  f"{'CROSSES' if cpsi_max > 0.25 else f'gap={gap:.4f}'}")


# ================================================================
# EXPLORATION 3: Echo effect - full characterization
# ================================================================
def exploration_3():
    print(f"\n{'=' * 80}")
    print("EXPLORATION 3: Echo effect - entanglement oscillation through mediator")
    print("=" * 80)
    
    # The echo: SA starts with Bell entanglement, SB starts empty.
    # Hamiltonian transfers entanglement SA → SB → SA → ...
    # Dephasing damps the oscillation.
    # Questions:
    # 1. What are the exact echo frequencies? (Bohr frequencies of H)
    # 2. How does damping scale? (palindromic decay rates)
    # 3. What happens at different J and γ?
    # 4. Does AB ever cross 1/4 from these echoes?
    # 5. N=4,5 star: does the echo pattern persist?
    
    # === Test 1: 3-qubit star, fine resolution, long time ===
    print(f"\n  --- Test 1: N=3 star, Bell_SA + |0⟩_B, γ=0.05 ---")
    N = 3; bonds = get_bonds('star', N)
    H = build_H(N, bonds)
    
    # Hamiltonian eigenvalues
    eig_H = np.linalg.eigvalsh(H)
    bohr_freqs = sorted(set([round(abs(eig_H[i]-eig_H[j]), 6) 
                             for i in range(len(eig_H)) for j in range(i+1, len(eig_H))
                             if abs(eig_H[i]-eig_H[j]) > 0.01]))
    bohr_periods = [2*np.pi/f for f in bohr_freqs]
    print(f"    H eigenvalues: {sorted(np.round(eig_H, 4))}")
    print(f"    Bohr frequencies: {bohr_freqs}")
    print(f"    Bohr periods: {[f'{p:.4f}' for p in bohr_periods]}")
    
    gamma = 0.05
    L = build_L(N, H, gamma)
    evolver = FastEvolver(L)
    
    # Long time, high resolution for Fourier
    times = np.linspace(0, 30, 6000)
    
    bell_SA = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
    psi0 = np.kron(bell_SA, up)
    rho0 = np.outer(psi0, psi0.conj())
    
    print(f"    Evolving {len(times)} timesteps...", flush=True)
    t0 = timer.time()
    states = evolver.evolve_many(rho0, times)
    print(f"    Done in {timer.time()-t0:.2f}s")
    
    pairs = {'SA': [0,1], 'SB': [0,2], 'AB': [1,2]}
    conc_data = {k: np.zeros(len(times)) for k in pairs}
    cpsi_data = {k: np.zeros(len(times)) for k in pairs}
    
    for i, rho in enumerate(states):
        for name, qubits in pairs.items():
            rho_p = ptrace(rho, N, qubits)
            C = concurrence_2q(rho_p)
            l1 = l1_coh(rho_p)
            conc_data[name][i] = C
            cpsi_data[name][i] = C * l1/3
    
    for name in pairs:
        c = conc_data[name]
        cp = cpsi_data[name]
        
        # Peaks
        peaks = []
        for i in range(2, len(c)-1):
            if c[i] > c[i-1] and c[i] > c[i+1] and c[i] > 0.003:
                peaks.append((times[i], c[i], cp[i]))
        
        print(f"\n    {name} pair:")
        print(f"      C(0)={c[0]:.4f}, CΨ(0)={cp[0]:.4f}")
        print(f"      C_max={np.max(c):.4f} at t={times[np.argmax(c)]:.4f}")
        print(f"      CΨ_max={np.max(cp):.6f}")
        
        if len(peaks) >= 2:
            # Extract intervals
            intervals = [peaks[i+1][0] - peaks[i][0] for i in range(len(peaks)-1)]
            print(f"      {len(peaks)} echo peaks found")
            print(f"      First 8 peaks: t = {[f'{p[0]:.3f}' for p in peaks[:8]]}")
            print(f"      First 8 amplitudes: C = {[f'{p[1]:.4f}' for p in peaks[:8]]}")
            print(f"      Intervals: {[f'{iv:.4f}' for iv in intervals[:8]]}")
            if intervals:
                print(f"      Mean interval: {np.mean(intervals[:10]):.4f}")
            
            # Envelope: fit exponential decay to peak amplitudes
            peak_times = np.array([p[0] for p in peaks])
            peak_amps = np.array([p[1] for p in peaks])
            if len(peak_amps) > 3 and np.all(peak_amps > 0):
                try:
                    log_amps = np.log(peak_amps)
                    coeffs = np.polyfit(peak_times, log_amps, 1)
                    decay_rate = -coeffs[0]
                    print(f"      Envelope decay rate: {decay_rate:.4f}")
                    print(f"      Compare: 2γ={2*gamma:.4f}, 8γ/3={8*gamma/3:.4f}, 10γ/3={10*gamma/3:.4f}")
                except:
                    pass
        
        # Fourier analysis of concurrence oscillation
        if np.max(c) > 0.01:
            dt = times[1] - times[0]
            c_detrended = c - np.mean(c)
            freqs = np.fft.rfftfreq(len(c), dt)
            spectrum = np.abs(np.fft.rfft(c_detrended))
            
            # Find dominant peaks
            spec_peaks = []
            for i in range(2, len(spectrum)-1):
                if spectrum[i] > spectrum[i-1] and spectrum[i] > spectrum[i+1]:
                    if spectrum[i] > 0.1 * np.max(spectrum):
                        spec_peaks.append((freqs[i], spectrum[i]))
            
            spec_peaks.sort(key=lambda x: -x[1])
            print(f"      Fourier peaks (freq, strength):")
            for f, s in spec_peaks[:5]:
                period = 1/f if f > 0 else np.inf
                # Match to Bohr frequencies
                matches = [bf for bf in bohr_freqs if abs(f*2*np.pi - bf) < 0.1]
                match_str = f" = ω/{2*np.pi:.3f} where ω={f*2*np.pi:.3f}" + (f" ≈ Bohr {matches[0]:.3f}" if matches else "")
                print(f"        f={f:.4f} (T={period:.4f}), strength={s:.2f}{match_str}")

    # === Test 2: γ dependence of echo ===
    print(f"\n  --- Test 2: Echo lifetime vs γ ---")
    print(f"  {'γ':>8} {'SA_peaks':>10} {'SB_C_max':>10} {'SB_peaks':>10} {'AB_CΨ_max':>12} {'AB_crosses':>12}")
    print(f"  {'-'*68}")
    
    for g in [0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50]:
        t_range = min(50, 10/max(g, 0.001))
        ts = np.linspace(0, t_range, 3000)
        L_g = build_L(N, H, g)
        ev_g = FastEvolver(L_g)
        states_g = ev_g.evolve_many(rho0, ts)
        
        data_g = {}
        for name, qubits in pairs.items():
            c_arr = np.zeros(len(ts))
            cp_arr = np.zeros(len(ts))
            for i, rho in enumerate(states_g):
                rho_p = ptrace(rho, N, qubits)
                C = concurrence_2q(rho_p)
                l1 = l1_coh(rho_p)
                c_arr[i] = C
                cp_arr[i] = C * l1/3
            
            peaks = sum(1 for i in range(2, len(c_arr)-1) 
                       if c_arr[i] > c_arr[i-1] and c_arr[i] > c_arr[i+1] and c_arr[i] > 0.003)
            data_g[name] = {'c_max': np.max(c_arr), 'peaks': peaks, 'cpsi_max': np.max(cp_arr)}
        
        ab_crosses = "YES" if data_g['AB']['cpsi_max'] > 0.25 else "no"
        print(f"  {g:>8.3f} {data_g['SA']['peaks']:>10} {data_g['SB']['c_max']:>10.4f} "
              f"{data_g['SB']['peaks']:>10} {data_g['AB']['cpsi_max']:>12.6f} {ab_crosses:>12}")
    
    # === Test 3: J dependence ===
    print(f"\n  --- Test 3: Echo vs J (γ=0.05) ---")
    print(f"  {'J':>8} {'SB_C_max':>10} {'AB_CΨ_max':>12} {'Echo_period':>12}")
    print(f"  {'-'*46}")
    
    for J in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
        H_J = build_H(N, bonds, J=J)
        L_J = build_L(N, H_J, 0.05)
        ev_J = FastEvolver(L_J)
        ts = np.linspace(0, max(30, 30/J), 3000)
        states_J = ev_J.evolve_many(rho0, ts)
        
        c_sb = np.array([concurrence_2q(ptrace(rho, N, [0,2])) for rho in states_J])
        cp_ab = np.array([concurrence_2q(ptrace(rho, N, [1,2])) * l1_coh(ptrace(rho, N, [1,2]))/3 
                         for rho in states_J])
        
        peaks_sb = []
        for i in range(2, len(c_sb)-1):
            if c_sb[i] > c_sb[i-1] and c_sb[i] > c_sb[i+1] and c_sb[i] > 0.003:
                peaks_sb.append(ts[i])
        
        period = np.mean(np.diff(peaks_sb[:5])) if len(peaks_sb) >= 2 else np.nan
        # Expected: period ~ 2π/(2J) = π/J (energy gap of Heisenberg bond)
        expected = np.pi / J
        
        print(f"  {J:>8.1f} {np.max(c_sb):>10.4f} {np.max(cp_ab):>12.6f} "
              f"{period:>12.4f}" + (f"  (π/J = {expected:.4f})" if not np.isnan(period) else ""))

    # === Test 4: N=4,5 star - does the echo persist? ===
    print(f"\n  --- Test 4: Echo in larger star topologies ---")
    
    for N_big in [4, 5]:
        bonds = get_bonds('star', N_big)
        H_big = build_H(N_big, bonds)
        L_big = build_L(N_big, H_big, 0.05)
        ev_big = FastEvolver(L_big)
        
        # Bell on S(0)-qubit1(1), |0> on rest
        bell_S1 = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
        rest = np.array([1.0], dtype=complex)
        for _ in range(N_big - 2):
            rest = np.kron(rest, up)
        psi0_big = np.kron(bell_S1, rest)
        rho0_big = np.outer(psi0_big, psi0_big.conj())
        
        ts = np.linspace(0, 20, 2000)
        states_big = ev_big.evolve_many(rho0_big, ts)
        
        print(f"\n    N={N_big} star, Bell_S0+|0⟩^{N_big-2}:")
        
        # Track S-qubit(0) with each leaf
        for leaf in range(1, N_big):
            c_arr = np.array([concurrence_2q(ptrace(rho, N_big, [0, leaf])) for rho in states_big])
            cp_arr = np.array([concurrence_2q(ptrace(rho, N_big, [0, leaf])) * 
                              l1_coh(ptrace(rho, N_big, [0, leaf]))/3 for rho in states_big])
            
            peaks = []
            for i in range(2, len(c_arr)-1):
                if c_arr[i] > c_arr[i-1] and c_arr[i] > c_arr[i+1] and c_arr[i] > 0.003:
                    peaks.append((ts[i], c_arr[i]))
            
            print(f"      S-{leaf}: C_max={np.max(c_arr):.4f}, CΨ_max={np.max(cp_arr):.6f}, "
                  f"peaks={len(peaks)}")
            if len(peaks) >= 2:
                intervals = [peaks[i+1][0] - peaks[i][0] for i in range(min(len(peaks)-1, 5))]
                print(f"             intervals: {[f'{iv:.3f}' for iv in intervals]}")
        
        # Track leaf-leaf pairs
        for i in range(1, min(N_big, 4)):
            for j in range(i+1, min(N_big, 4)):
                c_arr = np.array([concurrence_2q(ptrace(rho, N_big, [i, j])) for rho in states_big])
                cp_arr = np.array([c_arr[k] * l1_coh(ptrace(states_big[k], N_big, [i, j]))/3 
                                  for k in range(len(states_big))])
                print(f"      {i}-{j}: C_max={np.max(c_arr):.4f}, CΨ_max={np.max(cp_arr):.6f}")


# ================================================================
# MAIN
# ================================================================
if __name__ == '__main__':
    print(f"Start: {timer.strftime('%H:%M:%S')}")
    print(f"Machine: Python {sys.version.split()[0]}, NumPy {np.__version__}")
    print()
    
    t_total = timer.time()
    
    results_1 = exploration_1()
    exploration_2()
    exploration_3()
    
    print(f"\n{'=' * 80}")
    print(f"TOTAL TIME: {timer.time() - t_total:.1f}s")
    print(f"End: {timer.strftime('%H:%M:%S')}")
