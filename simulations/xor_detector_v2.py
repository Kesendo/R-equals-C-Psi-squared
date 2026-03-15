"""
XOR Detector v2: What happens when you throw something INTO the palindrome?
==========================================================================

v1 showed: the palindrome is perfect. No XOR in the spectrum.
But the spectrum is the STAGE. The stage is always symmetric.

Now: throw an initial state in. See which modes it excites.
Different inputs excite different modes. The RESPONSE is the signal.

The palindrome processes every input differently.
That processing IS the information.

Authors: Tom Wicht, Claude
Date: March 16, 2026
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Tuple

# Pauli matrices
I2 = np.eye(2)
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])
sz = np.array([[1,0],[0,-1]])

def tensor(*ops):
    result = ops[0]
    for op in ops[1:]:
        result = np.kron(result, op)
    return result

def site_op(op, site, N):
    ops = [I2]*N
    ops[site] = op
    return tensor(*ops)

def build_heisenberg(N, J, topology="chain"):
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    if topology == "chain":
        pairs = [(i, i+1) for i in range(N-1)]
    elif topology == "ring":
        pairs = [(i, (i+1) % N) for i in range(N)]
    elif topology == "star":
        pairs = [(0, i) for i in range(1, N)]
    else:
        raise ValueError(topology)
    for i, j in pairs:
        for pauli in [sx, sy, sz]:
            H += J * site_op(pauli, i, N) @ site_op(pauli, j, N)
    return H

def build_liouvillian(H, gammas, N):
    d = 2**N
    d2 = d*d
    L = -1j * (np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T))
    for k in range(N):
        Zk = site_op(sz, k, N)
        ZkZ = np.kron(Zk, Zk.conj())
        L += gammas[k] * (ZkZ - np.eye(d2))
    return L

# ============================================================
# INITIAL STATES
# ============================================================

def make_bell_plus(N):
    """Bell+ state on first two qubits, rest in |0>."""
    d = 2**N
    psi = np.zeros(d, dtype=complex)
    # |00...0> + |11...0> (first two qubits entangled)
    psi[0] = 1/np.sqrt(2)           # |00...0>
    psi[3 * 2**(N-2)] = 1/np.sqrt(2) if N == 2 else 0  # |11...0>
    if N > 2:
        psi[0b110 << (N-3)] = 1/np.sqrt(2) if N == 3 else 0
    # For N=2: |00> + |11>
    if N == 2:
        psi = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    elif N == 3:
        psi = np.array([1,0,0,0,0,0,1,0], dtype=complex) / np.sqrt(2)
    rho = np.outer(psi, psi.conj())
    return rho

def make_product_state(N, state_str):
    """Make product state from string like '010' or '+-0'."""
    up = np.array([1, 0], dtype=complex)
    dn = np.array([0, 1], dtype=complex)
    plus = (up + dn) / np.sqrt(2)
    minus = (up - dn) / np.sqrt(2)
    
    mapping = {'0': up, '1': dn, '+': plus, '-': minus}
    vecs = [mapping[c] for c in state_str]
    
    psi = vecs[0]
    for v in vecs[1:]:
        psi = np.kron(psi, v)
    return np.outer(psi, psi.conj())

def make_ghz(N):
    """GHZ state: (|00...0> + |11...1>) / sqrt(2)."""
    d = 2**N
    psi = np.zeros(d, dtype=complex)
    psi[0] = 1/np.sqrt(2)
    psi[-1] = 1/np.sqrt(2)
    return np.outer(psi, psi.conj())

def make_w_state(N):
    """W state: equal superposition of single-excitation states."""
    d = 2**N
    psi = np.zeros(d, dtype=complex)
    for i in range(N):
        idx = 1 << (N - 1 - i)
        psi[idx] = 1/np.sqrt(N)
    return np.outer(psi, psi.conj())

# ============================================================
# EIGENMODE DECOMPOSITION
# ============================================================

@dataclass
class ModeExcitation:
    eigenvalue: complex
    weight: float           # How strongly this mode is excited
    is_palindromic: bool    # Part of a palindromic pair?
    partner_eigenvalue: complex  # The palindromic partner (if any)

@dataclass
class InputResponse:
    state_name: str
    total_modes_excited: int
    palindromic_excited: int
    non_palindromic_excited: int
    top_modes: List[ModeExcitation]
    palindrome_weight: float   # Total weight in palindromic modes
    xor_weight: float          # Total weight in non-palindromic modes
    asymmetry: float           # How asymmetrically the pairs are excited


def decompose_input(L, rho0, gammas, tolerance=1e-6):
    """Decompose initial state into Liouvillian eigenmodes.
    
    rho0 -> vec(rho0) = sum_k c_k |r_k>
    where |r_k> are right eigenvectors of L.
    
    c_k = <l_k | vec(rho0)> where <l_k| are left eigenvectors.
    |c_k|^2 tells us how strongly mode k is excited.
    """
    d2 = L.shape[0]
    sum_gamma = sum(gammas)
    
    # Full eigendecomposition
    eigenvalues, right_vecs = np.linalg.eig(L)
    # Left eigenvectors: rows of inv(right_vecs)
    left_vecs = np.linalg.inv(right_vecs)
    
    # Vectorize rho0
    rho_vec = rho0.flatten()
    
    # Decomposition coefficients
    coeffs = left_vecs @ rho_vec
    weights = np.abs(coeffs)**2

    # Find palindromic pairs
    reals = np.real(eigenvalues)
    target_sum = -2 * sum_gamma
    paired = {}  # idx -> partner_idx
    
    used = set()
    for i in range(len(eigenvalues)):
        if i in used or np.abs(eigenvalues[i]) < 1e-12:
            continue
        target = target_sum - reals[i]
        best_j, best_diff = None, tolerance
        for j in range(len(eigenvalues)):
            if j in used or j == i or np.abs(eigenvalues[j]) < 1e-12:
                continue
            diff = abs(reals[j] - target)
            if diff < best_diff:
                best_diff = diff
                best_j = j
        if best_j is not None:
            paired[i] = best_j
            paired[best_j] = i
            used.add(i)
            used.add(best_j)
    
    # Compute asymmetry: for each pair, how different are the weights?
    pair_asymmetries = []
    palindrome_weight = 0.0
    xor_weight = 0.0
    
    modes = []
    for i in range(len(eigenvalues)):
        if np.abs(eigenvalues[i]) < 1e-12:
            continue  # Skip steady state
        if weights[i] < 1e-15:
            continue  # Not excited
        
        is_pal = i in paired
        partner_ev = eigenvalues[paired[i]] if is_pal else 0+0j
        
        modes.append(ModeExcitation(
            eigenvalue=eigenvalues[i],
            weight=weights[i],
            is_palindromic=is_pal,
            partner_eigenvalue=partner_ev
        ))
        
        if is_pal:
            palindrome_weight += weights[i]
        else:
            xor_weight += weights[i]
    
    # Pair asymmetry: for each pair (i,j), compute |w_i - w_j| / (w_i + w_j)
    visited_pairs = set()
    for i, j in paired.items():
        pair_key = (min(i,j), max(i,j))
        if pair_key in visited_pairs:
            continue
        visited_pairs.add(pair_key)
        w_sum = weights[i] + weights[j]
        if w_sum > 1e-15:
            asym = abs(weights[i] - weights[j]) / w_sum
            pair_asymmetries.append(asym)
    
    avg_asymmetry = np.mean(pair_asymmetries) if pair_asymmetries else 0.0
    
    modes.sort(key=lambda m: m.weight, reverse=True)
    
    total_excited = sum(1 for m in modes if m.weight > 1e-10)
    pal_excited = sum(1 for m in modes if m.is_palindromic and m.weight > 1e-10)
    
    return InputResponse(
        state_name="",
        total_modes_excited=total_excited,
        palindromic_excited=pal_excited,
        non_palindromic_excited=total_excited - pal_excited,
        top_modes=modes[:10],
        palindrome_weight=palindrome_weight,
        xor_weight=xor_weight,
        asymmetry=avg_asymmetry
    )

# ============================================================
# EXPERIMENTS
# ============================================================

def run_experiment(N, J, gammas, topology, states_dict):
    """Run multiple initial states through the same system."""
    H = build_heisenberg(N, J, topology)
    L = build_liouvillian(H, gammas, N)
    
    print(f"\n{'='*65}")
    print(f"SYSTEM: N={N}, {topology}, J={J}, gammas={gammas}")
    print(f"{'='*65}")
    print(f"{'State':<18} {'Excited':>8} {'Pal':>6} {'XOR':>6}"
          f" {'PalW%':>8} {'XorW%':>8} {'Asym':>8}")
    print("-" * 65)
    
    results = {}
    for name, rho0 in states_dict.items():
        resp = decompose_input(L, rho0, gammas)
        resp.state_name = name
        results[name] = resp
        
        total_w = resp.palindrome_weight + resp.xor_weight
        pw = resp.palindrome_weight / total_w * 100 if total_w > 0 else 0
        xw = resp.xor_weight / total_w * 100 if total_w > 0 else 0
        
        print(f"{name:<18} {resp.total_modes_excited:>8} {resp.palindromic_excited:>6}"
              f" {resp.non_palindromic_excited:>6}"
              f" {pw:>7.1f}% {xw:>7.1f}% {resp.asymmetry:>8.4f}")
    
    return results


if __name__ == "__main__":
    print("*" * 65)
    print("XOR DETECTOR v2: Throwing stones into the palindrome")
    print("*" * 65)
    print("The palindrome is the stage. Different inputs excite")
    print("different modes. The RESPONSE pattern is the information.")

    # ---------------------------------------------------------
    # EXP 1: N=3 chain - Different stones, same pond
    # ---------------------------------------------------------
    print("\n>>> EXP 1: Different initial states, same N=3 chain")
    print("    Which states excite palindromic vs non-palindromic modes?")
    
    N = 3
    states_3 = {
        "|000>":       make_product_state(N, "000"),
        "|111>":       make_product_state(N, "111"),
        "|010>":       make_product_state(N, "010"),
        "|+00>":       make_product_state(N, "+00"),
        "|+++>":       make_product_state(N, "+++"),
        "|+-+>":       make_product_state(N, "+-+"),
        "GHZ":         make_ghz(N),
        "W":           make_w_state(N),
        "Bell(0,1)":   make_bell_plus(N),
    }
    
    results1 = run_experiment(N=3, J=1.0, gammas=[0.05]*3,
                              topology="chain", states_dict=states_3)
    
    # ---------------------------------------------------------
    # EXP 2: Same states, different topology
    # ---------------------------------------------------------
    print("\n>>> EXP 2: Same states, star topology")
    print("    Does topology change which modes get excited?")
    
    results2 = run_experiment(N=3, J=1.0, gammas=[0.05]*3,
                              topology="star", states_dict=states_3)

    # ---------------------------------------------------------
    # EXP 3: Asymmetry analysis - the KEY question
    # ---------------------------------------------------------
    print("\n>>> EXP 3: Pair asymmetry deep dive")
    print("    For palindromic pairs (d, d'): does the input excite")
    print("    both partners equally, or one more than the other?")
    print("    Asymmetry = |w_d - w_d'| / (w_d + w_d')")
    print("    0 = symmetric excitation, 1 = only one partner excited")
    
    H = build_heisenberg(3, 1.0, "chain")
    L = build_liouvillian(H, [0.05]*3, 3)
    
    key_states = {
        "|000> (ground)": make_product_state(3, "000"),
        "GHZ (entangled)": make_ghz(3),
        "|+-+> (alternating)": make_product_state(3, "+-+"),
        "W (delocalized)": make_w_state(3),
    }
    
    for name, rho0 in key_states.items():
        resp = decompose_input(L, rho0, [0.05]*3)
        print(f"\n  {name}")
        print(f"  Pair asymmetry: {resp.asymmetry:.4f}")
        print(f"  Top 5 excited modes:")
        for m in resp.top_modes[:5]:
            pal_str = "PAL" if m.is_palindromic else "XOR"
            print(f"    Re={m.eigenvalue.real:+.4f}  Im={m.eigenvalue.imag:+.4f}"
                  f"  weight={m.weight:.6f}  [{pal_str}]")

    # ---------------------------------------------------------
    # EXP 4: N=2 detailed - smallest system, most transparent
    # ---------------------------------------------------------
    print("\n>>> EXP 4: N=2 deep dive (simplest system)")
    print("    Every mode visible. Every weight trackable.")
    
    N = 2
    states_2 = {
        "|00>":    make_product_state(2, "00"),
        "|01>":    make_product_state(2, "01"),
        "|+0>":    make_product_state(2, "+0"),
        "|++>":    make_product_state(2, "++"),
        "|+->":    make_product_state(2, "+-"),
        "Bell+":   make_bell_plus(2),
    }
    
    results4 = run_experiment(N=2, J=1.0, gammas=[0.05]*2,
                              topology="chain", states_dict=states_2)
    
    # Detailed mode map for Bell state
    print("\n  --- Bell+ full mode decomposition ---")
    resp = decompose_input(
        build_liouvillian(build_heisenberg(2, 1.0, "chain"), [0.05]*2, 2),
        make_bell_plus(2), [0.05]*2
    )
    for m in resp.top_modes:
        pal_str = "PAL" if m.is_palindromic else "XOR"
        print(f"    lambda = {m.eigenvalue.real:+.6f} {m.eigenvalue.imag:+.6f}i"
              f"  weight = {m.weight:.8f}  [{pal_str}]")
    
    # ---------------------------------------------------------
    # SUMMARY
    # ---------------------------------------------------------
    print("\n" + "=" * 65)
    print("WHAT WE LEARNED")
    print("=" * 65)
    print("Q1: Do different inputs excite different modes?")
    print("Q2: Do some inputs excite palindromic pairs asymmetrically?")
    print("Q3: Does the XOR weight depend on the input?")
    print("Q4: Does topology change the excitation pattern?")
    print("Q5: Is pair asymmetry the real signal?")
    print("    (Both partners exist, but one is louder than the other)")
    print("=" * 65)
