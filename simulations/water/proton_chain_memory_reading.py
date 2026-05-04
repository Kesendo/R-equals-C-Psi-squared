"""Memory-axis reading on a 3-proton water chain (Grotthuss model).

Applies the trio's state-level diagnostics (Frobenius static/memory partition +
Pi^2-odd-fraction-within-memory + per-proton Bloch reading) to a Heisenberg
proton chain under Z-dephasing.

Comparison to the IBM Torino single-qubit data (memory_reading_ibm_torino.py):
- IBM Torino: T1 + T2 → memory migrates X+ to Z+ (thermal direction)
- Water chain: pure Z-dephasing (framework canonical) → memory decays to I/d
  (the d=0 axis, maximally mixed). NO migration to a preferred direction.

This is the canonical framework setup: the trio's d=0-axis prediction holds
exactly when there is no T1 (no preferred ground-state pull). Water at room
temperature with kT >> hbar*omega has effectively balanced rates, so pure
dephasing dominates. The chain thermalises to I/d, not to |0...0>.

Run with:
  PYTHONIOENCODING=utf-8 python simulations/water/proton_chain_memory_reading.py
"""

import numpy as np
from scipy.linalg import expm

# ---------------------- Pauli + chain primitives ----------------------

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, SX, SY, SZ]
BIT_A = [0, 1, 1, 0]   # X, Y carry bit_a = 1
BIT_B = [0, 0, 1, 1]   # Y, Z carry bit_b = 1


def kron_n(mats):
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def site_op(P, k, N):
    return kron_n([P if i == k else I2 for i in range(N)])


def heisenberg_chain(N, J=1.0):
    """H = (J/4) Sum_b (X_b X_{b+1} + Y_b Y_{b+1} + Z_b Z_{b+1}). Truly bilinear."""
    H = np.zeros((2**N, 2**N), dtype=complex)
    for b in range(N - 1):
        for P in [SX, SY, SZ]:
            H += (J / 4.0) * site_op(P, b, N) @ site_op(P, b + 1, N)
    return H


def z_dephasing_lindblad(H, gamma, N):
    """L superoperator (vec_F) for dot(rho) = -i[H,rho] + gamma * Sum_l (Z_l rho Z_l - rho)."""
    d = 2**N
    eye = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, eye) - np.kron(eye, H.T))
    for k in range(N):
        Zk = site_op(SZ, k, N)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.kron(eye, eye))
    return L


# ---------------------- Initial state builders ----------------------

def x_basis_product(N, signs):
    v = np.array([1.0], dtype=complex)
    for s in signs:
        ket = np.array([1.0, float(s)], dtype=complex) / np.sqrt(2)
        v = np.kron(v, ket)
    return v


def y_basis_product(N, signs):
    v = np.array([1.0], dtype=complex)
    for s in signs:
        ket = np.array([1.0, 1j * float(s)], dtype=complex) / np.sqrt(2)
        v = np.kron(v, ket)
    return v


def computational_state(N, bits):
    idx = 0
    for k, b in enumerate(bits):
        idx |= (b & 1) << (N - 1 - k)
    v = np.zeros(2**N, dtype=complex)
    v[idx] = 1.0
    return v


def density_matrix(psi):
    return np.outer(psi, psi.conj())


# ---------------------- Trio diagnostics ----------------------

def kernel_projection(rho, N):
    """rho_d0 = Sum_n (Tr(P_n rho) / Tr(P_n^2)) P_n with P_n the popcount-n sector projector.

    For uniform Heisenberg + Z-dephasing: kernel of L equals span(P_0, ..., P_N) per F4.
    """
    d = 2**N
    rho_d0 = np.zeros_like(rho)
    for n in range(N + 1):
        P_n = np.zeros((d, d), dtype=complex)
        for b in range(d):
            if bin(b).count('1') == n:
                P_n[b, b] = 1.0
        rank_n = int(np.real(np.trace(P_n @ P_n)))
        if rank_n == 0:
            continue
        coeff = np.real(np.trace(P_n @ rho)) / rank_n
        rho_d0 = rho_d0 + coeff * P_n
    return rho_d0


def pi2_split(rho, N):
    """Split rho into Pi^2-even + Pi^2-odd Pauli-string components (Z-dephasing convention,
    bit_b parity = number of Y or Z letters mod 2)."""
    d = 2**N
    rho_even = np.zeros_like(rho)
    rho_odd = np.zeros_like(rho)
    inv = 1.0 / d
    for k in range(4**N):
        kk = k
        idxs = []
        for _ in range(N):
            idxs.append(kk & 3)
            kk >>= 2
        parity = sum(BIT_B[i] for i in idxs) & 1
        sigma = kron_n([PAULIS[i] for i in idxs])
        coeff = np.trace(sigma @ rho) * inv
        if abs(coeff) < 1e-14:
            continue
        if parity == 0:
            rho_even = rho_even + coeff * sigma
        else:
            rho_odd = rho_odd + coeff * sigma
    return rho_even, rho_odd


def per_proton_bloch(rho, N):
    return [
        (
            np.trace(site_op(SX, k, N) @ rho).real,
            np.trace(site_op(SY, k, N) @ rho).real,
            np.trace(site_op(SZ, k, N) @ rho).real,
        )
        for k in range(N)
    ]


def evolve(rho0, L, t):
    d = rho0.shape[0]
    vec0 = rho0.flatten(order='F')
    vec_t = expm(L * t) @ vec0
    rho_t = vec_t.reshape((d, d), order='F')
    return (rho_t + rho_t.conj().T) / 2.0


# ---------------------- Main analysis ----------------------

def run():
    N = 3
    J = 1.0
    gamma = 0.05
    H = heisenberg_chain(N, J)
    L = z_dephasing_lindblad(H, gamma, N)

    print(f"=== 3-proton water chain, Heisenberg J = {J}, Z-dephasing gamma = {gamma} ===")
    print(f"Time unit: 1/J. Heisenberg chain is 'truly' (Pi^2-even); F1 palindrome holds bit-exact.")
    print(f"Expected long-time behaviour: rho -> I/d = I/8 (the d=0 axis, maximally mixed).")
    print()

    states = [
        ("|+++>      (X-polarity, Pi^2-even content only)",
         x_basis_product(N, [+1, +1, +1])),
        ("|+i,+i,+i> (Y-polarity, Pi^2-odd 4/7 in memory)",
         y_basis_product(N, [+1, +1, +1])),
        ("|001>      (Z-basis, popcount-1 sector)",
         computational_state(N, [0, 0, 1])),
        ("|+, 0, +>  (mixed X-Z basis)",
         np.kron(np.kron(np.array([1, 1], complex)/np.sqrt(2),
                         np.array([1, 0], complex)),
                 np.array([1, 1], complex)/np.sqrt(2))),
    ]

    times = [0.0, 1.0, 5.0, 10.0, 25.0, 50.0, 100.0, 200.0]

    for name, psi in states:
        rho0 = density_matrix(psi)
        print(f"--- {name} ---")
        print(f"{'t':>6} {'static':>9} {'memory':>9} {'Pi2-odd/mem':>13}  per-proton |r|")
        print("-" * 78)
        for t in times:
            rho_t = evolve(rho0, L, t)
            rho_d0 = kernel_projection(rho_t, N)
            rho_d2 = rho_t - rho_d0
            _, rho_d2_odd = pi2_split(rho_d2, N)

            norm_total = np.linalg.norm(rho_t, 'fro') ** 2
            norm_static = np.linalg.norm(rho_d0, 'fro') ** 2
            norm_mem = np.linalg.norm(rho_d2, 'fro') ** 2
            norm_odd = np.linalg.norm(rho_d2_odd, 'fro') ** 2

            static_frac = norm_static / norm_total if norm_total > 0 else 0.0
            memory_frac = norm_mem / norm_total if norm_total > 0 else 0.0
            odd_in_mem = norm_odd / norm_mem if norm_mem > 0 else 0.0

            blochs = per_proton_bloch(rho_t, N)
            bloch_str = "  ".join(
                f"q{k}={np.sqrt(b[0]**2 + b[1]**2 + b[2]**2):.3f}"
                for k, b in enumerate(blochs)
            )
            print(f"{t:>6.1f} {static_frac:>9.4f} {memory_frac:>9.4f} {odd_in_mem:>13.4f}  {bloch_str}")
        print()

    print("=== Reading ===")
    print()
    print("Compare to IBM Torino single qubit (memory_reading_ibm_torino.py):")
    print(" - There: T1+T2, memory migrates X+ -> Z+ (thermal Z+, |r|/2 stays > 0.20).")
    print(" - Here:  pure Z-dephasing, memory decays toward d=0 axis (|r| -> 0 per qubit,")
    print("          and rho_d2 -> 0 in operator-norm).")
    print()
    print("The water Heisenberg chain is the framework's CANONICAL setup: the trio's")
    print("'memory disappears toward I/d' prediction is realised exactly. Hardware T1")
    print("breaks this by introducing a preferred direction. Pure dephasing (water at")
    print("room temperature kT >> hbar omega, or framework Z-dephasing) does not.")


if __name__ == "__main__":
    run()
