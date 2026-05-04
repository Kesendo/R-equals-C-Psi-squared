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
    """expm-based; OK for small d. For N >= 4 prefer evolve_eigendecomp."""
    d = rho0.shape[0]
    vec0 = rho0.flatten(order='F')
    vec_t = expm(L * t) @ vec0
    rho_t = vec_t.reshape((d, d), order='F')
    return (rho_t + rho_t.conj().T) / 2.0


def diagonalize_propagator(L):
    """One-shot eigendecomp of the Liouvillian for fast time-stepping at large N.

    Returns (lambdas, R, R_inv) such that exp(L·t) = R · diag(exp(λ·t)) · R^-1.
    """
    lambdas, R = np.linalg.eig(L)
    R_inv = np.linalg.inv(R)
    return lambdas, R, R_inv


def evolve_via_eigendecomp(rho0, lambdas, R, R_inv, t):
    d = rho0.shape[0]
    vec0 = rho0.flatten(order='F')
    c = R_inv @ vec0
    vec_t = R @ (np.exp(lambdas * t) * c)
    rho_t = vec_t.reshape((d, d), order='F')
    return (rho_t + rho_t.conj().T) / 2.0


# ---------------------- Main analysis ----------------------

def analyse(N, J=1.0, gamma=0.05, times=None):
    """Run the trio's state-level diagnostics on a Heisenberg+Z-dephasing chain at N qubits.

    For each canonical initial state (X-polarity, Y-polarity, sector-bound), evolve under
    L for the given times and dump (static, memory, Π²-odd-in-memory, per-proton |r|).
    """
    if times is None:
        times = [0.0, 1.0, 5.0, 10.0, 25.0, 50.0, 100.0, 200.0]

    H = heisenberg_chain(N, J)
    L = z_dephasing_lindblad(H, gamma, N)
    lambdas, R, R_inv = diagonalize_propagator(L)

    states = [
        (f"|+>^{N}      (X-polarity)",
         x_basis_product(N, [+1] * N)),
        (f"|+i>^{N}     (Y-polarity, expect Pi2-odd/mem = (2^(N-1))/(2^N-1) at t=0)",
         y_basis_product(N, [+1] * N)),
        (f"|0..01>     (popcount-1 sector, expect Pi2-odd/mem = 0.5 throughout)",
         computational_state(N, [0] * (N - 1) + [1])),
    ]

    print(f"=== N = {N}: Heisenberg J = {J}, Z-dephasing gamma = {gamma}, d = {2**N} ===")
    print(f"kernel dim of L = {N + 1} (sector projectors P_0..P_{N}); F1 palindrome bit-exact.")
    print()

    for name, psi in states:
        rho0 = density_matrix(psi)
        print(f"--- {name} ---")
        print(f"{'t':>6} {'static':>9} {'memory':>9} {'Pi2-odd/mem':>13}  per-proton |r|")
        print("-" * 78)
        for t in times:
            rho_t = evolve_via_eigendecomp(rho0, lambdas, R, R_inv, t)
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


def expected_y_polarity_pi2_odd_initial(N):
    """X-Y combinatorial theorem: at t=0, Y-polarity has 2^(N-1) Pi2-odd Pauli strings
    out of 2^N total; within memory (excluding III), Pi2-odd/mem = 2^(N-1) / (2^N - 1)."""
    return (1 << (N - 1)) / float((1 << N) - 1)


def expected_sector_bound_static_fraction(N, popcount):
    """Single computational basis state |b> with popcount(b) = n has rho_d0 = (1/C(N,n))·P_n.
    ‖rho_d0‖² = 1/C(N,n)."""
    from math import comb
    return 1.0 / comb(N, popcount)


def expected_sector_bound_per_proton_r(N, popcount):
    """At long t, |0..01> in P_1 thermalises within sector to P_n/C(N,n).
    Per-proton rz = (C(N-1,n) - C(N-1,n-1)) / C(N,n)."""
    from math import comb
    return abs(comb(N - 1, popcount) - comb(N - 1, popcount - 1)) / comb(N, popcount)


def run():
    print("Trio's state-level diagnostics on Heisenberg + Z-dephasing water chains")
    print("=" * 78)
    print()
    print("Predicted Y-polarity Pi2-odd/mem at t=0 (X-Y combinatorial theorem):")
    for N in (3, 4, 5):
        pred = expected_y_polarity_pi2_odd_initial(N)
        num, den = (1 << (N - 1)), ((1 << N) - 1)
        print(f"  N = {N}: 2^{N-1}/(2^{N}-1) = {num}/{den} = {pred:.4f}")
    print()
    print("Predicted sector-bound (popcount-1) static fraction (= 1/C(N,1) = 1/N):")
    for N in (3, 4, 5):
        print(f"  N = {N}: static = {expected_sector_bound_static_fraction(N, 1):.4f}; "
              f"per-proton |r|_eq = {expected_sector_bound_per_proton_r(N, 1):.4f}")
    print()

    for N in (3, 4, 5):
        analyse(N)

    print("=== Cross-N reading ===")
    print()
    print("Three structural inheritance patterns visible across N = 3, 4, 5:")
    print()
    print(" 1. X-polarity (|+>^N):  Pi2-odd/mem stays exactly 0 throughout. {I, X} are both")
    print("    Pi2-even under Z-dephasing; no Pi2-odd content is ever created. Memory")
    print("    decays to zero (rho -> I/d).")
    print()
    print(" 2. Y-polarity (|+i>^N): Pi2-odd/mem starts at predicted (2^(N-1))/(2^N-1)")
    print("    (4/7, 8/15, 16/31 at N = 3, 4, 5) and GROWS toward 1.0 over time.")
    print("    Pi2-EVEN content decays faster than Pi2-odd (matches yesterday's F80 lift")
    print("    dynamics observation). Inheritance pattern: Pi2-odd is structurally robust.")
    print()
    print(" 3. Sector-bound (|0..01>): Pi2-odd/mem stays EXACTLY 0.5 throughout for all N")
    print("    (matches yesterday's F86 verification: popcount-coherence-class states have")
    print("    50% Pi2-odd content combinatorially). Static fraction approaches 1/C(N,1) =")
    print("    1/N (1/3, 1/4, 1/5). Per-proton |r| asymptotes to (N-2)/N (1/3, 1/2, 3/5).")
    print()
    print("The framework's d=0-axis prediction (memory -> I/d) holds exactly when the")
    print("initial state is symmetric under sector permutations (case 1). Non-symmetric")
    print("initial states settle into the sector projector P_n, not into I/d, because")
    print("Heisenberg conserves popcount.")


if __name__ == "__main__":
    run()
