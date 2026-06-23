"""Compute both XOR (parity) and depth (full n_Y count / popcount distribution)
for each slow Liouvillian eigenmode on the N = 4 painter ring.

The Painter alternation showed a bit-exact Z₂ sectorization at the per-site
projection level: each slow mode is purely Y-axis or purely non-Y-axis at
each site. That's the parity reading (XOR is involutive; a single bit per
mode tells us which sector).

This script computes the second reading: the depth. Two complementary
depth diagnostics, one on the operator side, one on the state side.

  (A) Operator-side: decompose the full 4^N-dim eigenmode into the Pauli-
      string basis, group strings by n_Y (= count of Y letters in the
      string, 0..N), report |coeff|² per n_Y bin per mode. The parity
      reading is n_Y mod 2; the depth reading is the distribution within
      each parity sector.

  (B) State-side: treat the eigenmode as a d × d matrix ρ[i,j], bin |ρ[i,j]|²
      by popcount(i XOR j) (= Hamming distance between bra-ket indices in
      the computational basis, 0..N). This is the "Z-dephasing depth" each
      density-matrix element is exposed to: damping rate is 2γ·popcount(i
      XOR j) per LindbladPropagator.cs:44. Diagonal entries (popcount = 0)
      don't decay; off-diagonal entries decay faster at larger popcount.

Reading: for each slow mode, both decompositions give a histogram. XOR
tells us "which tower"; depth tells us "where in the tower". Both
together reconstruct the mode's full Pauli-string profile.
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))


PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def pauli_op(letters):
    op = PAULI[letters[0]]
    for L in letters[1:]:
        op = np.kron(op, PAULI[L])
    return op


def site_op(N, site, letter):
    letters = ["I"] * N
    letters[site] = letter
    return pauli_op(letters)


def two_site_op(N, a, b, la, lb):
    letters = ["I"] * N
    letters[a] = la
    letters[b] = lb
    return pauli_op(letters)


def hueckel_ring_H(N):
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for a in range(N):
        b = (a + 1) % N
        H = H + two_site_op(N, a, b, "X", "X") + two_site_op(N, a, b, "Y", "Y")
    return H


def zeeman_y_total(N):
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for l in range(N):
        H = H + site_op(N, l, "Y")
    return H


def commutator_superop_vec(H):
    d = H.shape[0]
    I = np.eye(d, dtype=complex)
    return -1j * (np.kron(H, I) - np.kron(I, H.T))


def dissipator_superop_vec(c):
    d = c.shape[0]
    I = np.eye(d, dtype=complex)
    c_dag_c = c.conj().T @ c
    return np.kron(c, c.conj()) - 0.5 * (np.kron(c_dag_c, I) + np.kron(I, c_dag_c.T))


def lindbladian_vec(H, c_list, gammas):
    L = commutator_superop_vec(H)
    for c, g in zip(c_list, gammas):
        L = L + g * dissipator_superop_vec(c)
    return L


def unvec_column(v, d):
    return v.reshape(d, d, order="F")


def pauli_string_decomp(rho_full, N):
    """Hilbert-Schmidt decomposition of a d × d operator in the 4^N Pauli
    string basis. Returns list of (string_label, n_Y, |coeff|²) tuples."""
    d = 2**N
    result = []
    for letters in itertools.product(["I", "X", "Y", "Z"], repeat=N):
        sigma = pauli_op(letters)
        coeff = np.trace(sigma.conj().T @ rho_full) / d
        n_Y = sum(1 for L in letters if L == "Y")
        result.append(("".join(letters), n_Y, abs(coeff) ** 2))
    return result


def n_Y_histogram(decomp, N):
    """Sum |coeff|² per n_Y bin (0..N)."""
    bins = [0.0] * (N + 1)
    for _, n_Y, w in decomp:
        bins[n_Y] += w
    return bins


def xor_popcount_histogram(rho_full, N):
    """Bin |ρ[i,j]|² by popcount(i XOR j), the Z-deph damping-depth per element."""
    d = 2**N
    bins = [0.0] * (N + 1)
    for i in range(d):
        for j in range(d):
            xor = i ^ j
            pc = bin(xor).count("1")
            bins[pc] += abs(rho_full[i, j]) ** 2
    return bins


def run(N=4, h_zeeman=0.5, gamma=1.0, n_slowest=8):
    print()
    print("=" * 100)
    print(f"N = {N} carbon ring,  h_y = {h_zeeman},  γ_Z-deph = {gamma}")
    print("=" * 100)
    print()

    H = hueckel_ring_H(N) + h_zeeman * zeeman_y_total(N)
    c_holstein = [site_op(N, l, "Z") for l in range(N)]
    L = lindbladian_vec(H, c_holstein, [gamma] * N)
    d = 2**N

    eigvals, eigvecs = np.linalg.eig(L)
    order = np.argsort(eigvals.real)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    # ---- (A) Operator-side: n_Y histogram per mode ----
    print("(A) Operator-side: n_Y depth distribution per slow mode")
    print("-" * 100)
    print(f"For each mode, |coeff|² summed over Pauli strings σ with n_Y(σ) = k.")
    print(f"XOR-reading = is the weight on even k or odd k? (parity of n_Y)")
    print(f"Depth-reading = within the parity, how distributed across k = 0, 2, 4 (even) or k = 1, 3 (odd)?")
    print()
    header = f"  {'k':>3}  {'Re(λ)':>10}  {'Im(λ)':>10}  " + "  ".join(
        [f"n_Y={i:>2}" for i in range(N + 1)]
    ) + "  parity"
    print(header)
    print("  " + "-" * (len(header) - 2))
    for k in range(n_slowest):
        v = eigvecs[:, k]
        rho = unvec_column(v, d)
        # normalize so we read fractions (not absolute weights)
        decomp = pauli_string_decomp(rho, N)
        bins = n_Y_histogram(decomp, N)
        total = sum(bins)
        if total < 1e-12:
            print(f"  {k:>3}  {eigvals[k].real:>+10.4f}  {eigvals[k].imag:>+10.4f}  (zero)")
            continue
        bins_n = [b / total for b in bins]
        even = sum(bins_n[i] for i in range(N + 1) if i % 2 == 0)
        odd = sum(bins_n[i] for i in range(N + 1) if i % 2 == 1)
        if odd < 1e-9:
            parity = "EVEN (n_Y mod 2 = 0)"
        elif even < 1e-9:
            parity = "ODD  (n_Y mod 2 = 1)"
        else:
            parity = f"mixed (even={even:.4f}, odd={odd:.4f})"
        bin_str = "  ".join([f"{b:>5.3f}" for b in bins_n])
        print(f"  {k:>3}  {eigvals[k].real:>+10.4f}  {eigvals[k].imag:>+10.4f}  {bin_str}  {parity}")
    print()

    # ---- (B) State-side: Z-deph XOR depth per mode ----
    print("(B) State-side: popcount(i XOR j) depth distribution per slow mode")
    print("-" * 100)
    print(f"For each mode, |ρ[i,j]|² summed over basis pairs with popcount(i XOR j) = k.")
    print(f"This is the Z-deph storage reading: each off-diagonal element decays at rate")
    print(f"2γ · popcount(i XOR j); diagonal (k = 0) does not decay.")
    print()
    header = f"  {'k':>3}  {'Re(λ)':>10}  {'Im(λ)':>10}  " + "  ".join(
        [f"pop={i:>2}" for i in range(N + 1)]
    ) + "  predicted Z-decay-mean"
    print(header)
    print("  " + "-" * (len(header) - 2))
    for k in range(n_slowest):
        v = eigvecs[:, k]
        rho = unvec_column(v, d)
        bins = xor_popcount_histogram(rho, N)
        total = sum(bins)
        if total < 1e-12:
            print(f"  {k:>3}  {eigvals[k].real:>+10.4f}  {eigvals[k].imag:>+10.4f}  (zero)")
            continue
        bins_n = [b / total for b in bins]
        # predicted Z-only-deph decay rate of this distribution: 2γ · weighted mean popcount
        mean_pc = sum(i * bins_n[i] for i in range(N + 1))
        predicted = 2.0 * gamma * mean_pc
        bin_str = "  ".join([f"{b:>5.3f}" for b in bins_n])
        print(f"  {k:>3}  {eigvals[k].real:>+10.4f}  {eigvals[k].imag:>+10.4f}  {bin_str}  "
              f"{predicted:>8.4f}  (actual |Re(λ)| = {abs(eigvals[k].real):.4f})")
    print()
    print("Note: the 'predicted Z-decay-mean' (2γ · mean popcount) equals the actual")
    print("|Re(λ)| exactly for every mode, Hamiltonian or not. This is the Absorption")
    print("Theorem read state-side: the Z-deph dissipator is diagonal in the coherence")
    print("basis with eigenvalue -2γ·popcount(i XOR j), and the Hamiltonian part (incl.")
    print("the y-Zeeman) is anti-Hermitian, so it contributes no real part to λ. The")
    print("y-Zeeman reshapes WHICH popcount profile each slow mode carries, but it")
    print("cannot break the mean-popcount = rate identity.")
    print()

    # ---- (C) Joint reading: XOR (A) + depth (B) together ----
    print("(C) Joint reading: parity bin (XOR) + depth shape per mode")
    print("-" * 100)
    print(f"This is what 'both compute' looks like as a single fingerprint per mode:")
    print(f"  parity = n_Y XOR (1 bit) → which tower (Y or non-Y)")
    print(f"  depth profile = the n_Y histogram inside the tower (the actual distribution)")
    print()
    for k in range(n_slowest):
        v = eigvecs[:, k]
        rho = unvec_column(v, d)
        decomp = pauli_string_decomp(rho, N)
        bins = n_Y_histogram(decomp, N)
        total = sum(bins)
        if total < 1e-12:
            continue
        bins_n = [b / total for b in bins]
        even = sum(bins_n[i] for i in range(N + 1) if i % 2 == 0)
        odd = sum(bins_n[i] for i in range(N + 1) if i % 2 == 1)
        if odd < 1e-9:
            parity = "EVEN"
            active = [(i, bins_n[i]) for i in range(N + 1) if i % 2 == 0 and bins_n[i] > 1e-9]
        elif even < 1e-9:
            parity = "ODD"
            active = [(i, bins_n[i]) for i in range(N + 1) if i % 2 == 1 and bins_n[i] > 1e-9]
        else:
            parity = "MIXED"
            active = [(i, bins_n[i]) for i in range(N + 1) if bins_n[i] > 1e-9]
        depth_str = ", ".join([f"n_Y={i}:{w:.3f}" for i, w in active])
        print(f"  Mode {k:>2} (Re(λ)={eigvals[k].real:+.4f}): parity={parity:>5}, depth = {{{depth_str}}}")
    print()


def main():
    print("=" * 100)
    print("Computing BOTH: XOR (parity, 1 bit) AND depth (full n_Y distribution)")
    print("=" * 100)
    print()
    print("XOR is the parity = n_Y mod 2 = which tower a slow mode lives in.")
    print("Depth is the full count = how the mode distributes across the n_Y values")
    print("inside its parity sector. Together they reconstruct the Pauli-string profile.")
    run(N=4, h_zeeman=0.5, gamma=1.0, n_slowest=8)


if __name__ == "__main__":
    main()
