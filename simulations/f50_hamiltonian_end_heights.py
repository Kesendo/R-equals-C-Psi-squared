"""The Hamiltonian end of block (p, p): the limiting heights and m_p, Heisenberg chain to N = 12.

Theorem C of docs/proofs/PROOF_WEIGHT1_DEGENERACY.md ("The count as a plane crossing") bounds the exceptional
count from below by C(N, p) - 1 - m_p, with m_p the number of non-stationary modes of block (p, p) whose height
k = <Hamming> = -Re lambda/(2 gamma) lies below 1 for every small gamma. Every mode with Im lambda -> E_a - E_b != 0
has height >= 1 by Theorem A, so only ker ad_H inside the block can hold such a mode. When the p-excitation
spectrum is simple that kernel is the span of the eigenstate populations |E_a><E_a|, and their limiting heights
are the eigenvalues of the Hamming operator compressed onto it (first-order perturbation of the normal operator
-i ad_H by gamma D, D = -2 Hamming on the cells). Since Hamming(x, y) = 2 (p - |x and y|) on two popcount-p
configurations, the compression is

    K = 2p I - 2 W_p,   W_p = sum_l n_l o n_l,   (n_l)_ab = <E_a| n_l |E_b>,   o = Hadamard product,

so the limiting heights are 2 (p - mu) over the eigenvalues mu of W_p. W_p is nonnegative with every row sum p
(sum_b (n_l)_ab^2 = (n_l)_aa and sum_l (n_l)_aa = p), so mu_1 = p is the stationary mode at height 0, and
m_p = 0 if the second-largest eigenvalue mu_2 lies below p - 1/2 (at mu_2 = p - 1/2 a limiting height sits on the
plane and first order does not decide). At p = 1, W_1 is the overlap matrix W
of the proof's "The count at p = 1, every N".

Rows (Pauli book, H = sum XX + YY + ZZ, J = 1, gamma/J the unit; every N from 4 to 12, every p from 1 to N/2):
  H1  direct compression sum_{x,y} E_a(x)E_a(y) Ham(x,y) E_b(x)E_b(y) against the closed form, N <= 6. Two float
      routes to one matrix: the difference is rounding, error model eps * d * ||K||, ratio printed and gated < 8.
  H2  the compressed heights against the full block eigensolver at J/gamma = 200 and 400 on the modes with
      |Im lambda| < 1e-3 (the kernel modes; their count must equal dim ker), on every block of N = 4, 5. Error
      model (gamma/J)^2: the deviation ratio 200/400 must lie in [3.8, 4.2] on every mode whose deviation is above
      rounding (1e-10); the modes at rounding are listed, not gated.
  H3  the sector spectrum is simple in every block, so the kernel is the populations: the smallest level gap is printed
      and read against the eigensolver's resolution (eps * ||H||, about 1e-15 here; the smallest gap found is 7.8e-6 at
      N = 12), a reading with its error model rather than a threshold that merely passes.
  H4  m_p read: #{heights < 1} - 1, and mu_2 against p - 1/2, every block N = 4..12 (--n16 adds N = 13..16, about
      fifteen minutes, the half-filling block at N = 16 having dimension 12870), with sigma_2(G)^2 beside it, G_xa = E_a(x)^2
      the doubly stochastic matrix of the eigenstates' weights: the diagonal fraction of a traceless commutant element is at
      most sigma_2(G)^2, and since every off-diagonal cell has Hamming >= 2 the height is at least 2 (1 - that fraction), so
      sigma_2(G)^2 < 1/2 is a sufficient condition one step stronger than mu_2 < p - 1/2. A reading, not a gate: it is
      what the proof cites as measured.
  H5  control: the ZZ term dropped (XY chain) at p = 1 must give heights 2N/(N+1) and 2 (the proof's sine-mode
      formula), and the N = 4 ring at p = 1 must show a height exactly 1 at the Hamiltonian end (the K = pi mode
      of gate R7 of f50_plane_crossing_count.py that never crosses).
  H6  the identity with F87's first-order omega = 0 block (experiments/F87_WINDOWED_CONVERSE_PER_BLOCK.md): on the
      populations M_0 = sum_l (Z_l o Z_l) - N I equals 4 W_p - 4p I (Z_l = 1 - 2 n_l), two routes to one matrix.
  H7  the uniform XY chain, every block N = 4..8, on the FULL kernel (its free-fermion levels are degenerate): m_p = 0
      and the smallest non-stationary height equals the l = 1 value 2N/(N + 1) (PROOF_FROZEN_BAND_SO4 Theorem 6.2,
      inherited to every rung by Corollary 4.2) and lies above the l >= 2 floor 2 l (N - l)/(N + 1) of Corollary 7.3.
  H8  negative control: on the complete graph K4 and the ring at N = 4, block (2,2), the first-order reader must register
      2 and 1 limiting heights strictly below the plane (0.845 on the ring), so a reader that cannot see a mode below the
      plane fails. K4's two are its whole Krein debt m_2 = 2; the ring's debt is m_2 = 3, two more modes approaching the
      plane from inside with limiting height exactly 1 (f50_inertia_identity.py I5), which first order cannot see.
Prints ALL GATES PASS when H1, H2, H3, H5, H6, H7 and H8 hold. About two minutes.
"""
import sys
import numpy as np
from math import comb

EPS = np.finfo(float).eps


def sector(N, p):
    confs = [c for c in range(2 ** N) if bin(c).count("1") == p]
    return confs, {c: i for i, c in enumerate(confs)}


def h_sector(N, p, bonds, delta=1.0):
    confs, idx = sector(N, p)
    d = len(confs)
    H = np.zeros((d, d))
    for i, c in enumerate(confs):
        for (l, m) in bonds:
            b1, b2 = (c >> l) & 1, (c >> m) & 1
            H[i, i] += delta * (1.0 if b1 == b2 else -1.0)
            if b1 != b2:
                H[idx[c ^ (1 << l) ^ (1 << m)], i] += 2.0
    return confs, H


def chain(N):
    return [(l, l + 1) for l in range(N - 1)]


def ring(N):
    return chain(N) + [(N - 1, 0)]


def hamming_matrix(confs):
    a = np.array(confs)
    x = a[:, None] ^ a[None, :]
    return np.vectorize(lambda v: bin(v).count("1"))(x).astype(float)


def min_gap(E):
    return float(np.min(np.diff(E)))


def k_direct(V, confs):
    Ham = hamming_matrix(confs)
    d = len(confs)
    K = np.zeros((d, d))
    for a in range(d):
        Ma = np.outer(V[:, a], V[:, a]) * Ham
        for b in range(d):
            K[a, b] = np.sum(Ma * np.outer(V[:, b], V[:, b]))
    return K


def w_closed(N, V, confs):
    d = len(confs)
    W = np.zeros((d, d))
    for l in range(N):
        bit = np.array([(c >> l) & 1 for c in confs], dtype=float)
        n_l = V.T @ (bit[:, None] * V)
        W += n_l * n_l
    return W


def heights(N, p, bonds, delta=1.0):
    confs, H = h_sector(N, p, bonds, delta)
    E, V = np.linalg.eigh(H)
    W = w_closed(N, V, confs)
    mu = np.sort(np.linalg.eigvalsh(W))[::-1]
    return E, V, confs, mu, np.sort(2 * (p - mu))


def full_block_kernel_heights(N, p, bonds, gamma_over_J, im_tol=1e-3, delta=1.0):
    confs, H = h_sector(N, p, bonds, delta)
    d = len(confs)
    Ham = hamming_matrix(confs)
    I = np.eye(d)
    B = -1j * (np.kron(H, I) - np.kron(I, H.T)) + gamma_over_J * np.diag((-2 * Ham).ravel())
    lam = np.linalg.eigvals(B)
    lam = lam[np.abs(lam.imag) < im_tol]
    return np.sort(-lam.real / (2 * gamma_over_J))


ok_all = True


def gate(cond, msg):
    global ok_all
    ok_all &= bool(cond)
    print(("  PASS  " if cond else "  FAIL  ") + msg)


print("H1  direct vs closed-form compression (two float routes; ratio of max |diff| to eps*d*||K||, gated < 8)")
for N in (4, 5, 6):
    for p in range(1, N // 2 + 1):
        confs, H = h_sector(N, p, chain(N))
        E, V = np.linalg.eigh(H)
        K1 = k_direct(V, confs)
        K2 = 2 * p * np.eye(len(confs)) - 2 * w_closed(N, V, confs)
        ratio = np.max(np.abs(K1 - K2)) / (EPS * len(confs) * np.linalg.norm(K2, 2))
        gate(ratio < 8, f"N={N} p={p}: max|diff| = {np.max(np.abs(K1 - K2)):.2e}, ratio to eps*d*||K|| = {ratio:.2f}")

print("H2  compressed heights vs full block at J/gamma = 200, 400 on the kernel modes (ratio of deviations, (gamma/J)^2 law)")
for N in (4, 5):
    for p in range(1, N // 2 + 1):
        E, V, confs, mu, k = heights(N, p, chain(N))
        h200 = full_block_kernel_heights(N, p, chain(N), 1 / 200.0)
        h400 = full_block_kernel_heights(N, p, chain(N), 1 / 400.0)
        gate(len(h200) == len(k) == len(h400), f"N={N} p={p}: kernel modes found {len(h200)} and {len(h400)}, dim ker = {len(k)}")
        dev200 = np.abs(h200 - k)
        dev400 = np.abs(h400 - k)
        live = dev200 > 1e-10
        ratios = dev200[live] / dev400[live]
        at_rounding = dev200[~live].max() if (~live).any() else 0.0
        gate(np.all((ratios > 3.8) & (ratios < 4.2)),
             f"N={N} p={p}: {live.sum()} modes above rounding, ratios {np.array2string(ratios, precision=3)}; "
             f"{(~live).sum()} at rounding (max dev {at_rounding:.1e})")

print("H3/H4  simple sector spectrum; m_p = #{heights < 1} - 1, mu_2 vs p - 1/2 and sigma_2(G)^2, chain N = 4..12 (--n16: ..16) (readings)")
Ns = [int(a) for a in sys.argv[1:] if a.isdigit()] or list(range(4, 17 if "--n16" in sys.argv else 13))
for N in Ns:
    for p in range(1, N // 2 + 1):
        E, V, confs, mu, k = heights(N, p, chain(N))
        gate(min_gap(E) > 1e-9, f"H3 N={N} p={p}: dim {comb(N, p):5d}, min level gap {min_gap(E):.3e}")
        below = int(np.sum(k < 1 - 1e-9))
        on = int(np.sum(np.abs(k - 1) <= 1e-9))
        G = V * V
        sigma2_sq = np.sort(np.linalg.eigvalsh(G.T @ G))[-2]
        print(f"        H4 N={N} p={p}: m_p = {below - 1}, on the plane {on}; mu_2 = {mu[1]:.6f} vs p - 1/2 = {p - 0.5}; "
              f"smallest nonzero height {k[1]:.6f}, largest {k[-1]:.6f} (= 2 min(p, N-p): {abs(k[-1] - 2 * min(p, N - p)) < 1e-9}); "
              f"sigma_2(G)^2 = {sigma2_sq:.6f}", flush=True)

print("H5  controls")
for N in (4, 5, 6, 8):
    E, V, confs, mu, k = heights(N, 1, chain(N), delta=0.0)
    m = (N - 1) // 2
    pred = np.sort(np.array([0.0] + [2 * N / (N + 1)] * m + [2.0] * (N - 1 - m)))
    gate(np.max(np.abs(k - pred)) < 1e-9,
         f"XY chain N={N} p=1: heights {np.array2string(k, precision=6)} = 0, 2N/(N+1) x{m}, 2 x{N - 1 - m}")


def kernel_heights_general(N, p, bonds, degtol=1e-9, delta=1.0):
    """The compression on the full ker ad_H of the block, |E_a><E_b| over exactly degenerate pairs (the degenerate case)."""
    confs, H = h_sector(N, p, bonds, delta)
    E, V = np.linalg.eigh(H)
    Ham = hamming_matrix(confs)
    groups, start = [], 0
    for i in range(1, len(E) + 1):
        if i == len(E) or E[i] - E[i - 1] > degtol:
            groups.append(list(range(start, i)))
            start = i
    basis = [(a, b) for g in groups for a in g for b in g]
    K = np.zeros((len(basis), len(basis)))
    for I, (a, b) in enumerate(basis):
        M = np.outer(V[:, a], V[:, b]) * Ham
        for J, (c, d) in enumerate(basis):
            K[I, J] = np.sum(M * np.outer(V[:, c], V[:, d]))
    return E, [len(g) for g in groups], np.sort(np.linalg.eigvalsh(K))


E, degs, k = kernel_heights_general(4, 1, ring(4))
print(f"        ring N=4 p=1 levels {np.round(E, 6)}, degeneracies {degs}, ker ad_H dim {len(k)} (the full kernel, not the populations alone)")
gate(np.sum(np.abs(k - 1) < 1e-9) >= 1, f"ring N=4 p=1 limiting heights {np.array2string(k, precision=6)}: a height exactly 1 present")
print("H6  F87's first-order omega = 0 block on the populations, M_0 = sum_l (Z_l)_aa' (Z_l)_bb' - N delta with Z_l = 1 - 2 n_l,")
print("    equals 4 W_p - 4p I (two routes to one matrix; ratio of max |diff| to eps*d*N gated < 8)")
for N in (4, 5, 6):
    for p in range(1, N // 2 + 1):
        confs, H = h_sector(N, p, chain(N))
        E, V = np.linalg.eigh(H)
        d = len(confs)
        M0 = -N * np.eye(d)
        for l in range(N):
            z = np.array([1.0 - 2.0 * ((c >> l) & 1) for c in confs])
            Z = V.T @ (z[:, None] * V)
            M0 += Z * Z
        W = w_closed(N, V, confs)
        diff = np.max(np.abs(M0 - (4 * W - 4 * p * np.eye(d))))
        gate(diff / (EPS * d * N) < 8, f"N={N} p={p}: max|M_0 - (4W_p - 4p I)| = {diff:.2e}, ratio {diff / (EPS * d * N):.2f}")

print("H7  XY chain, every block, the FULL kernel (degenerate free-fermion levels): m_p and the F144 floor")
print("    (PROOF_FROZEN_BAND_SO4 Cor 7.3: height = 2K >= 2 l (N - l)/(N + 1) on the rung-l multiplets, binding at l = 1)")
for N in (4, 5, 6, 7, 8):
    for p in range(1, N // 2 + 1):
        E, degs, k = kernel_heights_general(N, p, chain(N), delta=0.0)
        floor = 4 * (N - 2) / (N + 1)      # Corollary 7.3 at l = 2 in height units, the binding l >= 2 rung
        nonzero = k[k > 1e-9]
        below = int(np.sum(k < 1 - 1e-9))
        # the l = 1 heights are exactly 2N/(N+1) and 2; every other non-stationary height must sit at or above the l >= 2 floor
        others = nonzero[(np.abs(nonzero - 2 * N / (N + 1)) > 1e-9) & (np.abs(nonzero - 2.0) > 1e-9)]
        gate(below == 1 and abs(nonzero.min() - 2 * N / (N + 1)) < 1e-9 and (others.size == 0 or others.min() >= floor - 1e-9),
             f"N={N} p={p}: dim {comb(N, p):3d}, ker ad_H dim {len(k):4d} (degeneracies {sorted(set(degs))}), m_p = {below - 1}, "
             f"smallest nonzero height {nonzero.min():.6f} = 2N/(N+1) = {2 * N / (N + 1):.6f}; "
             f"{others.size} heights off the l = 1 values, smallest {others.min() if others.size else float('nan'):.6f} >= l >= 2 floor {floor:.6f}")

print("H8  negative control: the first-order reader must register the limiting heights strictly below the plane (N = 4, p = 2, full kernel)")
for label, bonds, expect in (("complete graph K4", [(a, b) for a in range(4) for b in range(a + 1, 4)], 2), ("ring", ring(4), 1)):
    E, degs, k = kernel_heights_general(4, 2, bonds)
    below = int(np.sum(k < 1 - 1e-9))
    gate(below - 1 == expect, f"{label} (2,2): limiting heights below 1 {np.array2string(np.sort(k)[:4], precision=6)}, strictly below the plane: {below - 1} (expected {expect})")

print("\nALL GATES PASS" if ok_all else "\nSOME GATE FAILED")
