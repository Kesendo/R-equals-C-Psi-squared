"""Gate-first PROOF of RING gap-dominance (the dihedral-lock ChiralK sibling of chain gap-dominance).

The chain result (PROVEN, chain_gap_dominance.py): max|Im| over the exact-(-2g) Liouvillian modes of the
OPEN XY chain = E1 = 2J cos(pi/(N+1)), the single-particle band edge. The RING is different: it is
translation-invariant (cyclic C_N), so its single-particle band TOP is the k=0 uniform mode at energy
2J = J*rho (rho = 2 = the ring adjacency spectral radius, matching TopologyBandEdgeClaim's "ring band edge").

RESULT (gate-verified below):
  * For N != 4:  max|Im| over the exact-(-2g) modes = 2J  (the periodic band top, the k=0 uniform free-
    fermion mode c_0^(dag).f(N_tot); reached exactly, gamma-independent). No symmetry enters: the (0,1)
    sector reaches J*rho on ANY connected graph, and "dihedral lock" is a historical label, not a mechanism
    (at even N the max |E_k| is attained TWICE, m=0 and m=N/2, so C_N singles out nothing there).
  * N = 4 is the SPECIAL CASE (the ring analogue of the chain's N=2, which also exceeds its
    band top for Q > 2/sqrt(3); the chain's N=3 case faces the other way and stays below E1):
    the half-filling (2,2) {0,2}-coherence
    sqrt-EP mode reaches Im = sqrt((2sqrt2 J)^2 - (2g)^2) -> 2sqrt2 J > 2J, EXCEEDING the band top. This is
    the "co-occupied floor mismatch for ring N=4" that TopologyBandEdgeClaim records; 2sqrt2 = sum of the two
    largest anti-periodic single-fermion energies (the even-parity / anti-periodic JW sector wraps the ring).
    It is the SAME half-filling (2,2) sector that makes K_4 special (StructuralCeiling 2-2/sqrt3).

  STAGE 0  the theorem (full Liouvillian): max|Im| at Re=-2g = 2J for N=3,5,6; = sqrt((2sqrt2)^2-(2g)^2) at N=4.
  STAGE 1  the band-top reach (general N): the k=0 uniform mode c_0^(dag).f(N_tot) is an exact -2g-/+i*2J
           L-eigenmode -- so 2J = J*rho is always achieved.
  STAGE 2  the N=4 exception: the (2,2) half-filling sqrt-EP at 2sqrt2 J (gamma-swept); N=6 (3,3) does NOT
           exceed 2J (the half-filling block closes on -2g only at N=4).
  STAGE 3  the exception's Q window: the floor count 26 -> 18 -> 16 across the two EPs (gamma = J and
           gamma = sqrt2 J) and the 16 + 2 + 8 split of the N=4 floor.
  STAGE 4  (--slow) the sector walk past the wall: the whole floor at N=7 and N=8 by (a,b) blocks.

Run: python simulations/ring_gap_dominance.py [--slow]
"""
import sys

import numpy as np
from math import cos, pi, sqrt

GAMMA = 0.05
J = 1.0
TOL = 1e-7
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_list(ops):
    out = np.array([[1.0 + 0j]])
    for o in ops:
        out = np.kron(out, o)
    return out


def site(op, l, N):
    return kron_list([op if k == l else I2 for k in range(N)])


def H_ring(N):
    d = 2 ** N
    H = np.zeros((d, d), complex)
    for i in range(N):
        j = (i + 1) % N
        H += (J / 2) * (site(X, i, N) @ site(X, j, N) + site(Y, i, N) @ site(Y, j, N))
    return H


def liouvillian(N, g=GAMMA):
    H = H_ring(N)
    d = 2 ** N
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in range(N):
        Zl = site(Z, l, N)
        L += g * (np.kron(Zl, Zl.conj()) - np.kron(Id, Id))
    return L


def vec(A):
    return A.flatten()


def ring_maximprime(N, g=GAMMA):
    ev = np.linalg.eigvals(liouvillian(N, g))
    at = np.abs(ev.real - (-2 * g)) < TOL
    return float(np.abs(ev.imag[at]).max()) if at.any() else float("nan")


# ====================================================================================================
# STAGE 0 -- THE THEOREM (full Liouvillian): max|Im| at Re=-2g
# ====================================================================================================
print("=" * 100)
print("STAGE 0 -- RING max|Im| at Re=-2g  vs  the band top 2J (=J*rho) and the N=4 (2,2) sqrt-EP")
print("=" * 100)
rho = 2.0                                                          # ring adjacency spectral radius (period., k=0)
band_top = J * rho                                                 # = 2J
print(f"{'N':>2} {'ring max|Im|':>13} {'prediction':>13} {'which':>28} {'|diff|':>9}")
for N in (3, 4, 5, 6):
    mx = ring_maximprime(N)
    if N == 4:
        pred = sqrt((2 * sqrt(2) * J) ** 2 - (2 * GAMMA) ** 2)     # the (2,2) sqrt-EP, -> 2sqrt2 J
        which = "(2,2) sqrt-EP -> 2sqrt2 J"
    else:
        pred = band_top                                            # 2J, the periodic band top, exact
        which = "2J = J*rho (band top)"
    d = abs(mx - pred)
    print(f"{N:>2} {mx:>13.6f} {pred:>13.6f} {which:>28} {d:>9.1e}")
    assert d < 1e-5, f"STAGE 0 GATE FIRED at N={N}: max|Im| {mx} != prediction {pred}"
print("STAGE 0 PASS: ring max|Im| = 2J (=J*rho, periodic band top) for N=3,5,6; N=4 = sqrt((2sqrt2)^2-(2g)^2)")
print("  -> 2sqrt2 J > 2J (the (2,2) half-filling sqrt-EP, the ring analogue of the chain's N=2")
print("     special, which also exceeds its band top; the chain N=3 one stays below E1).")

# ====================================================================================================
# STAGE 1 -- THE BAND-TOP REACH (general N): the (0,1) vacuum<->single-excitation sector reaches J*rho = 2J.
#   The naive JW c_0 is NOT a clean ring eigenoperator -- the wrap bond carries a parity factor (-1)^N_tot
#   (periodic for odd parity, anti-periodic for even), so c_0^(dag) flips parity and mixes the two BCs. But
#   the (0,1) sector is SINGLE-PARTICLE: |vac><psi_k| with psi_k an eigenvector of the ring adjacency A
#   (energy E_k = J*A-eigenvalue) is n_XY=1 (Re=-2g) and L_H gives [H,|vac><psi_k|] = -E_k|vac><psi_k|, so it
#   is an EXACT -2g + iE_k mode. max E_k = J*rho = 2J (the k=0 uniform psi_0), no JW parity needed.
# ====================================================================================================
print("\n" + "=" * 100)
print("STAGE 1 -- the dihedral lock: the (0,1) sector |vac><psi_k| are exact -2g + iE_k modes, max = J*rho = 2J")
print("=" * 100)


def ring_adjacency(N):
    A = np.zeros((N, N))
    for i in range(N):
        A[i, (i + 1) % N] = A[(i + 1) % N, i] = 1.0
    return A


for N in (3, 4, 5, 6):
    d = 2 ** N
    L = liouvillian(N)
    Ek, phi = np.linalg.eigh(J * ring_adjacency(N))                # single-excitation band = J * ring adjacency
    vac = np.zeros(d, complex); vac[0] = 1.0
    exc = [np.zeros(d, complex) for _ in range(N)]
    for i in range(N):
        e = np.zeros(d, complex); e[1 << (N - 1 - i)] = 1.0; exc[i] = e   # excitation at site i (bit i)
    allhit = True
    for k in range(N):
        psi_k = sum(phi[i, k] * exc[i] for i in range(N))
        rho = np.outer(vac, psi_k.conj())                          # |vac><psi_k|
        v = vec(rho)
        lam = -2 * GAMMA + 1j * Ek[k]                              # L_H[|vac><psi_k|] = +iE_k (H|vac>=0)
        if np.linalg.norm(L @ v - lam * v) > 1e-8 * max(1.0, np.linalg.norm(v)):
            allhit = False
    assert allhit, f"STAGE 1 GATE FIRED at N={N}: (0,1) coherences are not exact -2g+iE_k modes"
    top = Ek.max()
    assert abs(top - band_top) < 1e-9, f"STAGE 1 GATE FIRED at N={N}: single-excitation band top {top} != 2J"
    print(f"  N={N}: |vac><psi_k| exact at -2g + iE_k for all k; single-excitation band top = J*rho = {top:.6f} (=2J)")
print("STAGE 1 PASS: the (0,1) sector reaches J*rho = 2J exactly at every N (the C_N dihedral lock = the k=0")
print("  uniform single-excitation psi_0); the ring band top is always achieved on the exact-(-2g) subspace.")

# ====================================================================================================
# STAGE 2 -- THE N=4 EXCEPTION: the (2,2) half-filling sqrt-EP at 2sqrt2 J; N=6 (3,3) does NOT exceed 2J
# ====================================================================================================
print("\n" + "=" * 100)
print("STAGE 2 -- N=4 (2,2) half-filling sqrt-EP at 2sqrt2 J (gamma-swept); N=6 half-filling does NOT exceed")
print("=" * 100)
ap_two_top = 2 * (2 * J * cos(pi / 4))                             # sum of two largest anti-periodic E_k (N=4) = 2sqrt2
print(f"  2sqrt2 J = anti-periodic two-fermion top = {ap_two_top:.6f}")
for g in (0.05, 0.10, 0.20):
    mx = ring_maximprime(4, g)
    pred = sqrt(ap_two_top ** 2 - (2 * g) ** 2)
    print(f"  N=4, g={g:.2f}: max|Im| = {mx:.6f}   sqrt((2sqrt2)^2-(2g)^2) = {pred:.6f}   exceeds 2J = {mx > band_top + 1e-9}")
    assert abs(mx - pred) < 1e-5, f"STAGE 2 GATE FIRED: N=4 max|Im| {mx} != sqrt-EP {pred}"
    assert mx > band_top + 1e-9, "STAGE 2: N=4 should EXCEED the band top 2J"
# N=6: the half-filling (3,3) does NOT reach -2g with excess -> max|Im| stays 2J
assert abs(ring_maximprime(6) - band_top) < 1e-5, "STAGE 2: N=6 should NOT exceed 2J (max|Im|=2J)"
print("  N=6: max|Im| = 2J (the (3,3) half-filling does NOT close on -2g with excess) -- N=4 is the LONE")
print("       exception, the unique even half-filling where the (2,2) {0,2} block lands on -2g (cf K_4, ring-4).")
print("STAGE 2 PASS: N=4 is the sole exception (2sqrt2 J sqrt-EP); the same half-filling sector as K_4's ceiling.")

# ====================================================================================================
# STAGE 3 -- THE EXCEPTION'S Q WINDOW: a {0,2} block rides the floor only while its coupling beats the
#   watching. Its cells are 0 (population) and -4g (2-Hamming coherence), so the eigenvalues are
#   -2g +/- sqrt((2g)^2 - B^2): frequency sqrt(B^2-(2g)^2) above the EP at B = 2g, real split (OFF the
#   floor) below it. At N=4 two B's are in play -- B = 2J for the eight (1,1)/(3,3) modes, B = 2sqrt2 J
#   for the two (2,2) modes -- so the floor count falls 26 -> 18 -> 16 as g passes J and sqrt2*J, and
#   past the second crossing max|Im| is back at the band top 2J: the N=4 exception is a Q > Q* = 1/sqrt2
#   statement, not an all-Q one. (The 26 = 16 + 10 split is STAGE 3b.)
# ====================================================================================================
print("\n" + "=" * 100)
print("STAGE 3 -- the N=4 exception's Q window: the {0,2} blocks leave the floor past their own EP")
print("=" * 100)


def floor_dim_and_max(N, g):
    ev = np.linalg.eigvals(liouvillian(N, g))
    at = np.abs(ev.real - (-2 * g)) < TOL
    return int(at.sum()), float(np.abs(ev.imag[at]).max())


B_11, B_22 = 2 * J, 2 * sqrt(2) * J                                # the two {0,2} couplings at N=4
g_exceed = sqrt(B_22 ** 2 - band_top ** 2) / 2                     # sqrt(8-4)/2 = J: where the (2,2) freq = 2J
print(f"  the two EPs at N=4 (2g = B): g = {B_11/2:.6f} ((1,1)/(3,3)) and g = {B_22/2:.6f} ((2,2))")
print(f"  the EXCEEDING window ends earlier, at g = sqrt(B_22^2 - (J*rho)^2)/2 = {g_exceed:.6f} = J (Q* = 1),")
print("  where the (2,2) floor frequency has dropped back to the band top 2J while its block is still on the floor")
print(f"{'g':>8} {'floor dim':>10} {'max|Im|':>10} {'expected dim':>13} {'expected max':>13} {'exceeds 2J':>11}")
for g, exp_dim in ((0.05, 26), (0.50, 26), (1.20, 18), (1.35, 18), (1.50, 16), (1.60, 16)):
    dim, mx = floor_dim_and_max(4, g)
    freq_22 = sqrt(B_22 ** 2 - (2 * g) ** 2) if 2 * g < B_22 else 0.0     # 0 = gone from the floor
    exp_max = max(band_top, freq_22)
    print(f"{g:>8.2f} {dim:>10} {mx:>10.6f} {exp_dim:>13} {exp_max:>13.6f} {str(mx > band_top + 1e-9):>11}")
    assert dim == exp_dim, f"STAGE 3 GATE FIRED at g={g}: floor dim {dim} != {exp_dim}"
    assert abs(mx - exp_max) < 1e-5, f"STAGE 3 GATE FIRED at g={g}: max|Im| {mx} != {exp_max}"
    assert (mx > band_top + 1e-9) == (g < g_exceed - 1e-9), \
        f"STAGE 3 GATE FIRED at g={g}: exceeding the band top should hold exactly for g < J"
print("  -> the exception (max|Im| > 2J) is the g < J side, i.e. Q > 1; between J and sqrt2*J the (2,2) block")
print("     still rides the floor but below the band top; past sqrt2*J it is gone and only V_1 is left.")

# STAGE 3b -- the ten extras at N=4 are TWO (2,2) modes plus EIGHT (1,1)/(3,3) modes, all n_XY-mixed
w_ham = np.array([[bin(i ^ j).count("1") for j in range(16)] for i in range(16)]).reshape(-1)
ev4, vr4 = np.linalg.eig(liouvillian(4))
at4 = np.where(np.abs(ev4.real - (-2 * GAMMA)) < TOL)[0]
pure, ext_22, ext_11 = 0, 0, 0
for c in at4:
    v = np.abs(vr4[:, c]) ** 2
    v = v / v.sum()
    if v[w_ham != 1].sum() < 1e-9:
        pure += 1
    elif abs(abs(ev4[c].imag) - sqrt(B_22 ** 2 - (2 * GAMMA) ** 2)) < 1e-5:
        ext_22 += 1
    else:
        ext_11 += 1
        assert abs(abs(ev4[c].imag) - sqrt(B_11 ** 2 - (2 * GAMMA) ** 2)) < 1e-5, \
            "STAGE 3b GATE FIRED: an extra mode sits at neither {0,2} frequency"
print(f"  the 26 floor modes at N=4 split {pure} pure n_XY=1 (= V_1) + {ext_22} at sqrt(B_22^2-(2g)^2) "
      f"((2,2)) + {ext_11} at sqrt(B_11^2-(2g)^2) ((1,1)/(3,3))")
assert (pure, ext_22, ext_11) == (16, 2, 8), "STAGE 3b GATE FIRED: the 16 + 2 + 8 split does not hold"
print("STAGE 3 PASS: the exception holds on the Q > 1 side only, and the ten extras are 2 + 8, not ten (2,2)")
print("  modes: only the (2,2) pair exceeds the band top, the (1,1)/(3,3) eight sit below it.")

# ====================================================================================================
# STAGE 4 (opt-in: --slow) -- PAST THE 4^N WALL BY SECTORS. L is block-diagonal in the sector label
#   (a, b) = (excitation number of the ket, of the bra): H conserves excitation number and the Z-dephasing
#   is diagonal, so no matrix element leaves its (a, b) block. The largest block at N=8 is C(8,4)^2 = 4900
#   instead of 4^8 = 65536, which puts the FULL spectrum at N=7 and N=8 within reach: the completeness
#   section carries the n_XY=1 subspace to N=7 only, and nothing to N=8. If the exception really is the
#   smallest even half-filling, N=8 must read exactly 2J over ALL sectors, not just over V_1.
# ====================================================================================================
if "--slow" in sys.argv:
    from itertools import combinations

    def sector_floor(N, g=GAMMA):
        """max|Im| over the modes at Re = -2g, walking the (a,b) blocks of L one at a time."""
        states = {p: [sum(1 << i for i in c) for c in combinations(range(N), p)] for p in range(N + 1)}
        bonds = [(i, (i + 1) % N) for i in range(N)]

        def H_block(p):                                    # the p-excitation block of H, hopping J
            s = states[p]; idx = {v: k for k, v in enumerate(s)}
            M = np.zeros((len(s), len(s)), complex)
            for v in s:
                for (i, j) in bonds:
                    bi, bj = (v >> i) & 1, (v >> j) & 1
                    if bi != bj:
                        M[idx[v ^ (1 << i) ^ (1 << j)], idx[v]] += J
            return M

        Hb = {p: H_block(p) for p in range(N + 1)}
        best, count = 0.0, 0
        for a in range(N + 1):
            for b in range(N + 1):
                sa, sb = states[a], states[b]
                na, nb = len(sa), len(sb)
                L = np.kron(Hb[a], np.eye(nb)) - np.kron(np.eye(na), Hb[b].T)
                L = -1j * L
                deph = np.array([-2.0 * g * bin(sa[i] ^ sb[j]).count("1") for i in range(na) for j in range(nb)])
                ev = np.linalg.eigvals(L + np.diag(deph))
                at = np.abs(ev.real - (-2 * g)) < TOL
                count += int(at.sum())
                if at.any():
                    best = max(best, float(np.abs(ev.imag[at]).max()))
        return count, best

    print("\n" + "=" * 100)
    print("STAGE 4 (--slow) -- the FULL spectrum past the wall, by (a,b) sectors: N=7 and N=8")
    print("=" * 100)
    for N in (4, 5, 6, 7, 8):                              # 4..6 first: the sector walk must reproduce STAGE 0
        cnt, mx = sector_floor(N)
        pred = sqrt((2 * sqrt(2) * J) ** 2 - (2 * GAMMA) ** 2) if N == 4 else band_top
        print(f"  N={N}: floor dim {cnt:>4}  max|Im| = {mx:.6f}  (prediction {pred:.6f})")
        assert abs(mx - pred) < 1e-5, f"STAGE 4 GATE FIRED at N={N}: sector max|Im| {mx} != {pred}"
    print("STAGE 4 PASS: the sector walk reproduces the full-L answer at N=4,5,6 and then carries it past the")
    print("  wall: N=7 and N=8 read exactly 2J over EVERY (a,b) sector, so the exception really is N=4 alone.")

print("\n" + "=" * 100)
print("GATE-VERIFIED (N=3..6 full L, N=7..8 by sector walk; the all-N step stays open, see the proof):")
print("  max|Im| over the exact-(-2g) modes = 2J = J*rho (the periodic band top = ring adjacency radius,")
print("  reached by the (0,1) sector -- which reaches J*rho on any connected graph, so this half is not the")
print("  ring's), for all N EXCEPT N=4, where the half-filling (2,2) {0,2}")
print("  sqrt-EP reaches 2sqrt2 J > 2J -- on the Q > 1 side of that N, the exception's own window (STAGE 3).")
print("  Contrast the chain (open, no wrap): E1 = 2J cos(pi/(N+1)) < 2J, with")
print("  its OWN lone exception at N=3. Scope: ring (cyclic C_N); JW wrap = the periodic/anti-periodic split.")
print("DONE.")
