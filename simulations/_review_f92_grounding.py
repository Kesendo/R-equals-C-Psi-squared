"""Ad-hoc grounding check for the F92 review (throwaway, _-prefixed per CLAUDE.md).

Confirms FROM BELOW, not trusting the doc's own witness:
  (A) parameter-space facts about the table's "90 deg" map i: J_b -> 2*J_avg - J_{N-2-b}
        - is i an involution (i^2 = e) or order 4?
        - does i preserve the pair-DIFFERENCE and reflect the pair-SUM, or the reverse (line 39)?
  (B) load-bearing operator fact: R H(J_anti) R = -H(J_anti)  (R = spatial reflection / bit reversal)
  (C) F92's actual claim: the F71-refined diagonal-BLOCK spectrum is identical for uniform vs
        anti-palindromic J, while the FULL Liouvillian spectrum DIFFERS.
"""
import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def site_op(P, b, N):
    m = np.array([[1]], dtype=complex)
    for k in range(N):
        m = np.kron(m, P if k == b else I2)
    return m


def H_xy(J, N):
    """H = (1/2) sum_b J_b (X_b X_{b+1} + Y_b Y_{b+1})."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        XX = site_op(X, b, N) @ site_op(X, b + 1, N)
        YY = site_op(Y, b, N) @ site_op(Y, b + 1, N)
        H += 0.5 * J[b] * (XX + YY)
    return H


def bitrev(i, N):
    return int(format(i, f"0{N}b")[::-1], 2)


def reflection(N):
    d = 2 ** N
    R = np.zeros((d, d), dtype=complex)
    for i in range(d):
        R[bitrev(i, N), i] = 1.0
    return R


def liouvillian(J, N, gamma=0.3):
    d = 2 ** N
    H = H_xy(J, N)
    L = -1j * (np.kron(np.eye(d), H) - np.kron(H.T, np.eye(d)))
    for l in range(N):
        Zl = site_op(Z, l, N)
        L += gamma * (np.kron(Zl, Zl) - np.eye(d * d))  # Z real symmetric -> Z^T = Z
    return L


def revJ(J):
    return J[::-1].copy()


def sort_spec(vals, tol=9):
    v = np.round(vals, tol)
    return v[np.lexsort((v.imag, v.real))]


def block_spectra(L, Rsuper):
    w, V = np.linalg.eigh(Rsuper)  # eigenvalues +-1, real symmetric
    Up = V[:, w > 0]
    Um = V[:, w < 0]
    bp = Up.conj().T @ L @ Up
    bm = Um.conj().T @ L @ Um
    return np.linalg.eigvals(bp), np.linalg.eigvals(bm)


# ---------- (A) parameter-space facts ----------
print("=" * 70)
print("(A) The table's '90 deg' map  i: J_b -> 2*J_avg - J_{N-2-b}")
for N, J in [(4, np.array([2.0, 0.0, 0.0])), (5, np.array([1.7, 0.2, 0.9, -0.3]))]:
    Javg = J.mean()
    i1 = 2 * Javg - revJ(J)          # apply i once
    i2 = 2 * i1.mean() - revJ(i1)    # apply i twice (avg is i-invariant, but recompute honestly)
    T = J + revJ(J)                  # pair-sums
    B = J - revJ(J)                  # pair-differences
    Ti = i1 + revJ(i1)
    Bi = i1 - revJ(i1)
    print(f"  N={N}: J={J}  J_avg={Javg:.4f}")
    print(f"     i(J)      = {np.round(i1,4)}")
    print(f"     i^2(J)    = {np.round(i2,4)}   ->  i^2 == identity? {np.allclose(i2, J)}")
    print(f"     pair-SUM  T={np.round(T,4)} -> i -> {np.round(Ti,4)}   preserved? {np.allclose(Ti,T)}"
          f"   == 4*Javg - T? {np.allclose(Ti, 4*Javg - T)}")
    print(f"     pair-DIFF B={np.round(B,4)} -> i -> {np.round(Bi,4)}   preserved? {np.allclose(Bi,B)}")

# ---------- (B) operator anti-conjugation + (C) the actual F92 claim ----------
for N in (4, 5):
    print("=" * 70)
    print(f"(B)/(C)  N={N}")
    if N == 4:
        Jap = np.array([1.5, 1.0, 0.5])           # anti-palindromic: T=[2,2], J_avg=1
    else:
        Jap = np.array([1.3, 1.4, 0.6, 0.7])       # T0=2,T1=2, J_avg=1
    Jun = np.ones(N - 1)                            # uniform reference, same J_avg=1
    Janti = 0.5 * (Jap - revJ(Jap))
    Jsym = 0.5 * (Jap + revJ(Jap))
    print(f"   J_antipalin={Jap}  J_sym={np.round(Jsym,4)}  J_anti={np.round(Janti,4)}")

    R = reflection(N)
    Ha = H_xy(Janti, N)
    lhs = R @ Ha @ R
    print(f"   (B) R H(J_anti) R == -H(J_anti)?  {np.allclose(lhs, -Ha)}")
    Hs = H_xy(Jsym, N)
    print(f"       R H(J_sym)  R == +H(J_sym)?   {np.allclose(R @ Hs @ R, Hs)}")

    Lun = liouvillian(Jun, N)
    Lap = liouvillian(Jap, N)
    Rsuper = np.kron(R, R)
    print(f"       (RxR)^2 == I? {np.allclose(Rsuper @ Rsuper, np.eye((2**N)**2))}")

    bp_u, bm_u = block_spectra(Lun, Rsuper)
    bp_a, bm_a = block_spectra(Lap, Rsuper)
    match_plus = np.allclose(sort_spec(bp_u), sort_spec(bp_a))
    match_minus = np.allclose(sort_spec(bm_u), sort_spec(bm_a))
    print(f"   (C) diagonal-block(+) spectrum  uniform == anti-palindromic? {match_plus}")
    print(f"       diagonal-block(-) spectrum  uniform == anti-palindromic? {match_minus}")

    full_u = sort_spec(np.linalg.eigvals(Lun))
    full_a = sort_spec(np.linalg.eigvals(Lap))
    full_match = np.allclose(full_u, full_a)
    maxdiff = np.max(np.abs(full_u - full_a))
    print(f"       FULL-L spectrum uniform == anti-palindromic? {full_match}   max|diff|={maxdiff:.4e}")

# ---------- (D) multiset-vs-assignment: does F91 Step 6 'multiset' claim hold? ----------
# Both PALINDROMIC (so F71 is an exact symmetry: diagonal block = full structure, J_anti=0).
# Same pair-sum MULTISET, different ASSIGNMENT of which pair carries which sum.
print("=" * 70)
print("(D) F91 Step-6 'depends only on the MULTISET of pair-sums' -- gate check (N=6)")
N = 6
JA = np.array([5.0, 2.0, 3.0, 2.0, 5.0])  # palindromic; pair-sums {10,4,6}
JB = np.array([2.0, 5.0, 3.0, 5.0, 2.0])  # palindromic; pair-sums {4,10,6}  (same multiset)
print(f"   J_A={JA}  palindromic? {np.allclose(JA, revJ(JA))}  pair-sums={JA[:2]+JA[-2:][::-1]}, mid2={2*JA[2]}")
print(f"   J_B={JB}  palindromic? {np.allclose(JB, revJ(JB))}  pair-sums={JB[:2]+JB[-2:][::-1]}, mid2={2*JB[2]}")
R6 = reflection(N)
Rs6 = np.kron(R6, R6)
LA = liouvillian(JA, N)
LB = liouvillian(JB, N)
bpA, bmA = block_spectra(LA, Rs6)
bpB, bmB = block_spectra(LB, Rs6)
print(f"   diagonal-block(+) spectrum  J_A == J_B?  {np.allclose(sort_spec(bpA), sort_spec(bpB))}")
print(f"   diagonal-block(-) spectrum  J_A == J_B?  {np.allclose(sort_spec(bmA), sort_spec(bmB))}")
fullA, fullB = sort_spec(np.linalg.eigvals(LA)), sort_spec(np.linalg.eigvals(LB))
print(f"   FULL-L spectrum             J_A == J_B?  {np.allclose(fullA, fullB)}   max|diff|={np.max(np.abs(fullA-fullB)):.4e}")
print("   -> if NOT equal, the spectrum depends on the ASSIGNMENT l->S_l, not just the multiset.")
