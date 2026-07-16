# one_four_thesis.py -- the one-four thesis gate: one i, two Z4 gradings, tied at the square.
#
# Two candidate identities behind "every period-4 in the repo is the same quarter-turn of i":
#   (1) FACTORIZATION (holds, exact): V_g = Ad(U_g) o Pi^2 o K, with Pi^2 = Ad(X^N) (F1^2).
#       The F132 dead-set kill map contains Pi's SQUARE literally -- the Z2 center only.
#   (2) STRONG THESIS (refuted, every N): Pi's Z4 spectral grading == Majorana-degree-mod-4
#       grading. Multiplicities differ at EVERY N:
#         Pi sectors:      4^(N-1) each (uniform; orbit pairing, eigenvalues +-i^nYZ)
#         degree classes:  4^(N-1) + 2^(N-1) * Re(i^(N-r))   (roots-of-unity filter on C(2N,d))
#       The deviation +-2^(N-1) never vanishes, so no bijection exists at any N.
#
# Conventions pinned to sources:
#   Pi per-site: I<->X, Y->iZ, Z->iY (docs/proofs/MIRROR_SYMMETRY_PROOF.md, "Pi is order 4")
#   Majorana (left JW): a_{2l} = (prod_{m<l} Z_m) X_l, a_{2l+1} = (prod_{m<l} Z_m) Y_l
#   eps_odd = (-1)^(d(d-1)/2), eps_even = (-1)^(d(d+1)/2) (experiments/LATTICE_DEAD_SET_RULE.md)
#
# Checks: count check() calls before quoting a number in the doc.
import numpy as np
from itertools import product
from math import comb

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = [I2, X, Y, Z]
LET = "IXYZ"

_checks = [0, 0]
def check(name, ok, detail=""):
    _checks[0] += 1
    _checks[1] += ok
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    return ok

def kron_all(ms):
    out = np.array([[1.0 + 0j]])
    for m in ms:
        out = np.kron(out, m)
    return out

def string_matrix(s):
    return kron_all([PAULI[c] for c in s])

# ---------- combinatorial layer (any N) ----------

def majorana_mask(s):
    """Left-JW Majorana support of a Pauli string, as a 2N-bit mask (XOR of letter masks)."""
    mask = 0
    for l, c in enumerate(s):
        if c == 1:    # X_l = (prod_{m<l} Z_m) a_{2l}: bits 0..2l
            mask ^= (1 << (2 * l + 1)) - 1
        elif c == 2:  # Y_l = (prod_{m<l} Z_m) a_{2l+1}: bits 0..2l-1 and 2l+1
            mask ^= ((1 << (2 * l)) - 1) | (1 << (2 * l + 1))
        elif c == 3:  # Z_l = -i a_{2l} a_{2l+1}: bits 2l, 2l+1
            mask ^= 0b11 << (2 * l)
    return mask

def degree_classes(N):
    sizes = [0, 0, 0, 0]
    for s in product(range(4), repeat=N):
        sizes[bin(majorana_mask(s)).count("1") % 4] += 1
    return sizes

def degree_classes_closed_form(N):
    # sum_{d = r mod 4} C(2N, d) = 4^(N-1) + 2^(N-1) * Re(i^(N-r))
    re_i = [1, 0, -1, 0]  # Re(i^k) for k mod 4
    return [4 ** (N - 1) + 2 ** (N - 1) * re_i[(N - r) % 4] for r in range(4)]

def pi_multiplicities_combinatorial(N):
    # Pi pairs sigma with its I<->X, Y<->Z partner (never fixed); on each pair the
    # eigenvalues are +-i^nYZ. Even nYZ pairs feed {+1,-1}, odd feed {+i,-i}, evenly.
    E = sum(1 for s in product(range(4), repeat=N)
            if sum(1 for c in s if c >= 2) % 2 == 0)
    O = 4 ** N - E
    return {1: E // 2, -1: E // 2, 1j: O // 2, -1j: O // 2}

# ---------- matrix layer (N small) ----------

def build_pi(N, strings, idx):
    swap = {0: 1, 1: 0, 2: 3, 3: 2}
    P = np.zeros((4 ** N, 4 ** N), dtype=complex)
    for s in strings:
        nYZ = sum(1 for c in s if c >= 2)
        s2 = tuple(swap[c] for c in s)
        P[idx[s2], idx[s]] = 1j ** nYZ
    return P

def majoranas(N):
    out = []
    for l in range(N):
        tail = [Z] * l
        out.append(kron_all(tail + [X] + [I2] * (N - 1 - l)))
        out.append(kron_all(tail + [Y] + [I2] * (N - 1 - l)))
    return out

print("=== T1: Pi^4 = I and Pi^2 = Ad(X^N) on all strings (matrix, N=3,4) ===")
for N in (3, 4):
    strings = list(product(range(4), repeat=N))
    idx = {s: k for k, s in enumerate(strings)}
    P = build_pi(N, strings, idx)
    check(f"N={N}: Pi^4 = I", np.allclose(np.linalg.matrix_power(P, 4), np.eye(4 ** N)))
    P2 = P @ P
    ok = all(np.isclose(P2[idx[s], idx[s]], (-1) ** (sum(1 for c in s if c >= 2)))
             and np.isclose(np.abs(P2[:, idx[s]]).sum(), 1) for s in strings)
    check(f"N={N}: Pi^2 = (-1)^(nY+nZ) diagonal = Ad(X^N)", ok)

print("=== T2: Pi eigenvalue multiplicities uniform 4^(N-1) ===")
for N in (3, 4, 5):
    strings = list(product(range(4), repeat=N))
    idx = {s: k for k, s in enumerate(strings)}
    ev = np.linalg.eigvals(build_pi(N, strings, idx))
    mult = {k: int(np.isclose(ev, k).sum()) for k in (1, 1j, -1, -1j)}
    ok = all(v == 4 ** (N - 1) for v in mult.values())
    check(f"N={N}: dense eig multiplicities {list(mult.values())} == uniform {4 ** (N - 1)}", ok)
for N in range(2, 13):
    m = pi_multiplicities_combinatorial(N)
    check(f"N={N}: orbit-argument multiplicities uniform 4^(N-1)",
          all(v == 4 ** (N - 1) for v in m.values()))

print("=== T3: combinatorial Majorana degree == matrix-level monomial degree (N=3,4) ===")
for N in (3, 4):
    dim = 2 ** N
    maj = majoranas(N)
    fails = 0
    for s in product(range(4), repeat=N):
        mask = majorana_mask(s)
        m = np.eye(dim, dtype=complex)
        for k in range(2 * N):
            if (mask >> k) & 1:
                m = m @ maj[k]
        A = string_matrix(s)
        c = np.trace(m.conj().T @ A) / dim
        if abs(abs(c) - 1) > 1e-9:
            fails += 1
    check(f"N={N}: every string proportional to its combinatorial monomial", fails == 0)

print("=== T4: Hermitian phase c = +-i^(d(d-1)/2) (matrix, N=3,4,5) ===")
for N in (3, 4, 5):
    dim = 2 ** N
    maj = majoranas(N)
    fails = 0
    for s in product(range(4), repeat=N):
        mask = majorana_mask(s)
        d = bin(mask).count("1")
        m = np.eye(dim, dtype=complex)
        for k in range(2 * N):
            if (mask >> k) & 1:
                m = m @ maj[k]
        c = np.trace(m.conj().T @ string_matrix(s)) / dim
        expected = 1j ** (d * (d - 1) // 2)
        if not (np.isclose(c, expected) or np.isclose(c, -expected)):
            fails += 1
    check(f"N={N}: phase is +-i^(d(d-1)/2) on all {4 ** N} strings", fails == 0)

print("=== T5: degree-mod-4 class sizes == closed form (combinatorial, N=2..12) ===")
for N in range(2, 13):
    enum = degree_classes(N) if N <= 8 else None
    cf = degree_classes_closed_form(N)
    binom = [sum(comb(2 * N, d) for d in range(2 * N + 1) if d % 4 == r) for r in range(4)]
    ok = cf == binom and (enum is None or enum == cf)
    src = "enumeration+binomial" if enum is not None else "binomial"
    check(f"N={N}: classes {cf} == {src}", ok)

print("=== T6: the two gradings differ at every N (the strong thesis refuted) ===")
for N in range(2, 13):
    pi_m = sorted(pi_multiplicities_combinatorial(N).values())
    deg_m = sorted(degree_classes_closed_form(N))
    check(f"N={N}: Pi {pi_m[0]}x4 vs degree {deg_m} -> gradings differ", pi_m != deg_m,
          f"deviation +-{2 ** (N - 1)}")

print("=== T6b: even the mod-2 shadows are different partitions (N=2..8) ===")
# Clock 1's square is Pi^2 = Ad(X^N), grading strings by (nY+nZ) mod 2. The degree
# clock's mod-2 shadow grades by d mod 2. These are DIFFERENT partitions (witness
# ZI..I: d = 2 even, nY+nZ = 1 odd), so the two clocks do not even share their squares
# as gradings; what they share is X^N appearing as an operator FACTOR in V_g (T8).
for N in range(2, 9):
    mismatch = sum(1 for s in product(range(4), repeat=N)
                   if bin(majorana_mask(s)).count("1") % 2 != sum(1 for c in s if c >= 2) % 2)
    zii = tuple([3] + [0] * (N - 1))
    zii_witness = (bin(majorana_mask(zii)).count("1") % 2 == 0
                   and sum(1 for c in zii if c >= 2) % 2 == 1)
    check(f"N={N}: d mod 2 != (nY+nZ) mod 2 on {mismatch} strings; ZI..I is a witness",
          mismatch > 0 and zii_witness)

print("=== T7: eps_g == degree formula under direct matrix conjugation (N=3..6, both gauges) ===")
for N in (3, 4, 5, 6):
    XN = kron_all([X] * N)
    for gauge_name, parity in (("odd", 1), ("even", 0)):
        Ug = kron_all([Z if l % 2 == parity else I2 for l in range(N)])
        Wg = Ug @ XN
        fails = 0
        for s in product(range(4), repeat=N):
            A = string_matrix(s)
            d = bin(majorana_mask(s)).count("1")
            k = d * (d - 1) // 2 if gauge_name == "odd" else d * (d + 1) // 2
            if not np.allclose(Wg @ A.conj() @ Wg.conj().T, (-1) ** (k % 2) * A):
                fails += 1
        check(f"N={N} {gauge_name} gauge: eps matches on all {4 ** N} strings", fails == 0)

print("=== T8: factorization V_g = Ad(U_g) o Pi^2 o K (random matrix, N=3..6) ===")
rng = np.random.default_rng(7)
for N in (3, 4, 5, 6):
    dim = 2 ** N
    XN = kron_all([X] * N)
    Mrand = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
    for gauge_name, parity in (("odd", 1), ("even", 0)):
        Ug = kron_all([Z if l % 2 == parity else I2 for l in range(N)])
        Wg = Ug @ XN
        lhs = Wg @ Mrand.conj() @ Wg.conj().T
        rhs = Ug @ (XN @ Mrand.conj() @ XN) @ Ug
        check(f"N={N} {gauge_name} gauge: exact factorization",
              np.allclose(lhs, rhs, atol=1e-13), f"max dev {np.abs(lhs - rhs).max():.1e}")

print(f"\n{_checks[1]}/{_checks[0]} checks passed" + ("" if _checks[1] == _checks[0] else "  *** FAILURES ***"))
