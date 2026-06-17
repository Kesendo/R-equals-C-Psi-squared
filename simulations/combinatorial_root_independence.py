"""EXPLORATORY (gate-first): is the COMBINATORIAL root of the small-N specials an arithmetic consequence of
the Niven (cyclotomic) root, or genuinely INDEPENDENT?

The n3_special_cases arc found TWO roots: (A) NUMBER-THEORETIC = Niven's theorem on 2cos(π/(N+1)) (the band-
edge degree, quadratic surd √2/φ/√3 at N=3/4/5); (B) COMBINATORIAL = small-N filling maximality (the even-N
integer half-filling (N/2,N/2) sector, where the structural-ceiling surd LEAK and the (2,2) anomalies live).

Hypothesis to test (gate-first): A and B are INDEPENDENT properties of N — neither implies the other — and
N=4 is merely their COINCIDENCE (plus a third, rep-theoretic, thread: the ceiling 4/N reaching the floor).

THE GATE THAT CAN FIRE: if B were an arithmetic consequence of A, one N-set would contain the other (nested
or equal). Independence is PROVEN iff decoupling witnesses exist on BOTH sides:
   - some N is arithmetic-special but combinatorially-ordinary (A ∧ ¬B), and
   - some N is combinatorially-special but arithmetically-ordinary (B ∧ ¬A).
If only one side has a witness (one set ⊆ the other) the gate fires "not independent — reducible".

Axes (crisp, computable):
  A(N) = the band edge 2cos(π/(N+1)) is a QUADRATIC surd (algebraic degree φ_euler(2(N+1))/2 == 2)  [Niven]
  B(N) = N is EVEN, i.e. an integer half-filling (N/2,N/2) sector exists  [combinatorial]
  the physical LEAK = B(N) ∧ (the (N/2,N/2) commutant darkest < the (1,1) ceiling 4/N), the rep-theoretic
                      vacate condition 4/N ≥ 1 (computed directly, not assumed).

Run: python simulations/combinatorial_root_independence.py
"""
import sys
import numpy as np
import sympy as sp
from itertools import combinations

sys.stdout.reconfigure(encoding="utf-8")
NULLTOL, NZTOL = 1e-7, 1e-7
t = sp.Symbol("t")


def sector_states(n, p):
    return [sum(1 << i for i in c) for c in combinations(range(n), p)]


def sector_H_complete(n, states):
    bonds = [(i, j) for i in range(n) for j in range(i + 1, n)]
    idx = {s: a for a, s in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    for a, s in enumerate(states):
        for (i, j) in bonds:
            if ((s >> i) & 1) != ((s >> j) & 1):
                H[idx[s ^ (1 << i) ^ (1 << j)], a] += 1.0
    return H


def commutant_darkest(n, p, q):
    A, B = sector_states(n, p), sector_states(n, q)
    Hp, Hq = sector_H_complete(n, A), sector_H_complete(n, B)
    na, nb = len(A), len(B)
    adH = np.kron(Hp, np.eye(nb)) - np.kron(np.eye(na), Hq.T)
    w, V = np.linalg.eigh(adH)
    diag = np.array([bin(A[a] ^ B[b]).count("1") for a in range(na) for b in range(nb)], float)
    cols = [k for k in range(len(w)) if abs(w[k]) < NULLTOL]
    if not cols:
        return None
    U = V[:, cols]
    ev = np.linalg.eigvalsh(U.T @ (diag[:, None] * U))
    pos = ev[ev > NZTOL]
    return float(pos.min()) if len(pos) else None


def band_edge_degree(N):
    return int(sp.degree(sp.minimal_polynomial(2 * sp.cos(sp.pi / (N + 1)), t), t))


print("=" * 96)
print("COMBINATORIAL ROOT INDEPENDENCE — is the filling-maximality root an arithmetic consequence? (no)")
print("=" * 96)
print(f"{'N':>2} {'A: band-edge deg':>16} {'A quadratic?':>12} {'B: even?':>9} {'4/N>=1 (vacate)':>15} {'leak computed':>14}")

A_set, B_set, golden_set, vacate_set, leak_set = set(), set(), set(), set(), set()
for N in range(3, 9):
    deg = band_edge_degree(N)
    a = (deg == 2)
    b = (N % 2 == 0)
    golden = sp.simplify(2 * sp.cos(sp.pi / (N + 1)) - (1 + sp.sqrt(5)) / 2) == 0
    vacate = sp.Rational(4, N) >= 1
    leak = ""
    if b and N <= 8:
        c11 = commutant_darkest(N, 1, 1)
        chh = commutant_darkest(N, N // 2, N // 2)
        is_leak = chh is not None and c11 is not None and chh < c11 - 1e-7
        leak = f"{'LEAK' if is_leak else 'no'} ({chh:.3f} vs {c11:.3f})"
        if is_leak:
            leak_set.add(N)
    if a:
        A_set.add(N)
    if b:
        B_set.add(N)
    if golden:
        golden_set.add(N)
    if vacate:
        vacate_set.add(N)
    print(f"{N:>2} {deg:>16} {('YES' if a else 'no'):>12} {('YES' if b else 'no'):>9} "
          f"{('YES' if vacate else 'no'):>15} {leak:>14}")

print(f"\n  A (quadratic band edge, Niven)  = {sorted(A_set)}")
print(f"  B (even, integer half-filling)  = {sorted(B_set)}")
print(f"  golden (band edge = φ)          = {sorted(golden_set)}")
print(f"  vacate (4/N >= 1, rep-theoretic)= {sorted(vacate_set)} (within N>=3)")
print(f"  physical (N/2,N/2) leak         = {sorted(leak_set)}")

print("\n" + "=" * 96)
print("GATES (a firing gate is a finding: diagnose, don't loosen)")
print("=" * 96)

# G1: INDEPENDENCE — decoupling witnesses on BOTH sides (neither set contains the other)
a_not_b = A_set - B_set
b_not_a = B_set - A_set
assert a_not_b, f"G1 FIRED: no arithmetic-special-but-combinatorially-ordinary N (A ⊆ B?) -> reducible"
assert b_not_a, f"G1 FIRED: no combinatorially-special-but-arithmetically-ordinary N (B ⊆ A?) -> reducible"
print(f"G1 PASS (INDEPENDENT): A∖B = {sorted(a_not_b)} (e.g. N=5: quadratic √3 band edge, but ODD — no half-"
      f"filling) and B∖A = {sorted(b_not_a)} (e.g. N=6: even half-filling, but CUBIC heptagon band edge). "
      f"Neither set contains the other -> the combinatorial root is NOT an arithmetic consequence.")

# G2: N=4 is the UNIQUE triple coincidence (golden ∧ even ∧ vacate) among N>=3
triple = golden_set & B_set & vacate_set
assert triple == {4}, f"G2 FIRED: triple-coincidence set {sorted(triple)} != {{4}}"
print(f"G2 PASS (the coincidence): N=4 is the UNIQUE N≥3 where all three independent threads meet — golden "
      f"band edge (pentagon) ∧ even half-filling ∧ ceiling 4/N reaching the floor (4/4=1).")

# G3: the physical leak is exactly the even∧vacate event (rep-theoretic + combinatorial), NOT the arithmetic
assert leak_set == (B_set & vacate_set), f"G3 FIRED: leak {sorted(leak_set)} != even∧vacate {sorted(B_set & vacate_set)}"
assert leak_set == {4}, f"G3 FIRED: leak set {sorted(leak_set)} != {{4}}"
print(f"G3 PASS (the leak is rep-theoretic, not arithmetic): the (N/2,N/2) surd leak = even ∧ (4/N≥1) = "
      f"{{4}}; it needs the rep-theoretic ceiling to vacate AND an even half-filling — independent of the "
      f"cyclotomic golden, which only ALSO happens to single out N=4.")

print("\nROOT: the combinatorial root is INDEPENDENT of the Niven arithmetic root (decoupled at N=5 and N=6).")
print("There are really THREE independent threads of N — cyclotomic (Niven band-edge degree), rep-theoretic")
print("(S_N ceiling 4/N), and combinatorial (even half-filling) — and N=4 is their unique triple coincidence.")
print("That triple coincidence IS why N=4 keeps surfacing. The arc's 'two roots' sharpens to 'independent")
print("axes, N=4 = coincidence'. DONE.")
