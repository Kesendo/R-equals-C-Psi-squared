"""Gate-first probe: does the '4' in g2(K_N) = 4/N generalize as 2d, d², or stay 4 at qudit dimension d?

CONTEXT (the seam under test).
  * PolynomialDiscriminantAnchorClaim records 4 = discriminant of d²−2d=0 = d² for d=2 (the qubit).
  * The structural ceiling (F122, StructuralCeilingClaim) has g2(K_N) = 4/N. The "4" there is
    UNTYPED-linked to the discriminant 4: the two facts live in disjoint parts of the Claim graph.
  * Two-fours discipline (PROOF_ABSORPTION_THEOREM §2 remark): "same number, two genealogies; never
    interchangeably." At d=2, the rung-two 2d = 4 and the discriminant d² = 4 COLLAPSE. The qutrit
    (d=3) SPLITS them: 2d = 6 vs d² = 9 vs (stays 4).

HYPOTHESES (what the numerator of g2(K_N) becomes at d):
  H_2d : numerator -> 2d  (rung/cap genealogy; the (2d)^N F121 qudit cap trunk)  => g2(K_N,3) = 6/N
  H_d2 : numerator -> d²  (discriminant genealogy; "everything is one object")    => g2(K_N,3) = 9/N
  H_4  : numerator stays 4 (Hamming-2 between single excitations × S_N angle 2/N;  => g2(K_N,3) = 4/N
         both factors are d-INDEPENDENT)

THE RECIPE (faithful lift of topology_ceiling_rep_derivation.py / StructuralCeilingWitness).
  In the |a><b| computational-coherence basis the full-Cartan equidistant dephasing (F121) is DIAGONAL:
  N_diff(a,b) = #{ sites l : a_l != b_l } = the qudit Hamming distance (= n_XY at d=2). The high-Q
  structural ceiling = min nonzero eigenvalue of N_diff restricted to the H-commutant (ad_H kernel) of
  an excitation sector. The Hamiltonian is the FAITHFUL XY-lift = vacuum hopping (hard-core SU(d)
  bosons): a non-vacuum flavor f hops to an adjacent vacuum site, flavor preserved. At d=2 (one flavor)
  this IS the XY model, so the recipe reduces EXACTLY to the qubit one (Stage 0 gate).

  (Why vacuum-hopping and not SU(d) swap: swap at d=2 is Heisenberg/XXX, not XY, and would change even
   the d=2 answer. The reduction-to-4/N requirement forces the XY-lift.)

STAGE 0  -- port-fidelity gate: d=2 reproduces 4/N (complete) and 4/(N−1) (star). MUST pass.
STAGE 0b -- qutrit oracle gate: d=3 small-N degenerate-PT global g2 == full Liouvillian g2 at high Q.
STAGE 1  -- the probe: d=3 K_N (1,1) commutant darkest + global min; reads off the numerator. The
            primary gate asserts H_4 (the d-independent 4/N, my analytic prior). A FIRING gate IS the
            find: diagnose which hypothesis the data picks, do not loosen.

Run:  python simulations/qudit_g2_split.py
"""
import sys
import itertools
from math import sqrt
import numpy as np

GAMMA = 0.05
TOL_REL = 1e-6     # "same Omega" / "Omega=0" cluster tolerance, relative to the Omega scale
NZ = 1e-7          # nonzero-N_diff threshold (steady modes have eigenvalue 0 exactly)


# ----------------------------------------------------------------------------------------------------
# qudit basis, sectors, faithful XY-lift Hamiltonian, Hamming dephasing
# ----------------------------------------------------------------------------------------------------
def topo_bonds(topo, N):
    if topo == 'chain':
        return [(i, i + 1) for i in range(N - 1)]
    if topo == 'ring':
        return [(i, (i + 1) % N) for i in range(N)]
    if topo == 'star':
        return [(0, i) for i in range(1, N)]
    if topo == 'complete':
        return [(i, j) for i in range(N) for j in range(i + 1, N)]
    raise ValueError(topo)


def all_states(d, N):
    """All d^N computational strings as tuples in {0,..,d-1}^N (0 = vacuum)."""
    return list(itertools.product(range(d), repeat=N))


def excit(s):
    """Excitation number = number of non-vacuum digits."""
    return sum(1 for x in s if x != 0)


def sector_states(d, N, p):
    """Computational strings with exactly p non-vacuum digits."""
    return [s for s in all_states(d, N) if excit(s) == p]


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def sector_H(topo, d, N, states):
    """Faithful XY-lift (vacuum hopping) restricted to a fixed-excitation sector, real symmetric.
    A non-vacuum flavor f at site i hops to an adjacent VACUUM site j (and h.c.), flavor preserved.
    At d=2 this is exactly the XY hopping H = (1/2) sum (X_iX_j + Y_iY_j) on the popcount sector."""
    m = len(states)
    idx = {s: a for a, s in enumerate(states)}
    H = np.zeros((m, m))
    bonds = topo_bonds(topo, N)
    for a, s in enumerate(states):
        for (i, j) in bonds:
            si, sj = s[i], s[j]
            if (si == 0) != (sj == 0):          # exactly one site is vacuum -> the flavor can hop
                s2 = list(s)
                s2[i], s2[j] = sj, si           # move the flavor across the bond (swap with the vacuum)
                s2 = tuple(s2)
                H[idx[s2], a] += 1.0            # symmetric: the reverse hop fills the transpose entry
    return H


def commutant_darkest(topo, d, N, p, q):
    """min nonzero eigenvalue of N_diff restricted to the ad_H kernel of the (p,q) coherence sector.
    Faithful qudit port of StructuralCeilingWitness.CommutantDarkest. None if empty / no non-steady mode."""
    A = sector_states(d, N, p)
    B = sector_states(d, N, q)
    na, nb = len(A), len(B)
    if na == 0 or nb == 0:
        return None
    Hp = sector_H(topo, d, N, A)
    Hq = sector_H(topo, d, N, B)
    adH = np.kron(Hp, np.eye(nb)) - np.kron(np.eye(na), Hq)      # real symmetric; eigenvalues = Omega
    Omega, V = np.linalg.eigh(adH)
    diag = np.array([hamming(A[a], B[b]) for a in range(na) for b in range(nb)], dtype=float)
    oscale = 1.0 + (np.abs(Omega).max() if Omega.size else 0.0)
    z = np.abs(Omega) < TOL_REL * oscale                        # the ad_H kernel (commutant, [H,A]=0)
    if not z.any():
        return None
    U = V[:, z]
    Ntil = U.T @ (diag[:, None] * U)                            # N_diff projected onto the commutant
    w = np.linalg.eigvalsh(Ntil)
    nz = w[w > NZ]
    return float(nz.min()) if nz.size else None


def global_g2(topo, d, N, pmax=None):
    """min nonzero commutant-darkest over all (p,q) sectors (the high-Q structural ceiling).
    pmax caps the excitation number scanned (sectors get large for middle p at big N/d)."""
    if pmax is None:
        pmax = N
    best, where = np.inf, None
    table = {}
    for p in range(pmax + 1):
        for q in range(pmax + 1):
            v = commutant_darkest(topo, d, N, p, q)
            if v is None:
                continue
            table[(p, q)] = v
            if v < best:
                best, where = v, (p, q)
    return (None if best == np.inf else best), where, table


# ----------------------------------------------------------------------------------------------------
# full Liouvillian oracle (Stage 0b): build L on the full d^N Hilbert space, read g2 = strict_gap/2γ
# ----------------------------------------------------------------------------------------------------
def full_H(topo, d, N, J):
    """Vacuum-hopping H on the full d^N space (all excitation sectors at once)."""
    states = all_states(d, N)
    m = len(states)
    idx = {s: a for a, s in enumerate(states)}
    H = np.zeros((m, m))
    bonds = topo_bonds(topo, N)
    for a, s in enumerate(states):
        for (i, j) in bonds:
            si, sj = s[i], s[j]
            if (si == 0) != (sj == 0):
                s2 = list(s); s2[i], s2[j] = sj, si; s2 = tuple(s2)
                H[idx[s2], a] += J
    return H, states


def full_L_g2(topo, d, N, Q=1500.0):
    """Full Liouvillian g2 = smallest nonzero decay rate / (2γ) at high Q. Dephasing = full-Cartan
    equidistant: D[|a><b|] = -2γ·hamming(a,b)·|a><b| (diagonal in the coherence basis; = Z-dephasing at d=2)."""
    J = Q * GAMMA
    H, states = full_H(topo, d, N, J)
    m = len(states)
    I = np.eye(m)
    Lham = -1j * (np.kron(H, I) - np.kron(I, H.T))              # -i[H, .]  (H real symmetric)
    deph = np.array([-2.0 * GAMMA * hamming(states[a], states[b])
                     for a in range(m) for b in range(m)])      # diagonal dissipator, coherence basis
    L = Lham + np.diag(deph)
    ev = np.linalg.eigvals(L)
    rates = -ev.real
    dec = rates > 1e-9
    return float(rates[dec].min()) / (2 * GAMMA)


def fmt(x):
    return "—" if x is None else f"{x:.9f}"


# ====================================================================================================
# STAGE 0 -- PORT-FIDELITY GATE: d=2 reproduces 4/N (complete) and 4/(N−1) (star)
# ====================================================================================================
print("=" * 100)
print("STAGE 0 -- PORT FIDELITY: the qudit recipe at d=2 reproduces the known g2 closed forms")
print("=" * 100)
print(f"{'topo':9} {'N':>2} {'(1,1) darkest':>15} {'closed form':>14} {'==exact?':>9}")
worst0 = 0.0
for N in (4, 5, 6):
    v = commutant_darkest('complete', 2, N, 1, 1)
    cf = 4.0 / N
    worst0 = max(worst0, abs(v - cf))
    print(f"{'complete':9} {N:>2} {fmt(v):>15} {f'4/{N}':>14} {('YES' if abs(v-cf)<1e-9 else 'no'):>9}")
for N in (6, 7):
    v = commutant_darkest('star', 2, N, 1, 1)
    cf = 4.0 / (N - 1)
    worst0 = max(worst0, abs(v - cf))
    print(f"{'star':9} {N:>2} {fmt(v):>15} {f'4/{N-1}':>14} {('YES' if abs(v-cf)<1e-9 else 'no'):>9}")
assert worst0 < 1e-9, f"STAGE 0 GATE FIRED: qudit recipe does not reduce to the qubit g2 (worst {worst0:.2e})"
print(f"\nSTAGE 0 PASS: the qudit recipe reduces EXACTLY to the qubit g2 at d=2 (worst |diff| = {worst0:.1e}).")

# ====================================================================================================
# STAGE 0b -- QUTRIT ORACLE GATE: degenerate-PT global g2 == full Liouvillian g2 at high Q
# ====================================================================================================
print("\n" + "=" * 100)
print("STAGE 0b -- QUTRIT ORACLE: degenerate-PT global g2  vs  full Liouvillian g2 (Q=1500), d=3")
print("=" * 100)
print(f"{'topo':9} {'N':>2} {'g2_PT':>12} {'g2_fullL':>12} {'win sector':>11} {'|diff|':>9}")
worst0b = 0.0
for topo in ('chain', 'complete'):
    for N in (2, 3):
        gpt, where, _ = global_g2(topo, 3, N)
        gfl = full_L_g2(topo, 3, N)
        dd = abs(gpt - gfl)
        worst0b = max(worst0b, dd)
        print(f"{topo:9} {N:>2} {fmt(gpt):>12} {fmt(gfl):>12} {str(where):>11} {dd:>9.1e}")
assert worst0b < 3e-3, f"STAGE 0b GATE FIRED: qutrit degenerate-PT disagrees with full L by {worst0b:.2e}"
print(f"\nSTAGE 0b PASS: the qutrit dephasing+H model and the sector recipe agree with the full Liouvillian "
      f"(worst |diff| = {worst0b:.1e}, the O(1/Q) residual).")

# ====================================================================================================
# STAGE 1 -- THE PROBE: d=3 K_N (1,1) commutant darkest + global min; read off the numerator
# ====================================================================================================
print("\n" + "=" * 100)
print("STAGE 1 -- THE QUTRIT SPLIT: g2(K_N) at d=3 -- is the numerator 2d=6, d^2=9, or still 4?")
print("=" * 100)
print(f"{'N':>2} {'(1,1) darkest':>15} {'g2*N':>8} || {'4/N':>9} {'6/N':>9} {'9/N':>9} || {'picks':>7}")
for N in (4, 5, 6, 7):
    v = commutant_darkest('complete', 3, N, 1, 1)
    num = v * N
    cands = {'4 (H_4)': 4.0 / N, '6 (H_2d)': 6.0 / N, '9 (H_d2)': 9.0 / N}
    pick = min(cands, key=lambda k: abs(cands[k] - v))
    pick = pick if abs(cands[pick] - v) < 1e-6 else f"OTHER({num:.4f})"
    print(f"{N:>2} {fmt(v):>15} {num:>8.4f} || {4.0/N:>9.5f} {6.0/N:>9.5f} {9.0/N:>9.5f} || {pick:>7}")

# global-min cross-check at the N where every sector is tractable: is 4/N still the GLOBAL ceiling
# (or does a flavor-enriched middle sector at d=3 undercut it)?  N=5 takes ~75s (the (3,3)/(4,4) sectors).
print("\n  global-min cross-check (d=3, all sectors), confirming where the ceiling lives:")
print(f"  {'N':>2} {'global g2':>12} {'win sector':>11} {'(1,1) value':>13} {'note':>26}")
gres = {}
for N in (3, 4, 5):
    g, where, table = global_g2('complete', 3, N)
    gres[N] = (g, where, table)
    note = ("band edge 1.0 (4/N>1)" if abs(g - 1.0) < 1e-9 and 4.0 / N >= 1.0 else
            "(2,2) outlier 2-2/sqrt3" if N == 4 else
            "= 4/N global ceiling" if abs(g - 4.0 / N) < 1e-9 else "(1,1) UNDERCUT")
    print(f"  {N:>2} {fmt(g):>12} {str(where):>11} {fmt(table.get((1,1))):>13} {note:>26}")
# the N=4 (2,2) outlier value is itself d-independent: 2-2/sqrt(3) at BOTH d=2 and d=3;
# and at N=5 the GLOBAL ceiling (not just the (1,1) sector) is 4/5 -- no flavor sector undercuts it.
assert abs(gres[4][0] - (2.0 - 2.0 / sqrt(3.0))) < 1e-7, f"K_4 d=3 outlier expected 2-2/sqrt3, got {gres[4][0]:.9f}"
assert abs(gres[5][0] - 4.0 / 5) < 1e-9, f"K_5 d=3 GLOBAL ceiling expected 4/5, got {gres[5][0]:.9f} in {gres[5][1]}"

# PRIMARY GATE: H_4 (the d-independent 4/N). A firing gate is the find -> diagnose.
print("\n" + "-" * 100)
val5 = commutant_darkest('complete', 3, 5, 1, 1)
val6 = commutant_darkest('complete', 3, 6, 1, 1)
assert abs(val5 - 4.0 / 5) < 1e-9 and abs(val6 - 4.0 / 6) < 1e-9, (
    f"STAGE 1 GATE FIRED (H_4 refuted): g2(K_N,d=3) (1,1) darkest is NOT 4/N -- "
    f"N=5 gave {val5:.9f} (4/5={4/5:.9f}), N=6 gave {val6:.9f} (4/6={4/6:.9f}). "
    f"·N numerators: N=5 -> {val5*5:.6f}, N=6 -> {val6*6:.6f}. Diagnose: 2d=6? d²=9? other?")
print(f"STAGE 1 PASS (H_4 holds): g2(K_N, d=3) = 4/N to machine precision -- the numerator 4 is "
      f"d-INDEPENDENT,\n  confirmed at the GLOBAL ceiling (N=5: 4/5, not just the (1,1) sector). It is "
      f"Hamming-2 (between two\n  single-excitation strings) * S_N angle 2/N, NOT 2d and NOT d^2 -- both "
      f"factors are d-independent.\n  The seam to the discriminant four (d^2 -> 9 at the qutrit) is "
      f"REFUTED for g2: the two-fours discipline\n  holds and EXTENDS -- the qutrit is the prism that "
      f"splits the d=2 coincidence d^2=2d=4 into {{9, 6, 4}},\n  and g2 rides the 4-ray, a third "
      f"genealogy distinct from both the discriminant (d^2) and the cap (2d).")

# ====================================================================================================
# STAGE 2 -- THE CAP LEG (2d) vs F121: the discriminant four and the cap four ARE the two terms of the
#            trunk polynomial d^2 - 2d = 0, equal only at its root d=2, split by the qutrit.
# ====================================================================================================
# The foundational polynomial d^2 - 2d has TWO terms; the qutrit splits them:
#   * the d^2 term = the per-site coherence-space dimension = the squared dimension = the discriminant /
#                    a_-1 reading (PolynomialDiscriminantAnchorClaim).            -> d^2 = 9 at d=3.
#   * the 2d term  = the F121 product-mirror cap base: each per-site mirror q_l is a strict class swap
#                    between the d DARK letters {(x,x)} and the d^2-d LIT letters {(i,j), i!=j}, so
#                    rank(q_l) <= min(d, d^2-d) + min(d^2-d, d) = 2d, and the cap = (2d)^N.  -> 2d = 6.
# They coincide at d=2 (= 4, the trunk root d^2-2d=0; F121's "third appearance" of the trunk). We
# RECOMPUTE the cap base directly (generic class-swap rank, F121 PROOF_QUDIT_PARTIAL_PALINDROME BLOCK C)
# and cross-check the full cap (2d)^N against F121's measured 36/216/64. The g2 ceiling four (Stage 1)
# is the THIRD reading, EXTERNAL to the trunk (2 x 2/N, stays 4).
print("\n" + "=" * 100)
print("STAGE 2 -- THE CAP LEG (2d) vs F121: the two terms of the trunk d^2-2d=0, split by the qutrit")
print("=" * 100)


def random_class_swap(d, seed):
    """A generic per-site dark<->lit class swap on the d^2-dim single-site coherence space (the F121
    per-site mirror q_l). Block-anti-diagonal (dark->lit and lit->dark only), so rank = 2d for d>=2."""
    dark = [i * d + i for i in range(d)]                          # (x,x) coherence indices, d of them
    lit = [i * d + j for i in range(d) for j in range(d) if i != j]   # (i,j) i!=j, d^2-d of them
    q = np.zeros((d * d, d * d))
    rng = np.random.default_rng(seed)
    for c in dark:
        for r in lit:
            q[r, c] = rng.standard_normal()                      # dark -> lit block
    for c in lit:
        for r in dark:
            q[r, c] = rng.standard_normal()                      # lit  -> dark block
    return q


print(f"{'d':>2} {'cap base (rank)':>16} {'2d':>5} {'d^2':>5} {'trunk d^2-2d':>13} {'rides':>8}")
for d in (2, 3, 4, 5):
    base = int(np.linalg.matrix_rank(random_class_swap(d, 20260617 + d), tol=1e-9))
    twod, dsq = 2 * d, d * d
    rides = "2d=d^2" if twod == dsq else ("2d" if base == twod else "?")
    print(f"{d:>2} {base:>16} {twod:>5} {dsq:>5} {dsq - 2 * d:>13} {rides:>8}")
    assert base == twod, f"STAGE 2 GATE FIRED: per-site cap base at d={d} is {base}, not 2d={twod}"

# direct (3,2) product rank = 36 (reproduces F121 BLOCK C): the cap is 6^2, not 9^2 = 81
q3 = random_class_swap(3, 99)
rank32 = int(np.linalg.matrix_rank(np.kron(q3, q3), tol=1e-9))
assert rank32 == 36, f"STAGE 2 GATE FIRED: (3,2) class-swap product rank {rank32} != 36 = (2d)^N"
print(f"\n  direct (d=3,N=2) class-swap product rank = {rank32} = (2d)^N = 6^2  (the d^2-base would give "
      f"d^(2N) = 81; F121 measures 36)")

# cross-check the closed form (2d)^N against F121's measured values (qudit_product_mirror_cap.py)
F121_CAP = {(3, 2): 36, (3, 3): 216, (4, 2): 64}                 # = 6^2, 6^3, 8^2 ; NOT 81/729/256 = d^(2N)
for (d, N), measured in F121_CAP.items():
    cap, dsq_pow = (2 * d) ** N, (d * d) ** N
    assert cap == measured and cap != dsq_pow, f"STAGE 2 GATE FIRED: (2d)^N at (d={d},N={N}) != F121 {measured}"
    print(f"  F121 cap({d},{N}) = (2d)^N = {cap}  (the d^2-base would give d^(2N) = {dsq_pow})")
assert (3 * 3 - 2 * 3) != 0 and (2 * 2 - 2 * 2) == 0, "STAGE 2: trunk d^2-2d=0 holds only at d=2"

# the three fours, fanned by the qutrit prism (the ceiling four pulled live from Stage 1)
ceil4 = round(commutant_darkest('complete', 3, 5, 1, 1) * 5)    # = 4 (d-independent, Stage 1)
assert ceil4 == 4, "STAGE 2: the ceiling four should stay 4 at d=3"
print(f"\n  THE THREE FOURS (all = 4 at d=2, fanned by the qutrit prism):")
print(f"    {'reading':32} {'d=2':>4} {'d=3':>4}  {'origin':<36}")
print(f"    {'squared dimension  d^2':32} {4:>4} {9:>4}  {'d^2 term of the trunk (a_-1, dyadic)':<36}")
print(f"    {'product-mirror cap  2d':32} {4:>4} {6:>4}  {'2d term of the trunk (F121 cap base)':<36}")
print(f"    {'structural ceiling  2 x (2/N)':32} {4:>4} {ceil4:>4}  {'trunk-EXTERNAL (Hamming x S_N angle)':<36}")
print(f"\nSTAGE 2 PASS: the cap leg rides 2d (= 6 at d=3, F121 cap 36 = 6^2, NOT 81 = 9^2) and the "
      f"discriminant\n  leg rides d^2 (= 9). They are the TWO TERMS of the trunk d^2 - 2d = 0 -- equal only "
      f"at its root\n  d=2, fanned to {{9, 6}} by the qutrit. The g2 ceiling four (stays 4) is the THIRD, "
      f"trunk-EXTERNAL\n  reading. The one d=2 four fans into {{9, 6, 4}}.")
print("\nDONE.")
