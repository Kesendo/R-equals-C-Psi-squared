"""qd_letter_gates.py: pre-registered gates for the record letter law (F136).

Proof: docs/proofs/PROOF_RECORD_LETTER_LAW.md. Committed output copy:
simulations/results/qd_letter/qd_letter_gates_out.txt. The numerics never use
the channel formulas under test: each configuration is built as the full
2^N x 2^N closed-form state (Absorption substrate) and partial-traced to the
pair, so the algebra and the measurement are independent computational paths.

Derivation under test (from F135 Prop 1, all by hand BEFORE this run):
pair (S,j), all S-bonds Delta_S, t* = pi/(4 Delta_S). D = shared dressers
(ratio r = Delta_jk/Delta_S), P = private watchers of j, Q = private
neighbors of S (besides j).
  H-hinge: S has a neighbor besides j (D u Q nonempty <=> V_S(t*) = 0);
         a pure pendant S instead ROLE-SWAPS into j's witness (G15).
  H-Q:  Q nonempty -> both double coherences 0 (kills the Bell record).
  H-par: BELL record alive iff Q empty, m >= 1, every D-ratio ODD and every
         P-ratio EVEN. POINTER record alive iff the write bond is present,
         deg(S) >= 2, and EVERY watcher of j (D and P alike) is EVEN --
         the corrected classification (a review round caught the first
         draft's D-empty requirement; G14 gates the corrected region;
         a second round caught the missing hinge; G15 gates the pendant).
  Then c1 = (1/4)*prod_D(-1)^((1+r)/2)*prod_P(-1)^(r/2),
       c2 = (1/4)*prod_D(-1)^((1-r)/2)*prod_P(-1)^(r/2),
  letter: m=|D| odd -> YY = 2(c2-c1) = +-1, XX = 0;
          m even    -> XX = 2(c1+c2) = +-1, YY = 0.
  Bell mixture: eigenvalues {1/2,1/2,0,0}, marginals I/2, TD(Z-cond)=0.
  gamma: c1,c2 *= exp(-2(gamma_S+gamma_j)t*); I = 1 - h2((1+kappa)/2).
         gamma on D sites: exactly invisible.
  Pointer side (Law A + signed corollary): unwatched leaf ZY=+1; private
  even watcher r=2 -> ZY=-1 (pi rotation), r=4 -> ZY=+1.
Every gate prints PRED vs GOT and OK/FAIL.
"""

import numpy as np
from itertools import product

t_star = np.pi / 4.0
FAILS = []

PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def state(N, bonds, t, gammas=None):
    """bonds: list of (a, b, Delta). gammas: per-site rates or None."""
    z = np.array(list(product([1, -1], repeat=N)))
    E = np.zeros(len(z))
    for (a, b, d) in bonds:
        E += d * z[:, a] * z[:, b]
    rho = np.exp(-1j * t * (E[:, None] - E[None, :])) / 2**N
    if gammas is not None:
        g = np.asarray(gammas, dtype=float)
        ham = ((z[:, None, :] != z[None, :, :]) * g).sum(axis=2)
        rho = rho * np.exp(-2.0 * t * ham)
    return rho


def ptrace_pair(rho, N, a, b):
    keep = [a, b]
    dims = [2] * N
    r = rho.reshape(dims + dims)
    order = keep + [s for s in range(N) if s not in keep]
    r = np.transpose(r, order + [N + s for s in order])
    r = r.reshape(4, 2 ** (N - 2), 4, 2 ** (N - 2))
    return np.einsum("ikjk->ij", r)


def vn(rho):
    w = np.linalg.eigvalsh(rho)
    w = w[w > 1e-12]
    return float(-(w * np.log2(w)).sum())


def h2(p):
    if p <= 0 or p >= 1:
        return 0.0
    return float(-p * np.log2(p) - (1 - p) * np.log2(1 - p))


def measure(N, bonds, S, j, gammas=None):
    rho = state(N, bonds, t_star, gammas)
    pair = ptrace_pair(rho, N, S, j)
    rS = np.einsum("ikjk->ij", pair.reshape(2, 2, 2, 2))
    rj = np.einsum("kikj->ij", pair.reshape(2, 2, 2, 2))
    I = vn(rS) + vn(rj) - vn(pair)
    corr = {A + B: float(np.real(np.trace(pair @ np.kron(PAULI[A], PAULI[B]))))
            for A in "XYZ" for B in "XYZ"}
    p4 = pair.reshape(2, 2, 2, 2)
    c0, c1 = p4[0, :, 0, :], p4[1, :, 1, :]
    n0, n1 = np.trace(c0).real, np.trace(c1).real
    td = float("nan")
    if min(n0, n1) > 1e-12:
        d = c0 / n0 - c1 / n1
        td = float(0.5 * np.abs(np.linalg.eigvalsh(d)).sum())
    ev = np.sort(np.linalg.eigvalsh(pair))[::-1]
    return I, corr, td, ev


COUNT = [0]


def gate(name, pred, got, tol=1e-9):
    COUNT[0] += 1
    ok = abs(pred - got) <= tol
    if not ok:
        FAILS.append(name)
    print(f"  [{'OK' if ok else 'FAIL'}] {name}: pred {pred:+.6f} got {got:+.6f}")


print("record-letter GATES (predictions fixed in code before running)")

# --- G1-G3: single dresser K2,1 at odd ratios: letter Y, sign (-1)^((1-r)/2)
for r, sgn in [(1, +1.0), (3, -1.0), (5, +1.0)]:
    N, bonds = 3, [(0, 2, 1.0), (1, 2, float(r))]
    I, corr, td, ev = measure(N, bonds, 0, 1)
    print(f"G K2,1 r={r}:")
    gate(f"I=1 (r={r})", 1.0, I)
    gate(f"YY sign (r={r})", sgn, corr["YY"])
    gate(f"XX=0 (r={r})", 0.0, corr["XX"])
    gate(f"TD=0 (r={r})", 0.0, td)

# --- G4: K2,2 ratios (1,1) -> XX=+1 ; (1,3) -> XX=-1
for rr, sgn in [((1, 1), +1.0), ((1, 3), -1.0)]:
    N = 4
    bonds = [(0, 2, 1.0), (0, 3, 1.0), (1, 2, float(rr[0])), (1, 3, float(rr[1]))]
    I, corr, td, ev = measure(N, bonds, 0, 1)
    print(f"G K2,2 r={rr}:")
    gate(f"I=1 {rr}", 1.0, I)
    gate(f"XX sign {rr}", sgn, corr["XX"])
    gate(f"YY=0 {rr}", 0.0, corr["YY"])

# --- G5: triangle + private watcher on j: r_P=2 -> YY=-1 ; r_P=4 -> YY=+1
for rp, sgn in [(2, -1.0), (4, +1.0)]:
    N = 4
    bonds = [(0, 1, 1.0), (0, 2, 1.0), (1, 2, 1.0), (1, 3, float(rp))]
    I, corr, td, ev = measure(N, bonds, 0, 1)
    print(f"G triangle + private r={rp} on j:")
    gate(f"I=1 (rp={rp})", 1.0, I)
    gate(f"YY sign (rp={rp})", sgn, corr["YY"])

# --- G6: kill conditions -> exactly dark
print("G kill conditions:")
I, corr, td, ev = measure(3, [(0, 2, 1.0), (1, 2, 2.0)], 0, 1)     # even dresser
gate("even dresser r=2 -> I=0", 0.0, I)
I, corr, td, ev = measure(4, [(0, 1, 1.0), (0, 2, 1.0), (1, 2, 1.0),
                              (1, 3, 1.0)], 0, 1)                   # odd private
gate("odd private r=1 -> I=0", 0.0, I)
I, corr, td, ev = measure(4, [(0, 2, 1.0), (1, 2, 1.0), (0, 3, 1.0)], 0, 1)  # Q nonempty
gate("S-private neighbor -> I=0", 0.0, I)

# --- G7: structure at the uniform square
I, corr, td, ev = measure(4, [(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0),
                              (3, 0, 1.0)], 0, 2)
print("G square structure:")
gate("eig1=1/2", 0.5, float(ev[0]))
gate("eig2=1/2", 0.5, float(ev[1]))
gate("eig3=0", 0.0, float(ev[2]))
gate("XY=0", 0.0, corr["XY"])
gate("YX=0", 0.0, corr["YX"])

# --- G8: gamma. Square, gamma_S=gamma_j=0.05 -> I = 1-h2((1+kappa)/2)
g = 0.05
kappa = np.exp(-2 * (g + g) * t_star)
gam = [g, 0.0, g, 0.0]  # S=0, j=2 watched; dressers 1,3 unwatched
I, corr, td, ev = measure(4, [(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0),
                              (3, 0, 1.0)], 0, 2, gammas=gam)
print("G gamma race (Bell record pays BOTH sites):")
gate("I = 1-h2((1+kappa)/2)", 1 - h2((1 + kappa) / 2), I, tol=1e-7)
gate("XX = kappa", float(kappa), corr["XX"], tol=1e-7)
gam = [0.0, 0.3, 0.0, 0.3]  # gamma only on the dressers
I, corr, td, ev = measure(4, [(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0),
                              (3, 0, 1.0)], 0, 2, gammas=gam)
gate("gamma on dressers invisible: I=1", 1.0, I, tol=1e-9)

# --- G9: pointer side, signed corollary. chain3 S=1, j=0 leaf + private watcher.
print("G pointer side (Law A + sign):")
I, corr, td, ev = measure(3, [(0, 1, 1.0), (1, 2, 1.0)], 1, 0)
gate("leaf ZY=+1", 1.0, corr["ZY"])
gate("leaf I=1", 1.0, I)
for rp, sgn in [(2, -1.0), (4, +1.0)]:
    bonds = [(0, 1, 1.0), (1, 2, 1.0), (0, 3, float(rp))]
    I, corr, td, ev = measure(4, bonds, 1, 0)
    gate(f"watched-leaf rp={rp}: ZY sign", sgn, corr["ZY"])
    gate(f"watched-leaf rp={rp}: I=1", 1.0, I)

# --- G10: the asymmetry: pointer record ignores gamma_S, Bell record does not.
print("G gamma asymmetry pointer vs Bell:")
I, corr, td, ev = measure(3, [(0, 1, 1.0), (1, 2, 1.0)], 1, 0,
                          gammas=[0.0, 0.3, 0.0])  # gamma_S only
gate("pointer record: gamma_S invisible, I=1", 1.0, I, tol=1e-9)

# --- G11: the uniform sightings that opened the arc, as gates. Triangle and
# square baselines; pentagon and hexagon exactly dark at t* (kill law on
# 5- and 6-cycles: odd private watchers and nonempty Q, per the theorem).
print("G uniform-cycle catalogue:")
N = 5
k23 = [(0, 2 + k, 1.0) for k in range(3)] + [(1, 2 + k, 1.0) for k in range(3)]
I, corr, td, ev = measure(N, k23, 0, 1)
gate("K2,3 physical frame: YY=+1 (m=3, the ladder refuter)", 1.0, corr["YY"])
gate("K2,3 physical frame: XX=0", 0.0, corr["XX"])
I, corr, td, ev = measure(3, [(0, 1, 1.0), (1, 2, 1.0), (0, 2, 1.0)], 0, 1)
gate("triangle YY=+1", 1.0, corr["YY"])
I, corr, td, ev = measure(4, [(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0),
                              (3, 0, 1.0)], 0, 2)
gate("square XX=+1", 1.0, corr["XX"])
pent = [(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 4, 1.0), (4, 0, 1.0)]
for j in (1, 2):
    I, corr, td, ev = measure(5, pent, 0, j)
    gate(f"pentagon (0,{j}) dark", 0.0, I)
hexg = [(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 4, 1.0), (4, 5, 1.0),
        (5, 0, 1.0)]
for j in (1, 2, 3):
    I, corr, td, ev = measure(6, hexg, 0, j)
    gate(f"hexagon (0,{j}) dark", 0.0, I)

# --- G12: mutual exclusivity at full strength: a perfect pointer record has
# dead Bell channels, a perfect Bell record has a dead pointer channel.
print("G exclusivity pointer vs Bell:")
I, corr, td, ev = measure(3, [(0, 1, 1.0), (1, 2, 1.0)], 1, 0)
gate("leaf: XX=0", 0.0, corr["XX"])
gate("leaf: YY=0", 0.0, corr["YY"])
I, corr, td, ev = measure(4, [(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0),
                              (3, 0, 1.0)], 0, 2)
gate("square: ZY=0", 0.0, corr["ZY"])
gate("square: ZZ=0", 0.0, corr["ZZ"])

# --- G14: the pointer family with SHARED even watchers (reviewer-found
# region: the theorem's first draft wrongly required D = empty). An even
# dresser kills the doubles but forgives the single-j record; with the
# direct bond present this is a perfect pointer record, sign
# prod_{D u P}(-1)^{r/2}. And without the direct bond there is no writer:
# D = empty + P even + no S-j bond is dark.
print("G pointer family with even dressers:")
I, corr, td, ev = measure(3, [(0, 1, 1.0), (0, 2, 1.0), (1, 2, 2.0)], 0, 1)
gate("even-dresser triangle: I=1", 1.0, I)
gate("even-dresser triangle: ZY=-1", -1.0, corr["ZY"])
gate("even-dresser triangle: TD=1", 1.0, td)
gate("even-dresser triangle: XX=0", 0.0, corr["XX"])
gate("even-dresser triangle: YY=0", 0.0, corr["YY"])
I, corr, td, ev = measure(4, [(0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0),
                              (1, 2, 2.0), (1, 3, 2.0)], 0, 1)
gate("direct + two even dressers: I=1", 1.0, I)
gate("direct + two even dressers: ZY=+1", 1.0, corr["ZY"])
I, corr, td, ev = measure(5, [(0, 2, 1.0), (0, 4, 1.0), (1, 3, 2.0)], 0, 1)
gate("no write bond, D empty, P even: dark", 0.0, I)

# --- G15: the hinge and the pendant role-swap (second-review-found region:
# without a second S-neighbor, V_S survives and the classification's "two
# ways" does not apply). Pendant S with an ODD watcher on j: exactly 1
# separable bit, but it is j's pointer recorded in S's equator (YZ channel,
# the role swap; MI is symmetric, the system label is the reading's choice).
# EVEN watcher: 2 bits, entangled. Restoring the hinge (any second
# S-neighbor) closes the channel.
print("G hinge and pendant role-swap:")
I, corr, td, ev = measure(3, [(0, 1, 1.0), (1, 2, 1.0)], 0, 1)
gate("pendant S, odd watcher: I=1", 1.0, I)
gate("pendant S, odd watcher: YZ=+1 (records Z_j)", 1.0, corr["YZ"])
gate("pendant S, odd watcher: ZY=0 (no Z_S record)", 0.0, corr["ZY"])
I, corr, td, ev = measure(3, [(0, 1, 1.0), (1, 2, 2.0)], 0, 1)
gate("pendant S, even watcher: I=2 (entangled)", 2.0, I)
I, corr, td, ev = measure(4, [(0, 1, 1.0), (1, 2, 1.0), (0, 3, 1.0)], 0, 1)
gate("hinge restored (Q nonempty), odd watcher: dark", 0.0, I)

# --- G16: the fan-out corollary. Law B bounds POINTER redundancy by deg(S);
# Bell witnesses need not be neighbors, so the anti-pointer redundancy is
# not deg-bounded: K_{4,2} (S = 0 and corners 1,2,3 sharing dressers 4,5)
# gives deg(S) = 2 but THREE perfect XX bits, and the corners Bell-record
# each other (the clique).
print("G fan-out (anti-pointer redundancy beyond deg(S)):")
k42 = [(i, k, 1.0) for i in range(4) for k in (4, 5)]
for (a, b) in [(0, 1), (0, 2), (0, 3), (1, 2)]:
    I, corr, td, ev = measure(6, k42, a, b)
    gate(f"K4,2 ({a},{b}): I=1", 1.0, I)
    gate(f"K4,2 ({a},{b}): XX=+1", 1.0, corr["XX"])

# --- G17: two review-found refinements of the pendant/edge scope. (a) The
# role-swap is EXISTENTIAL: one odd-integer watcher on j suffices (it zeroes
# beta_j and both double coherences at once); other watchers, even or
# non-integer, do not spoil it. (b) The fully-shared corner (write bond +
# ONE dresser + nothing else) holds a FORCED full bit at EVERY far-bond
# ratio: the ratio rotates which channel carries it, (YY, ZY) =
# (cos, sin) of pi(1-r)/2, but I = 1 always; the blanket "non-integer ->
# generic" has this one exception.
print("G pendant existential + the forced triangle bit:")
I, corr, td, ev = measure(4, [(0, 1, 1.0), (1, 2, 1.0), (1, 3, 2.0)], 0, 1)
gate("pendant, watchers (1,2): I=1 (one odd suffices)", 1.0, I)
gate("pendant, watchers (1,2): YZ=+1", 1.0, corr["YZ"])
I, corr, td, ev = measure(4, [(0, 1, 1.0), (1, 2, 1.0), (1, 3, 1.5)], 0, 1)
gate("pendant, watchers (1,1.5): I=1", 1.0, I)
for r in (1.5, 2.6):
    ang = np.pi * (1.0 - r) / 2.0
    I, corr, td, ev = measure(3, [(0, 1, 1.0), (0, 2, 1.0), (1, 2, r)], 0, 1)
    gate(f"bare triangle r={r}: I=1 (forced)", 1.0, I)
    gate(f"bare triangle r={r}: YY=cos(pi(1-r)/2)", np.cos(ang), corr["YY"], tol=1e-9)
    gate(f"bare triangle r={r}: ZY=sin(pi(1-r)/2)", np.sin(ang), corr["ZY"], tol=1e-9)

# --- G13: frame reconciliation with the graph-state literature (HEB 2004,
# quant-ph/0307130). U(t*) = prod_bonds e^{-i(pi/4)ZZ} equals prod CZ times
# the degree rotation e^{-i(pi/4)deg_j Z_j} per site (up to global phase), so
# HEB's adjacency-keyed letters (leaf X(x)Z; disconnected matched X(x)X;
# connected matched Y(x)Y) must map to ours under that rotation. Gates: undo
# the degree rotation on the measured pair and recover HEB's letters.
print("G frame reconciliation (HEB letters after undoing degree rotation):")


def undo_degree(pair, dS, dj):
    RS = np.diag([np.exp(1j * np.pi / 4 * dS), np.exp(-1j * np.pi / 4 * dS)])
    Rj = np.diag([np.exp(1j * np.pi / 4 * dj), np.exp(-1j * np.pi / 4 * dj)])
    R = np.kron(RS, Rj)
    return R @ pair @ R.conj().T


def corr_of(pair, A, B):
    return float(np.real(np.trace(pair @ np.kron(PAULI[A], PAULI[B]))))


rho = state(5, [(0, 2, 1.0), (0, 3, 1.0), (0, 4, 1.0),
                (1, 2, 1.0), (1, 3, 1.0), (1, 4, 1.0)], t_star)  # K2,3
pair = undo_degree(ptrace_pair(rho, 5, 0, 1), 3, 3)
gate("K2,3 graph frame: XX=+1 (HEB disconnected)", 1.0, corr_of(pair, "X", "X"))
gate("K2,3 graph frame: YY=0", 0.0, corr_of(pair, "Y", "Y"))
rho = state(3, [(0, 1, 1.0), (1, 2, 1.0), (0, 2, 1.0)], t_star)  # triangle
pair = undo_degree(ptrace_pair(rho, 3, 0, 1), 2, 2)
gate("triangle graph frame: YY=+-1 (HEB connected)", 1.0,
     abs(corr_of(pair, "Y", "Y")))
rho = state(3, [(0, 1, 1.0), (1, 2, 1.0)], t_star)  # chain3, S=1, leaf j=0
pair = undo_degree(ptrace_pair(rho, 3, 1, 0), 2, 1)
gate("leaf graph frame: |ZX|=1 (HEB leaf X(x)Z)", 1.0,
     abs(corr_of(pair, "Z", "X")))

print()
print(f"RESULT: {COUNT[0]} gates, "
      f"{'ALL OK' if not FAILS else 'FAILURES: ' + ', '.join(FAILS)}")
