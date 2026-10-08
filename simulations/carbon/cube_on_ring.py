"""The letter cube laid on the selected carbon ring: which dephasing letter the pairing can see, and the rule.

Selected C_4 spin ring, Hueckel hopping XX + YY on every bond (coefficient 1), gamma = 0.3 on every site, one
dephasing letter P in {Z, X, Y} on every site ("reading along P"). The cube places a term at (k_Z, k_X, k_Y), k_Q
the number of its letters anticommuting with Q. Reading along P, the lit strings (anticommuting with every jump,
F158) are the strings in the two other letters Q, R; the uniform string Q^N commutes with H exactly when NO term of
H is odd in Q (conjugation by Q^N is the sign (-1)^{k_Q} on every term), and a commuting lit string gives the
reflection S L S^-1 = -L^dagger - 2 sigma entry for entry, so the spectrum pairs.

THE RULE, at the level of the term SET: the uniform colouring along P survives iff (no term odd in Q) or (no term
odd in R). It fails iff some term is odd in Q AND some term, the same or another, is odd in R. A single term odd in
both is ONE way to fail (a one-site field along P, (k_Q, k_R) = (1, 1), fails P reading alone; a two-letter bond
term PQ fails the reading along the third letter); two terms each odd in one of them is the other way (a Y field
with transverse DM YZ - ZY fails Z reading, an X field with a Z field fails Y reading), and no term is odd in both
there. Terms with one letter twice (XX, YY, ZZ) are even in every letter and never fail anything, so the bare
Hueckel ring and Hueckel + ZZ keep the uniform colouring under all three letters: on that model the PAIRING cannot
tell which letter the environment reads (the Z spectrum differs from the X and Y spectra, which coincide by the
quarter-turn about Z; the pairing does not differ).

Where the uniform colouring fails the spectrum decides, and it can still pair, in two ways seen here. (i) A rotated
colouring: an X field and a Y field together (common axis n) fail the uniform colouring under Z reading the second
way, but (n . sigma)^N, a lit string of the circle letters, commutes with H and carries the reflection, exactly in
rationals at n = (4/5, 3/5). (ii) A frame: on the open chain the axial DM term XY - YX is carried by the one-site
rotations U = prod_l exp(i l phi Z_l / 2), phi = arctan(D/J), to sqrt(J^2 + D^2)/J Hueckel, and
F = U^dagger X^N U = prod_l (cos(l phi) X_l + sin(l phi) Y_l), a lit element that is not itself a uniform string,
commutes with H and carries the reflection, the colours-on-the-circle colouring of THE_PALINDROME_AS_A_COLOURING
with angle l phi at site l; at D/J = 3/4 every cos(l phi), sin(l phi) is rational and the row is exact. On the ring
the same term is a flux N phi through the cycle, the pairing breaks at D = 0.37 J, and at D = J (flux pi at N = 4,
real hopping with one bond's sign flipped, which X^N keeps) it pairs again (read).

Gates: the set rule on the parities, with the separation from the per-term reading asserted on the two mixed rows;
every commuting uniform string's shift compared to 0.0 as a reflection; every failed cell gated on the spectrum
against the eigensolver's error model (a reflection of any kind would force the pairing); the rotated colouring and
the chain frame's F compared to 0 exactly in rationals (with a wrong angle as control), the frame U H U^dagger at
D = 0.37 read against the model; the centre-line row [H, Z^N] = 0 and the U(1) row [H, sum Z] = 0 gated per model.
"""
import sys
import numpy as np
from itertools import product
from scipy.optimize import linear_sum_assignment
from scipy.linalg import expm
sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2); X = np.array([[0, 1], [1, 0]], dtype=complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1.0, -1.0]).astype(complex)
LET = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
def kron_all(ms):
    r = np.eye(1, dtype=complex)
    for m in ms: r = np.kron(r, m)
    return r
def op(N, l, P): return kron_all([P if j == l else I2 for j in range(N)])
def string(s): return kron_all([LET[c] for c in s])
def liou(N, H, g, P):
    d = 2 ** N; L = -1j * (np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T))
    for l in range(N):
        Q = op(N, l, P); L += g * (np.kron(Q, Q.T) - np.eye(d * d))
    return L
def right_mult(N, F): return np.kron(np.eye(2 ** N), F.T)
def maxabs(A): return float(np.max(np.abs(A)))
def pairing_distance(L, sigma):
    ev = np.linalg.eigvals(L); cost = np.abs(ev[:, None] - (-2 * sigma - ev)[None, :])
    r, c = linear_sum_assignment(cost); return float(cost[r, c].max())
def scale(L): return float(np.max(np.abs(np.linalg.eigvals(L))))
EPS = np.finfo(float).eps
FAILS = []
def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}  {detail}")
    if not ok: FAILS.append(name)

def k_coords(term_letters):
    """(k_Z, k_X, k_Y) of a Pauli string given as its letters."""
    return tuple(sum(1 for c in term_letters if c not in ('I', P)) for P in ('Z', 'X', 'Y'))

N = 4; GAM = 0.3
bonds = [(l, (l + 1) % N) for l in range(N)]
chain_bonds = [(l, l + 1) for l in range(N - 1)]
def two(P, Q, a, b): return op(N, a, P) @ op(N, b, Q)
def hueckel(bs): return sum(two(X, X, a, b) + two(Y, Y, a, b) for a, b in bs)
def field(P, h): return h * sum(op(N, l, P) for l in range(N))
def dm_axial(bs, d): return d * sum(two(X, Y, a, b) - two(Y, X, a, b) for a, b in bs)
def dm_trans(bs, d): return d * sum(two(Y, Z, a, b) - two(Z, Y, a, b) for a, b in bs)
def zz(bs, k): return k * sum(two(Z, Z, a, b) for a, b in bs)
huck = hueckel(bonds)
# each model: name, H, its term letters (for the set rule), readings whose uniform colouring fails
models = [
    ("Hueckel ring XX+YY", huck, ["XX", "YY"], []),
    ("+ ZZ 0.5", huck + zz(bonds, 0.5), ["XX", "YY", "ZZ"], []),
    ("+ Y field 0.4", huck + field(Y, 0.4), ["XX", "YY", "Y"], ['Y']),
    ("+ X field 0.4 (water's tunnelling shape)", huck + field(X, 0.4), ["XX", "YY", "X"], ['X']),
    ("+ Z field 0.4 (water's bias shape)", huck + field(Z, 0.4), ["XX", "YY", "Z"], ['Z']),
    ("+ axial DM 0.37 (XY-YX, the bond current)", huck + dm_axial(bonds, 0.37), ["XX", "YY", "XY"], ['Z']),
    ("+ transverse DM 0.37 (YZ-ZY)", huck + dm_trans(bonds, 0.37), ["XX", "YY", "YZ"], ['X']),
    ("+ Y field 0.4 + transverse DM 0.37 (mixed: no term odd in both X and Y)", huck + field(Y, 0.4) + dm_trans(bonds, 0.37), ["XX", "YY", "Y", "YZ"], ['Z', 'X', 'Y']),
    ("+ X field 0.4 + Z field 0.3 (mixed: no term odd in both X and Z)", huck + field(X, 0.4) + field(Z, 0.3), ["XX", "YY", "X", "Z"], ['X', 'Z', 'Y']),
]
letters = {'Z': Z, 'X': X, 'Y': Y}
others = {'Z': ('X', 'Y'), 'X': ('Y', 'Z'), 'Y': ('X', 'Z')}
idx = {'Z': 0, 'X': 1, 'Y': 2}

print("=== the set rule on the parities: the uniform colouring along P fails iff some term is odd in Q and some term is odd in R ===")
for name, H, terms, failing in models:
    coords = {t: k_coords(t) for t in terms}
    set_rule = [P for P in 'ZXY' if all(any(coords[t][idx[Q]] % 2 == 1 for t in terms) for Q in others[P])]
    per_term = [P for P in 'ZXY' if any(all(coords[t][idx[Q]] % 2 == 1 for Q in others[P]) for t in terms)]
    check(f"{name}: coordinates {coords}; set rule fails {set_rule or 'none'}", set(set_rule) == set(failing))
    if "mixed" in name:
        check(f"{name}: the per-term reading ('a term odd in both') says {per_term or 'none'}, the set rule says {set_rule}: they differ", set(set_rule) != set(per_term))

print("\n=== commutators, reflections, and the spectrum where the uniform colouring fails ===")
table = {}
for name, H, terms, failing in models:
    row = {}
    for P in 'ZXY':
        L = liou(N, H, GAM, letters[P]); sigma = N * GAM
        commuting = [Q for Q in others[P] if maxabs(string(Q * N) @ H - H @ string(Q * N)) == 0.0]
        pd = pairing_distance(L, sigma); model = EPS * scale(L); row[P] = (pd, pd / model, commuting)
        if P not in failing:
            Q = commuting[0]; S = right_mult(N, string(Q * N))
            res = maxabs(S @ L @ np.linalg.inv(S) + L.conj().T + 2 * sigma * np.eye(L.shape[0]))
            check(f"{name}, {P} reading: {Q}^N commutes with H and its shift is a reflection, residual == 0.0", bool(commuting) and res == 0.0, f"commuting {commuting}, residual {res}")
        else:
            check(f"{name}, {P} reading: no uniform lit string commutes with H and the spectrum does not pair", commuting == [] and pd > 1e6 * model, f"{pd:.2e} = {pd / model:.1e} x model")
    table[name] = row

print("\n=== the table (pairing distance about -N gamma; small entries are the eigensolver's floor, ratio to eps x spectral scale in brackets) ===")
print(f"  {'model':74s} {'Z reading':>22s} {'X reading':>22s} {'Y reading':>22s}")
for name, row in table.items():
    print(f"  {name:74s} " + " ".join(f"{row[P][0]:9.2e} ({row[P][1]:8.1e})" for P in 'ZXY'))

print("\n=== where the uniform colouring fails the spectrum can still pair: (i) a rotated colouring on the ring, exact in rationals ===")
import sympy as sp
sX = sp.Matrix([[0, 1], [1, 0]]); sY = sp.Matrix([[0, -sp.I], [sp.I, 0]]); sZ = sp.diag(1, -1); sI = sp.eye(2)
def skron(ms):
    r = sp.Matrix([[1]])
    for m in ms: r = sp.kronecker_product(r, m)
    return r
def sop(l, P): return skron([P if j == l else sI for j in range(N)])
def stwo(P, Q, a, b): return sop(a, P) * sop(b, Q)
def sliou(H, g):
    d = 2 ** N; Id = sp.eye(d); L = -sp.I * (sp.kronecker_product(H, Id) - sp.kronecker_product(Id, H.T))
    for l in range(N):
        Zl = sop(l, sZ); L += g * (sp.kronecker_product(Zl, Zl.T) - sp.eye(d * d))
    return L
def sright(F): return sp.kronecker_product(sp.eye(2 ** N), F.T)
def zero_exact(M):
    """True only when every entry expands to exactly 0 (sympy's is_zero_matrix is None on unexpanded entries)."""
    return M.applyfunc(sp.expand).is_zero_matrix is True
def is_reflection_exact(L, F, sigma):
    S = sright(F); return zero_exact(S * L * S.inv() + L.H + 2 * sigma * sp.eye(L.shape[0]))
R = sp.Rational
sH_huck = sum((stwo(sX, sX, a, b) + stwo(sY, sY, a, b) for a, b in bonds), sp.zeros(16, 16))
hx, hy = R(2, 5), R(3, 10)                     # X field 2/5, Y field 3/10: common axis n = (4/5, 3/5)
sH_xy = sH_huck + sum((hx * sop(l, sX) + hy * sop(l, sY) for l in range(N)), sp.zeros(16, 16))
G = skron([R(4, 5) * sX + R(3, 5) * sY] * N)
check("ring + X field 2/5 + Y field 3/10, Z reading: no uniform lit string commutes (set rule, second way)",
      all(not zero_exact(skron([P] * N) * sH_xy - sH_xy * skron([P] * N)) for P in (sX, sY)))
check("ring + X field + Y field: the circle string G = ((4/5) X + (3/5) Y)^N anticommutes with every Z jump and commutes with H, exactly",
      all(zero_exact(G * sop(l, sZ) + sop(l, sZ) * G) for l in range(N)) and zero_exact(G * sH_xy - sH_xy * G))
sL_xy = sliou(sH_xy, R(3, 10))
check("ring + X field + Y field, Z reading: the shift by G is a reflection, exactly in rationals", is_reflection_exact(sL_xy, G, N * R(3, 10)))
Gw = skron([R(3, 5) * sX + R(4, 5) * sY] * N)
check("control: the circle string at the wrong angle ((3/5) X + (4/5) Y)^N does not commute with H", not zero_exact(Gw * sH_xy - sH_xy * Gw))
Lxy = liou(N, huck + field(X, 0.4) + field(Y, 0.3), GAM, Z); pdxy = pairing_distance(Lxy, N * GAM)
print(f"       ring + X field 0.4 + Y field 0.3, Z reading: pairing distance {pdxy:.2e} = {pdxy / (EPS * scale(Lxy)):.1f} x model (read)")

print("\n=== (ii) the chain frame: axial DM on the open chain, F = prod (cos(l phi) X + sin(l phi) Y), exact at D/J = 3/4 ===")
Dq, Jq = R(3, 4), 1                              # cos phi = 4/5, sin phi = 3/5
sH_chain = sum((stwo(sX, sX, a, b) + stwo(sY, sY, a, b) + Dq * (stwo(sX, sY, a, b) - stwo(sY, sX, a, b)) for a, b in chain_bonds), sp.zeros(16, 16))
cos_l = [sp.cos(l * sp.atan(R(3, 4))) for l in range(N)]; sin_l = [sp.sin(l * sp.atan(R(3, 4))) for l in range(N)]
cos_l = [sp.nsimplify(sp.expand_trig(c)) for c in cos_l]; sin_l = [sp.nsimplify(sp.expand_trig(s_)) for s_ in sin_l]
check("D/J = 3/4: every cos(l phi), sin(l phi) is rational", all(c.is_Rational for c in cos_l) and all(s_.is_Rational for s_ in sin_l), f"{cos_l} {sin_l}")
Fq = skron([cos_l[l] * sX + sin_l[l] * sY for l in range(N)])
check("chain + axial DM 3/4: F anticommutes with every Z jump and commutes with H, exactly; no uniform lit string does",
      all(zero_exact(Fq * sop(l, sZ) + sop(l, sZ) * Fq) for l in range(N)) and zero_exact(Fq * sH_chain - sH_chain * Fq)
      and all(not zero_exact(skron([P] * N) * sH_chain - sH_chain * skron([P] * N)) for P in (sX, sY)))
check("chain + axial DM 3/4, Z reading: the shift by F is a reflection, exactly in rationals", is_reflection_exact(sliou(sH_chain, R(3, 10)), Fq, N * R(3, 10)))
Fw = skron([cos_l[l] * sX - sin_l[l] * sY for l in range(N)])
check("control: the frame at the wrong angle sign does not commute with H", not zero_exact(Fw * sH_chain - sH_chain * Fw))
print(f"       the X^N coefficient of F, prod cos(l phi) = {sp.prod(cos_l)}: F is not a uniform string but contains one")
D = 0.37; J = 1.0; phi = np.arctan(D / J)
Hc = hueckel(chain_bonds) + dm_axial(chain_bonds, D)
U = kron_all([expm(1j * l * phi * Z / 2) for l in range(N)])
res_frame = maxabs(U @ Hc @ U.conj().T - np.sqrt(J * J + D * D) / J * hueckel(chain_bonds)); model_H = EPS * maxabs(Hc) * N
print(f"       D = 0.37 (irrational angle): U H U^dagger - sqrt(J^2 + D^2)/J Hueckel = {res_frame:.1e} = {res_frame / model_H:.1f} x (eps N |H|) (read)")
Lc = liou(N, Hc, GAM, Z); model_L = EPS * scale(Lc); pdc = pairing_distance(Lc, N * GAM)
print(f"       chain + axial DM 0.37, Z reading: pairing distance {pdc:.2e} = {pdc / model_L:.1f} x model (read)")
Lr = liou(N, huck + dm_axial(bonds, 1.0), GAM, Z); pdr = pairing_distance(Lr, N * GAM)
print(f"       ring + axial DM at D = J (flux pi at N = 4), Z reading: pairing distance {pdr:.2e} = {pdr / (EPS * scale(Lr)):.1f} x model (read: the break depends on D through the flux, not on the term's presence)")

print("\n=== the Z-side rows of the table of moves: the centre line needs [H, Z^N] = 0 (even k_Z everywhere); the quarter-turn holds under U(1), [H, sum Z] = 0 ===")
ZN = string('Z' * N); SZ = sum(op(N, l, Z) for l in range(N))
u1_expected = {"Hueckel ring XX+YY": True, "+ ZZ 0.5": True, "+ Y field 0.4": False, "+ X field 0.4 (water's tunnelling shape)": False,
               "+ Z field 0.4 (water's bias shape)": True, "+ axial DM 0.37 (XY-YX, the bond current)": True, "+ transverse DM 0.37 (YZ-ZY)": False}
for name, H, terms, failing in models[:7]:
    cl = maxabs(H @ ZN - ZN @ H) == 0.0; u1 = maxabs(H @ SZ - SZ @ H) == 0.0
    even_kz = all(k_coords(t)[0] % 2 == 0 for t in terms)
    check(f"{name}: centre line [H, Z^N] = 0 is {cl} (= every term even in k_Z: {even_kz}); U(1) [H, sum Z] = 0 is {u1} (expected {u1_expected[name]})", cl == even_kz and u1 == u1_expected[name])

print(f"\n{'ALL OK' if not FAILS else 'FAILED: ' + ', '.join(FAILS)}")
sys.exit(1 if FAILS else 0)
