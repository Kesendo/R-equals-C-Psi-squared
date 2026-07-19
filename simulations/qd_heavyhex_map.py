"""qd_heavyhex_map.py: the record map of the heavy-hex geometry.

Who records whom on IBM's Heron-style lattice? Engine = F135 Proposition 1
DIRECTLY (the pair page as a cosine product; exact at any N, no 2^N state).
Cross-checked against the full-state path on small graphs before use.

Patch: one heavy-hex cell = 12-ring (corners at even positions, edge qubits
at odd) + three dangling bridge qubits on alternating corners 0, 4, 8
(sites 12, 13, 14). N = 15, 105 pairs. gamma = 0, t* = pi/4, base Delta = 1.

Findings (gated below; the dark-lattice corollary is Corollary 8 of
docs/proofs/PROOF_RECORD_LETTER_LAW.md):

PLAY A  The uniform cell is dark except its leaves: 3/105 luminous pairs,
        exactly the three corner-bridge pointer records (ZY = +1). The
        cell's bulk instantiates the corollary: girth >= 5 + leafless
        => every pair exactly dark at uniform coupling, so the infinite
        heavy-hex bulk (girth 12, degrees 2 and 3, no leaves) holds ZERO
        luminous pairs; only the finite patch's dangling bridges record.
PLAY B  Aiming adds without erasing: pushing the two private-watcher bonds
        of corner 0's edge witnesses to ratio 2 (even: forgiven, F135 Law A)
        gives corner 0 a full pointer broadcast R_perfect = deg = 3, and NO
        pair goes dark. But the edit is not confined to the aim: corner 0's
        three witnesses now also record each other (a weave around the aim),
        and each edited bond's far corner (c2, c10) picks up one rotated
        pointer record of its own edge qubit. 10/105 luminous pairs in
        total; the full engineered map is gated below.
PLAY C  Girth >= 5 is sufficient for darkness, not necessary: the uniform
        5x5 square-lattice TORUS (girth 4, leafless) also reads 0 luminous
        pairs -- an interior plaquette's diagonal pair keeps private
        watchers in the surrounding lattice (Q != 0 kills both records).
        The X (x) X plaquette weave belongs to the ISOLATED C4 only.

Run: python simulations/qd_heavyhex_map.py   [~5 s]
Output: simulations/results/qd_witness/qd_heavyhex_map_out.txt
"""

import numpy as np
from itertools import combinations

T = np.pi / 4.0
PAULI = {"I": np.eye(2, dtype=complex),
         "X": np.array([[0, 1], [1, 0]], dtype=complex),
         "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
         "Z": np.array([[1, 0], [0, -1]], dtype=complex)}


def pair_page(bonds, N, a, b, t=T):
    """F135 Prop 1: the exact two-site reduced state, gamma = 0."""
    dab = 0.0
    third = {}                                  # site k -> (Delta_ak, Delta_bk)
    for (x, y, d) in bonds:
        if {x, y} == {a, b}: dab = d
        elif x in (a, b) or y in (a, b):
            k = y if x in (a, b) else x
            s = x if x in (a, b) else y
            da, db = third.get(k, (0.0, 0.0))
            third[k] = (da + d, db) if s == a else (da, db + d)
    rho = np.zeros((4, 4), dtype=complex)
    zs = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
    for i, (za, zb) in enumerate(zs):
        for j, (zA, zB) in enumerate(zs):
            v = 0.25 * np.exp(-1j * t * dab * (za * zb - zA * zB))
            for (dak, dbk) in third.values():
                v *= np.cos(t * (dak * (za - zA) + dbk * (zb - zB)))
            rho[i, j] = v
    return rho


def read(bonds, N, a, b):
    rho = pair_page(bonds, N, a, b)
    rS = np.einsum("ikjk->ij", rho.reshape(2, 2, 2, 2))
    rj = np.einsum("kikj->ij", rho.reshape(2, 2, 2, 2))
    def vn(m):
        w = np.linalg.eigvalsh(m); w = w[w > 1e-12]
        return float(-(w * np.log2(w)).sum())
    I = vn(rS) + vn(rj) - vn(rho)
    corr = {A + B: float(np.real(np.trace(rho @ np.kron(PAULI[A], PAULI[B]))))
            for A in "XYZ" for B in "XYZ"}
    return I, corr


# ---- cross-check the Prop-1 engine against the full-state path (small graphs)
exec(open(__file__.replace("qd_heavyhex_map.py", "qd_letter_gates.py"))
     .read().split('print("record-letter GATES')[0])
print("engine cross-check (Prop 1 vs full state):")
checks = [("triangle", 3, [(0, 1, 1.0), (0, 2, 1.0), (1, 2, 1.0)], 0, 1),
          ("square", 4, [(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 0, 1.0)], 0, 2),
          ("K2,2 (1,3)", 4, [(0, 2, 1.0), (0, 3, 1.0), (1, 2, 1.0), (1, 3, 3.0)], 0, 1),
          ("chain3 pendant", 3, [(0, 1, 1.0), (1, 2, 1.0)], 0, 1)]
ok = True
for name, N, bonds, a, b in checks:
    I1, c1 = read(bonds, N, a, b)
    I2, c2, _, _ = measure(N, bonds, a, b)
    d = max(abs(I1 - I2), max(abs(c1[k] - c2[k]) for k in c1))
    ok &= d < 1e-9
    print(f"  {name}: max dev {d:.2e}")
print(f"  -> {'ENGINE OK' if ok else 'ENGINE MISMATCH'}")
assert ok

# ---- the heavy-hex cell
N = 15
ring = [(i, (i + 1) % 12, 1.0) for i in range(12)]
bridges = [(0, 12, 1.0), (4, 13, 1.0), (8, 14, 1.0)]
hh = ring + bridges
names = {**{i: f"c{i}" if i % 2 == 0 else f"e{i}" for i in range(12)},
         12: "b0", 13: "b4", 14: "b8"}

print()
print("PLAY A: uniform heavy-hex cell -- the full 105-pair map")
lum = []
for a, b in combinations(range(N), 2):
    I, corr = read(hh, N, a, b)
    if I > 0.9999:
        ch = max(("YY", "XX", "ZY", "YZ"), key=lambda c: abs(corr[c]))
        lum.append((names[a], names[b], ch, corr[ch]))
print(f"  luminous pairs: {len(lum)} of {N * (N - 1) // 2}")
for a, b, ch, v in lum:
    print(f"    ({a},{b}): {ch} = {v:+.3f}")
assert sorted((a, b, ch) for a, b, ch, v in lum) == \
    [("c0", "b0", "ZY"), ("c4", "b4", "ZY"), ("c8", "b8", "ZY")], \
    f"luminous set is not exactly the three corner-bridge pointer records: {lum}"
assert all(abs(v - 1.0) < 1e-9 for _, _, _, v in lum)

print()
print("PLAY B: aim the light -- engineer pointer records for corner 0")
# S = corner 0 keeps its three bonds at 1; its edge witnesses e1, e11 have
# private watchers c2, c10 -- push those bonds to ratio 2 (even: forgiven).
eng = [(x, y, (2.0 if {x, y} in ({1, 2}, {10, 11}) else d)) for (x, y, d) in hh]
for j in (1, 11, 12):
    I, corr = read(eng, N, 0, j)
    ch = max(("YY", "XX", "ZY", "YZ"), key=lambda c: abs(corr[c]))
    print(f"  S=c0, j={names[j]}: I={I:.4f}  {ch}={corr[ch]:+.3f}")
    assert abs(I - 1.0) < 1e-9 and ch == "ZY", f"aimed record at j={j} failed: I={I}, {ch}"
print("  the full engineered map -- what the edit really does:")
lum_eng = sorted((a, b) for a, b in combinations(range(N), 2)
                 if read(eng, N, a, b)[0] > 0.9999)
lum_uni = sorted((a, b) for a, b in combinations(range(N), 2)
                 if read(hh, N, a, b)[0] > 0.9999)
newly = sorted(set(lum_eng) - set(lum_uni))
gone = sorted(set(lum_uni) - set(lum_eng))
print(f"    engineered luminous: {len(lum_eng)}/105: "
      f"{[(names[a], names[b]) for a, b in lum_eng]}")
print(f"    newly lit: {[(names[a], names[b]) for a, b in newly]}")
print(f"    went dark: {gone}")
assert gone == [], f"aiming erased existing records: {gone}"
assert lum_eng == sorted([(0, 1), (0, 11), (0, 12), (1, 11), (1, 12), (2, 3),
                          (4, 13), (8, 14), (9, 10), (11, 12)]), \
    f"engineered luminous set moved: {lum_eng}"

print()
print("PLAY C: the square-lattice bulk is dark too (girth >= 5 sufficient, not necessary)")
L = 5
def _idx(i, j): return (i % L) * L + (j % L)
torus = []
for i in range(L):
    for j in range(L):
        torus.append((_idx(i, j), _idx(i + 1, j), 1.0))
        torus.append((_idx(i, j), _idx(i, j + 1), 1.0))
lum_torus = [(a, b) for a, b in combinations(range(L * L), 2)
             if read(torus, L * L, a, b)[0] > 0.9999]
I_diag = read(torus, L * L, _idx(0, 0), _idx(1, 1))[0]
I_c4 = read([(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 0, 1.0)], 4, 0, 2)[0]
print(f"  5x5 torus (uniform): luminous pairs {len(lum_torus)}/300")
print(f"  plaquette diagonal in the lattice: I={I_diag:.4f}   isolated C4 diagonal: I={I_c4:.4f}")
assert lum_torus == [] and I_diag < 1e-9 and abs(I_c4 - 1.0) < 1e-9

print()
print("heavy-hex GATES all OK (engine, PLAY A dark map, PLAY B aiming: "
      "broadcast + side-light, nothing erased; PLAY C square-lattice bulk dark)")
