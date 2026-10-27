"""The spin-3/2 nucleus as a key: the sixteen string rates as a linear code of the three spectral densities,
read back, and the densities read back to the motion (experiments/THE_NUCLEUS_AS_A_KEY.md).

At zero splitting and offset the sixteen strings of experiments/THE_SODIUM_NUCLEUS_TWO_MIRRORS.md are the
dissipator's eigenoperators, each decaying at twice the summed rate of the jumps ZZ, ZX, ZY, XI, YI (rates
3J0/2, 3J1/2, 3J1/2, 3J2/2, 3J2/2) it anticommutes with. This script
  1. builds the 16 x 3 integer generator G (rate = G @ (J0, J1, J2)), checks it against the page's rate table,
     its rank and how its 13 linear relations split into trivial and genuine ones (exact, integers);
  2. checks that with a splitting each bright pair's two rates still sum to the two cells' rates, all four pairs
     to the one number 6(J0 + J1 + J2), and reads the densities back from the dark cells and that number;
  3. reads the densities back from noisy rates, three readings against the nine decays a measurement can separate
     (each weighted by its size, the textbook weighting for relative noise), averaged over many draws (a
     statistic, printed, not gated as a number);
  4. under one Lorentzian density J(w) = A*tau/(1 + (w tau)^2) checks the single-motion relation
     (J0/J2 - 1) = 4 (J0/J1 - 1) exactly in rationals, reads w0*tau and A*tau back, and shows that at w0*tau << 1
     only the product A*tau survives.
"""
from fractions import Fraction as F
import numpy as np

def anti(a, b):
    return sum(1 for x, y in zip(a, b) if x != "I" and y != "I" and x != y) % 2 == 1

W = [a + b for a in "IXYZ" for b in "IXYZ"]
JUMPS, WHICH = ["ZZ", "ZX", "ZY", "XI", "YI"], [0, 1, 1, 2, 2]
G = np.array([[sum(3 for j, k in zip(JUMPS, WHICH) if anti(w, j) and k == c) for c in range(3)] for w in W])

# 1. the code. The page's table, times (J0, J1, J2) coefficients, in Souza's labelling
PAGE = {"II": (0, 0, 0), "ZZ": (0, 6, 6), "IZ": (0, 6, 0), "ZI": (0, 0, 6),
        "XX": (0, 3, 3), "XY": (0, 3, 3), "YX": (0, 3, 3), "YY": (0, 3, 3),
        "ZX": (3, 3, 6), "ZY": (3, 3, 6), "IX": (3, 3, 0), "IY": (3, 3, 0),
        "XI": (3, 6, 3), "YI": (3, 6, 3), "XZ": (3, 0, 3), "YZ": (3, 0, 3)}
for w, row in zip(W, G):
    assert tuple(int(x) for x in row) == PAGE[w], (w, row, PAGE[w])
assert np.linalg.matrix_rank(G) == 3
rows = [tuple(int(x) for x in r) for r in G]
distinct = sorted(set(rows))
nonzero = [r for r in distinct if any(r)]
zero_cells = sum(1 for r in rows if not any(r))
duplicates = len(rows) - len(distinct)
genuine = len(nonzero) - 3
assert zero_cells + duplicates + genuine == 16 - 3 == 13
print(f"1. the code: 16 cells, rank 3, 13 relations = {zero_cells} cell fixed at 0 (II) + {duplicates} cells repeating "
      f"another cell's rate + {genuine} genuine relations among the {len(nonzero)} distinct non-zero rates")
assert (zero_cells, duplicates, genuine) == (1, 7, 5)
idx = {w: i for i, w in enumerate(W)}
assert all(G[idx["ZZ"]] == G[idx["IZ"]] + G[idx["ZI"]])          # one genuine relation, by name

# 2. with a splitting: the bright pairs and their sums
PAIRS = [("IX", "ZX"), ("IY", "ZY"), ("XI", "XZ"), ("YI", "YZ")]
for a_, b_ in PAIRS:
    assert tuple(G[idx[a_]] + G[idx[b_]]) == (6, 6, 6), (a_, b_)
J = np.array([1.0, 0.6, 0.3]); key = G @ J
pair_sum = key[idx["IX"]] + key[idx["ZX"]]
J1s, J2s = key[idx["IZ"]] / 6, key[idx["ZI"]] / 6
assert np.allclose([pair_sum / 6 - J1s - J2s, J1s, J2s], J, rtol=0, atol=1e-12)
print("2. with a splitting every bright pair sums to 6(J0 + J1 + J2), one number for all four pairs; the dark cells "
      "IZ, ZI give J1, J2 and the pair sum gives J0")

# 3. reading back with noise: the nine decays a measurement separates (a coherence and its conjugate are one decay;
# the central transition and the triple quantum share a rate but are measured apart)
READ = {"ZZ": G[idx["ZZ"]], "IZ": G[idx["IZ"]], "ZI": G[idx["ZI"]], "CT": G[idx["XX"]], "TQ": G[idx["XX"]],
        "sat1": G[idx["ZX"]], "sat2": G[idx["IX"]], "dq1": G[idx["XI"]], "dq2": G[idx["XZ"]]}
R = np.array(list(READ.values()), dtype=float)
assert np.linalg.matrix_rank(R) == 3
def three(r):
    J1, J2 = r[1] / 6, r[2] / 6
    return np.array([r[6] / 3 - J1, J1, J2])
rng = np.random.default_rng(20261027)
true = R @ J
e3, e9 = [], []
for _ in range(20000):
    noisy = true * (1 + 0.02 * rng.standard_normal(len(true)))
    e3.append(three(noisy) / J - 1)
    wgt = 1 / noisy
    e9.append(np.linalg.lstsq(R * wgt[:, None], noisy * wgt, rcond=None)[0] / J - 1)
rms = lambda e: np.sqrt(np.mean(np.square(e), axis=0))
r3, r9 = rms(e3), rms(e9)
print("3. 2% relative noise on each of nine decays, J = (1, 0.6, 0.3), 20000 draws, rms relative error of (J0, J1, J2):")
print("   three readings:", np.round(r3, 4), "  nine readings, weighted by size:", np.round(r9, 4))
assert all(r9 < r3)

# 4. back to the motion, exactly in rationals
def lorentz(A, tau, w0):
    return [A * tau / (1 + (n * w0 * tau) ** 2) for n in (0, 1, 2)]

for a in [F(1, 100), F(1, 2), F(3), F(17, 5)]:          # a = (w0 tau)^2
    J0, J1, J2 = [F(1) / (1 + n * n * a) for n in (0, 1, 2)]   # A = tau = 1, w0 tau = sqrt(a)
    assert J0 / J2 - 1 == 4 * (J0 / J1 - 1)               # the single-motion relation, exact
    assert J0 / J1 - 1 == a                                 # w0*tau read back from one ratio
J0, J1, J2 = [x + y for x, y in zip([F(1) / (1 + n * n * F(1, 4)) for n in (0, 1, 2)], [F(1, 3) / (1 + n * n * 9) for n in (0, 1, 2)])]
assert J0 / J2 - 1 != 4 * (J0 / J1 - 1)
print("4. single-motion relation exact at four values of (w0 tau)^2; fails for a sum of two motions, as it should")

w0, A = 2 * np.pi * 105.8e6, 1.0e10                       # 23Na at 9.4 T
for name, tau, w in [("w0 tau ~ 1", 2e-9, w0), ("w0 tau << 1 (water, a few ps)", 5e-12, w0)]:
    Jt = np.array(lorentz(A, tau, w))
    a = Jt[0] / Jt[1] - 1
    tau_back = np.sqrt(a) / w
    print(f"   {name}: tau back {tau_back:.3e} (true {tau:.1e}), A*tau back {Jt[0]:.4e} (true {A * tau:.4e})")
    print(f"      (w0 tau)^2 = {(w * tau) ** 2:.2e}")
print("   at w0 tau << 1 the ratio J0/J1 - 1 is far below any measured precision: only A*tau can be read")
