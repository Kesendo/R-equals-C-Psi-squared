"""Gate for PROOF_F108_KLEIN_V4_EQUIVALENCE: which move carries which F108 mirror, and which carries the Lindbladian.

Everything is exact: per-site 4x4 maps on the Pauli basis (I, X, Y, Z) with entries in {0, +-1, +-i}, and N = 3
Liouvillians with integer couplings, integer rates and the unnormalised Cliffords X + Z (Hadamard) and I + iX (quarter
turn about X), whose inverses are (X + Z)/2 and (I - iX)/2, so every comparison is ==.

T1  Per site: M_Y = -M_Z^-1, and the Klein-V4 permutes the four oriented mirrors {M_Z, M_X, M_Y, -M_X^-1} simply
    transitively: D (Y -> -Y) as (M_Z, M_Y)(M_X, -M_X^-1), H (X <-> Z, Y fixed) as (M_Z, M_X)(M_Y, -M_X^-1), Q_zx (X <-> Z,
    Y -> -Y, the Hadamard's action) as (M_X, M_Y)(M_Z, -M_X^-1); D.H = Q_zx. By tensor powers these hold at every N.
    The label maps D, H, Q_zx are the same operators as the transpose, the Hadamard conjugation and their product
    used in T2 (checked at N = 3).
T2  On Liouvillians (N = 3, a Hamiltonian with every two-site string on some bond and two fields, integer couplings):
    Q_zx carries L_Z(H) to L_X(U H U^-1); H carries L_Z(H) to L_X(-U H^T U^-1); D carries L_Z(H) to L_Z(-H^T), the same
    letter; the quarter turn about X carries L_Z(H) to L_Y(R H R^-1). All four images are Lindbladians.
T3  The mirrors at N = 3, on a Hamiltonian of each Part's bilinear set with integer couplings and rates 1, 2, 3:
    Pi_5b(Z) palindromizes Part 1, Pi_5b(X) Part 2, Pi_5b(Y) Part 3, and the wrong letter's mirror fails on each (the
    controls); the Hadamard transport Q_zx Pi_5b(Z) Q_zx^-1 is -Pi_5b(X)^-1 at N = 3 and palindromizes Part 2; the
    quarter turn fixes Pi_5b(Z), and Pi_5b(Z) palindromizes Part 3 too.
Run: python simulations/f108_klein_transport_gate.py"""
import itertools, sys
import numpy as np

if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass
FAILS = []
def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not ok: FAILS.append(name)

L4 = "IXYZ"
def site_map(images):  # images: letter -> (coefficient, letter); column = input, row = output
    M = np.zeros((4, 4), complex)
    for a, (c, b) in images.items(): M[L4.index(b), L4.index(a)] = c
    return M
PI5B = {"Z": site_map({"I": (1, "X"), "X": (-1, "I"), "Y": (1j, "Z"), "Z": (-1j, "Y")}),
        "Y": site_map({"I": (1, "X"), "X": (-1, "I"), "Y": (-1j, "Z"), "Z": (1j, "Y")}),
        "X": site_map({"I": (1, "Z"), "Z": (-1, "I"), "X": (-1j, "Y"), "Y": (1j, "X")})}
D1 = site_map({"I": (1, "I"), "X": (1, "X"), "Y": (-1, "Y"), "Z": (1, "Z")})
H1 = site_map({"I": (1, "I"), "X": (1, "Z"), "Y": (1, "Y"), "Z": (1, "X")})
Q1 = site_map({"I": (1, "I"), "X": (1, "Z"), "Y": (-1, "Y"), "Z": (1, "X")})
conj = lambda K, M: K @ M @ K                                  # all three are involutions
inv1 = lambda M: M.conj().T                                    # the per-site maps are unitary
check("T1 M_Y = -M_Z^-1 per site", np.array_equal(PI5B["Y"], -inv1(PI5B["Z"])))
FOUR = {"M_Z": PI5B["Z"], "M_X": PI5B["X"], "M_Y": PI5B["Y"], "-M_X^-1": -inv1(PI5B["X"])}
WANT = {"D": {"M_Z": "M_Y", "M_Y": "M_Z", "M_X": "-M_X^-1", "-M_X^-1": "M_X"},
        "H": {"M_Z": "M_X", "M_X": "M_Z", "M_Y": "-M_X^-1", "-M_X^-1": "M_Y"},
        "Q_zx": {"M_X": "M_Y", "M_Y": "M_X", "M_Z": "-M_X^-1", "-M_X^-1": "M_Z"}}
for name, K in (("D", D1), ("H", H1), ("Q_zx", Q1)):
    got = {a: [b for b, B in FOUR.items() if np.array_equal(conj(K, A), B)] for a, A in FOUR.items()}
    check(f"T1 {name} permutes the four oriented mirrors as {WANT[name]}", all(got[a] == [WANT[name][a]] for a in FOUR), f"{got}")
check("T1 D.H = Q_zx per site, so any two of the three generate the third", np.array_equal(D1 @ H1, Q1) and np.array_equal(H1 @ D1, Q1))
check("T1 the per-site mirrors are unitary and D, H, Q_zx are involutions (so inverse = adjoint, conjugation = K M K)",
      all(np.array_equal(M @ M.conj().T, np.eye(4)) for M in PI5B.values()) and all(np.array_equal(K @ K, np.eye(4)) for K in (D1, H1, Q1)))
def kpow(M, n):
    out = np.array([[1]], complex)
    for _ in range(n): out = np.kron(out, M)
    return out
for n_ in (2, 3, 4):
    sgn = (-1) ** n_
    ok_y = np.array_equal(kpow(PI5B["Y"], n_), sgn * kpow(PI5B["Z"], n_).conj().T)
    ok_t = np.array_equal(kpow(Q1, n_) @ kpow(PI5B["Z"], n_) @ kpow(Q1, n_), sgn * kpow(PI5B["X"], n_).conj().T)
    check(f"T1 N={n_} on the tensor powers: Pi_5b(Y) = (-1)^N Pi_5b(Z)^-1 and Q_zx Pi_5b(Z) Q_zx = (-1)^N Pi_5b(X)^-1 (sign {sgn:+d})",
          ok_y and ok_t)

# ---------------------------------------------------------------- T2, T3 on N = 3 Liouvillians
P = {"I": np.eye(2, dtype=complex), "X": np.array([[0, 1], [1, 0]], complex), "Y": np.array([[0, -1j], [1j, 0]]),
     "Z": np.diag([1, -1]).astype(complex)}
N = 3; Dm = 2 ** N; I = np.eye(Dm)
def op(d):
    m = np.array([[1]], complex)
    for i in range(N): m = np.kron(m, d.get(i, np.eye(2)))
    return m
def lr(A, B): return np.kron(A, B.T)                           # rho -> A rho B, row-major vec
def L(H, letter, rates=(1, 2, 3)):
    return -1j * (lr(H, I) - lr(I, H)) + sum(rates[l] * (lr(op({l: P[letter]}), op({l: P[letter]})) - np.eye(Dm * Dm)) for l in range(N))
T = np.zeros((Dm * Dm, Dm * Dm))
for a in range(Dm):
    for b in range(Dm): T[b * Dm + a, a * Dm + b] = 1               # rho -> rho^T, the operator D
rng = np.random.default_rng(7)
Hgen = sum(int(rng.integers(-3, 4)) * op({b: P[x], b + 1: P[y]}) for b in range(N - 1) for x, y in itertools.product("XYZ", repeat=2))
Hgen = Hgen + 2 * op({0: P["X"]}) + 3 * op({2: P["Y"]})
U = op({l: P["X"] + P["Z"] for l in range(N)}); Ui = op({l: (P["X"] + P["Z"]) / 2 for l in range(N)})
R = op({l: np.eye(2) + 1j * P["X"] for l in range(N)}); Ri = op({l: (np.eye(2) - 1j * P["X"]) / 2 for l in range(N)})
Qzx, Qzx_i = lr(U, Ui), lr(Ui, U)
Hop, Hop_i = Qzx @ T, T @ Qzx_i
Rop, Rop_i = lr(R, Ri), lr(Ri, R)
check("T2 Q_zx carries L_Z(H) to L_X(U H U^-1)", np.array_equal(Qzx @ L(Hgen, "Z") @ Qzx_i, L(U @ Hgen @ Ui, "X")))
check("T2 H = Q_zx.D carries L_Z(H) to L_X(-U H^T U^-1)", np.array_equal(Hop @ L(Hgen, "Z") @ Hop_i, L(U @ (-Hgen.T) @ Ui, "X")))
check("T2 D carries L_Z(H) to L_Z(-H^T): a Lindbladian of the same letter", np.array_equal(T @ L(Hgen, "Z") @ T, L(-Hgen.T, "Z")))
r1, r1i = np.eye(2) + 1j * P["X"], (np.eye(2) - 1j * P["X"]) / 2
check("T2 the quarter turn's per-site letter map is X -> X, Y -> -Z, Z -> Y",
      np.array_equal(r1 @ P["X"] @ r1i, P["X"]) and np.array_equal(r1 @ P["Y"] @ r1i, -P["Z"]) and np.array_equal(r1 @ P["Z"] @ r1i, P["Y"]))
check("T2 the quarter turn about X carries L_Z(H) to L_Y(R H R^-1)", np.array_equal(Rop @ L(Hgen, "Z") @ Rop_i, L(R @ Hgen @ Ri, "Y")))

def pauli_superop(site):                                       # N-site Pauli-basis map -> superoperator on row-major vec
    S = np.zeros((Dm * Dm, Dm * Dm), complex)
    for word in itertools.product(L4, repeat=N):
        src = op({i: P[w] for i, w in enumerate(word)})
        img = np.array([[1]], complex)
        for w in word:
            col = site[:, L4.index(w)]; k = int(np.nonzero(col)[0][0])
            img = np.kron(img, col[k] * P[L4[k]])
        S += np.outer(img.reshape(-1), src.conj().reshape(-1)) / Dm
    return S
PI = {d: pauli_superop(PI5B[d]) for d in "XYZ"}
check("T1-T2 the label maps are the T2 superoperators at N = 3: D = transpose, Q_zx = Hadamard conjugation, H = their product",
      np.array_equal(pauli_superop(D1), T) and np.array_equal(pauli_superop(Q1), Qzx) and np.array_equal(pauli_superop(H1), Hop))
PARTS = {"Part 1": ("Z", ["XX", "YY", "YZ", "ZY", "ZZ"]), "Part 2": ("X", ["ZZ", "XX", "XY", "YX", "YY"]),
         "Part 3": ("Y", ["XX", "YY", "YZ", "ZY", "ZZ"])}
def part_H(strings):
    return sum((k % 3 + 1) * (1 if k % 2 == 0 else -1) * op({b: P[s[0]], b + 1: P[s[1]]})
               for b in range(N - 1) for k, s in enumerate(strings))
sigma = 6
def palin(S, Lm):                                             # the mirrors are unitary: S^-1 = S^dagger, checked exactly
    Si = S.conj().T
    return np.array_equal(S @ Si, np.eye(Dm * Dm)) and np.array_equal(S @ Lm @ Si, -Lm - 2 * sigma * np.eye(Dm * Dm))
WRONG = {"Part 1": "X", "Part 2": "Z", "Part 3": "X"}
for part, (d, strings) in PARTS.items():
    Lp = L(part_H(strings), d)
    check(f"T3 Pi_5b({d}) palindromizes {part} (jumps {d}) exactly; the control Pi_5b({WRONG[part]}) does not",
          palin(PI[d], Lp) and not palin(PI[WRONG[part]], Lp))
Pt = Qzx @ PI["Z"] @ Qzx_i
check("T3 the Hadamard transport Q_zx Pi_5b(Z) Q_zx^-1 is -Pi_5b(X)^-1 at N = 3 and palindromizes Part 2",
      np.array_equal(Pt, -PI["X"].conj().T) and palin(Pt, L(part_H(PARTS["Part 2"][1]), "X")))
check("T3 the quarter turn about X fixes Pi_5b(Z), which palindromizes Part 3 as well",
      np.array_equal(Rop @ PI["Z"] @ Rop_i, PI["Z"]) and palin(PI["Z"], L(part_H(PARTS["Part 3"][1]), "Y")))
print("\nALL PASS" if not FAILS else f"\n{len(FAILS)} FAIL")
sys.exit(1 if FAILS else 0)
