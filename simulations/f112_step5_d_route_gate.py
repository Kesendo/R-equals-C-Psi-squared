"""Gate: F114's transpose gives a second route to F112 Step 5.

F112 Step 5 (PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE) needs ||L_{H,+i}||^2 = ||L_{H,-i}||^2 for the commutator superoperator
L_H = -i[H, .], with +-i the eigenvalues of conjugation by Pi = Pi_Z. The parent proof uses the dagger (Lemma A) and
L_H^dagger = -L_H for Hermitian H (Lemma B); the non-Hermitian extension uses Lemma N-A (||L_{sigma,-i}||^2 = 4^N for bit_b-odd sigma) and
Lemma N-B (distinct bit_b-odd strings have orthogonal -i parts). The transpose route: D = diag((-1)^{n_Y}) inverts Pi
(F118: D Pi D = Pi^-1, the dihedral relation), so conjugation by D is a unitary map sending the Pi-eigenvalue lambda to
conj(lambda), the +i part of A onto the -i part of D A D; and D L_sigma D = eps(sigma) L_sigma (F114). Hence
  ||L_{H,+i}|| = ||(D L_H D)_{-i}|| = ||L_{H,-i}||  whenever the bit_b-odd strings of H share one n_Y parity (bit_b-even
strings have no +-i part), any coefficients; and for any H, since the norms differ by -4 Re<L_{H_e,-i}, L_{H_o,-i}>
(H_e, H_o the parts of even and odd n_Y), a sum over distinct strings that N-B kills.
Everything is exact: Pi and D are signed permutations with phases in {+-1, +-i}, L_sigma has Gaussian-integer entries in
the Pauli basis, coefficients are Gaussian integers, and 4 A_lambda = sum_m lambda^-m Pi^m A Pi^-m, so every compared
quantity is a Gaussian integer held exactly in complex128 (entries far below 2^53), compared with ==.
G1  D Pi D = Pi^-1 at N = 1, 2, 3, with Pi built from F118's rho^T X^N and checked against build_pi_full; control:
    diag((-1)^{n_X}) does not invert Pi at odd N.
G2  D L_sigma D = eps(sigma) L_sigma, eps = (-1)^{n_Y + 1}, for every string at N = 2 and 3.
G3  Non-Hermitian H whose bit_b-odd strings share one n_Y parity, plus bit_b-even strings of the other parity
    (Gaussian-integer coefficients), N = 2 and 3: D L_H D is not +-L_H, yet (D L_H D)_-i = eps (L_H)_-i and
    ||L_{H,+i}||^2 = ||L_{H,-i}||^2; control: with bit_b-odd strings of both parities (D L_H D)_-i is not +-(L_H)_-i.
G4  N-B: <L_{sigma,lambda}, L_{tau,lambda}> = 0 for distinct strings, lambda = +-i, N = 2 and 3 (sampled at 3), and
    ||L_{H,+i}||^2 = ||L_{H,-i}||^2 for H with bit_b-odd strings of both n_Y parities (a confirmation of the theorem).
G5  The cross-term identity on generic superoperators: for integer A with D A D = -A and B with D B D = B,
    ||(A+B)_{+i}||^2 - ||(A+B)_{-i}||^2 = -4 Re<A_{-i}, B_{-i}>, with a nonzero right side.
Output: simulations/results/f112_step5_d_route_gate.txt. Run from the repo root: python simulations/f112_step5_d_route_gate.py"""
import itertools, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from framework.symmetry import build_pi_full

OUT = []
FAILS = []
def say(s): print(s); OUT.append(s)
def check(name, ok, detail=""):
    say(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not ok: FAILS.append(name)

P1 = {"I": np.eye(2), "X": np.array([[0, 1], [1, 0]]), "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.diag([1, -1])}
def string(s):
    out = np.array([[1]], dtype=complex)
    for c in s: out = np.kron(out, P1[c])
    return out
def labels(N): return ["".join(t) for t in itertools.product("IXYZ", repeat=N)]
def setup(N):
    lab = labels(N)
    S = [string(s) for s in lab]
    # Pi from F118's rule Pi(rho) = rho^T X^N in this file's own labels, checked against the framework's Pi
    F = string("X" * N); d = 2 ** N
    Pi = np.array([[np.trace(S[a].conj().T @ (S[b].T @ F)) / d for b in range(len(S))] for a in range(len(S))])
    assert np.array_equal(Pi, build_pi_full(N, "Z")), "F118's rho^T X^N differs from build_pi_full in this labelling"
    Pinv = Pi.conj().T
    assert np.array_equal(Pi @ Pinv, np.eye(4 ** N)), "Pi not unitary"
    D = np.diag([(-1) ** s.count("Y") for s in lab]).astype(complex)
    return lab, Pi, Pinv, D, S
def L_of(H, S, N):
    """Pauli-basis matrix of -i[H, .]: entry (a, b) = Tr(S_a^dagger (-i[H, S_b])) / 2^N (Gaussian integers here)"""
    d = 2 ** N
    M = np.array([[np.trace(S[a].conj().T @ (-1j * (H @ S[b] - S[b] @ H))) / d for b in range(len(S))] for a in range(len(S))])
    R = np.round(M.real) + 1j * np.round(M.imag)
    assert np.array_equal(R, M), "L entries are not Gaussian integers"
    return R
def proj4(A, lam, Pi, Pinv):
    """4 A_lambda = sum_m lambda^-m Pi^m A Pi^-m (exact for Gaussian-integer A)"""
    out = np.zeros_like(A); Pm = np.eye(len(A), dtype=complex); Pmi = Pm.copy()
    for m in range(4):
        out = out + (lam ** (-m)) * (Pm @ A @ Pmi)
        Pm, Pmi = Pi @ Pm, Pmi @ Pinv
    return out
def n2(A):
    v = np.sum(A.real ** 2 + A.imag ** 2)
    assert v == int(v), "a norm that should be an integer is not"
    return int(v)
def ip(A, B): return complex(np.vdot(A, B))

for N in (1, 2, 3):
    lab, Pi, Pinv, D, S = setup(N)
    check(f"G1 N={N}: D Pi D = Pi^-1 exactly", np.array_equal(D @ Pi @ D, Pinv))
    if N in (1, 3):
        DX = np.diag([(-1) ** s.count("X") for s in lab]).astype(complex)
        check(f"G1 control N={N}: the sign diag((-1)^n_X) does not invert Pi at odd N", not np.array_equal(DX @ Pi @ DX, Pinv))

rng = np.random.default_rng(20261007)
for N in (2, 3):
    lab, Pi, Pinv, D, S = setup(N)
    Ls = {s: L_of(S[i], S, N) for i, s in enumerate(lab) if s != "I" * N}
    ok2 = all(np.array_equal(D @ L @ D, (-1) ** (s.count("Y") + 1) * L) for s, L in Ls.items())
    check(f"G2 N={N}: D L_sigma D = (-1)^(n_Y+1) L_sigma for all {len(Ls)} strings", ok2)
    def gauss(): return complex(int(rng.integers(-3, 4)), int(rng.integers(-3, 4)))
    def bitb_odd(s): return (s.count("Y") + s.count("Z")) % 2 == 1
    ok3 = True; det = []
    for par in (0, 1):
        odd_pool = [s for s in Ls if s.count("Y") % 2 == par and bitb_odd(s)]
        even_other = [s for s in Ls if s.count("Y") % 2 != par and not bitb_odd(s)]
        for _ in range(3):
            terms = list(rng.choice(odd_pool, 3, replace=False)) + list(rng.choice(even_other, 2, replace=False))
            L = sum(gauss() * Ls[t] for t in terms)
            eps = (-1) ** (par + 1)
            a, b = n2(proj4(L, 1j, Pi, Pinv)), n2(proj4(L, -1j, Pi, Pinv))
            ok3 &= (not np.array_equal(D @ L @ D, eps * L)) and a == b and a > 0 \
                and np.array_equal(D @ proj4(L, 1j, Pi, Pinv) @ D, eps * proj4(L, -1j, Pi, Pinv))
            det.append(a)
    check(f"G3 N={N}: bit_b-odd strings of one n_Y parity plus bit_b-even strings of the other: D L_H D != eps L_H, yet "
          "D (L_H)_+i D = eps (L_H)_-i and ||L_+i||^2 = ||L_-i||^2 > 0", ok3, f"16 ||L_+i||^2 = {det[:3]}...")
    oe = [s for s in Ls if s.count("Y") % 2 == 0 and bitb_odd(s)]; oo = [s for s in Ls if s.count("Y") % 2 == 1 and bitb_odd(s)]
    Lmix = Ls[oe[0]] + 1j * Ls[oo[0]]
    Pm, Pp = proj4(Lmix, -1j, Pi, Pinv), proj4(Lmix, 1j, Pi, Pinv)
    check(f"G3 N={N} control: with bit_b-odd strings of both n_Y parities ({oe[0]}, {oo[0]}), D (L_H)_+i D is neither +(L_H)_-i nor -(L_H)_-i",
          not np.array_equal(D @ Pp @ D, Pm) and not np.array_equal(D @ Pp @ D, -Pm))
    keys = list(Ls)
    if N == 3: keys = list(rng.choice(keys, 24, replace=False))
    P = {s: (proj4(Ls[s], 1j, Pi, Pinv), proj4(Ls[s], -1j, Pi, Pinv)) for s in keys}
    ok4 = all(ip(P[s][0], P[t][0]) == 0 and ip(P[s][1], P[t][1]) == 0 for s, t in itertools.combinations(keys, 2))
    check(f"G4 N={N}: <L_sigma,lambda, L_tau,lambda> = 0 for distinct strings, lambda = +-i ({len(keys)} strings{', sampled' if N == 3 else ''})", ok4)
    ok5 = True
    for _ in range(4):
        terms = list(rng.choice(oe, 2, replace=False)) + list(rng.choice(oo, 2, replace=False))
        L = sum(gauss() * Ls[t] for t in terms)
        ok5 &= n2(proj4(L, 1j, Pi, Pinv)) == n2(proj4(L, -1j, Pi, Pinv))
    check(f"G4 N={N}: bit_b-odd strings of both n_Y parities, non-Hermitian: ||L_+i||^2 = ||L_-i||^2 (a confirmation of the theorem)", ok5)
    dim = 4 ** N
    A0 = rng.integers(-3, 4, (dim, dim)) + 1j * rng.integers(-3, 4, (dim, dim))
    B0 = rng.integers(-3, 4, (dim, dim)) + 1j * rng.integers(-3, 4, (dim, dim))
    A, B = A0 - D @ A0 @ D, B0 + D @ B0 @ D            # D A D = -A, D B D = B
    lhs = n2(proj4(A + B, 1j, Pi, Pinv)) - n2(proj4(A + B, -1j, Pi, Pinv))
    rhs = -4 * ip(proj4(A, -1j, Pi, Pinv), proj4(B, -1j, Pi, Pinv)).real
    check(f"G5 N={N}: ||(A+B)_+i||^2 - ||(A+B)_-i||^2 = -4 Re<A_-i, B_-i> on generic superoperators, right side nonzero",
          lhs == rhs and rhs != 0, f"{lhs} = {rhs}")

say("\nALL PASS" if not FAILS else f"\n{len(FAILS)} FAIL")
os.makedirs("simulations/results", exist_ok=True)
open("simulations/results/f112_step5_d_route_gate.txt", "w", encoding="utf-8", newline="\n").write("\n".join(OUT) + "\n")
sys.exit(1 if FAILS else 0)
