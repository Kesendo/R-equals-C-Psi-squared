"""Which soft mirrors are Π in a turned frame: the gate behind the two-term Klein routing's frame section.

A one-site unitary that keeps Z-dephasing (a turn about Z, or X times one) leaves I in place,
sends Z to ±Z and mixes X with Y. A product of them keeps the number of Z letters in every
Pauli string; F85's truly rule asks every string for an even number of Y AND of Z. So:

  G1  the three alternating soft pairs (XY+YX, XY+ZZ, YX+ZZ) carry no string with an odd
      number of Z, and with every odd site turned by (X+Y)/√2, or by the quarter turn S, they
      are truly: checked on the letters (exact, N = 3, 4, 5) and on the dense residual
      (== 0.0 exactly at γ = 0.5, N = 3, 4). XY+YX turned by (X+Y)/√2 is XX+YY itself.
  G2  the other sixteen soft pairs each carry an odd-Z string; F85's residual is additive over
      strings with weight 4·4^N per unit coefficient for a Π²-odd string and 8·4^N for a
      Π²-even non-truly one, so ‖M‖² ≥ 4·4^N·c² with c² the summed squared coefficients of the
      odd-Z strings, which the turn keeps. The additivity holds because Π acts on the letters as a
      translation: each string P of H sends a basis string S only to P·S, so different strings of
      H fill different entries of M; checked: the weights and the additivity at N = 3 (exact), the bound on random local
      frames (relative tolerance 1e-12, rounding of a sum of O(4^N·4^N) unit-size terms), the
      bound reached exactly by the frame as given where every odd-Z term is X-lit and by S on
      every site where it is Y-lit, reached by R_z(−π/4) for the two continuous pairs, and a
      Nelder-Mead minimisation landing on it (reported).
  G3  the hardware consequence. Exact: the superoperator of V = (X+Y)/√2 on the middle qubit
      (built as (X+Y)⊗conj(X+Y)/2, dyadic) carries L(XX+YY) onto L(XY+YX) entry by entry at the
      dyadic γ = 0.5, so the equality holds at every t. Then soft from |+−+⟩ equals truly from
      |+,−i,+⟩ on all 16 two-qubit Paulis of (q0, q2), continuous and in the runner's Trotter
      order, and the reverse; exact for the Trotter circuit too, since V turns truly's ordered step
      into soft's up to one swap of two commuting terms (checked on the letters). The 1e-12 bound on
      the propagated tables is a reading of two different matrix routes, not the evidence. With
      amplitude damping V fails, and the quarter turn S on the middle qubit, which keeps σ⁻,
      carries XY+YX onto YY − XX (truly strings): soft from |+−+⟩ equals YY − XX from |+,−i,+⟩
      at γ_T1 = 0.05 and 0.3 (rounding gated at 1e-12), while the V equality itself breaks there
      (control, difference > 1e-3).
  G4  wider frames go further, the two examples the page names: the CZ star CZ(0,2)·CZ(1,2), diagonal
      and so keeping the Z-dephasing, carries XZ+ZZ at N = 3 onto a truly Hamiltonian; and Π also
      palindromizes Y-dephasing, so R_x(π/2) on the middle site, which turns the light there from Z
      to Y, carries YZ+ZY at N = 3 onto Σ(ZZ − YY) with Π an exact mirror of the turned generator
      (letters and residual exact, unit rates). G4b: one-site Clifford frames of period two with the
      light kept in {Z, Y} carry exactly six of the ten (XX+YZ, XX+ZY, XZ+YY, XZ+ZX, YY+ZX, YZ+ZY) onto
      truly letters at N = 3, 4, 5; the truly letters with Z/Y light make Π an exact mirror.
Control: XX+XY (hard) is not made truly by the odd-site turn, and the bound of G2 is not 0
for any of the sixteen, so G1's zero cannot pass by construction.
"""
import itertools
import sys

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

sys.path.insert(0, __file__.rsplit("simulations", 1)[0] + "simulations")
import framework as fw  # noqa: E402
from framework.lindblad import lindbladian_z_dephasing, palindrome_residual  # noqa: E402
from framework.symmetry import _pauli_tuple_is_truly  # noqa: E402
from framework.workflows.predict_signature_table import predict_signature_table  # noqa: E402

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8")

P = {'I': np.eye(2, dtype=complex), 'X': np.array([[0, 1], [1, 0]], complex),
     'Y': np.array([[0, -1j], [1j, 0]], complex), 'Z': np.diag([1, -1]).astype(complex)}
ALTERNATING = ['XY+YX', 'XY+ZZ', 'YX+ZZ']
UNIFORM = ['XX+XZ', 'XX+YZ', 'XX+ZX', 'XX+ZY', 'XZ+YY', 'XZ+ZX', 'XZ+ZZ', 'YY+YZ', 'YY+ZX',
           'YY+ZY', 'YZ+ZY', 'YZ+ZZ', 'ZX+ZZ', 'ZY+ZZ']
CONTINUOUS = ['XZ+YZ', 'ZX+ZY']
F_MAP = {'I': [('I', 1)], 'X': [('Y', 1)], 'Y': [('X', 1)], 'Z': [('Z', -1)]}   # (X+Y)/√2
S_MAP = {'I': [('I', 1)], 'X': [('Y', 1)], 'Y': [('X', -1)], 'Z': [('Z', 1)]}   # quarter turn about Z
GAMMA = 0.5
FAILS = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not ok:
        FAILS.append(name)


def chain_letters(pair, N):
    H = {}
    for b in range(N - 1):
        for t in pair.split('+'):
            s = ['I'] * N
            s[b], s[b + 1] = t[0], t[1]
            H[''.join(s)] = H.get(''.join(s), 0) + 1
    return H


def turn(H, maps):
    out = {}
    for s, c in H.items():
        acc = [('', c)]
        for i, L in enumerate(s):
            imgs = maps[i][L] if i in maps else [(L, 1)]
            acc = [(pre + l2, cc * c2) for (pre, cc) in acc for (l2, c2) in imgs]
        for k, cc in acc:
            out[k] = out.get(k, 0) + cc
    return {k: v for k, v in out.items() if v != 0}


def dense(H, N):
    M = np.zeros((2 ** N, 2 ** N), complex)
    for s, c in H.items():
        op = P[s[0]]
        for L in s[1:]:
            op = np.kron(op, P[L])
        M += c * op
    return M


def resid2(H, N, gamma=GAMMA):
    """‖M‖² of the dense residual. M is γ-independent (Master Lemma); at γ = 0 only the
    Hamiltonian part is built, which cancels exactly at any coupling, and at the dyadic γ = 0.5
    the dissipator part cancels exactly too. At other γ that part rounds at ~1e-30, which is
    why the exact comparisons below use γ = 0 and 0.5 only."""
    L = lindbladian_z_dephasing(dense(H, N), [gamma] * N)
    return float(np.sum(np.abs(palindrome_residual(L, gamma * N, N)) ** 2))


def letters_truly(H):
    return all(_pauli_tuple_is_truly(tuple(k)) for k in H)


def frame(N, thetas, flips):
    maps = {}
    for l in range(N):
        c, s, f = np.cos(thetas[l]), np.sin(thetas[l]), (-1 if flips[l] else 1)
        maps[l] = {'I': [('I', 1)], 'X': [('X', c), ('Y', s * f)], 'Y': [('X', -s), ('Y', c * f)],
                   'Z': [('Z', f)]}
    return maps


# G1
for N in (3, 4, 5):
    for name, M in (('(X+Y)/√2', F_MAP), ('S', S_MAP)):
        maps = {i: M for i in range(1, N, 2)}
        for pair in ALTERNATING:
            H = chain_letters(pair, N)
            Ht = turn(H, maps)
            ok = letters_truly(Ht) and not letters_truly(H)
            if N <= 4:
                ok = ok and resid2(Ht, N, 0.0) == 0.0 and resid2(Ht, N) == 0.0 and resid2(H, N) > 0.0
            check(f"G1 N={N} {pair} odd sites turned by {name} is truly", ok)
    check(f"G1 N={N} XY+YX turned by (X+Y)/√2 is XX+YY",
          turn(chain_letters('XY+YX', N), {i: F_MAP for i in range(1, N, 2)}) == chain_letters('XX+YY', N))
    check(f"G1 control N={N} XX+XY stays non-truly under the odd-site turn",
          not letters_truly(turn(chain_letters('XX+XY', N), {i: F_MAP for i in range(1, N, 2)})))

# G2
N = 3
w = {s: resid2({s: 1}, N) for s in ['XZI', 'ZXI', 'YZI', 'ZYI', 'XYI', 'XXI', 'YYI', 'ZZI']}
check("G2 F85 weights at N=3: Π²-odd 4·4^N, Π²-even non-truly 8·4^N, truly 0",
      w['XZI'] == w['ZXI'] == w['XYI'] == 4 * 4 ** N and w['YZI'] == w['ZYI'] == 8 * 4 ** N
      and w['XXI'] == w['YYI'] == w['ZZI'] == 0.0, str(w))
rng = np.random.default_rng(20261123)
for pair in UNIFORM + CONTINUOUS:
    H = chain_letters(pair, N)
    oddZ = sum(v ** 2 for k, v in H.items() if k.count('Z') % 2 == 1)
    floor = 4 * 4 ** N * oddZ
    additive = resid2(H, N) == sum(resid2({k: v}, N) for k, v in H.items())
    samples = [resid2(local := turn(H, frame(N, rng.uniform(0, 2 * np.pi, N), rng.integers(0, 2, N))), N)
               for _ in range(40)]
    best = np.inf
    for flips in itertools.product([0, 1], repeat=N):
        for _ in range(3):
            r = minimize(lambda th: resid2(turn(H, frame(N, th, flips)), N), rng.uniform(0, 2 * np.pi, N),
                         method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-12, 'maxiter': 4000})
            best = min(best, r.fun)
    check(f"G2 {pair}: odd-Z content {oddZ}, additive, every sampled frame ≥ floor {floor}",
          oddZ > 0 and additive and min(samples) >= floor * (1 - 1e-12),
          f"min sampled {min(samples):.6f}, minimised {best:.9f}, ratio {best / floor:.9f}")
X_LIT = ['XX+XZ', 'XX+ZX', 'XZ+YY', 'XZ+ZX', 'XZ+ZZ', 'YY+ZX', 'ZX+ZZ']
Y_LIT = ['XX+YZ', 'XX+ZY', 'YY+YZ', 'YY+ZY', 'YZ+ZY', 'YZ+ZZ', 'ZY+ZZ']
check("G2 the sixteen split into X-lit, Y-lit and continuous", sorted(X_LIT + Y_LIT + CONTINUOUS) == sorted(UNIFORM + CONTINUOUS))
for pair in X_LIT + Y_LIT:
    H = chain_letters(pair, N)
    floor = 4 * 4 ** N * sum(v ** 2 for k, v in H.items() if k.count('Z') % 2 == 1)
    Ht = H if pair in X_LIT else turn(H, {i: S_MAP for i in range(N)})
    check(f"G2 {pair}: the floor {floor} is reached exactly ({'as given' if pair in X_LIT else 'S on every site'})",
          resid2(Ht, N) == floor and resid2(Ht, N, 0.0) == floor)
for pair in CONTINUOUS:
    H = chain_letters(pair, N)
    floor = 4 * 4 ** N * sum(v ** 2 for k, v in H.items() if k.count('Z') % 2 == 1)
    r = resid2(turn(H, frame(N, [-np.pi / 4] * N, [0] * N)), N)
    # R_z(−π/4) has irrational entries: relative rounding of order 1e-15 on a sum of O(4^N·4^N) terms
    check(f"G2 {pair}: the floor {floor} is reached by R_z(−π/4) on every site", abs(r / floor - 1) < 1e-12, f"{r!r}")
for pair in ALTERNATING:
    check(f"G2 control {pair} has no odd-Z string",
          all(k.count('Z') % 2 == 0 for k in chain_letters(pair, N)))

# G3
ket = {'+': np.array([1, 1]) / np.sqrt(2), '-': np.array([1, -1]) / np.sqrt(2),
       '-i': np.array([1, -1j]) / np.sqrt(2)}


def rho(labels):
    psi = ket[labels[0]]
    for l in labels[1:]:
        psi = np.kron(psi, ket[l])
    return np.outer(psi, psi.conj())


A, B = rho(['+', '-', '+']), rho(['+', '-i', '+'])
terms = {'truly': [('X', 'X'), ('Y', 'Y')], 'soft': [('X', 'Y'), ('Y', 'X')]}
chain = fw.ChainSystem(N=3)
tA = predict_signature_table(chain, terms, t=0.8, n_trotter=3, gamma_z=0.1, rho_0=A)
tB = predict_signature_table(chain, terms, t=0.8, n_trotter=3, gamma_z=0.1, rho_0=B)
d_trot = max(max(abs(tA['soft'][k] - tB['truly'][k]), abs(tA['truly'][k] - tB['soft'][k])) for k in tA['soft'])


def cont(pair, r0):
    L = lindbladian_z_dephasing(dense(chain_letters(pair, 3), 3), [0.1] * 3)
    return (expm(L * 0.8) @ r0.flatten()).reshape(8, 8)


obs = {(a, b): np.kron(np.kron(P[a], P['I']), P[b]) for a in 'IXYZ' for b in 'IXYZ'}
cs, ct = cont('XY+YX', A), cont('XX+YY', B)
cs2, ct2 = cont('XY+YX', B), cont('XX+YY', A)
d_cont = max(max(abs(np.trace(cs @ o) - np.trace(ct @ o)), abs(np.trace(cs2 @ o) - np.trace(ct2 @ o))) for o in obs.values())
check("G3 soft(|+−+⟩) = truly(|+,−i,+⟩) and the reverse, 16 Paulis, Trotter (runner order)", d_trot < 1e-12, f"{d_trot:.2e}")
check("G3 the same, continuous Lindblad", d_cont < 1e-12, f"{d_cont:.2e}")
print(f"  ⟨X₀Z₂⟩ Trotter: soft(+−+) {tA['soft'][('X', 'Z')]:+.4f}, truly(+,−i,+) {tB['truly'][('X', 'Z')]:+.4f}")



# G3, exact generator identity and the amplitude-damping route
from framework.lindblad import lindbladian_z_plus_t1  # noqa: E402

XpY = P['X'] + P['Y']
Vmid = np.kron(np.kron(P['I'], XpY), P['I'])
Wsup = np.kron(Vmid, Vmid.conj()) / 2
Lt = lindbladian_z_dephasing(dense(chain_letters('XX+YY', 3), 3), [0.5] * 3)
Ls = lindbladian_z_dephasing(dense(chain_letters('XY+YX', 3), 3), [0.5] * 3)
check("G3 exact: the superoperator of V carries L(XX+YY) onto L(XY+YX) entry by entry",
      np.array_equal(Wsup @ Lt @ Wsup.conj().T, Ls))
Smid = np.kron(np.kron(P['I'], np.diag([1, 1j])), P['I'])
Hss = Smid @ dense(chain_letters('XY+YX', 3), 3) @ Smid.conj().T
TRUL_ORDER = ['XXI', 'YYI', 'IXX', 'IYY']      # the runner's order: bonds outside, terms inside
SOFT_ORDER = ['XYI', 'YXI', 'IXY', 'IYX']
imgs = [turn({t: 1}, {1: F_MAP}) for t in TRUL_ORDER]
check("G3 exact Trotter route: V turns truly's ordered step into soft's, up to one swap of two commuting terms",
      imgs == [{'XYI': 1}, {'YXI': 1}, {'IYX': 1}, {'IXY': 1}]
      and np.array_equal(dense({'IXY': 1}, 3) @ dense({'IYX': 1}, 3), dense({'IYX': 1}, 3) @ dense({'IXY': 1}, 3)))
check("G3 S on the middle qubit carries XY+YX onto YY − XX",
      np.array_equal(Hss, dense({'YYI': 1, 'IYY': 1, 'XXI': -1, 'IXX': -1}, 3)))
for g1 in (0.05, 0.3):
    La = lindbladian_z_plus_t1(dense(chain_letters('XY+YX', 3), 3), [0.1] * 3, [g1] * 3)
    Lb = lindbladian_z_plus_t1(Hss, [0.1] * 3, [g1] * 3)
    ra = (expm(La * 0.8) @ A.flatten()).reshape(8, 8)
    rb = (expm(Lb * 0.8) @ B.flatten()).reshape(8, 8)
    d = max(abs(np.trace(ra @ o) - np.trace(rb @ o)) for o in obs.values())
    check(f"G3 with amplitude damping γ_T1={g1}: soft(|+−+⟩) = (YY − XX)(|+,−i,+⟩)", d < 1e-12, f"{d:.2e}")
    Lv = lindbladian_z_plus_t1(dense(chain_letters('XX+YY', 3), 3), [0.1] * 3, [g1] * 3)
    rv = (expm(Lv * 0.8) @ B.flatten()).reshape(8, 8)
    dv = max(abs(np.trace(ra @ o) - np.trace(rv @ o)) for o in obs.values())
    check(f"G3 control γ_T1={g1}: V fails, soft(|+−+⟩) differs from truly(|+,−i,+⟩)", dv > 1e-3, f"{dv:.3f}")


# G4, the two wider frames the page names
def resid2_dense(Hm, N, cops=None, sigma=0.0):
    from framework.lindblad import lindbladian_general
    L = lindbladian_general(Hm, cops or [])
    return float(np.sum(np.abs(palindrome_residual(L, sigma, N)) ** 2))


def site(N, l, L):
    ops = [P['I']] * N
    ops[l] = P[L]
    out = ops[0]
    for o in ops[1:]:
        out = np.kron(out, o)
    return out


star = np.diag([(-1) ** ((((x >> 2) & 1) + ((x >> 1) & 1)) * (x & 1)) for x in range(8)]).astype(complex)  # CZ(0,2)·CZ(1,2)
Hm = dense(chain_letters('XZ+ZZ', 3), 3)
check("G4 the CZ star carries XZ+ZZ at N=3 onto X₀Z₁Z₂ + X₁ + Z₀Z₁ + Z₁Z₂",
      np.array_equal(star @ Hm @ star.conj().T, dense({'XZZ': 1, 'IXI': 1, 'ZZI': 1, 'IZZ': 1}, 3)))
check("G4 ... which is truly (Hamiltonian part of the residual exactly 0, before > 0)",
      resid2_dense(star @ Hm @ star.conj().T, 3) == 0.0 and resid2_dense(Hm, 3) > 0.0)
RX = {'I': [('I', 1)], 'X': [('X', 1)], 'Y': [('Z', 1)], 'Z': [('Y', -1)]}   # R_x(π/2): Y → Z, Z → −Y
Hy = turn(chain_letters('YZ+ZY', 3), {1: RX})
check("G4 R_x(π/2) on the middle site carries YZ+ZY at N=3 onto Σ(ZZ − YY)",
      Hy == {'ZZI': 1, 'IZZ': 1, 'YYI': -1, 'IYY': -1})
cops = [site(3, 0, 'Z'), site(3, 1, 'Y'), site(3, 2, 'Z')]   # the Z-dephasing turned: Y on the middle site, unit rates
check("G4 ... and Π is an exact mirror there: residual 0.0 with Y-dephasing on the middle site",
      resid2_dense(dense(Hy, 3), 3, cops, 3.0) == 0.0
      and resid2_dense(dense(chain_letters('YZ+ZY', 3), 3), 3, [site(3, l, 'Z') for l in range(3)], 3.0) > 0.0)

TEN_ALL = ['XX+XZ', 'XX+YZ', 'XX+ZX', 'XX+ZY', 'XZ+YY', 'XZ+ZX', 'YY+YZ', 'YY+ZX', 'YY+ZY', 'YZ+ZY']


# G4b, the six by period-two one-site Clifford frames with the light in {Z, Y}
def clifford_maps():
    out = []
    for lx, lz in itertools.permutations('XYZ', 2):
        for sx in (1, -1):
            for sz in (1, -1):
                M = 1j * sx * sz * P[lx] @ P[lz]
                for ly in 'XYZ':
                    for sy in (1, -1):
                        if np.array_equal(M, sy * P[ly]):
                            out.append({'I': [('I', 1)], 'X': [(lx, sx)], 'Y': [(ly, sy)], 'Z': [(lz, sz)]})
    return out


CL = clifford_maps()
SIX_LIGHT = {'XX+YZ', 'XX+ZY', 'XZ+YY', 'XZ+ZX', 'YY+ZX', 'YZ+ZY'}
check("G4b the 24 one-site Cliffords", len(CL) == 24)
for n in (3, 4, 5):
    hit = set()
    for pair in TEN_ALL:
        H = chain_letters(pair, n)
        for ce, co in itertools.product(range(24), repeat=2):
            maps = {i: (CL[ce] if i % 2 == 0 else CL[co]) for i in range(n)}
            if any(maps[i]['Z'][0][0] not in 'ZY' for i in range(n)):
                continue          # the turned light must stay in {Z, Y}, which Π palindromizes
            if letters_truly(turn(H, maps)):
                hit.add(pair)
                break
    check(f"G4b N={n}: such frames carry exactly the six onto chains Π mirrors", hit == SIX_LIGHT, f"{sorted(hit)}")

print()
print(f"{len(FAILS)} failure(s)" if FAILS else "ALL PASS")
sys.exit(1 if FAILS else 0)
