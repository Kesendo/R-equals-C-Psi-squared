"""Gate for the readings of the F88b lens and of the Lebensader skeleton on the 2026 flights.

The flights: the framework_snapshots runs on Marrakesh and the soft_break runs on Marrakesh and
Kingston, read through the F88b lens, and the Lebensader runs on Marrakesh 2026-04-26 (EQ-030).
All at N = 3 from the X-basis product |+−+⟩, every term on both bonds.

  G1  The F88b lens reads the Π² parity of H. Π² parity of a Pauli string is (−1)^(n_Y + n_Z),
      a character of the letters (checked on the 16 products); so a Π²-even string sends each
      parity sector to itself, Z-dephasing is diagonal on strings, and |+−+⟩ carries only I and
      X: a Π²-even H gives exactly zero Π²-odd content at every time and in any order of its
      gates. Ten of the 36 pairs are Π²-even, the three truly ones and seven soft ones (YZ+ZY
      among them). The converse is one-way: a Π²-odd pair can read small (the table). The
      hardware values are recomputed from the two soft_break files in the repo.
  G2  The census of silent strings under T1 keeps the silence complex conjugation carries.
      σ⁻ is real, Z σ⁻ Z = −σ⁻, X σ⁻ X = σ⁺: conjugation and diagonal phases keep the
      cooling, a flip turns it into heating. For a purely imaginary H (H* = −H, every string
      with an odd number of Y) from a real start under real dissipators, ρ(t) stays real, so the 28
      strings with odd n_Y are silent with T1 on; the cockpit's y_parity panel reports the same
      set. The census (the cockpit's instrument, at the cockpit test's γ = 0.05, γ_T1 = 0.005;
      EQ-030's 0.1 and 0.01 give the same counts; an integer Krylov census, with no eigensolver,
      gives the same silent sets) keeps 28 or 29 strings for the six imaginary
      pairs and at most 2 for the other thirty. YZ+ZY annihilates |+−+⟩, so under Z dephasing
      ρ(t) is the dephased start and its 56 silent strings are those with a Y or Z; T1 moves the
      state off the start and leaves the 28 of conjugation. The split belongs to the pair of H and start:
      XX+YY from |+, −i, +⟩ keeps 29, 28 of them by Z₁·K.
  G3  The flown observable ⟨X₀Z₂⟩: for every Π²-even H from |+−+⟩ under Z-dephasing it is
      silent by X^⊗3 (YZ+ZY, the truly XX+YY and the Heisenberg chain alike), a flip symmetry, so
      T1 lifts it; at the run's simulator rate γ_T1 = 0.01 Z + T1 keeps YZ+ZY's at most 0.0051 at
      the four flown times, and at t = 0.8 it gives −0.001 (−0.004 at the calibration-derived
      0.046) in the continuous model, against the first hardware reading, +0.1338 in EQ-030's
      table, 8.6σ. Where X^⊗3 is no symmetry, other mirrors can hold the same string at zero:
      X₂ alone on XX+YX, X^⊗3 followed by complex conjugation on XZ+ZX; on XZ+ZZ no Pauli
      string, alone or followed by complex conjugation, does, and the string is silent anyway.
  G4  Record: in the two soft_break tables in the repo, the (q0, q2) cells of the truly chain
      that T1 should lift read no larger than the cells of XY+YX that stay silent under T1, so at
      these rates the tables do not show the asymmetry; Z + T1 (continuous, γ_Z = 0.1, t = 0.8)
      lifts the truly chain's four cells to at most 0.014 for γ_T1 up to γ_Z, below the shot
      noise 1/√4096.
  G5  The frames and flips the prose names: S on every site carries YZ+ZY onto −(XZ+ZX), a
      Π²-odd pair, keeping the dephasing and the decay, so the F83 anti-fraction of one dynamics reads 0 in one frame
      and 1/2 in the other; XZ+ZX from |+−+⟩ is silent by X^⊗3·K, a flip, which T1 lifts.

Run: python simulations/pi2_lens_and_t1_silence_gate.py   (under a minute)
"""
import itertools
import json
import sys
from fractions import Fraction
from math import lcm
from pathlib import Path

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).resolve().parent))
if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import framework as fw
from framework.pauli import _build_bilinear, _vec_to_pauli_basis_transform, pauli_basis_vector, _pauli_label
from framework.lindblad import lindbladian_z_dephasing, lindbladian_z_plus_t1
from framework.lebensader import cockpit_panel
from framework.symmetry import y_parity_panel

REPO = Path(__file__).resolve().parents[1]
FAILS = []


def check(name, ok, detail=''):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ''))
    if not ok:
        FAILS.append(name)


N = 3
P = {'I': np.eye(2, dtype=complex), 'X': np.array([[0, 1], [1, 0]], dtype=complex),
     'Y': np.array([[0, -1j], [1j, 0]], dtype=complex), 'Z': np.diag([1, -1]).astype(complex)}
SM = np.array([[0, 1], [0, 0]], dtype=complex)            # σ⁻, the cockpit's convention


def kron(*ms):
    out = ms[0]
    for m in ms[1:]:
        out = np.kron(out, m)
    return out


def pi2(s):
    return (s.count('Y') + s.count('Z')) % 2


STRINGS = [''.join(s) for s in itertools.product('IXYZ', repeat=N)]
MAT = {s: kron(*[P[c] for c in s]) for s in STRINGS}
TERMS9 = [a + b for a in 'XYZ' for b in 'XYZ']
PAIRS = list(itertools.combinations(TERMS9, 2))
chain = fw.ChainSystem(N=N, gamma_0=0.05)


def ket(*names):
    v = {'+': [1, 1], '-': [1, -1], '+i': [1, 1j], '-i': [1, -1j]}
    return kron(*[np.array(v[n], dtype=complex) / np.sqrt(2) for n in names])


S0 = ket('+', '-', '+'); RHO0 = np.outer(S0, S0.conj())
ST = ket('+', '-i', '+'); RHOT = np.outer(ST, ST.conj())


def H_of(terms):
    return _build_bilinear(N, chain.bonds, [(a, b, c * chain.J) for a, b, c in terms])


def H_pair(pair):
    return H_of([(t[0], t[1], 1) for t in pair])


def strings_of_pair(pair):
    out = set()
    for (i, j) in chain.bonds:
        for t in pair:
            s = ['I'] * N
            s[i], s[j] = t[0], t[1]
            out.add(''.join(s))
    return out


print("G1  the F88b lens reads the Π² parity of H")
ok = True
for a, b in itertools.product('IXYZ', repeat=2):
    prod = P[a] @ P[b]
    c = [l for l in 'IXYZ' if any(np.array_equal(prod, ph * P[l]) for ph in (1, -1, 1j, -1j))]
    ok &= len(c) == 1 and pi2(c[0]) == (pi2(a) + pi2(b)) % 2
check("G1a Π² parity (−1)^(n_Y + n_Z) is a character of the letters (16 products)", ok)
even_pairs = [p for p in PAIRS if all(pi2(t) == 0 for t in p)]
classes = {p: fw.classify_pauli_pair(chain, [tuple(t) for t in p]) for p in PAIRS}
check("G1b ten of the 36 pairs are Π²-even: the three truly ones and seven soft ones",
      len(even_pairs) == 10 and sorted(classes[p] for p in even_pairs).count('truly') == 3
      and sorted(classes[p] for p in even_pairs).count('soft') == 7,
      ', '.join('+'.join(p) for p in even_pairs))
check("G1b ... YZ+ZY among them", ('YZ', 'ZY') in even_pairs)
check("G1c |+−+⟩ carries only I and X strings (start value exactly 0 on every string with a Y or Z)",
      all(np.trace(RHO0 @ MAT[s]) == 0.0 for s in STRINGS if any(c in 'YZ' for c in s)))
# the ideal lens at t = 0.8, γ = 0.1, continuous (a reading; the zero rests on G1a to G1c)
M_basis = _vec_to_pauli_basis_transform(N, order='C')


def lens_2q(r02):
    P0 = np.diag([1, 0, 0, 0]).astype(complex); Pm = np.diag([0, 1, 1, 0]).astype(complex)
    P2 = np.diag([0, 0, 0, 1]).astype(complex)
    d0 = np.trace(P0 @ r02).real * P0 + np.trace(Pm @ r02).real / 2 * Pm + np.trace(P2 @ r02).real * P2
    d2 = r02 - d0
    odd = sum(abs(np.trace(kron(P[a], P[b]) @ d2) / 4) ** 2 * 4
              for a in 'IXYZ' for b in 'IXYZ' if (pi2(a) + pi2(b)) % 2)
    mem = np.linalg.norm(d2) ** 2
    return odd / mem if mem > 1e-14 else 0.0


def reduce02(rho):
    return np.einsum('abcdbf->acdf', rho.reshape(2, 2, 2, 2, 2, 2)).reshape(4, 4)


ideal = {}
for p in PAIRS:
    L = lindbladian_z_dephasing(H_pair(p), [0.1] * N)
    rho = (expm(L * 0.8) @ RHO0.reshape(-1)).reshape(8, 8)
    ideal[p] = lens_2q(reduce02(rho))
odd_vals = sorted((ideal[p], '+'.join(p)) for p in PAIRS if p not in even_pairs)
check("G1d reading: ideal lens at t = 0.8 is at rounding for the ten Π²-even pairs",
      max(ideal[p] for p in even_pairs) < 1e-20, f"max {max(ideal[p] for p in even_pairs):.1e}")
check("G1d ... and above 1e-4 for every Π²-odd or mixed pair, the smallest small (one-way rule)",
      odd_vals[0][0] > 1e-4, f"smallest {odd_vals[0][1]} {odd_vals[0][0]:.4f}, largest {odd_vals[-1][1]} {odd_vals[-1][0]:.4f}")


def reconstruct(exp):
    r = np.zeros((4, 4), dtype=complex)
    for key, v in exp.items():
        a, b = key.split(',')
        r += v * kron(P[a], P[b]) / 4
    return r


kj = json.loads((REPO / 'data/ibm_soft_break_april2026/soft_break_ibm_kingston_20260505_102806.json').read_text(encoding='utf-8'))
mj = json.loads((REPO / 'data/ibm_soft_break_april2026/soft_break_ibm_marrakesh_20260426_001101.json').read_text(encoding='utf-8'))
kl = {c: lens_2q(reconstruct(e)) for c, e in kj['expectations'].items()}
ml = {c: lens_2q(reconstruct(e)) for c, e in mj['expectations'].items()}
check("G1e Kingston 2026-05-05: truly 0.0022, XY+YX 0.7409, YZ+ZY (Π²-even, soft) 0.0046, XY+YZ 0.6220",
      [round(kl[c], 4) for c in ('truly_unbroken', 'pi2_odd_pure', 'pi2_even_nontruly', 'mixed_anti_one_sixth')]
      == [0.0022, 0.7409, 0.0046, 0.6220], str({c: round(v, 4) for c, v in kl.items()}))
check("G1e Marrakesh 2026-04-26 [48, 49, 50]: truly 0.0013, XY+YX 0.7646",
      round(ml['truly_unbroken'], 4) == 0.0013 and round(ml['soft_broken'], 4) == 0.7646, str({c: round(v, 4) for c, v in ml.items()}))

print()
print("G2  under T1 the census keeps the silence complex conjugation carries")
check("G2a σ⁻ is real, Z σ⁻ Z = −σ⁻, X σ⁻ X = σ⁺",
      np.array_equal(SM.conj(), SM) and np.array_equal(P['Z'] @ SM @ P['Z'], -SM)
      and np.array_equal(P['X'] @ SM @ P['X'], SM.conj().T))


def silent(L, rho0, threshold=1e-9, cluster_tol=1e-8):
    """The cockpit's census (cluster sums of lebensader.py), returning the strings."""
    Lp = (M_basis.conj().T @ L @ M_basis) / (2 ** N)
    ev, V = np.linalg.eig(Lp)
    c = np.linalg.inv(V) @ pauli_basis_vector(rho0, N)
    used = np.zeros(len(ev), bool); clusters = []
    for i in range(len(ev)):
        if used[i]:
            continue
        cl = [i]; used[i] = True
        for j in range(i + 1, len(ev)):
            if not used[j] and abs(ev[j] - ev[i]) < cluster_tol:
                cl.append(j); used[j] = True
        clusters.append(cl)
    return {_pauli_label(a, N) for a in range(1, 4 ** N)
            if max(abs(sum(V[a, k] * c[k] for k in cl)) for cl in clusters) < threshold}


def to_int(arr):
    """An array of rationals as integers over one common denominator (each entry checked)."""
    fr = [[Fraction(float(x)).limit_denominator(10 ** 6) for x in row] for row in arr]
    assert all(abs(float(fr[i][j]) - float(arr[i][j])) < 1e-12 for i in range(len(arr)) for j in range(len(arr[0])))
    den = 1
    for row in fr:
        for x in row:
            den = lcm(den, x.denominator)
    return [[int(x * den) for x in row] for row in fr]


def exact_silent(L, rho0):
    """The census without an eigensolver: ⟨P⟩(t) vanishes for all t exactly when coordinate P of
    L^k v₀ vanishes for k = 0 .. 4^N − 1 (Cayley-Hamilton), computed in integers."""
    Lp = (M_basis.conj().T @ L @ M_basis) / (2 ** N)
    assert np.max(np.abs(Lp.imag)) < 1e-12
    A = to_int(Lp.real)
    v0 = pauli_basis_vector(rho0, N)
    assert np.max(np.abs(v0.imag)) < 1e-12
    v = [r[0] for r in to_int(v0.real.reshape(-1, 1))]
    seen = [x != 0 for x in v]
    for _ in range(4 ** N):
        v = [sum(a * b for a, b in zip(row, v) if b) for row in A]
        seen = [s or (x != 0) for s, x in zip(seen, v)]
    return {_pauli_label(a, N) for a in range(1, 4 ** N) if not seen[a]}


def census(H, rho0):
    return (silent(lindbladian_z_dephasing(H, [0.05] * N), rho0),
            silent(lindbladian_z_plus_t1(H, [0.05] * N, [0.005] * N), rho0))


ODD_Y = {s for s in STRINGS if s.count('Y') % 2}
imag_pairs = [p for p in PAIRS if np.array_equal(H_pair(p).conj(), -H_pair(p))]
check("G2b six of the 36 pairs are purely imaginary (H* = −H): XY+YX, YZ+ZY and four hard ones",
      sorted('+'.join(p) for p in imag_pairs) == ['XY+YX', 'XY+YZ', 'XY+ZY', 'YX+YZ', 'YX+ZY', 'YZ+ZY']
      and sorted(classes[p] for p in imag_pairs) == ['hard'] * 4 + ['soft'] * 2)
check("G2b ... every string of those six has an odd number of Y, and 28 strings at N = 3 do",
      all(s.count('Y') % 2 for p in imag_pairs for s in strings_of_pair(p)) and len(ODD_Y) == 28)
rows = {}
for p in PAIRS:
    H = H_pair(p)
    sp, st = census(H, RHO0)
    sk = cockpit_panel(H, [0.05] * N, RHO0, N, gamma_t1_l=[0.005] * N, t_max=1.0, dt=0.01)['lebensader']['skeleton']
    yp = set(y_parity_panel(H, [0.05] * N, RHO0, N, gamma_t1_l=[0.005] * N)['Y_parity_protected'])
    rows[p] = (sp, st, sk, yp)
check("G2c the census here counts as the cockpit's own (pure and T1, all 36)",
      all(len(r[0]) == r[2]['n_protected_pure'] and len(r[1]) == r[2]['n_protected_t1'] for r in rows.values()))
check("G2c ... and it is exact: integer Krylov (no eigensolver, no threshold) gives the same silent sets, pure and T1, all 36",
      all(exact_silent(lindbladian_z_dephasing(H_pair(p), [0.05] * N), RHO0) == rows[p][0]
          and exact_silent(lindbladian_z_plus_t1(H_pair(p), [0.05] * N, [0.005] * N), RHO0) == rows[p][1] for p in PAIRS))
check("G2c under T1 the six imaginary pairs keep 28 or 29 silent strings, containing the 28 odd-Y ones",
      all(len(rows[p][1]) in (28, 29) and ODD_Y <= rows[p][1] for p in imag_pairs),
      str({'+'.join(p): (len(rows[p][0]), len(rows[p][1])) for p in imag_pairs}))
others = [p for p in PAIRS if p not in imag_pairs]
check("G2c ... and the other thirty keep at most 2",
      max(len(rows[p][1]) for p in others) <= 2, f"max {max(len(rows[p][1]) for p in others)}")
Hyz = H_pair(('YZ', 'ZY'))
yz_strings = {s for s in STRINGS if any(c in 'YZ' for c in s)}
check("G2c YZ+ZY annihilates |+−+⟩ (exactly), so its pure-Z census is the dephased start: every string with a Y or Z (56)",
      np.array_equal(Hyz @ S0, np.zeros(8, dtype=complex)) and rows[('YZ', 'ZY')][0] == yz_strings)
check("G2c ... and under T1 it keeps exactly the 28 strings conjugation silences", rows[('YZ', 'ZY')][1] == ODD_Y)
eq030 = {}
for name, H0, r0 in (('XY+YX', H_pair(('XY', 'YX')), RHO0), ('YZ+ZY', H_pair(('YZ', 'ZY')), RHO0),
                     ('XZ+ZX', H_pair(('XZ', 'ZX')), RHO0), ('XX+YY', H_pair(('XX', 'YY')), RHO0),
                     ('IY+YI', H_of([('I', 'Y', 1), ('Y', 'I', 1)]), RHO0), ('XX+YY turned', H_pair(('XX', 'YY')), RHOT)):
    sp1 = silent(lindbladian_z_dephasing(H0, [0.05] * N), r0); st1 = silent(lindbladian_z_plus_t1(H0, [0.05] * N, [0.005] * N), r0)
    sp2 = silent(lindbladian_z_dephasing(H0, [0.1] * N), r0); st2 = silent(lindbladian_z_plus_t1(H0, [0.1] * N, [0.01] * N), r0)
    eq030[name] = ((len(sp1), len(st1)), (len(sp2), len(st2)))
check("G2c EQ-030's rates (γ = 0.1, γ_T1 = 0.01) give the same counts as the cockpit test's",
      all(a == b for a, b in eq030.values()), str({k: v[1] for k, v in eq030.items()}))
check("G2d y_parity_panel reports exactly the 28 odd-Y strings for the six and nothing for the thirty",
      all(rows[p][3] == ODD_Y for p in imag_pairs) and all(not rows[p][3] for p in others))
drops = {'+'.join(p): rows[p][2]['drop'] for p in PAIRS}
check("G2d reading: drop ≤ 1 exactly for XY+YX and the four imaginary hard pairs; YZ+ZY 28, XX+YY 31",
      sorted(k for k, v in drops.items() if v <= 1) == ['XY+YX', 'XY+YZ', 'XY+ZY', 'YX+YZ', 'YX+ZY']
      and drops['YZ+ZY'] == 28 and drops['XX+YY'] == 31)
Htr = H_of([('X', 'X', 1), ('Y', 'Y', 1)])
spt, stt = census(Htr, RHOT)
Z1 = MAT['IZI']
check("G2e XX+YY from |+, −i, +⟩: Z₁·K is a symmetry of H, start, dephasing and decay (Z₁ H* Z₁ = −H, Z₁ ρ₀* Z₁ = ρ₀, Z σ⁻* Z = −σ⁻)",
      np.array_equal(Z1 @ Htr.conj() @ Z1, -Htr) and np.array_equal(Z1 @ RHOT.conj() @ Z1, RHOT)
      and np.array_equal(P['Z'] @ SM.conj() @ P['Z'], -SM))
check("G2e ... and its census keeps its skeleton under T1, 30 → 29, against 32 → 1 from |+−+⟩",
      (len(spt), len(stt)) == (30, 29) and (len(rows[('XX', 'YY')][0]), len(rows[('XX', 'YY')][1])) == (32, 1))
z1k = {Q for Q in STRINGS[1:] if np.array_equal(Z1 @ MAT[Q].conj() @ Z1, -MAT[Q])}
check("G2e ... 28 of the 29 are the strings Z₁·K reverses", len(z1k) == 28 and z1k <= stt, f"the rest {sorted(stt - z1k)}")
Hiy = H_of([('I', 'Y', 1), ('Y', 'I', 1)])
spi, sti = census(Hiy, RHO0)
check("G2e ... exact as well (integer Krylov)", exact_silent(lindbladian_z_plus_t1(Htr, [0.05] * N, [0.005] * N), RHOT) == stt)
check("G2f IY+YI (a Y field, imaginary) from |+−+⟩: 37 → 37, every string with a Y (each site stays real)",
      (len(spi), len(sti)) == (37, 37) and sti == {s for s in STRINGS if 'Y' in s})

print()
print("G3  the flown observable ⟨X₀Z₂⟩")
X3 = MAT['XXX']
for p in [('YZ', 'ZY'), ('XX', 'YY')]:
    H = H_pair(p)
    check(f"G3 {'+'.join(p)}: X^⊗3 keeps H and |+−+⟩ and reverses X₀Z₂, so ⟨X₀Z₂⟩ = 0 under Z-dephasing",
          np.array_equal(X3 @ H @ X3, H) and np.array_equal(X3 @ RHO0 @ X3, RHO0)
          and np.array_equal(X3 @ MAT['XIZ'] @ X3, -MAT['XIZ']) and 'XIZ' in rows[p][0])
    check(f"G3 ... and T1 lifts it: XIZ is not in {'+'.join(p)}'s T1 census", 'XIZ' not in rows[p][1])
check("G3 XY+YX: ⟨X₀Z₂⟩ is active (not silent) under Z-dephasing", 'XIZ' not in rows[('XY', 'YX')][0])
Hheis = H_of([('X', 'X', 1), ('Y', 'Y', 1), ('Z', 'Z', 1)])
check("G3 the Heisenberg chain: X^⊗3 keeps H and |+−+⟩, so its ⟨X₀Z₂⟩ zero is the same parity",
      np.array_equal(X3 @ Hheis @ X3, Hheis) and 'XIZ' in silent(lindbladian_z_dephasing(Hheis, [0.05] * N), RHO0))
Hxyx = H_pair(('XX', 'YX'))
X2 = MAT['IIX']
check("G3 XX+YX, where X^⊗3 is no symmetry: X₂ alone keeps H and |+−+⟩ and reverses X₀Z₂, so ⟨X₀Z₂⟩ is silent by another mirror",
      not np.array_equal(X3 @ Hxyx @ X3, Hxyx) and np.array_equal(X2 @ Hxyx @ X2, Hxyx)
      and np.array_equal(X2 @ RHO0 @ X2, RHO0) and np.array_equal(X2 @ MAT['XIZ'] @ X2, -MAT['XIZ'])
      and 'XIZ' in rows[('XX', 'YX')][0])
Hxzx = H_pair(('XZ', 'ZX'))
check("G3 XZ+ZX: X^⊗3 reverses H, which is real, so X^⊗3 followed by complex conjugation keeps the generator; with a real start and a real X^⊗3-odd readout ⟨X₀Z₂⟩ is silent",
      np.array_equal(X3 @ Hxzx @ X3, -Hxzx) and not np.any(Hxzx.imag) and not np.any(RHO0.imag)
      and not np.any(MAT['XIZ'].imag) and 'XIZ' in rows[('XZ', 'ZX')][0])
Hxzz = H_pair(('XZ', 'ZZ'))
lin = [s for s in STRINGS if np.array_equal(MAT[s] @ Hxzz @ MAT[s], Hxzz) and np.array_equal(MAT[s] @ RHO0 @ MAT[s], RHO0)
       and np.array_equal(MAT[s] @ MAT['XIZ'] @ MAT[s], -MAT['XIZ'])]
anti = [s for s in STRINGS if np.array_equal(MAT[s] @ Hxzz.conj() @ MAT[s], -Hxzz)
        and np.array_equal(MAT[s] @ RHO0.conj() @ MAT[s], RHO0) and np.array_equal(MAT[s] @ MAT['XIZ'].conj() @ MAT[s], -MAT['XIZ'])]
check("G3 XZ+ZZ: ⟨X₀Z₂⟩ is silent (exact census), with no Pauli string, alone or followed by complex conjugation, as its mirror",
      'XIZ' in rows[('XZ', 'ZZ')][0] and not lin and not anti, f"linear {lin}, antilinear {anti}")
Lyz = lindbladian_z_plus_t1(H_pair(('YZ', 'ZY')), [0.1] * N, [0.01] * N)
model = max(abs(np.trace(MAT['XIZ'] @ (expm(Lyz * t) @ RHO0.reshape(-1)).reshape(8, 8)).real) for t in (0.8, 1.4, 1.7, 2.2))
at08 = []
for gt1 in (0.01, 0.046):
    Lg = lindbladian_z_plus_t1(H_pair(('YZ', 'ZY')), [0.1] * N, [gt1] * N)
    at08.append(np.trace(MAT['XIZ'] @ (expm(Lg * 0.8) @ RHO0.reshape(-1)).reshape(8, 8)).real)
check("G3 YZ+ZY at the run's simulator rate γ_T1 = 0.01: Z + T1 keeps ⟨X₀Z₂⟩ below a tenth of the hardware's +0.1338 at t = 0.8, 1.4, 1.7, 2.2",
      model < 0.01338, f"model max {model:.4f}")
check("G3 ... at t = 0.8 it stays below 0.005 in magnitude (continuous) at γ_T1 = 0.01 and 0.046, against +0.1338 in EQ-030's table",
      all(abs(v) < 0.005 for v in at08), f"{[round(v, 4) for v in at08]}")
print("      record arithmetic: 0.1338 at 4096 shots is", round(0.1338 * 64, 1), "σ")

print()
print("G4  the two soft_break tables in the repo against the census classes (record)")
CELLS = [a + 'I' + b for a in 'XYZ' for b in 'XYZ']
for tag, js, tk, sk_ in (('Marrakesh [48, 49, 50]', mj, 'truly_unbroken', 'soft_broken'),
                         ('Kingston [43, 56, 63]', kj, 'truly_unbroken', 'pi2_odd_pure')):
    lifted = [c for c in CELLS if c in rows[('XX', 'YY')][0] and c not in rows[('XX', 'YY')][1]]
    kept = [c for c in CELLS if c in rows[('XY', 'YX')][1]]
    lv = max(abs(js['expectations'][tk][c[0] + ',' + c[2]]) for c in lifted)
    kv = max(abs(js['expectations'][sk_][c[0] + ',' + c[2]]) for c in kept)
    print(f"      {tag}: truly's cells T1 lifts {lifted} read at most {lv:.3f}; "
          f"XY+YX's cells silent under both {kept} read up to {kv:.3f}")
    check(f"G4 {tag}: the cells T1 should lift read no larger than the cells that stay silent", lv < kv)
# the comparator: how large Z + T1 makes those four cells (continuous, γ_Z = 0.1, t = 0.8, the runs' |+−+⟩)
sigma_shot = 1 / np.sqrt(4096)
lifts = []
for gt1 in np.linspace(0.0, 0.1, 11):
    L = lindbladian_z_plus_t1(Htr, [0.1] * N, [gt1] * N)
    r = (expm(L * 0.8) @ RHO0.reshape(-1)).reshape(8, 8)
    lifts.append(max(abs(np.trace(MAT[c] @ r).real) for c in ('XIY', 'XIZ', 'YIX', 'ZIX')))
check(f"G4 for γ_T1 from 0 to γ_Z = 0.1, Z + T1 lifts those cells monotonically to at most {max(lifts):.4f}, below the shot noise 1/√4096 = {sigma_shot:.4f}",
      all(b >= a for a, b in zip(lifts, lifts[1:])) and max(lifts) < sigma_shot, f"at 0.05: {lifts[5]:.4f}")

print()
print("G5  the frames and flips the prose names")
S1 = np.diag([1, 1j])
S3 = kron(S1, S1, S1)
Hyz = H_pair(('YZ', 'ZY')); Hxz = H_pair(('XZ', 'ZX'))
check("G5 S on every site carries YZ+ZY onto −(XZ+ZX), a Π²-odd pair (letters, exact)",
      np.array_equal(S3 @ Hyz @ S3.conj().T, -Hxz) and all(pi2(t) == 1 for t in ('XZ', 'ZX')))
check("G5 ... and S keeps the dephasing and the decay (it is diagonal: S Z S† = Z, S σ⁻ S† = −i·σ⁻)",
      np.array_equal(S1 @ P['Z'] @ S1.conj().T, P['Z']) and np.array_equal(S1 @ SM @ S1.conj().T, -1j * SM))
from framework.diagnostics.f83_anti_fraction import predict_pi_decomposition_anti_fraction as anti_fraction
check("G5 ... so the F83 anti-fraction moves with that frame: 0 for YZ+ZY, 1/2 for XZ+ZX (fw closed form)",
      anti_fraction(chain, [('Y', 'Z'), ('Z', 'Y')]) == 0.0 and anti_fraction(chain, [('X', 'Z'), ('Z', 'X')]) == 0.5)
oddz = {s for s in STRINGS if s.count('Z') % 2}
check("G5 XZ+ZX: the 28 strings with an odd number of Z (which X^⊗3·K reverses) are 28 of its 30 silent ones",
      len(oddz) == 28 and oddz <= rows[('XZ', 'ZX')][0], f"the rest {sorted(rows[('XZ', 'ZX')][0] - oddz)}")
check("G5 XZ+ZX from |+−+⟩: X^⊗3·K is a symmetry (X^⊗3 H* X^⊗3 = −H, start kept), a flip, and T1 leaves at most 1 string",
      np.array_equal(X3 @ Hxz.conj() @ X3, -Hxz) and np.array_equal(X3 @ RHO0.conj() @ X3, RHO0)
      and len(rows[('XZ', 'ZX')][0]) == 30 and len(rows[('XZ', 'ZX')][1]) <= 1)

print()
print(f"{len(FAILS)} failure(s)" + (': ' + '; '.join(FAILS) if FAILS else ', ALL PASS'))
sys.exit(1 if FAILS else 0)
