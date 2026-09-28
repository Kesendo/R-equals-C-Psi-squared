#!/usr/bin/env python3
"""F5's depolarizing error is a lower bound for every H, and exact when ad_H has an eigenvector
among the operators traceless on every site.

Setting (F5's class): L(rho) = -i[H, rho] + sum_l eta_l sum_{P in X,Y,Z} (P_l rho P_l - rho),
isotropic per site with total budget Gamma_l = 3 eta_l > 0 at EVERY site, sigma = sum_l Gamma_l.
The palindrome pairs the steady state (lambda = 0) with a partner at -2 sigma; F5's error is the
RATE shortfall 2 sigma - r_max, r_max = max over the spectrum of -Re(lambda).

The dissipator is diagonal on Pauli strings, -4 * (sum of eta_l over the sites where the string
is not the identity), and the Hamiltonian part is skew-adjoint in the Hilbert-Schmidt product. So
Re(lambda) = <v, L_D v>/<v, v> for an eigenvector v (the absorption theorem's Rayleigh reading),
r_max <= 4 sum eta_l = (4/3) sigma, and the error is at least (2/3) sigma for EVERY H. Equality
needs an eigenvector at the extreme of the Hermitian part: by the window-edge lemma of
PROOF_CODIM1_BY_ADDITIVITY section 6, v then lies in the bottom eigenspace of L_D, which with every
Gamma_l > 0 is E, the span of the strings with no identity letter (the operators traceless on every
site), and -i[H, v] is a multiple of v. Conversely such a v is an eigenvector of L with real part
-(4/3) sigma. Hence

    error = (2/3) sigma   iff   ad_H has an eigenvector in E,

and error > (2/3) sigma otherwise. A global Pauli string with no identity letter that commutes
with H (Z^N for every chain that conserves the excitation number or its parity, X^N for Ising in an
X field) is such an eigenvector; so is |0...0><1...1| whenever |0...0> and |1...1> are eigenvectors
of H, at ad_H eigenvalue E_0 - E_1, with no Pauli symmetry at all. The criterion governs the rate
shortfall only; the repo's other depolarizing metrics (a complex pairing distance, a best-pairing
error) are different numbers and are not decided here.

Stages (all must pass; prints "ALL STAGES PASS"):
  A  EXACT. For the Hamiltonians in the repo's depolarizing measurements (Heisenberg, XY, XX, Ising,
     XXZ at Delta = 2, DM, Heisenberg + DM) and Ising in an X field, at N = 2, 3, 4, a global string
     P with no identity letter commutes with H, entry by entry (Pauli strings are monomial matrices,
     so both products are the same few numbers placed in the same cells). At one non-uniform rational
     profile L(X^3) = -(4/3) sigma X^3 is checked by applying L. Control: a Z field on one site of the
     Ising + X chain leaves no commuting string.
  B  The criterion against the spectrum, float. Random 2-local chains in three classes (generic,
     even under X^N, even under Z^N) and a LADDER class (H diagonal on |0...0> and |1...1> at
     distinct energies, a generic block on the rest, so no global string without an identity
     letter commutes and the witnesses |0...0><1...1| and its adjoint sit at mu = +-(E_0 - E_1) != 0; at N = 3, 4, since at N = 2 the ladder commutes with ZZ), the chains
     at N = 2, 3, random positive site budgets. The largest dimension of an
     ad_H eigenspace intersected with E is read by singular values, r_max / ((4/3) sigma) by an
     eigensolver; the criterion must say "attained" exactly on the three non-generic classes, and
     the spectrum must agree with the class, which the decade check below reads. The two float readings are asserted to separate by at least
     six decades; the measured separations are printed.
  C  The size off the locus, read, not gated: error / ((2/3) sigma) on the generic rows.

Run:  python simulations/f5_depolarizing_attainment.py
   >  simulations/results/f5_depolarizing_attainment.txt
"""
import itertools
import sys
from fractions import Fraction

import numpy as np

FAIL = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


I2 = np.array([[1, 0], [0, 1]], dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
LET = {"I": I2, "X": X, "Y": Y, "Z": Z}


def kron_all(mats):
    out = np.array([[1]], dtype=complex)
    for m in mats:
        out = np.kron(out, m)
    return out


def site(N, s, P):
    return kron_all([P if k == s else I2 for k in range(N)])


def bond(N, i, P, Q):
    return site(N, i, P) @ site(N, i + 1, Q)


def chain(N, words, fields=()):
    """words: (coefficient, letter pair) on every bond; fields: (coefficient, letter, site)."""
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for i in range(N - 1):
        for c, (a, b) in words:
            H = H + float(c) * bond(N, i, LET[a], LET[b])
    for c, a, s in fields:
        H = H + float(c) * site(N, s, LET[a])
    return H


def measured_chains(N):
    return {
        "Heisenberg": chain(N, [(1, "XX"), (1, "YY"), (1, "ZZ")]),
        "XY": chain(N, [(1, "XX"), (1, "YY")]),
        "XX": chain(N, [(1, "XX")]),
        "Ising": chain(N, [(1, "ZZ")]),
        "XXZ Delta=2": chain(N, [(1, "XX"), (1, "YY"), (2, "ZZ")]),
        "DM": chain(N, [(1, "XY"), (-1, "YX")]),
        "Heisenberg+DM": chain(N, [(1, "XX"), (1, "YY"), (1, "ZZ"), (0.3, "XY"), (-0.3, "YX")]),
        "Ising + X field": chain(N, [(1, "ZZ")], [(0.7, "X", s) for s in range(N)]),
    }


def commuting_string(N, H):
    for t in itertools.product("XYZ", repeat=N):
        P = kron_all([LET[c] for c in t])
        if np.array_equal(H @ P, P @ H):
            return "".join(t)
    return None


def stage_a():
    print("## Stage A: exact, a commuting global string with no identity letter")
    for N in (2, 3, 4):
        for name, H in measured_chains(N).items():
            s = commuting_string(N, H)
            check(f"N={N} {name}: a global string with no identity letter commutes with H", s is not None,
                  f"P = {s}")
    # apply L to X^3 at one non-uniform rational profile, scaled by 30 so every entry is an integer
    N = 3
    eta30 = [3, 7, 11]                     # 30 * eta_l for Gamma = (3, 7, 11)/10
    H = chain(N, [(1, "ZZ")])
    P = kron_all([X, X, X])
    LP30 = -30j * (H @ P - P @ H) + sum(eta30[l] * (site(N, l, LET[q]) @ P @ site(N, l, LET[q]) - P)
                                         for l in range(N) for q in "XYZ")
    field_commutes = not sum(site(N, l, X) @ P - P @ site(N, l, X) for l in range(N)).any()
    check("N=3 Ising + X field at Gamma = (3, 7, 11)/10: the field commutes with X^3 and, applying L with the "
          "ZZ bonds, 30 L(X^3) = -84 X^3 entry by entry, i.e. L(X^3) = -(4/3) sigma X^3", np.array_equal(LP30, -84 * P) and field_commutes,
          "-(4/3) sigma = -14/5")
    Hc = measured_chains(N)["Ising + X field"] + site(N, 0, Z)
    check("control: with a Z field on site 0 added, no global string with no identity letter commutes",
          commuting_string(N, Hc) is None)


def liouvillian(N, H, Gam):
    d = 2 ** N
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for s in range(N):
        for P in (X, Y, Z):
            A = site(N, s, P)
            L = L + (Gam[s] / 3) * (np.kron(A, A.conj()) - np.kron(Id, Id))
    return L


def criterion(N, H):
    """largest dim of (an eigenspace of ad_H) intersected with E, and all singular values read."""
    d = 2 ** N
    E = np.array([kron_all([LET[c] for c in t]).reshape(-1) for t in itertools.product("XYZ", repeat=N)]).T / np.sqrt(d)
    w, V = np.linalg.eigh(H)
    groups = []
    for a in range(d):
        for b in range(d):
            mu = w[a] - w[b]
            for g in groups:
                if abs(g[0] - mu) < 1e-8:
                    g[1].append((a, b))
                    break
            else:
                groups.append([mu, [(a, b)]])
    best, svals = 0, []
    for mu, pairs in groups:
        K = np.array([np.outer(V[:, a], V[:, b].conj()).reshape(-1) for a, b in pairs]).T
        s = np.linalg.svd(K.conj().T @ E, compute_uv=False)
        svals.extend(s.tolist())
        best = max(best, int(np.sum(s > 1 - 1e-6)))
    return best, svals


def random_chain(rng, N, allowed):
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for i in range(N - 1):
        for a in "XYZ":
            for b in "XYZ":
                if allowed(a + b):
                    H = H + rng.normal() * bond(N, i, LET[a], LET[b])
    for i in range(N):
        for a in "XYZ":
            if allowed(a):
                H = H + rng.normal() * site(N, i, LET[a])
    return H


def ladder(rng, N):
    """|0...0> and |1...1> eigenvectors at distinct energies, a generic Hermitian block elsewhere."""
    d = 2 ** N
    A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    H = (A + A.conj().T) / 2
    for k in (0, d - 1):
        H[k, :] = 0
        H[:, k] = 0
    H[0, 0], H[d - 1, d - 1] = 0.7, -1.3
    return H


def stage_b_c():
    print()
    print("## Stage B: the criterion against the spectrum; Stage C: the size off the locus")
    rng = np.random.default_rng(20260928)
    classes = {
        "generic": lambda N: random_chain(rng, N, lambda w: True),
        "even under X^N": lambda N: random_chain(rng, N, lambda w: sum(c in "YZ" for c in w) % 2 == 0),
        "even under Z^N": lambda N: random_chain(rng, N, lambda w: sum(c in "XY" for c in w) % 2 == 0),
        "ladder |0..0><1..1|": lambda N: ladder(rng, N),
    }
    attained_dev, generic_gap, sv_one, sv_rest = [], [], [], []
    agree, ladder_no_string = True, True
    # the chains at N = 2, 3; the ladder at N = 3, 4 (at N = 2 it commutes with ZZ, so it tests nothing)
    plan = [(N, name) for N in (2, 3) for name in classes if not name.startswith("ladder")]
    plan += [(N, name) for N in (3, 4) for name in classes if name.startswith("ladder")]
    for N, name in plan:
            make = classes[name]
            for _ in range(4):
                H = make(N)
                Gam = rng.uniform(0.3, 2.0, N)
                sigma = Gam.sum()
                ev = np.linalg.eigvals(liouvillian(N, H, Gam))
                ratio = -ev.real.min() / (4 / 3 * sigma)
                dim, svals = criterion(N, H)
                sv_one += [s for s in svals if abs(1 - s) < 1e-6]
                sv_rest += [s for s in svals if not abs(1 - s) < 1e-6]
                if name == "generic":
                    agree &= dim == 0
                    generic_gap.append(1 - ratio)
                    print(f"    C  N={N} generic: error / ((2/3) sigma) = {(2 * sigma + ev.real.min()) / (2 / 3 * sigma):.6f}")
                else:
                    agree &= dim > 0
                    attained_dev.append(abs(1 - ratio))
                if name.startswith("ladder"):
                    ladder_no_string &= commuting_string(N, H) is None
    check("the criterion says 'not attained' on the generic class and 'attained' on the two symmetric "
          "classes and on the ladder (the spectrum's side is the decade check below)", agree)
    check("the ladder rows commute with no global string without an identity letter, so their witness is "
          "not such a symmetry",
          ladder_no_string)
    worst_att, best_gen = max(attained_dev), min(generic_gap)
    dec_spec = np.log10(best_gen / max(worst_att, 1e-300))
    check("spectral readings separate by at least six decades", dec_spec >= 6,
          f"worst attained |1 - ratio| {worst_att:.1e}, smallest generic shortfall {best_gen:.1e}, "
          f"{dec_spec:.1f} decades measured")
    gap_sv = min(abs(1 - s) for s in sv_rest)
    worst_one = max(abs(1 - s) for s in sv_one)
    dec_sv = np.log10(gap_sv / max(worst_one, 1e-300))
    check("singular-value readings separate by at least six decades", dec_sv >= 6,
          f"worst |1 - s| counted {worst_one:.1e}, nearest uncounted {gap_sv:.1e}, {dec_sv:.1f} decades measured")


if __name__ == "__main__":
    stage_a()
    stage_b_c()
    print()
    if FAIL:
        print(f"{len(FAIL)} FAILURE(S): {FAIL}")
        sys.exit(1)
    print("ALL STAGES PASS")
