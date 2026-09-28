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
  D  Off the isotropic class: an anisotropic Pauli channel, dephasing along X, Y, Z at rates
     gamma_P^l of their own on every site. A letter P at site l decays at 2 (sigma_l - gamma_P^l),
     sigma_l the site's rate sum, so the fastest string decays at 2 sum_l (sigma_l - min_P gamma_P^l)
     and the shortfall is at least 2 sum_l min_P gamma_P^l. The same window-edge argument gives
     equality exactly when ad_H has an eigenvector in the span E_fast of the strings whose letter at
     each site is one of that site's fastest letters (the identity letter among them at a site with
     every rate zero, so no positivity is needed). The fastest rate is also read off the dissipator
     exactly (Q S Q = +-S entry by entry, rates as Fractions). Rational profiles, N = 2, 3:
     Heisenberg and XXZ at Delta = 2 with one smallest letter P at every site (P^N commutes, attained),
     Heisenberg with a unique smallest letter at each site, differing between sites (the fast span is
     one mixed string, and Heisenberg commutes with a string without identity letters only when its
     letters agree: not attained), generic 2-local chains, Heisenberg with X and Y tied smallest at one
     site (attained through X^N), zero minima (Z dephasing alone, fast span {X, Y} per site, and X plus Z
     dephasing with Y dark, fast span Y^N: both attained at bound 0 on Heisenberg, Z dephasing alone not
     attained on generic chains), and Heisenberg with one dark site (attained through P^N). The criterion
     must agree with the predicted rows, and the criterion and the spectrum must agree on every row,
     read by the same six-decade separations.

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


def criterion(N, H, letters=None):
    """largest dim of (an eigenspace of ad_H) intersected with E, and all singular values read.
    E is spanned by the strings whose letter at site l lies in letters[l] (default "XYZ" at every site)."""
    d = 2 ** N
    letters = letters or ["XYZ"] * N
    E = np.array([kron_all([LET[c] for c in t]).reshape(-1) for t in itertools.product(*letters)]).T / np.sqrt(d)
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


def pauli_liouvillian(N, H, rates):
    """rates[l] = {"X": g, "Y": g, "Z": g}: dephasing along each letter at its own rate, per site."""
    d = 2 ** N
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for s in range(N):
        for q in "XYZ":
            A = site(N, s, LET[q])
            L = L + rates[s][q] * (np.kron(A, A.conj()) - np.kron(Id, Id))
    return L


def fastest_letters(r):
    """the letters (identity included) at the site's largest decay rate 2 * (sum of the rates of the
    letters anticommuting with it), and that rate; r holds exact Fractions so ties are exact."""
    tot = r["X"] + r["Y"] + r["Z"]
    rate = {"I": Fraction(0), **{q: 2 * (tot - r[q]) for q in "XYZ"}}
    top = max(rate.values())
    return "".join(q for q in "IXYZ" if rate[q] == top), top


def stage_d():
    print()
    print("## Stage D: an anisotropic Pauli channel, the bound 2 sum_l min_P gamma_P^l and its criterion")
    rng = np.random.default_rng(20260929)

    def profile(N, argmin):
        """rational rates in tenths; argmin[l] is the site's strictly smallest letter, "" for a dark site,
        or an explicit {letter: Fraction} site profile (ties and zero rates)."""
        out = []
        for a in argmin:
            if isinstance(a, dict):
                out.append(dict(a))
                continue
            if a == "":
                out.append({q: Fraction(0) for q in "XYZ"})
                continue
            hi = rng.choice(np.arange(3, 12), size=2, replace=False)
            others = [q for q in "XYZ" if q != a]
            r = {a: Fraction(int(rng.integers(1, 3)), 10)}
            r.update({others[0]: Fraction(int(hi[0]), 10), others[1]: Fraction(int(hi[1]), 10)})
            out.append(r)
        return out

    heis = lambda N: chain(N, [(1, "XX"), (1, "YY"), (1, "ZZ")])
    classes = {
        # predicted attained: P^N commutes with H and is the fast span's only string
        "Heisenberg, one smallest letter everywhere": (heis, lambda N: [str(rng.choice(list("XYZ")))] * N, True),
        "XXZ Delta=2, one smallest letter everywhere": (lambda N: chain(N, [(1, "XX"), (1, "YY"), (2, "ZZ")]),
                                                       lambda N: [str(rng.choice(list("XYZ")))] * N, True),
        # predicted not attained: the fast span is one string mixing letters, and Heisenberg commutes with
        # a string without identity letters only when all its letters agree
        "Heisenberg, smallest letters differ": (heis, lambda N: ["X", "Z"] + [str(rng.choice(list("XYZ"))) for _ in range(N - 2)], False),
        "generic 2-local": (lambda N: random_chain(rng, N, lambda w: True),
                            lambda N: [str(rng.choice(list("XYZ"))) for _ in range(N)], False),
        # a site with every rate zero puts the identity letter into the fast span (no positivity needed)
        # P^N lies in the span (P is among the dark site's four letters) and commutes with Heisenberg
        "Heisenberg, site 0 dark, one smallest letter elsewhere": (heis, lambda N: [""] + [str(rng.choice(list("XYZ")))] * (N - 1), True),
        # a tie: X and Y smallest at site 0, X smallest elsewhere; the span {X, Y} (x) X... holds X^N
        "Heisenberg, X and Y tied at site 0, X elsewhere": (
            heis, lambda N: [{"X": Fraction(1, 10), "Y": Fraction(1, 10), "Z": Fraction(int(rng.integers(3, 12)), 10)}]
            + ["X"] * (N - 1), True),
        # zero minima, the bound is 0: Z dephasing alone (fast span {X, Y} per site, X^N commutes with
        # Heisenberg) and X plus Z dephasing with Y dark (fast span Y^N), against a generic chain
        "Heisenberg, Z dephasing only": (heis, lambda N: [{"X": Fraction(0), "Y": Fraction(0),
                                                            "Z": Fraction(int(rng.integers(1, 12)), 10)} for _ in range(N)], True),
        "Heisenberg, X and Z dephasing, Y dark": (heis, lambda N: [{"X": Fraction(int(rng.integers(1, 12)), 10), "Y": Fraction(0),
                                                                     "Z": Fraction(int(rng.integers(1, 12)), 10)} for _ in range(N)], True),
        "generic 2-local, Z dephasing only": (lambda N: random_chain(rng, N, lambda w: True),
                                              lambda N: [{"X": Fraction(0), "Y": Fraction(0),
                                                          "Z": Fraction(int(rng.integers(1, 12)), 10)} for _ in range(N)], False),
    }
    agree, predicted, att_dev, gap, sv_one, sv_rest = True, True, [], [], [], []
    for N in (2, 3):
        for name, (make_h, make_arg, expect) in classes.items():
            for _ in range(4):
                H = make_h(N)
                rates = profile(N, make_arg(N))
                fast = [fastest_letters(r) for r in rates]
                sigma = float(sum(r["X"] + r["Y"] + r["Z"] for r in rates))
                bound = float(sum(2 * min(r.values()) for r in rates))
                ev = np.linalg.eigvals(pauli_liouvillian(N, H, {l: {q: float(v) for q, v in r.items()}
                                                                  for l, r in enumerate(rates)}))
                shortfall = 2 * sigma + ev.real.min()
                dim, svals = criterion(N, H, [f[0] for f in fast])
                sv_one += [s for s in svals if abs(1 - s) < 1e-6]
                sv_rest += [s for s in svals if not abs(1 - s) < 1e-6]
                rel = (shortfall - bound) / sigma
                if dim > 0:
                    att_dev.append(abs(rel))
                else:
                    gap.append(rel)
                # the fastest string read off the dissipator itself, exactly: Q_l S Q_l = +-S entry by entry,
                # so L_D(S) = -(rate) S with rate = sum of 2 gamma_q^l over the (l, q) where the sign is -1
                read = Fraction(-1)
                for t in itertools.product("IXYZ", repeat=N):
                    S = kron_all([LET[c] for c in t])
                    rate = Fraction(0)
                    for l, r in enumerate(rates):
                        for q in "XYZ":
                            QSQ = site(N, l, LET[q]) @ S @ site(N, l, LET[q])
                            assert np.array_equal(QSQ, S) or np.array_equal(QSQ, -S)
                            if np.array_equal(QSQ, -S):
                                rate += 2 * r[q]
                    read = max(read, rate)
                top_exact = sum(f[1] for f in fast)
                bound_exact = sum(2 * min(r.values()) for r in rates)
                sigma_exact = sum(r["X"] + r["Y"] + r["Z"] for r in rates)
                agree &= read == top_exact and 2 * sigma_exact - read == bound_exact
                if expect is not None:
                    predicted &= (dim > 0) == expect
                print(f"    D  N={N} {name}: fast letters {'/'.join(f[0] for f in fast)}, "
                      f"(shortfall - bound) / sigma = {rel:+.3e}, criterion dim = {dim}")
    check("the fastest string read off the dissipator's diagonal decays at 2 sum_l (sigma_l - min_P gamma_P^l), "
          "so the bound is 2 sum_l min_P gamma_P^l", agree)
    check("the criterion says 'attained' on the one-letter Heisenberg and XXZ rows, the tied row, the "
          "zero-minimum Heisenberg rows and the dark-site rows, and 'not attained' on the mixed-letter Heisenberg and generic rows", predicted)
    worst_att = max(att_dev) if att_dev else 0.0
    dec = np.log10(min(gap) / max(worst_att, 1e-300))
    check("where the criterion says attained the shortfall is the bound, where it says not the shortfall "
          "exceeds it, separated by at least six decades", min(gap) > 0 and dec >= 6,
          f"worst attained |shortfall - bound| / sigma {worst_att:.1e}, smallest excess {min(gap):.1e}, {dec:.1f} decades")
    dec_sv = np.log10(min(abs(1 - s) for s in sv_rest) / max(max(abs(1 - s) for s in sv_one), 1e-300))
    check("singular-value readings separate by at least six decades", dec_sv >= 6, f"{dec_sv:.1f} decades measured")


if __name__ == "__main__":
    stage_a()
    stage_b_c()
    stage_d()
    print()
    if FAIL:
        print(f"{len(FAIL)} FAILURE(S): {FAIL}")
        sys.exit(1)
    print("ALL STAGES PASS")
