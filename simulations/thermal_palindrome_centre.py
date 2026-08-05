"""The palindrome centre under a thermal bath: F137 extended to n_bar > 0.

WHAT THE REPO ALREADY HOLDS (the Stage-0 sweep, before anything here was run)

  F137 (docs/ANALYTICAL_FORMULAS.md, 2026-07-21, Tier1Candidate): pure
  amplitude damping, jump sigma^- per site at rates gamma_i, leaves the
  Liouvillian spectrum an exact palindrome centred at Re lambda = -Sum(gamma_i)/2,
  HALF the dephasing centre -Sum(gamma_i). Its own lesson: "Pi fails" and "the
  palindrome fails" are different claims. Derived for H = 0 by tensor sum,
  measured for H != 0 (18/18 configurations, N=2-5).

  docs/proofs/MIRROR_SYMMETRY_PROOF.md sec. "T1 alone": the same, plus the
  composition rule. Transverse dephasing composes with T1 (X 64/64, Y 64/64);
  CO-AXIAL Z-dephasing breaks it, down to 8/64.

  docs/KMS_DETAILED_BALANCE.md: the per-site Pauli rates of a THERMAL bath are
  [0, r/2, r/2, r] with r = gamma(2*n_bar + 1), and "(0, r) and (r/2, r/2) both
  sum to r. So rate-pairing IS possible in principle."

  F84: the Pi-conjugation violation of M is governed by the NET rate
  dgamma = gamma_down - gamma_up, and vanishes at detailed balance.

WHAT WAS MISSING, and what this script measures

  F137 covers sigma^- only, i.e. n_bar = 0. KMS has the thermal per-site rates
  but never writes them as a centre. Neither cites the other, and
  experiments/THERMAL_BREAKING.md has called the question open since
  2026-03-30. Joining them predicts, for a bath with per-site cooling rate
  gamma_down_l and heating rate gamma_up_l:

      centre = -Sum_l (gamma_down_l + gamma_up_l) / 2

  because KMS's r is exactly gamma_down + gamma_up, and F137's tensor-sum
  argument then settles H = 0. KMS asked for exactly this: "Numerical test
  recommended: compute the Liouvillian spectrum for XXZ + thermal bath at
  various temperatures and check whether any palindromic structure survives."
  This is that test, and it goes further than numerics where it can.

TWO THINGS THIS SCRIPT IS CAREFUL NOT TO CONFUSE

  THE CENTRE IS AN IDENTITY, NOT EVIDENCE. If a multiset is closed under
  lambda -> 2c - lambda then each pair sums to 2c, so sum(lambda) = n*c and
  c = mean(lambda) exactly. The mean is trace(L)/dim, and the commutator part
  of L is traceless, so the centre is fixed by the DISSIPATOR alone and does
  not depend on H at all. It therefore needs no proof at H != 0 and carries no
  information about whether a palindrome exists: a broken spectrum has the same
  well-defined centre (section 6). What the identity does buy is that there is
  exactly ONE candidate centre, so a large distance at it means NO centre works.
  That is what turns "broken" into a measurement rather than a failed search,
  and it retires THERMAL_BREAKING's caution that the thermal centre was only
  "estimated".

  THE PAIRING IS THE CLAIM, and where the rates are rational it is settled
  EXACTLY (section 2): the characteristic polynomial over Q(i) satisfies
  p(2c - x) == p(x) identically, with no eigensolver and no tolerance. Section 3
  then reaches the N the exact route cannot afford, using the canonical F1
  distance (fw.max_f1_pairing_distance) in units of eps * spectral radius. That
  floor grows with N and the greedy matcher adds an order-dependence spread on
  top, so read a row against its own N and do not rank rows against each other.

THE SHARP PART is which rate appears. The palindrome centre follows the SUM
gamma_down + gamma_up; F84's Pi-conjugation violation follows the NET
difference. Section 4 shows why, and it is not a coincidence of two formulas:
in the (I, X, Y, Z) Pauli basis the one-site thermal dissipator is TRIANGULAR,
with the sum on the diagonal and the net rate as its only off-diagonal entry,
where it cannot move an eigenvalue. The repo's name for that general fact is
the SOFT BREAK (reflections/ON_THE_SOFT_BREAK.md, 2026-04-25); the entry itself
was tabulated in PROOF_F82 Step 3 (2026-04-30) and the sum/difference split is
the single-qubit Bloch equation in experiments/F81_VIOLATION_HARDWARE_BRIDGE.md
(2026-07-05). Only the join is new. That also resolves KMS's summary table
against its own body: the obstruction it names is against Pi, not the pairing.

Run: python simulations/thermal_palindrome_centre.py   (~2 min, the exact
section dominates; writes simulations/results/thermal_palindrome_centre.txt)

Depends on sympy for section 2 (exact rational characteristic polynomials).
"""

import sys
from pathlib import Path

import numpy as np
from fractions import Fraction

sys.path.insert(0, str(Path(__file__).parent))
from framework import max_f1_pairing_distance  # noqa: E402

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
# sigma^- as the framework defines it: |1> -> |0>, the ground-state convention.
# The other sign convention flips gamma_down and gamma_up, which swaps cooling
# and heating; every result here reads only their SUM or the magnitude of their
# difference, so it is convention-independent.
SM = np.array([[0, 1], [0, 0]], dtype=complex)
SP = SM.conj().T


def at(op, site, n):
    """Single-site operator placed at `site` of an n-qubit register."""
    m = np.eye(1, dtype=complex)
    for s in range(n):
        m = np.kron(m, op if s == site else I2)
    return m


def heisenberg(n, J=1.0):
    d = 2 ** n
    H = np.zeros((d, d), dtype=complex)
    for b in range(n - 1):
        for P in (X, Y, Z):
            H += J * at(P, b, n) @ at(P, b + 1, n)
    return H


def liouvillian(H, n, gz=None, g_down=None, g_up=None, gx=None, gy=None):
    """Column-stacked Liouvillian with any subset of the five local channels."""
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))

    def add(L, Lk):
        A = Lk.conj().T @ Lk
        return (L + np.kron(Lk.conj(), Lk)
                - 0.5 * np.kron(Id, A) - 0.5 * np.kron(A.T, Id))

    for k in range(n):
        for rates, op in ((gz, Z), (g_down, SM), (g_up, SP), (gx, X), (gy, Y)):
            if rates is not None and rates[k] > 0:
                L = add(L, np.sqrt(rates[k]) * at(op, k, n))
    return L


def trace_centre(evals):
    """The ONLY possible palindrome centre of this spectrum (see module docstring)."""
    return float(np.mean(evals).real)


def distance_at(evals, centre):
    """Canonical F1 distance about `centre`, in eps * spectral radius."""
    eps = np.finfo(float).eps
    rho = float(np.max(np.abs(evals)))
    return max_f1_pairing_distance(evals, -centre) / (eps * rho)


def exact_liouvillian(N, J, g_down, g_up, gz=None):
    """The same Liouvillian over the RATIONALS, for an exact palindrome test.

    Every rate enters bilinearly, so with rational rates and rational J the
    whole matrix is rational over Q(i) and its characteristic polynomial is
    exact. That turns "is the spectrum palindromic about c?" into the exact
    identity p(2c - x) == p(x): no eigensolver, no tolerance. Case 1 of the
    house rule, an exact route exists, so compare exactly.
    """
    import sympy as sp

    Id = sp.eye(2 ** N)
    Xs = sp.Matrix([[0, 1], [1, 0]])
    Ys = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    Zs = sp.Matrix([[1, 0], [0, -1]])
    SMs = sp.Matrix([[0, 1], [0, 0]])
    SPs = SMs.T

    def site(op, k):
        m = sp.Matrix([[1]])
        for t in range(N):
            m = sp.Matrix(sp.kronecker_product(m, op if t == k else sp.eye(2)))
        return m

    H = sp.zeros(2 ** N)
    for b in range(N - 1):
        for P in (Xs, Ys, Zs):
            H += sp.Rational(J) * site(P, b) * site(P, b + 1)
    L = -sp.I * (sp.Matrix(sp.kronecker_product(Id, H))
                 - sp.Matrix(sp.kronecker_product(H.T, Id)))
    for k in range(N):
        for rates, op in ((gz, Zs), (g_down, SMs), (g_up, SPs)):
            if rates is None or rates[k] == 0:
                continue
            g = sp.Rational(rates[k])
            A = site(op, k)
            AdA = A.H * A * g
            L = (L + g * sp.Matrix(sp.kronecker_product(A.conjugate(), A))
                 - sp.Rational(1, 2) * sp.Matrix(sp.kronecker_product(Id, AdA))
                 - sp.Rational(1, 2) * sp.Matrix(sp.kronecker_product(AdA.T, Id)))
    return L


def exact_palindrome_verdict(N, J, g_down, g_up, gz=None):
    """EXACT: is p(2c - x) equal to p(x)? Returns (verdict, centre)."""
    import sympy as sp

    x = sp.symbols('x')
    L = exact_liouvillian(N, J, g_down, g_up, gz)
    p = sp.expand(L.charpoly(x).as_expr())
    c = (-(sum(g_down) + (sum(g_up) if g_up else 0)) / 2
         - (sum(gz) if gz else 0))
    q = sp.expand(p.subs(x, 2 * c - x))
    if sp.simplify(sp.expand(q - p)) == 0:
        return "EXACT palindrome, p(2c-x) == p(x)", c
    if sp.simplify(sp.expand(q + p)) == 0:
        return "EXACT anti-palindrome, p(2c-x) == -p(x)", c
    return "BROKEN, neither identity holds", c


def single_site_dissipator_pauli(g_down, g_up):
    """The one-site thermal dissipator as a 4x4 matrix on (I, X, Y, Z)."""
    basis = [I2 / np.sqrt(2), X / np.sqrt(2), Y / np.sqrt(2), Z / np.sqrt(2)]
    D = np.zeros((4, 4), dtype=complex)
    for j, B in enumerate(basis):
        img = np.zeros((2, 2), dtype=complex)
        for g, A in ((g_down, SM), (g_up, SP)):
            img += g * (A @ B @ A.conj().T
                        - 0.5 * (A.conj().T @ A @ B + B @ A.conj().T @ A))
        for i, C in enumerate(basis):
            D[i, j] = np.trace(C.conj().T @ img)
    return D


def main():
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    out_path = results_dir / "thermal_palindrome_centre.txt"

    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    rng = np.random.default_rng(20260805)

    out("=== THE PALINDROME CENTRE UNDER A THERMAL BATH ===")
    out("F137 extended from sigma^- to sigma^- + sigma^+.")
    out()

    out("--- 1. the centre is an IDENTITY, not evidence ---")
    out("    trace(L)/dim is the mean eigenvalue, hence the only possible centre,")
    out("    and the commutator part of L is traceless: the centre is fixed by the")
    out("    dissipator alone and does not depend on H. Shown by varying J at fixed")
    out("    rates. The residual below is float arithmetic on a sum, NOT a reading")
    out("    about any symmetry; a broken spectrum has the same centre (section 6).")
    out(f"    {'case':30s} {'predicted':>16s} {'trace/dim':>16s} {'residual':>11s}")
    N = 3
    gd = np.array([0.30, 0.20, 0.10])
    gu = np.array([0.10, 0.25, 0.20])
    gz = np.array([0.10, 0.10, 0.10])
    for J in (0.0, 1.0, 5.0):
        H = heisenberg(N, J)
        for tag, kw, extra in (("thermal", dict(g_down=gd, g_up=gu), 0.0),
                               ("thermal + Z", dict(g_down=gd, g_up=gu, gz=gz),
                                float(sum(gz)))):
            L = liouvillian(H, N, **kw)
            pred = -float(sum(gd) + sum(gu)) / 2 - extra
            meas = float(np.trace(L).real / L.shape[0])
            out(f"    J={J:<4g} {tag:24s} {pred:+16.10f} {meas:+16.10f}"
                f" {abs(pred - meas):11.2e}")

    out()
    out("--- 2. EXACT, no eigensolver: p(2c - x) == p(x) over the rationals ---")
    out("    Rational rates and rational J make L rational over Q(i), so the")
    out("    characteristic polynomial is exact and the palindrome becomes a")
    out("    polynomial identity instead of a small number.")
    exact_cases = [
        ("N=2 thermal", 2, 1,
         [Fraction(3, 10), Fraction(1, 5)], [Fraction(1, 10), Fraction(1, 4)], None),
        ("N=2 cooling only (F137)", 2, 1,
         [Fraction(3, 10), Fraction(1, 5)], None, None),
        ("N=2 detailed balance", 2, 1,
         [Fraction(1, 5), Fraction(1, 5)], [Fraction(1, 5), Fraction(1, 5)], None),
        ("N=3 thermal", 3, 1,
         [Fraction(3, 10), Fraction(1, 5), Fraction(1, 10)],
         [Fraction(1, 10), Fraction(1, 4), Fraction(1, 5)], None),
        ("N=2 thermal + Z", 2, 1,
         [Fraction(3, 10), Fraction(1, 5)], [Fraction(1, 10), Fraction(1, 4)],
         [Fraction(1, 10), Fraction(1, 10)]),
    ]
    for tag, n, J, a, b, c_z in exact_cases:
        verdict, centre = exact_palindrome_verdict(n, J, a, b, c_z)
        out(f"    {tag:26s} centre {str(centre):>8s}   {verdict}")
    out("    At these N the thermal palindrome is PROVEN for the Heisenberg H,")
    out("    not measured, and the Z + amplitude break is proven too. The last")
    out("    row is the case that must fail; without it the section could not.")

    out()
    out("--- 3. larger N, where only the eigensolver reaches ---")
    out(f"    {'case':30s} {'centre':>16s} {'F1 dist':>12s}")
    for N in (2, 3, 4, 5):
        H = heisenberg(N)
        for trial in range(2):
            gd = rng.uniform(0.02, 0.4, N)
            gu = rng.uniform(0.00, 0.4, N)
            ev = np.linalg.eigvals(liouvillian(H, N, g_down=gd, g_up=gu))
            pred = -float(sum(gd) + sum(gu)) / 2
            out(f"    N={N} thermal, trial {trial:<12d} {pred:+16.10f}"
                f" {distance_at(ev, pred):12.2f}")
    out("    Read each row against its own N: the eigensolver floor grows with N,")
    out("    and greedy matching adds an order-dependence spread on top, so these")
    out("    are all 'at the floor' and are not comparable with each other.")

    out()
    out("--- 4. WHY the spectrum sees the sum: the net rate is triangular ---")
    D = single_site_dissipator_pauli(0.37, 0.13)
    out("    one-site thermal dissipator on the (I, X, Y, Z) Pauli basis,")
    out("    gamma_down = 0.37, gamma_up = 0.13:")
    for row in D.real:
        out("      [" + "  ".join(f"{v:+8.4f}" for v in row) + "]")
    out(f"    diagonal = [0, -r/2, -r/2, -r] with r = gd + gu = {0.37 + 0.13:.2f}")
    out(f"    the single off-diagonal entry is (Z, I) = {D[3, 0].real:+.6f}"
        f" = gd - gu, the NET rate")
    out("    In this ordering that entry sits strictly below the diagonal, so the")
    out("    matrix is TRIANGULAR and the net rate cannot move an eigenvalue.")
    out("    That is the mechanism: the spectrum sees only the sum while F84's")
    out("    Pi-conjugation violation sees only the difference. Not two readings")
    out("    of one quantity; two different entries of one matrix.")

    out()
    out("--- 5. the sum/net separation, against F84's own closed form ---")
    N = 5
    H = heisenberg(N)
    gd = np.array([0.30, 0.10, 0.20, 0.05, 0.15])
    gu = np.array([0.10, 0.10, 0.05, 0.05, 0.15])
    ev = np.linalg.eigvals(liouvillian(H, N, g_down=gd, g_up=gu))
    f84 = float(np.sqrt(sum((gd - gu) ** 2)) * 2 ** (N - 1))
    centre = -float(sum(gd) + sum(gu)) / 2
    out(f"    F84 violation ||D_odd||_F = sqrt(sum (gd-gu)^2)*2^(N-1) = {f84:.6f}")
    out(f"    palindrome centre                                      = {centre:+.6f}")
    out(f"    F1 distance there                                      = "
        f"{distance_at(ev, centre):.2f}")
    g = np.array([0.10, 0.25, 0.05, 0.20, 0.15])
    ev_b = np.linalg.eigvals(liouvillian(H, N, g_down=g, g_up=g))
    f84_b = float(np.sqrt(sum((g - g) ** 2)) * 2 ** (N - 1))
    out("    at detailed balance, gamma_down = gamma_up:")
    out(f"      F84 violation     = {f84_b:.6f}   exactly zero, its rate is the net")
    out(f"      palindrome centre = {-float(sum(g)):+.6f}   NOT zero, its rate is the sum")
    out(f"      F1 distance there = {distance_at(ev_b, -float(sum(g))):.2f}")
    out("    Section 2 settles that balanced case exactly, so this is not a")
    out("    tolerance question. F84's number here is its published closed form")
    out("    evaluated on these rates, not a re-derivation of it.")

    out()
    out("--- 6. the centre alone decides nothing ---")
    N = 4
    H = heisenberg(N)
    gd = rng.uniform(0.02, 0.3, N)
    gu = rng.uniform(0.00, 0.3, N)
    g2 = rng.uniform(0.02, 0.3, N)
    for tag, kw in (("thermal alone", {}),
                    ("thermal + Z  (co-axial)", {'gz': g2}),
                    ("thermal + X  (transverse)", {'gx': g2}),
                    ("thermal + Y  (transverse)", {'gy': g2})):
        ev = np.linalg.eigvals(liouvillian(H, N, g_down=gd, g_up=gu, **kw))
        c = trace_centre(ev)
        out(f"    {tag:26s} centre {c:+.8f}   F1 dist {distance_at(ev, c):14.4g}")
    out("    The three added-channel rows share a centre to every digit, because")
    out("    the centre is the trace and all three add the same total rate. Only")
    out("    the co-axial one breaks the pairing. So the centre must be computed")
    out("    AND the distance checked; neither decides alone. WHICH axis survives")
    out("    is MIRROR_SYMMETRY_PROOF's composition rule, measured there on a")
    out("    pairing-fraction metric with sigma^- at N=3; this is the same")
    out("    statement read with the F1 distance on a thermal bath at N=4, not a")
    out("    reproduction of those counts.")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\n>>> Results saved to: {out_path}")


if __name__ == "__main__":
    main()
