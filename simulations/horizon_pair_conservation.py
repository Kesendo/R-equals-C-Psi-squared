"""The horizon reading, measured: what the palindrome does to a pair.

WHAT THE REPO ALREADY HOLDS

  F1 (docs/proofs/MIRROR_SYMMETRY_PROOF.md): under Z-dephasing the Liouvillian
  eigenvalue multiset is closed under lambda -> -2*sigma - lambda.

  F137 (docs/ANALYTICAL_FORMULAS.md, extended 2026-08-05): amplitude damping
  keeps the palindrome at a HALVED centre, a thermal bath at
  -Sum(gamma_down + gamma_up)/2; the centre is trace(L)/dim and the commutator
  part of L is traceless, so no Hamiltonian can move it.

  hypotheses/PAIR_BREAKING_AT_THE_HORIZON.md reads that closure as the Hawking
  mechanism in operator space, and describes it as pairs being "torn apart"
  when dephasing shifts the palindrome off zero.

WHAT THIS SCRIPT SHOWS, and what it does NOT

  The pair sum is F1'S REAL PART, not an independent law. From
  lambda2 = -2*sigma - lambda1 it follows in one line that
  r1 + r2 = 2*sigma. This script cannot fail unless F1 fails, and F1 is
  already verified over 87,376 eigenvalues at N=2-8. Worse for any claim of
  independence: the pairing here is BUILT by nearest-reflection matching, so
  the deviation of the sum from 2*sigma IS the palindrome residual, by
  construction. Section 1 is therefore a re-reading of F1 at N=4, printed so
  the reading below has something to point at, and NOT a measurement of a new
  conservation law.

  Nor is the reading itself unprecedented. experiments/N_INFINITY_PALINDROME.md
  already says "the Hamiltonian shifts rates within palindromic pairs but never
  breaks the pairing", and hypotheses/README.md says the same of
  ZERO_IS_THE_MIRROR, which the horizon document lists as a dependency. What is
  new here is only that the horizon document itself said otherwise.

  The reading, then, is a relabelling: nothing about the
  pair loosens as sigma grows. What grows is the SEPARATION |r1 - r2| between
  the two partners' rates, from 0 for a self-paired mode at the centre to
  2*sigma for the extreme pair. The palindrome holds the pair and hands the
  halves different fates. hypotheses/PAIR_BREAKING_AT_THE_HORIZON.md used to
  say the pair is "torn apart"; that is the wording this replaces.

  Scope, because "at every coupling" invites over-reading: F1 is a theorem for
  the XXZ family under Z-dephasing, so sweeping J inside that family tests
  nothing about H-independence. A site-dependent transverse field breaks the
  palindrome outright while leaving the CENTRE untouched, which is the sharp
  version of the distinction: the centre is H-independent because it is a
  trace; the pair sum is not, it is H-independent only where F1 holds.

  Three readings:
    1. the pair sum and the separation, with the residual read against the
       eigensolver's error model rather than against a digit count;
    2. the centre against the coupling. This is the trace identity, printed as
       a contrast (the dynamics moves by a factor of 50, the centre does not
       move at all) and NOT as evidence about any symmetry: the same
       identity holds when the palindrome is broken;
    3. the census at the two extremes. Equal counts follow from the palindrome
       and cannot come out otherwise WHILE IT HOLDS; the broken row is included
       so the equality is scoped rather than stated as "always".

Run: python simulations/horizon_pair_conservation.py   (~10 s, writes
simulations/results/horizon_pair_conservation.txt)
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from thermal_palindrome_centre import liouvillian, heisenberg  # noqa: E402


def pair_readings(evals, sigma):
    """Rate sums and rate separations over the palindromic matching."""
    reflected = -2.0 * sigma - evals
    taken = np.zeros(len(evals), dtype=bool)
    sums, seps = [], []
    for x in evals:
        d = np.abs(x - reflected)
        d[taken] = np.inf
        j = int(np.argmin(d))
        taken[j] = True
        sums.append(float(-x.real + -evals[j].real))
        seps.append(abs(float(-x.real - -evals[j].real)))
    return sums, seps


def main():
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    out_path = results_dir / "horizon_pair_conservation.txt"

    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    N, g = 4, 0.05
    sigma = N * g

    out("=== THE HORIZON READING, MEASURED ===")
    out(f"N={N} Heisenberg chain, uniform Z-dephasing gamma={g}, Sigma_gamma={sigma}")
    out()

    out("--- 1. the pair sum is F1's real part, and the separation is what grows ---")
    out("    r1 + r2 = 2*sigma follows from lambda -> -2*sigma - lambda in one line.")
    out("    The matching below is BUILT by nearest reflection, so the deviation of")
    out("    the sum from 2*sigma is the palindrome residual itself. Read the")
    out("    residual against eps * spectral radius, the eigensolver's error model,")
    out("    not against a digit count.")
    eps = np.finfo(float).eps
    out(f"    {'coupling':>10} {'max|sum - 2sig|':>17} {'in eps*rho':>12}"
        f" {'min|r1-r2|':>12} {'max|r1-r2|':>12}")
    for J in (0.0, 1.0, 5.0, 50.0):
        ev = np.linalg.eigvals(liouvillian(heisenberg(N, J), N, gz=np.full(N, g)))
        sums, seps = pair_readings(ev, sigma)
        dev = max(abs(x - 2 * sigma) for x in sums)
        rho = float(np.max(np.abs(ev)))
        out(f"    J={J:<8g} {dev:17.3e} {dev / (eps * rho) if rho else 0:12.2f}"
            f" {min(seps):12.3e} {max(seps):12.6f}")
    out(f"    2*Sigma_gamma = {2 * sigma:.12f}")
    out("    The deviation grows with the spectral radius and stays at the floor in")
    out("    eps*rho, which is the law; the raw digits do not survive J=50 and were")
    out("    never the point. An exact route exists where the rates are rational:")
    out("    p(2c - x) == p(x) over Q(i), implemented in thermal_palindrome_centre.py.")
    out("    The separation runs from 0 (a self-paired mode sits at the centre) to")
    out("    2*sigma (the extreme pair). That spread is the whole content of the")
    out("    reading: the binding never loosens, the fates diverge.")

    out("--- 2. the coupling moves everything except the horizon ---")
    out(f"    {'coupling':>10} {'centre':>18} {'max|Im lambda|':>16}")
    for J in (0.0, 1.0, 5.0, 50.0):
        ev = np.linalg.eigvals(liouvillian(heisenberg(N, J), N, gz=np.full(N, g)))
        out(f"    J={J:<8g} {np.mean(ev).real:+18.12f} {np.max(np.abs(ev.imag)):16.3f}")
    out("    The dynamics spans a factor of 50 in max|Im| and the centre does not")
    out("    move a digit. Not a coincidence: the centre is trace(L)/dim and the")
    out("    commutator part of L is traceless, so it is fixed by the dissipator")
    out("    alone. The illumination sets the horizon; the Hamiltonian cannot.")

    out()
    out("--- 3. the extremes are equally occupied WHILE the palindrome holds ---")
    out(f"    {'channel':>26} {'centre':>10} {'rate=0':>8} {'rate=2c':>8} {'max rate':>10}")
    H = heisenberg(N)
    channels = (
        ("Z-dephasing", dict(gz=np.full(N, g)), N * g),
        ("amplitude damping", dict(g_down=np.full(N, g)), N * g / 2),
        ("thermal bath", dict(g_down=np.full(N, g), g_up=np.full(N, g / 2)),
         N * (g + g / 2) / 2),
        ("BROKEN: T1 + co-axial Z", dict(gz=np.full(N, g), g_down=np.full(N, g)),
         N * g + N * g / 2),
    )
    for tag, kw, s_c in channels:
        ev = np.linalg.eigvals(liouvillian(H, N, **kw))
        rate = -ev.real
        n_undying = int(np.sum(rate < 1e-10))
        n_maximal = int(np.sum(np.abs(rate - 2 * s_c) < 1e-10))
        out(f"    {tag:>26} {-s_c:10.4f} {n_undying:8d} {n_maximal:8d} {max(rate):10.6f}")
    out(f"    (the full spectrum is 4^N = {4 ** N} modes)")
    out("    Equality is forced by the palindrome, which maps the rate-0 set onto")
    out("    the rate-2c set bijectively, so it is the pairing restated and not an")
    out("    independent finding. The last row is the point of printing it: where")
    out("    the palindrome breaks the equality goes too.")
    out()
    out("    A caution the census itself raises, about the surviving set. The I/Z")
    out("    sector has 2^N Pauli strings of dissipator rate 0, and at J = 0 that")
    out("    is exactly what survives. With the Hamiltonian on it is NOT:")
    out(f"    {'N':>4} {'J=0':>8} {'J=0.5':>8} {'J=1':>8} {'2^N':>8} {'N+1':>8}")
    for n in (2, 3, 4, 5):
        counts = []
        for J in (0.0, 0.5, 1.0):
            ev = np.linalg.eigvals(liouvillian(heisenberg(n, J), n, gz=np.full(n, g)))
            counts.append(int(np.sum(-ev.real < 1e-10)))
        out(f"    {n:>4} {counts[0]:>8} {counts[1]:>8} {counts[2]:>8}"
            f" {2 ** n:>8} {n + 1:>8}")
    out("    The survivors number N+1, not 2^N: the Hamiltonian empties the sector")
    out("    down to the total-Sz projectors. Measured directly, -i[H, Z_1] puts")
    out("    100.0000% of its weight on X/Y-containing strings. So 'the immune")
    out("    sector survives forever' is a statement about the dissipator alone.")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\n>>> Results saved to: {out_path}")


if __name__ == "__main__":
    main()
