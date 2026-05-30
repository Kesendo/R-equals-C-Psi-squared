"""_only_the_relative_survives.py - the absolute is gauge, only the relative is physical.

A research probe, built to be REFUTABLE. It tests one claim on two axes - the
tick (gamma) and the rotation (omega) - and lets the algebra judge, not the
framing. Every check prints its raw numbers so a human reads the verdict from
the deviations, not from a green PASS.

THE CLAIM: in a closed-plus-dephasing quantum system there is no absolute tick
and no absolute rotation. What you can read from inside is only relative: ratios
of rates, phase angles of modes, and the differences between rotation rates. The
absolute second and the absolute energy zero are gauge.

  PART 1  the gamma axis (the tick)
    A. Absolute gamma-scale is gauge. Rescaling the WHOLE clock (J and gamma
       together, L -> g*L) is a pure change of time unit: the dimensionless
       intra-system observables (per-mode angle theta_k, and rate_k/rate_min)
       are bit-exact unchanged, while the absolute spectrum scales by exactly g.
       Refutation: any dimensionless observable that moves under pure scale.
       (Honest caveat learned while building this: scaling gamma ALONE, with J
       fixed, is NOT a pure time rescale and does move the angles - see the note
       printed in Part A. The gauge freedom is the rescaling of the time UNIT,
       which rescales every rate, the Hamiltonian one included.)
    B. The gamma-RATIO is physical. Hold Sigma gamma fixed but change the ratio
       across sites (uniform vs skewed). The dimensionless observables DIFFER.
       Refutation: if the ratio were invisible.

  PART 2  the omega axis (rotation, the relative)
    C. Absolute rotation/energy is gauge. Shift H -> H + c*I. The Liouvillian L
       and every observable are bit-exact unchanged, because [c*I, .] = 0.
       Refutation: if a global energy shift changed anything.
    D. The relative rotation is physical, and the framework sorts it by Pi^2
       parity (F80). The observable rotation is the set of Bohr differences
       {Im lambda_a - Im lambda_b}; show they are nonzero. Then the
       framework-specific structure: the F1 palindrome residual M for a single
       chain bond has a purely imaginary spectrum whose |Im| clusters are
         - Pi^2-ODD bond (XY):                2*|single H eigenvalues|  (absolute-like)
         - Pi^2-EVEN non-truly bond (YZ+ZY):  2*|H eigenvalue DIFFERENCES| (Bohr, relative)
       i.e. the framework has its OWN absolute-vs-relative-rotation sorter, the
       Pi^2 parity deciding single-eigenvalues vs differences.

HONESTY: this construction is only a spin-chain Hamiltonian plus local dephasing,
nothing imported from outside. "Only the relative is physical" is a READING of the
algebra, not a derivation. A and C are general QM / dimensional gauge freedom
(A: dimensionless = scale-free; C: [c*I, .] = 0 is automatic); only D is
framework-specific (the F80 Pi^2-parity single-vs-difference sort).

Tom + Claude, 2026-05-30. Run: python "simulations/_only_the_relative_survives.py"
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from framework.pauli import _build_bilinear
from framework.lindblad import lindbladian_pauli_dephasing, palindrome_residual

ATOL = 1e-10                                   # "bit-exact" floor for this probe
HEIS = [("X", "X", 1.0), ("Y", "Y", 1.0), ("Z", "Z", 1.0)]
XY = [("X", "X", 0.5), ("Y", "Y", 0.5)]        # the XY model: (XX + YY)/2

_results = []


def report(name, ok, extra=""):
    _results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{extra}")
    return ok


def fmt(xs, p=6):
    return "[" + ", ".join(f"{x:.{p}f}" for x in xs) + "]"


# ----------------------------------------------------------------------
# builders + dimensionless observables of the Liouvillian spectrum
# ----------------------------------------------------------------------

def build_L(N, terms, J, gamma_l):
    """Open chain: H = J * sum_bond(terms), uniform-letter Z-dephasing gamma_l."""
    bonds = [(i, i + 1) for i in range(N - 1)]
    H = J * _build_bilinear(N, bonds, terms)
    return lindbladian_pauli_dephasing(H, gamma_l, dephase_letter="Z"), H


def mode_angles(L, rel_floor=1e-9):
    """Per-mode angle theta_k = atan2(|Im lambda|, |Re lambda|) in degrees, as a
    sorted multiset. Defined ONLY on modes whose |lambda| is above a tiny
    relative floor: the phase of a (near-)zero eigenvalue is numerically
    meaningless (Re, Im both ~1e-16), so including it would inject noise, not
    physics. This floor is the honest scope of the angle observable."""
    e = np.linalg.eigvals(L)
    mag = np.abs(e)
    keep = mag > rel_floor * mag.max()
    th = np.degrees(np.arctan2(np.abs(e[keep].imag), np.abs(e[keep].real)))
    return np.sort(np.round(th, 8))


def rate_ratios(L, floor=1e-9):
    """Dimensionless decay-rate ratios rate_k / rate_min over the decaying modes
    (rate_k = -Re lambda_k > 0), as a sorted multiset."""
    r = -np.linalg.eigvals(L).real
    r = r[r > floor]
    return np.sort(np.round(r / r.min(), 8))


def sorted_spec(L):
    s = np.linalg.eigvals(L)
    return np.array(sorted(s, key=lambda z: (round(z.real, 9), round(z.imag, 9))))


# ----------------------------------------------------------------------
# PART 1A - absolute gamma-scale is gauge (a pure rescale of the time unit)
# ----------------------------------------------------------------------

def part_1a(g=2.7):
    print("PART 1A - absolute clock-scale is gauge: rescale the WHOLE clock L -> g*L")
    print(f"          (general QM / dimensional gauge freedom - rescaling the time unit), g={g}")
    for N, terms, gam, name in [
        (2, XY,   [0.1, 0.1],        "XY  N=2"),
        (3, HEIS, [0.1, 0.2, 0.15],  "Heis N=3"),
    ]:
        # World B = World A run g times faster: BOTH the rotation J and the tick
        # gamma scale by g. That, not "gamma alone", is the pure time rescale.
        LA, _ = build_L(N, terms, 1.0, gam)
        LB, _ = build_L(N, terms, g, [g * x for x in gam])
        sA, sB = sorted_spec(LA), sorted_spec(LB)

        # The absolute spectrum is NOT invariant - it scales by exactly g.
        scale_dev = float(np.max(np.abs(sB - g * sA)))
        abs_move = float(np.max(np.abs(sB - sA)))
        report(f"{name}: spec(L_B) = g * spec(L_A)  (absolute spectrum scales, not invariant)",
               scale_dev < ATOL, f"   max|s_B - g*s_A| = {scale_dev:.2e},  |s_B - s_A| = {abs_move:.3f} (moved)")

        # The dimensionless observables ARE invariant.
        ang_dev = float(np.max(np.abs(mode_angles(LA) - mode_angles(LB))))
        report(f"{name}: per-mode angle theta_k multiset unchanged under the rescale",
               ang_dev < ATOL, f"   max|d theta| = {ang_dev:.2e} deg")
        rr_dev = float(np.max(np.abs(rate_ratios(LA) - rate_ratios(LB))))
        report(f"{name}: rate_k / rate_min multiset unchanged under the rescale",
               rr_dev < ATOL, f"   max|d ratio| = {rr_dev:.2e}")

    # The honest caveat, shown with numbers: gamma ALONE (J fixed) is NOT a pure
    # time rescale, and it DOES move the angles. The gauge is the time unit.
    N = 3
    L1, _ = build_L(N, HEIS, 1.0, [0.1, 0.2, 0.15])
    L2, _ = build_L(N, HEIS, 1.0, [g * 0.1, g * 0.2, g * 0.15])   # J fixed, gamma*g
    ang_only = float(np.max(np.abs(mode_angles(L1) - mode_angles(L2))))
    print(f"          caveat (gamma alone, J fixed, is NOT the gauge): angle moves by "
          f"max|d theta| = {ang_only:.3f} deg  <- expected to MOVE, not a failure")
    print()


# ----------------------------------------------------------------------
# PART 1B - the gamma-ratio is physical
# ----------------------------------------------------------------------

def part_1b():
    print("PART 1B - the gamma-RATIO is physical: same Sigma gamma, different ratio -> observables MOVE")
    print("          (claim: the ratio is visible in AT LEAST ONE dimensionless observable.")
    print("           NOTE the honest surprise below: at N=2 the XY decay RATES are ratio-blind.)")
    for N, terms, name in [(2, XY, "XY  N=2"), (3, HEIS, "Heis N=3")]:
        # Same total dephasing Sigma gamma; only the per-site RATIO differs.
        if N == 2:
            gam_unif, gam_skew = [0.1, 0.1], [0.05, 0.15]
        else:
            gam_unif, gam_skew = [0.1, 0.1, 0.1], [0.02, 0.10, 0.18]
        assert abs(sum(gam_unif) - sum(gam_skew)) < 1e-15, "Sigma gamma must match"
        Lu, _ = build_L(N, terms, 1.0, gam_unif)
        Ls, _ = build_L(N, terms, 1.0, gam_skew)

        au, as_ = mode_angles(Lu), mode_angles(Ls)
        ru, rs = rate_ratios(Lu), rate_ratios(Ls)
        ang_move = float(np.max(np.abs(au - as_))) if len(au) == len(as_) else float("inf")
        rr_move = float(np.max(np.abs(ru - rs))) if len(ru) == len(rs) else float("inf")

        # The CLAIM under test: the ratio is physical, i.e. it moves SOME dimensionless
        # observable. Refutation would be: it moves NONE (ratio fully invisible).
        report(f"{name}: gamma-ratio is visible in some dimensionless observable "
               f"(Sigma gamma = {sum(gam_unif)})",
               max(ang_move, rr_move) > 1e-6,
               f"   angle moves {ang_move:.4f} deg, rate-ratio moves {rr_move:.4f}")
        # Print BOTH channels so the per-observable behaviour is on the table, not hidden.
        print(f"          -> angle observable:      {'MOVES' if ang_move > 1e-6 else 'BLIND'} "
              f"(max|d theta| = {ang_move:.4f} deg)")
        print(f"          -> rate-ratio observable: {'MOVES' if rr_move > 1e-6 else 'BLIND'} "
              f"(max|d ratio| = {rr_move:.4f})")
        if rr_move <= 1e-6:
            r = -np.linalg.eigvals(Lu).real
            r = np.sort(np.round(r[r > 1e-9], 6))
            print(f"          !! honest finding: {name} decay-rate multiset is RATIO-BLIND. The")
            print(f"             rates sit at {fmt(sorted(set(r.tolist())), 4)} = multiples of Sigma gamma")
            print(f"             only (popcount-locked), so the ratio hides from the RATES here but")
            print(f"             still shows in the ANGLES. The ratio is physical; not every channel")
            print(f"             sees it. (A clean partial-refutation of 'every observable moves'.)")
    print()


# ----------------------------------------------------------------------
# PART 2C - absolute rotation/energy is gauge ([c*I, .] = 0, automatic)
# ----------------------------------------------------------------------

def part_2c(c=3.0):
    print(f"PART 2C - absolute rotation/energy is gauge: H -> H + c*I leaves L bit-exact unchanged")
    print(f"          (general QM - the commutator [c*I, .] = 0 is automatic), c={c}")
    for N, terms, gam, name in [
        (2, XY,   [0.1, 0.1],       "XY  N=2"),
        (3, HEIS, [0.1, 0.1, 0.1],  "Heis N=3"),
    ]:
        L, H = build_L(N, terms, 1.0, gam)
        d = 2 ** N
        L_shift = lindbladian_pauli_dephasing(H + c * np.eye(d), gam, dephase_letter="Z")
        dev = float(np.linalg.norm(L_shift - L))
        report(f"{name}: ||L(H + c*I) - L(H)||_F = 0  (global energy shift invisible)",
               dev < ATOL, f"   ||dL||_F = {dev:.2e}")
        # and therefore every spectral observable is identical
        sdev = float(np.max(np.abs(sorted_spec(L_shift) - sorted_spec(L))))
        report(f"{name}: full Liouvillian spectrum identical after the shift",
               sdev < ATOL, f"   max|d spec| = {sdev:.2e}")
    print()


# ----------------------------------------------------------------------
# PART 2D - the relative rotation is physical; F80 Pi^2-parity sorts it
# ----------------------------------------------------------------------

def m_imag_clusters(M, tol=1e-6, dp=6):
    """Distinct |Im| clusters of M's spectrum (rounded), plus max|Re| to test
    that M is purely imaginary (anti-Hermitian)."""
    e = np.linalg.eigvals(M)
    max_re = float(np.max(np.abs(e.real)))
    cl = sorted({round(abs(z.imag), dp) for z in e if abs(z.imag) > tol})
    return max_re, cl


def h_single_and_diff(H, dp=6, tol=1e-6):
    """2*|single eigenvalues| and 2*|eigenvalue differences| of Hermitian H,
    each as a sorted distinct cluster set (matching M's |Im| rounding)."""
    ev = np.linalg.eigvalsh(H)
    single = sorted({round(2 * abs(l), dp) for l in ev if abs(l) > tol})
    diff = sorted({round(2 * abs(a - b), dp) for a in ev for b in ev if abs(a - b) > tol})
    return single, diff


def part_2d():
    print("PART 2D - the relative rotation is physical; F80 Pi^2-parity sorts single vs differences")
    print("          (THIS is the framework-specific structure - not generic QM)")
    N = 3
    gam = [0.1, 0.1, 0.1]
    sigma = sum(gam)
    bonds = [(i, i + 1) for i in range(N - 1)]

    # (i) The observable rotation = Bohr differences {Im lambda_a - Im lambda_b}.
    #     Show they are nonzero (the system genuinely rotates, relatively).
    L_heis, _ = build_L(N, HEIS, 1.0, gam)
    ims = np.sort(np.unique(np.round(np.linalg.eigvals(L_heis).imag, 6)))
    bohr = sorted({round(abs(a - b), 6) for a in ims for b in ims if abs(a - b) > 1e-6})
    report("Heisenberg N=3: Bohr differences {Im l_a - Im l_b} are nonzero (it rotates)",
           len(bohr) > 0, f"   #distinct |Bohr| = {len(bohr)}, max = {max(bohr):.4f}")

    # (ii) The F80 Pi^2-parity sorter on the F1 palindrome residual M.
    #      Pi^2-ODD bond (XY)            -> M sees 2*|single eigenvalues|  (absolute-like)
    #      Pi^2-EVEN non-truly (YZ+ZY)   -> M sees 2*|eigenvalue DIFFERENCES| (Bohr, relative)
    cases = [
        ("Pi^2-ODD  bond  (X,Y)",      [("X", "Y", 1.0)],                 "single"),
        ("Pi^2-EVEN bond  (Y,Z)+(Z,Y)", [("Y", "Z", 1.0), ("Z", "Y", 1.0)], "diff"),
    ]
    for label, terms, expect in cases:
        H = _build_bilinear(N, bonds, terms)
        L = lindbladian_pauli_dephasing(H, gam, dephase_letter="Z")
        M = palindrome_residual(L, sigma, N)
        max_re, mcl = m_imag_clusters(M)
        single, diff = h_single_and_diff(H)
        m_single = (mcl == single)
        m_diff = (mcl == diff)

        print(f"    {label}:  ||M||_F = {np.linalg.norm(M):.4f}   max|Re spec(M)| = {max_re:.2e}")
        print(f"        M |Im| clusters        = {fmt(mcl)}")
        print(f"        2*|single H eigenvalues| = {fmt(single)}   (match: {m_single})")
        print(f"        2*|H eigenvalue diffs|   = {fmt(diff)}   (match: {m_diff})")

        report(f"{label}: M is purely imaginary (anti-Hermitian)", max_re < 1e-9)
        if expect == "single":
            report(f"{label}: M clusters = 2*|single eigenvalues| (absolute-like), NOT differences",
                   m_single and not m_diff)
        else:
            report(f"{label}: M clusters = 2*|eigenvalue DIFFERENCES| (Bohr, relative), NOT single",
                   m_diff and not m_single)
    print()


# ----------------------------------------------------------------------

def main():
    print("=" * 80)
    print("ONLY THE RELATIVE SURVIVES - the absolute is gauge, the relative is physical")
    print("=" * 80)
    print("Only a spin-chain Hamiltonian + local Z-dephasing here, nothing imported.")
    print("\"Only the relative is physical\" is a READING of the algebra, not a derivation:")
    print("A and C are general QM gauge freedom, only D (F80) is framework-specific.")
    print("-" * 80)
    part_1a()
    part_1b()
    part_2c()
    part_2d()

    n_ok, n_tot = sum(_results), len(_results)
    print("=" * 80)
    print(f"RESULT: {n_ok}/{n_tot} bit-exact ({'ALL PASS' if n_ok == n_tot else 'SOME FAILED - read the numbers'})")
    print("=" * 80)
    print("""
Honest reading - what the algebra said:

  CONFIRMED BUT BUILT-IN (general QM / dimensional gauge freedom, not our structure):
    1A  Rescaling the whole clock (J and gamma together, L -> g*L) is a pure change
        of time unit. Every dimensionless observable - per-mode angle, rate ratios -
        is bit-exact invariant; the absolute spectrum scales by exactly g. This is
        just "you cannot read the second's length from inside". CAVEAT, printed
        above: gamma ALONE (J fixed) is NOT this gauge and DOES move the angles - the
        gauge is the time unit, which rescales the Hamiltonian rate too. The angle
        observable is also only defined on modes with nonzero |lambda| (a zero mode
        has no phase); that scope is honest, not a rig.
    1B  The gamma-RATIO across sites is physical: at fixed Sigma gamma, uniform vs
        skewed moves a dimensionless observable, so the relative profile survives
        while the absolute scale washes. HONEST SURPRISE the probe surfaced: WHICH
        observable sees it is not universal. At N=2 the XY decay-RATE multiset is
        ratio-blind - the rates lock to {Sigma gamma, 2 Sigma gamma} (popcount-set)
        and ignore the split - so the ratio hides from the rates there and shows
        only in the angles. At Heisenberg N=3 both channels see it. A clean partial
        refutation of the stronger claim "every dimensionless observable moves": the
        ratio is physical, but it is not visible in every channel.
    2C  H -> H + c*I changes nothing, bit-exact, because [c*I, .] = 0. The absolute
        energy zero is gauge. This is automatic in any Lindblad/commutator dynamics.

  FRAMEWORK-SPECIFIC (this is ours - the F80 Pi^2-parity sorter):
    2D  The observable rotation is the set of Bohr differences {Im l_a - Im l_b},
        and they are nonzero. The framework then has its OWN absolute-vs-relative
        sorter on the F1 palindrome residual M: Pi^2-ODD bonds (XY) give M a spectrum
        of 2*|single eigenvalues| (absolute-like), while Pi^2-EVEN non-truly bonds
        (YZ+ZY) give M a spectrum of 2*|eigenvalue DIFFERENCES| (the Bohr frequencies,
        relative). The Klein parity of the bond decides single-energies vs
        differences, bit-exact. That sorting is not generic QM - it is the structure
        F80 names.

  In one line: A and C are the price of having no master clock and no energy origin
  (everybody's QM pays it); only D is the framework saying something the bare
  Liouvillian does not - that Pi^2-parity is itself an absolute-vs-relative switch.
""")


if __name__ == "__main__":
    main()
