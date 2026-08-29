"""The rate is a unit, gamma-alone scaling is a Q-move, and a Q-move is
invisible on a sector that does not see J.

NOT an F-number. The scale law itself is owned, Tier 1, by
experiments/Q_SCALE_THREE_BANDS.md:83 ("L(lambda*J, lambda*gamma_0) =
lambda*L(J, gamma_0)"), its unit reading by docs/THE_GENESIS_OF_AN_OSCILLATION.md,
and the statement that scaling gamma alone moves Q by
docs/proofs/INCOMPLETENESS_PROOF.md. This file measures them, and measures the
one thing none of those says: that a gamma-only sweep at fixed J is invisible
on a state whose sector cannot see J, so the SAME sweep gives opposite answers
on two states of one chain (part F). Plus two scope traps (part G).

NOT claimed here, and withdrawn during the work: any statement about F14's
K = gamma*t_cross. Testing that needs the repository's CPsi (purity times the
l1 off-diagonal over 3), and the candidate counterexample |01> has CPsi(0) = 0,
i.e. it starts below the threshold and crosses it 10, 4 and 2 times at the three
rates, so it has no first crossing to compare. Bell+ alone does give
K = 0.03735 at every gamma, spread exactly 0.000%, matching F25's printed value,
but one arm is not a discriminator.

WHAT CARRIES CONTENT AND WHAT DOES NOT, because the count of gates is not
the size of the evidence. A, B and D are one line of algebra read three ways,
and at dyadic a their exact 0.0 is an IEEE fact rather than a result: scaling
by a power of two is exact in radix 2, so that input class CANNOT break them.
They are kept because a construction can still be miswired, not because they
test the law. The independent readings are A at non-dyadic a, C, E1, E3b, F
and G; E2 and E3a are the law restated and are labelled as such in place.

  A  The JOINT scaling:  L(a*H, a*gamma) = a * L(H, gamma).
     Exact (== 0.0) at dyadic a, by radix-2 arithmetic and not by physics.
     At non-dyadic a there is no exact route, so the tolerance is an error
     law (eps * a * max|L|) and the ratio to it is printed.
  B  gamma-ALONE scaling is NOT a symmetry, and the leftover has a closed
     form: exactly (a-1) times the commutator part, the part carrying no rate.
  D  What gamma-alone scaling IS:  L(H, a*gamma) = a * L(H/a, gamma), exactly.
     At H homogeneous of degree 1 in J this reads J -> J/a, i.e. Q -> Q/a.
     A move along the Q axis wearing the costume of a broken symmetry.
  C  The trajectory face: (a*H, a*gamma) at time t is (H, gamma) at time a*t.
  E  The negative control, rehabilitated. GAMMA_TIME_DISTINCTION's Part 3
     measured tau = gamma*t failing to collapse irreversible observables and
     read it as "gamma provides the arrow, J provides the content". E1
     reproduces that run's published NUMBERS (0.790 for S(rho_A), 0.057 for
     purity) at its own grid and its own all-pairs statistic, then E3b shows
     the dynamics does move along Q, so the break has somewhere to come from:
     that sweep holds J = 1.0 while gamma runs over 20x, i.e. it sweeps Q from
     100 to 5 and compares five different systems.
  F  The same sweep, two states, opposite answers, on the TRAJECTORY (purity
     at matched tau) rather than on any crossing time: Bell+ is flat because
     its sector is H-dead, so Q has nothing to act on; |01> is not. The
     mechanism is measured ([H, rho_Bell+] = 0 exactly), not asserted. F2 then
     fences the phrase "a different system": that is unconditional about the
     GENERATOR and conditional about what is OBSERVED.
  G  The two fences as counterexamples rather than as prose, both of which a
     reader can walk into: the c-operator parameterization (c -> sqrt(a)*c,
     NOT a*c) and a bath scale the group does not touch (a fixed temperature).

Scope: H homogeneous of degree 1 in J; rates entering L linearly; and (H,
gamma) exhausting the dimensionful parameters. The framework's Heisenberg/XY
chains under Pauli dephasing meet all three. G is what happens when the last
two are read carelessly.

Run: python simulations/gamma_unit_scaling_gate.py
Out: simulations/results/gamma_unit_scaling_gate.txt
"""
import sys
import os
import math
import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import framework as fw
from framework.lindblad import lindbladian_z_dephasing

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "gamma_unit_scaling_gate.txt")

_lines = []


def log(msg=""):
    print(msg)
    _lines.append(str(msg))


PASS = 0
FAIL = 0


def gate(name, ok, detail=""):
    global PASS, FAIL
    if ok:
        PASS += 1
    else:
        FAIL += 1
    log("  [{}] {}{}".format("PASS" if ok else "FAIL", name,
                             "   " + detail if detail else ""))


rng = np.random.default_rng(20260829)

log("=" * 78)
log("The rate is a unit; a gamma-only sweep is a Q-move; a Q-move needs a")
log("sector that sees J. Plus the two scope traps.")
log("=" * 78)

# ---------------------------------------------------------------- A, B, D, C
for N in (2, 3, 4):
    chain = fw.ChainSystem(N=N, gamma_0=0.05, J=1.0)
    H = chain.H
    for label, gam in (("uniform", np.full(N, 0.05)),
                       ("random ", rng.uniform(0.01, 0.9, size=N))):
        L = lindbladian_z_dephasing(H, gam)
        nrmL = np.abs(L).max()
        log("\nN={}  gamma profile {}  max|L| = {:.6g}".format(N, label, nrmL))

        # A: an exact route exists at dyadic a, so compare exactly.
        for a in (0.5, 2.0, 4.0):
            r = np.abs(lindbladian_z_dephasing(a * H, a * gam) - a * L).max()
            gate("A joint scaling, exact at dyadic a={}".format(a),
                 r == 0.0, "residual = {!r}".format(r))

        # A': no exact route at non-dyadic a. State the error law, print the ratio.
        for a in (0.1, 1.0 / 3.0, math.pi):
            r = np.abs(lindbladian_z_dephasing(a * H, a * gam) - a * L).max()
            bound = 8 * np.finfo(float).eps * a * nrmL
            gate("A joint scaling, eps-law at a={:.6g}".format(a), r <= bound,
                 "residual = {:.3e}  bound = {:.3e}  ratio = {:.3f}".format(
                     r, bound, r / bound if bound else 0.0))

        # B: gamma alone is not a symmetry, and the leftover has a closed form.
        comm = np.abs(lindbladian_z_dephasing(H, np.zeros(N))).max()
        for a in (2.0, 4.0):
            broken = np.abs(lindbladian_z_dephasing(H, a * gam) - a * L).max()
            gate("B gamma-alone BREAKS at a={}".format(a), broken > 0.1 * nrmL,
                 "break = {:.6g}".format(broken))
            gate("B  break == (a-1)*max|commutator part| at a={}".format(a),
                 abs(broken - (a - 1) * comm) <= 1e-12 * nrmL,
                 "delta = {:.3e}".format(abs(broken - (a - 1) * comm)))

        # D: what it actually is.
        for a in (2.0, 4.0, 0.5):
            r = np.abs(lindbladian_z_dephasing(H, a * gam)
                       - a * lindbladian_z_dephasing(H / a, gam)).max()
            gate("D gamma-alone IS a Q-move: L(H,a.g) == a*L(H/a,g), a={}".format(a),
                 r == 0.0, "residual = {!r}".format(r))

        # C: the trajectory face.
        d = 2 ** N
        psi = rng.normal(size=d) + 1j * rng.normal(size=d)
        psi /= np.linalg.norm(psi)
        rho0 = np.outer(psi, psi.conj())
        t, a = 1.7, math.pi / 2.0
        rho_fast = (expm(lindbladian_z_dephasing(a * H, a * gam) * t)
                    @ rho0.flatten()).reshape(d, d)
        rho_slow = (expm(L * (a * t)) @ rho0.flatten()).reshape(d, d)
        r = np.abs(rho_fast - rho_slow).max()
        gate("C (aH,a.g) at t == (H,g) at a*t", r < 40 * np.finfo(float).eps * abs(a * t) * nrmL,
             "max|drho| = {:.3e}".format(r))

# ------------------------------------------------------------------------ E
log("\n" + "=" * 78)
log("E  The negative control, rehabilitated")
log("   GAMMA_TIME_DISTINCTION Part 3 settings, reproduced then explained")
log("=" * 78)


def _H2(J):
    """The two-site Heisenberg Hamiltonian, homogeneous of degree 1 in J."""
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    return J * (np.kron(X, X) + np.kron(Y, Y) + np.kron(Z, Z))


def _observables(L, rho0, t):
    """The two irreversible observables the source's Part 3 tabulates."""
    d = rho0.shape[0]
    rho = (expm(L * t) @ rho0.flatten()).reshape(d, d)
    rho = (rho + rho.conj().T) / 2
    rA = np.trace(rho.reshape(2, 2, 2, 2), axis1=1, axis2=3)
    rA = (rA + rA.conj().T) / 2
    ev = np.clip(np.linalg.eigvalsh(rA), 1e-15, None)
    S = float(-np.sum(ev * np.log2(ev)))
    P = float(np.real(np.trace(rho @ rho)))
    return S, P


rho0 = np.zeros((4, 4), dtype=complex)
rho0[1, 1] = 1.0                     # |01>, the source's own initial state
gammas_test = [0.01, 0.02, 0.05, 0.10, 0.20]
# The source's own grid and its own statistic, both read off the committed
# script simulations/gamma_is_time_proof.py: linspace(0, 1.0, 100) at :236, and
# the maximum over ALL PAIRS of curves at :270-272, not against the first. Both
# matter: against-the-first gives 0.729 and reproduces nothing. Taken from the
# script rather than chosen here, which is the difference between provenance
# and a statistic tuned to hit its target.
tau_points = np.linspace(0.0, 1.0, 100)


def _allpairs(curves, col):
    ks = list(curves.keys())
    return max(np.abs(curves[ks[i]][:, col] - curves[ks[j]][:, col]).max()
               for i in range(len(ks)) for j in range(i + 1, len(ks)))


# E1: the break at fixed J. This is the published reading's setting, and the
# gate reproduces its published NUMBER rather than merely its verdict.
curves_fixedJ = {}
for g in gammas_test:
    L = lindbladian_z_dephasing(_H2(1.0), np.array([g, g]))
    curves_fixedJ[g] = np.array([_observables(L, rho0, tau / g)
                                 for tau in tau_points])

dS = _allpairs(curves_fixedJ, 0)
dP = _allpairs(curves_fixedJ, 1)
log("\n  fixed J = 1.0, gamma over {}  =>  Q from {:.0f} down to {:.0f}".format(
    gammas_test, 1.0 / gammas_test[0], 1.0 / gammas_test[-1]))
gate("E1 tau-scaling BREAKS at fixed J, reproducing the source's 0.790",
     abs(dS - 0.790) < 0.0005,
     "max delta S(rho_A) = {:.6f}   max delta purity = {:.6f}".format(dS, dP))

# E2: the same observables scaled JOINTLY (Q held fixed) collapse to one curve.
#
# READ THIS BEFORE CITING E2 AS EVIDENCE. It is not independent evidence for
# the law, and saying so is the honest version. Building L at (a*J, a*gamma)
# produces a matrix that IS a*L_ref exactly, by part A, and the matched time is
# t_ref/a, so E2 evaluates expm(a*L_ref * t_ref/a) against expm(L_ref * t_ref):
# the same product, and its residual is a report on expm's rounding. Its role
# is to stand beside E1 at the SAME grid and statistic, so that the 0.790 has
# something to be 0.790 against; the independent content of this file is A at
# non-dyadic alpha, C, E1, E3b, F, and the two fence gates G.
curves_joint = {}
for g in gammas_test:
    a = g / gammas_test[0]
    L = lindbladian_z_dephasing(_H2(1.0 * a), np.array([g, g]))
    curves_joint[g] = np.array([_observables(L, rho0, tau / g)
                                for tau in tau_points])

dS = _allpairs(curves_joint, 0)
dP = _allpairs(curves_joint, 1)
log("\n  joint scaling, Q = J/gamma held at {:.0f} throughout".format(
    1.0 / gammas_test[0]))
gate("E2 same grid and statistic, Q held: the 0.790 goes away (restatement,"
     " not evidence)", dS < 1e-9 and dP < 1e-9,
     "max delta S(rho_A) = {:.3e}   max delta purity = {:.3e}".format(dS, dP))

# E3: the curve is a function of Q alone, not of gamma. E3b is the half that
#     carries content: it shows the dynamics moves along Q at all, so E1's
#     break has somewhere to come from. E3a is joint scaling by exactly 4 and
#     is therefore the same restatement as E2, kept only as its control.
L_a = lindbladian_z_dephasing(_H2(1.0), np.array([0.05, 0.05]))       # Q = 20
L_b = lindbladian_z_dephasing(_H2(4.0), np.array([0.20, 0.20]))       # Q = 20
same_Q = max(abs(_observables(L_a, rho0, tau / 0.05)[0]
                 - _observables(L_b, rho0, tau / 0.20)[0])
             for tau in tau_points)
L_c = lindbladian_z_dephasing(_H2(0.25), np.array([0.05, 0.05]))      # Q = 5
diff_Q = max(abs(_observables(L_a, rho0, tau / 0.05)[0]
                 - _observables(L_c, rho0, tau / 0.05)[0])
             for tau in tau_points)
gate("E3a same Q, gamma apart by 4x: curves agree (restatement)", same_Q < 1e-9,
     "max delta S = {:.3e}".format(same_Q))
gate("E3b same gamma, Q apart by 4x: curves differ", diff_Q > 0.1,
     "max delta S = {:.3f}".format(diff_Q))

# ------------------------------------------------------------------------ F
log("\n" + "=" * 78)
log("F  Why F14 holds where GAMMA_TIME_DISTINCTION's Part 3 fails")
log("   Same fixed J, same gamma range, two initial states")
log("=" * 78)

# F14 (K = gamma * t_cross = constant) and GAMMA_TIME_DISTINCTION Part 3
# (tau-scaling breaks) both sweep gamma at fixed J and reach opposite verdicts.
# The corollary says what separates them: gamma-alone scaling is a Q-move, so it
# is invisible exactly where Q has nothing to act on. F14's state is Bell+, whose
# sector is H-dead (LATTICE_OPENING_LAW: the cat pair is J-free); Part 3's state
# is |01>, which is not. Neither result needs correcting.
bell = np.zeros((4, 4), dtype=complex)
for _i in (0, 3):
    for _j in (0, 3):
        bell[_i, _j] = 0.5

for _name, _rho0, _expect_flat in (("Bell+, H-dead sector", bell, True),
                                   ("|01>, H-live sector ", rho0, False)):
    _curves = []
    for g in (0.01, 0.05, 0.20):
        L = lindbladian_z_dephasing(_H2(1.0), np.array([g, g]))
        _curves.append(np.array([_observables(L, _rho0, tau / g)[1]
                                 for tau in tau_points]))
    _d = max(np.abs(c - _curves[0]).max() for c in _curves)
    if _expect_flat:
        gate("F  " + _name + ": tau-scaling HOLDS at fixed J", _d < 1e-12,
             "max delta purity = {:.3e}".format(_d))
    else:
        gate("F  " + _name + ": tau-scaling BREAKS at fixed J", _d > 1e-3,
             "max delta purity = {:.3e}".format(_d))

# The mechanism, checked rather than asserted: the cat coherence is annihilated
# by the commutator, so no J can enter that sector's dynamics.
_comm = _H2(1.0) @ bell - bell @ _H2(1.0)
gate("F  the mechanism: [H, rho_Bell+] = 0, so J cannot enter",
     np.abs(_comm).max() < 1e-12, "max|[H,rho]| = {:.3e}".format(np.abs(_comm).max()))
_comm01 = _H2(1.0) @ rho0 - rho0 @ _H2(1.0)
gate("F  and the control: [H, rho_|01>] != 0",
     np.abs(_comm01).max() > 0.1, "max|[H,rho]| = {:.3f}".format(np.abs(_comm01).max()))

# F2: the sharp form, and it FENCES the phrase "gamma-alone scaling gives a
# different system". That phrase is about the GENERATOR and is unconditional
# (part D). As a statement about what is OBSERVED it is conditional: on a sector
# that does not see J, one-sided gamma scaling is an exact symmetry up to the
# time rescaling, because there is no Q for the move to move.
for _name, _r, _flat in (("Bell+, H-dead ", bell, True),
                         ("|01>, H-live  ", rho0, False)):
    _g, _a = 0.3, 4.0
    _L1 = lindbladian_z_dephasing(_H2(1.0), np.array([_g, _g]))
    _L2 = lindbladian_z_dephasing(_H2(1.0), np.array([_a * _g, _a * _g]))
    _d = max(abs(_observables(_L1, _r, _t)[1] - _observables(_L2, _r, _t / _a)[1])
             for _t in np.linspace(0.1, 8.0, 80))
    if _flat:
        gate("F2 " + _name + ": ONE-SIDED gamma scaling IS a symmetry here",
             _d < 1e-14, "max|dPurity| = {:.3e}".format(_d))
    else:
        gate("F2 " + _name + ": ONE-SIDED gamma scaling is NOT a symmetry",
             _d > 0.1, "max|dPurity| = {:.3f}".format(_d))

# ------------------------------------------------------------------------ F3
# F14's own number, under a gamma sweep, in F14's own book. This is ONE ARM and
# is not a discriminator: it confirms that the K-invariance F14 reports is real
# for Bell+ and reproduces the exact closed form, which is what the scope note
# in ANALYTICAL_FORMULAS cites. WATCH THE BOOK: CPsi = C * l1/(d-1) with C the
# WOOTTERS CONCURRENCE gives ln(4/3)/8; with C = purity it gives F25's 0.03735.
# Quoting one book's number beside the other book's measurement is the trap
# CROSSING_TAXONOMY's own "which C, before which book" section warns about, and
# an earlier version of the scope note walked into it.
log("\n" + "=" * 78)
log("F3 F14's K = gamma*t_cross for Bell+, in the Wootters-concurrence book")
log("=" * 78)


def _concurrence(r):
    sy = np.kron(np.array([[0, -1j], [1j, 0]], dtype=complex),
                 np.array([[0, -1j], [1j, 0]], dtype=complex))
    ev = np.sqrt(np.clip(np.real(np.linalg.eigvals(r @ sy @ r.conj() @ sy)),
                         0.0, None))
    ev = np.sort(ev)[::-1]
    return float(max(0.0, ev[0] - ev[1] - ev[2] - ev[3]))


def _cpsi_wootters(r):
    l1 = float(np.sum(np.abs(r)) - np.sum(np.abs(np.diag(r))))
    return _concurrence(r) * l1 / 3.0


_Ks, _ncross = [], []
for _g in (0.02, 0.05, 0.10):
    _L = lindbladian_z_dephasing(_H2(1.0), np.array([_g, _g]))
    _ts = np.linspace(0.0, 300.0, 15001)
    _c = np.empty(_ts.size)
    for _i, _t in enumerate(_ts):
        _x = (expm(_L * _t) @ bell.flatten()).reshape(4, 4)
        _c[_i] = _cpsi_wootters((_x + _x.conj().T) / 2.0)
    _ab = _c > 0.25
    _ncross.append(int(np.sum(_ab[1:] != _ab[:-1])))
    _idx = np.where(_ab[:-1] & ~_ab[1:])[0]
    _lo, _hi = _ts[_idx[0]], _ts[_idx[0] + 1]
    for _ in range(60):                      # bisect, so the grid is not the answer
        _m = (_lo + _hi) / 2.0
        _x = (expm(_L * _m) @ bell.flatten()).reshape(4, 4)
        if _cpsi_wootters((_x + _x.conj().T) / 2.0) > 0.25:
            _lo = _m
        else:
            _hi = _m
    _Ks.append(_g * (_lo + _hi) / 2.0)

gate("F3a Bell+ crosses the 1/4 threshold exactly once at every rate",
     _ncross == [1, 1, 1], "crossings = {}".format(_ncross))
_spread = (max(_Ks) - min(_Ks)) / (sum(_Ks) / len(_Ks))
gate("F3b K = gamma*t_cross is rate-independent for Bell+", _spread < 1e-4,
     "K = {}  spread = {:.4%}".format([round(float(k), 6) for k in _Ks], _spread))
gate("F3c and it is F14's exact closed form ln(4/3)/8",
     all(abs(k - math.log(4.0 / 3.0) / 8.0) < 1e-5 for k in _Ks),
     "ln(4/3)/8 = {:.6f}".format(math.log(4.0 / 3.0) / 8.0))

# ------------------------------------------------------------------------ G
log("\n" + "=" * 78)
log("G  The two fences, as counterexamples rather than as prose")
log("=" * 78)

from framework.lindblad import lindbladian_general   # noqa: E402

# G1: the law is stated in the RATE parameterization. In the c-operator
# parameterization the rate sits under a square root, c = sqrt(gamma)*A, and the
# dissipator is QUADRATIC in c. The framework's own general builder takes c_ops,
# so a reader applying the law there with c -> a*c gets it wrong by 100%.
_H1 = np.diag([1.0, -1.0]).astype(complex)
_A = np.diag([1.0, -1.0]).astype(complex)
_c = [math.sqrt(0.3) * _A]
_Lg = lindbladian_general(_H1, _c)
_wrong = np.abs(lindbladian_general(3.0 * _H1, [3.0 * x for x in _c])
                - 3.0 * _Lg).max()
_right = np.abs(lindbladian_general(3.0 * _H1, [math.sqrt(3.0) * x for x in _c])
                - 3.0 * _Lg).max()
gate("G1a c-parameterization: c -> a*c is WRONG, and grossly so", _wrong > 1.0,
     "residual = {:.3f}   (this is the trap, not a failure)".format(_wrong))
gate("G1b c-parameterization: c -> sqrt(a)*c is the law", _right < 1e-14,
     "residual = {:.3e}".format(_right))

# G2: (H, gamma) must EXHAUST the dimensionful parameters. A thermal bath at a
# FIXED temperature holds a third scale, and the law breaks; scale T with them
# and it comes back exactly. Without this fence the law would be false as
# advertised for every finite-temperature channel in the literature.
_sm = np.array([[0.0, 1.0], [0.0, 0.0]], dtype=complex)


def _thermal(w, Gam, T):
    n = 1.0 / (math.exp(w / T) - 1.0)
    return lindbladian_general((w / 2.0) * np.diag([1.0, -1.0]).astype(complex),
                               [math.sqrt(Gam * (1.0 + n)) * _sm,
                                math.sqrt(Gam * n) * _sm.conj().T])


_r0 = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)


def _run(L, t):
    return (expm(L * t) @ _r0.flatten()).reshape(2, 2)


_w, _G, _T, _tt = 1.0, 0.1, 0.5, 4.0
# Swept over several a rather than one, and NOT asserted as exactly 0.0. An
# earlier version did assert that, and it held only for one (a, t, literal)
# spelling: writing the scaled rate as 3*0.1 instead of 0.3 already moves it to
# 7e-18. There is no exact route through expm, so the tolerance is stated as a
# law (a few eps against the state's own scale) and the residual is printed.
for _aa in (2.0, 3.0, 5.0, 7.0):
    _base = _run(_thermal(_w, _G, _T), _tt)
    _Tfixed = _run(_thermal(_aa * _w, _aa * _G, _T), _tt / _aa)
    _Tscaled = _run(_thermal(_aa * _w, _aa * _G, _aa * _T), _tt / _aa)
    _broke = np.abs(_base - _Tfixed).max()
    _kept = np.abs(_base - _Tscaled).max()
    gate("G2a thermal bath at a={:g}, T held fixed: the law BREAKS".format(_aa),
         _broke > 1e-3, "|drho| = {:.6f}".format(_broke))
    gate("G2b thermal bath at a={:g}, T scaled too: the law holds".format(_aa),
         _kept <= 16 * np.finfo(float).eps,
         "|drho| = {:.3e}  (bound {:.3e})".format(
             _kept, 16 * np.finfo(float).eps))

log("\nGATES: {} passed, {} failed".format(PASS, FAIL))

os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
with open(OUT_PATH, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(_lines) + "\n")
print("\nwritten: " + OUT_PATH)
sys.exit(0 if FAIL == 0 else 1)
