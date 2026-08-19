"""The theta_D refreeze harness (corner beat).

WHY THIS FILE EXISTS. Section 8a of the pre-registration registers the
refreeze as ONE unit, and its defining property is the CHAIN it runs
on: "re-freeze theta_D through the ANALYZE-SIDE chain (the runner's
fit_rate, whose dead-cell NaN differs from the gate's silent 0; pooled
inversion + single clip; fixed R_MAX_FIT; NESTED tables; the
REALIZED-DOSE time axis)". The committed gate cannot do that by
itself: its counts path applies readout confusion as a per-shot
keep/reroute and never inverts, so the numbers it froze were measured
on a chain the flight does not run.

WHAT IS REUSED AND WHAT IS NEW, stated because a second physics path
that silently differs from the gate would be worse than no harness:

  PHYSICS, verbatim from the committed gate (imported, never copied):
    strang, floquet_modes, phase_tables, run_bindings, and the
    gate-error bracket built from p2/f_leak/gates-per-block.
  COUNTS + ESTIMATOR, the runner's analyze chain (imported):
    per-binding RAW 64-dim counts -> pool -> ONE confusion inversion
    -> ONE clip -> kept-count floor -> arm_rates on fit_time_axis().

  The ONE new thing here is the bridge between them: turning the
  gate's per-binding one-magnon occupations into per-binding RAW
  64-dim count vectors, which is what the runner's chain consumes and
  what hardware actually returns. That bridge forces a modelling
  choice the gate never had to make, and it is a SWITCH here rather
  than a silent decision (see OUT_OF_SECTOR below), because the point
  of the refreeze is to measure what the chain does, not to pick the
  variant that flatters it.

Run:  python corner_beat_refreeze.py --reps 500 --heldout 500
      python corner_beat_refreeze.py --reps 20 --quick   (a smoke run)
      python corner_beat_refreeze.py --offlocus --reps 500
                                     (the specificity diagnostic)
The fit axis DEFAULTS to realized_dose, the axis section 5 registers
for this refreeze; --axis nominal is available only to measure the
difference, never to freeze on.
Every RATE (P(detect), the H0 false rates, the floor trip) is measured
on HELD-OUT campaigns and re-read at the REGISTERED max-over-bracket
thresholds, which are the two things section 7 and section 8a ask for
and the two things this file's first version did not do.

SIX CAMPAIGNS PER CORNER since 2026-08-18, not five: H0-C, H0-both and
H1 build the estimated objects, and held-out H0-C, H0-both and H1 read
them. H1 joined the held-out side when the s2(C) floor became a
FUNCTION whose coefficients are chosen against H1 draws (see the block
above main()); an estimated object read on the sample it was fitted to
is exactly the defect two reviewers found in theta_D on 2026-08-17.
That is ~2 h for the four-corner freeze rather than ~1.7 h.
THE FREEZE NEEDS >= 299 HELD-OUT H1 DRAWS PER CORNER, and this is a
hard floor, not a preference: below it the reachability ceiling of 0.01
cannot be demonstrated by any campaign, however clean (min_n_for_ceiling
below). A smoke run reports that it hit the sample size, never that the
floor failed.
"""

import argparse
import importlib.util as _ilu
import json
import math
import time
from datetime import datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent


def _load(name, path):
    spec = _ilu.spec_from_file_location(name, path)
    mod = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


rcb = _load("rcb_refreeze", HERE / "run_corner_beat.py")
gate = _load("cbgate_refreeze", rcb.GATE_PATH)

N = rcb.N
D6 = rcb.D6
DY = rcb.DY

# ------------------------------------------------------ frozen config
# The refreeze runs the FROZEN configuration only. It is not a search.
Q = rcb.Q
JDT = rcb.JDT
STEPS = rcb.STEPS_MAX
K = rcb.K_GRID
SHOTS = rcb.SHOTS_SCIENCE
M_BIND = rcb.M_BIND
GPB = 2                      # fractional: 2 two-qubit gates per XY block
P2 = gate.P2_RECORDED        # 0.005, and since Amendment 1.8 the guard
F_LEAK_ENDS = gate.F_LEAK_BRACKET          # (8/15, 0.9)
HOP_FRACTIONS = (0.5, 0.0)   # the within-sector split; the gate pins
                             # 0.5 and section 8a registers the MAXIMUM
                             # over this bracket (a 4.7x lever on a
                             # theta_D-sized quantity)

# ------------------------------------------------- the counted curve
# THE TWO BRACKETS ARE NOT INDEPENDENT (2026-08-18). Section 10 derives
# f_leak from a Pauli count: section 5's "7 of 15 two-qubit Pauli errors
# preserve the excitation count ON a differing pair: IZ/ZI/ZZ outright,
# XX/XY/YX/YY on the one-magnon pair subspace only". That sentence has
# already SORTED the survivors into the two classes the hop fraction is
# about, and the document then treats the split as a free modelling
# choice. It is not free: the same count fixes it.
#
# Recomputed exactly on the six-site one-magnon sector (every two-qubit
# Pauli on every bond on every basis state, classified as left-sector /
# moved / stayed), which reproduces BOTH numbers section 10 states:
#   excitation ON the acted bond   : 80/150 leave  -> f_leak = 8/15
#                                    of the 70 left, 40 move, 30 do not
#                                                   -> hop    = 4/7
#   excitation OFF the acted bond  : 240/300 leave -> f_leak = 4/5
#                                    of the 60 left, 0 move  -> hop = 0
# So with q = the fraction of gate errors landing on a bond that
# currently carries the excitation, BOTH parameters are functions of
# that ONE unknown:
#       f_leak(q) = (12 - 4q)/15        hop(q) = 4q/(3 + 4q)
# q = 1 gives (0.5333, 0.5714), q = 0 gives (0.8000, 0.0000).
#
# WHY THIS MATTERS FOR EVERY THRESHOLD HERE. The flown bracket is a
# PRODUCT BOX over the two, so two of its four corners lie off the
# curve, and one of them, (f_leak = 8/15, hop = 0), is as far off it as
# the box allows: at f_leak = 8/15 the count says hop = 4/7, and hop = 0
# requires f_leak = 4/5. That corner is also the one at which the
# floor's existence condition fails (P(s2C < theta_D | H1) = 0.0140
# against the registered 0.01), and it is a corner the design's own
# arithmetic excludes. The corner that SETS theta_D, (8/15, 0.5), sits
# near the q = 1 end and is physical.
# Honest fence: this is the uniform two-qubit depolarizing model. A
# device with extra Z-flavoured error (coherent ZZ, 1/f) sits at lower
# hop AND lower f_leak, i.e. it moves ALONG the curve's spirit and off
# its low-f_leak end, never into the anti-correlated corner.
def leak_hop_at(q):
    """(f_leak, hop) at bond-occupancy fraction q, from the count."""
    if not 0.0 <= q <= 1.0:
        raise ValueError(f"q must lie in [0, 1], got {q}")
    return (12.0 - 4.0 * q) / 15.0, 4.0 * q / (3.0 + 4.0 * q)

# Readout confusion of the synthetic device, per qubit and ASYMMETRIC.
# The gate's counts path is symmetric at 1.5%; the analyze side already
# implements asymmetry and section 8a lists "asymmetric readout
# confusion" as outstanding gate work, so the harness carries it.
EPS0 = np.array([0.017, 0.015, 0.015, 0.016, 0.014, 0.016])  # P(1|true 0)
EPS1 = np.array([0.015, 0.017, 0.015, 0.016, 0.014, 0.016])  # P(0|true 1)

# THE READOUT LADDER, and the CAL that inverts it (2026-08-18).
# Section 8a lists "the readout-confusion ladder to the 2% guard bound
# (the G1 model used symmetric 1.5%)" as OWED. EPS_SCALE walks it: the
# Class-1 line rule admits readout <= 2%, so the pinned ~1.5% above is
# not the admissible worst case, and every other unknown-at-flight-time
# parameter in this design (f_leak, hop) is bracketed with the MAXIMUM
# registered.
EPS_SCALE = 1.0              # 1.0 = the pinned ~1.5%; 4/3 reaches 2%
# CAL_SHOTS closes a second gap, and this one made the refreeze blind
# to the risk it exists for. The harness inverted with Minv =
# inv(Mfwd), the EXACT inverse of the generating matrix, so readout
# mitigation was perfect by construction and contributed zero bias.
# The flight builds Minv from CAL PUBs: Amendment 1.5 flies one head
# pair at 16384 plus per-job pairs at 4096, pooled, so each eps is
# estimated on ~57k shots and carries a shot-noise SE of ~5e-4. That
# matters here far more than its size suggests, because at the deep end
# 92-99% of the mass the inversion operates on is OUT of sector, so a
# mis-specified eps injects contamination of order eps_err*(1-s)
# against an in-sector mass s: site-dependent and depth-growing, which
# is the shape of the signal itself. One CAL estimate is drawn per REP,
# since one flight buys one calibration. 0 = the old perfect inversion.
# WHICH END OF AMENDMENT 1.5's BRACKET (2026-08-18, review finding).
# The job count is payload-determined at submit, "8 by size alone and
# up to ~11", giving 7 to 10 CAL PAIRS after the head pair. The first
# version of this constant used 10 pairs = 57344 shots, which is the
# end section 10 registers FOR COST, where more PUBs is the
# conservative direction. In a PRECISION calculation the direction
# reverses: fewer CAL shots means a noisier inversion, a wider H0 and a
# LARGER theta_D. Taking the cost end here would have flattered a
# frozen number, in the repair introduced to stop this harness
# flattering itself. 7 pairs is the precision-worst end, and it is
# bracketed and maximised over exactly as f_leak and hop are.
CAL_SHOTS = 45056            # 16384 head + 7 pairs x 4096, per state

# WHERE A SECTOR-LEAVING SHOT LANDS. The gate DISCARDS it: the shot is
# never counted, so the cell simply has fewer shots. Hardware MEASURES
# it: the shot lands on a bitstring of the wrong magnon number, is
# removed by post-selection, and still passes through the confusion
# inversion, where out-of-sector mass leaks back into the one-magnon
# entries. Neither is obviously right for a synthetic refreeze, so both
# run and the campaign REPORTS the difference instead of hiding it.
OUT_OF_SECTOR = "measured"   # "measured" | "dropped"

# WHICH AXIS CARRIES THE BEAT PHASE (2026-08-18). MEASURED, REJECTED,
# and kept switchable so the rejection stays reproducible. The envelope
# counts the dephasing DOSE, n-1 injection layers at grid depth n
# (section 5); the cosine counts the Trotter STEPS, n. The committed
# estimator uses one axis for both, so throwing the section-5 axis
# switch also shifts the cosine. Every depth n >= 1 shifts by exactly
# one step, and a uniform shift is absorbed by the fit's free amplitude
# and phase, so the whole effect sits on depth 0: the switch repairs a
# ~5% dose error there and introduces a 1 - cos(om*dt) = 4.5% phase
# error in its place. That is section 8a's "the registered axis
# correction COSTS detection", explained.
# "dose" is the committed single-axis behaviour; "nominal" gives the
# cosine its own counter. The verdict, from the gate's deterministic
# NOISELESS construction against the closed-form rate triples: the
# two-axis fit has the LARGER residual on all three arms and inflates
# s2(C) further from the bare ideal. Its +15% on d and +8% on power
# (200-rep paired) is section 6.3's inflation trap, not accuracy. See
# run_corner_beat.fit_rate's ts_phase note for the numbers.
PHASE_AXIS = "dose"          # "dose" | "nominal"

# THE NULL THE FLIGHT IS ACTUALLY PAYING TO EXCLUDE (2026-08-18, from
# the physics review). The two registered nulls both fly the C arm
# UNIFORM, so they measure the estimator-basis artifact and nothing
# else. But the only arm-dependent quantity anywhere in this model is
# the per-site rate profile: the compression of diag(gamma) onto the
# three-dimensional room is scalar for the uniform profile ONLY because
# sum_l psi_i(l)^2 psi_j(l)^2 = 1/6 on all three room dyads, and ANY
# departure from uniform breaks that generically. Nothing in the chain
# from profile to d uses the mirror-transversal structure, the R90
# locus, or the frozen divisor. So D-sign, on its own, tests that a
# NON-UNIFORM dephasing profile makes the compressed 3x3 non-scalar,
# which is weaker than it reads. Section 1 already concedes the
# registered pair's scope in words; this constant makes the concession
# MEASURABLE: an equal-dose profile that is neither a mirror
# transversal nor on the locus (3*gbar on sites 0 and 1, dark
# elsewhere: total dose 6*gbar as every engineered arm, but
# gamma_0 + gamma_5 = 3*gbar != 2*gbar, so off the locus). Flown as C
# and analyzed with the C estimator, its d says how much of the
# detection the class physics is actually responsible for. If it clears
# theta_D, that belongs on the page BEFORE the flight, not in the
# record afterwards.
# MATCH THE SECOND MOMENT, NOT THE DOSE (2026-08-18, review finding).
# The first version of this control was [3,3,0,0,0,0]: equal dose and
# off the locus, but profile VARIANCE 2.0 against the corner arms' 1.0.
# Dose is the wrong thing to match. What breaks the scalar compression
# is the profile's second moment on the room dyads, so a control with
# twice the inhomogeneity over-states the share generic inhomogeneity
# carries, and that share is the whole point of the measurement.
# [2,2,0,0,2,0]: dose 6, variance 1.000 (equal to CORN and CORNP),
# pair sums (2,4,0) so OFF the locus, and not a transversal since the
# pair {1,4} is doubly lit. Both are kept and both are run; the pair is
# the honest reading, since the first also brackets a stronger
# perturbation.
#
# ... AND MATCHING THE VARIANCE WAS ALSO THE WRONG QUANTITY. Computed
# 2026-08-18: the mirror-odd part of CORN is (1,-1,-1,1,1,-1), and the
# projection of each control onto it is
#     CORN                1.000   (the flown C, by definition)
#     [3,3,0,0,0,0]       0.000   <- the dose-matched control
#     [2,2,0,0,2,0]       0.333   <- the variance-matched control
# so the first control was ORTHOGONAL to the very direction the C
# estimator reads as a split, and the second carries a third of it.
# G4 registered this law before either was built: "the first-order
# background sensitivity is the PROJECTION of the background profile
# onto C's mirror-odd part, NOT its scale". Dose and variance are both
# the wrong match; the projection is the coordinate.
# THE LADDER, therefore, instead of a single control: p(t) = UNIF +
# t*odd(CORN), t in {0, 1/3, 2/3, 1}, every one at dose 6 and
# non-negative, t = 1 being CORN itself and t = 0 the uniform arm. If d
# is a function of t alone, that IS what D-sign tests, stated exactly,
# and the concession section 1 makes in words acquires its curve. The
# fifth entry carries the same t = 2/3 OFF the locus (an even
# perturbation, which is what leaving the locus means: the locus is
# exactly the profiles whose even part is constant), to separate the
# locus from the projection.
def _profile_at(t, even=None):
    odd_c = (gate.CORN - gate.CORN[::-1]) / 2
    p = gate.UNIF + t * odd_c + (0.0 if even is None else even)
    assert abs(p.sum() - 6.0) < 1e-12 and p.min() >= -1e-12, p
    return p


_EVEN = np.array([1 / 3, -1 / 3, 0.0, 0.0, -1 / 3, 1 / 3])
OFF_LOCUS_PROFILES = {
    "t=0 (uniform)": _profile_at(0.0),
    "t=1/3": _profile_at(1 / 3),
    "t=2/3": _profile_at(2 / 3),
    "t=1 (CORN itself)": _profile_at(1.0),
    "t=2/3 OFF locus": _profile_at(2 / 3, _EVEN),
}
OFF_LOCUS = OFF_LOCUS_PROFILES["t=2/3 OFF locus"]


def nominal_time_axis():
    """t = n*dt on the frozen depth grid: the axis the beat phase runs
    on, whatever the dose axis is doing."""
    return np.arange(0, STEPS + 1, K) * JDT


def confusion_forward(eps0=EPS0, eps1=EPS1):
    """The 64x64 true -> measured map, as the kron the runner inverts."""
    M = np.array([[1.0]])
    for q in range(N - 1, -1, -1):
        Cq = np.array([[1 - eps0[q], eps1[q]],
                       [eps0[q], 1 - eps1[q]]])
        M = np.kron(M, Cq)
    return M


def popcount_table():
    return np.array([bin(i).count("1") for i in range(D6)])


_POP = popcount_table()
_ONE_MAGNON = np.array(rcb.CFG_1M)
_OFF_SECTOR = np.array([i for i in range(D6) if _POP[i] != 1])


# ------------------------------------------------- RE-ENTRANT SHOTS
# WHAT THE SURVIVAL MODEL LEAVES OUT (2026-08-18). Everywhere in this
# design, sector loss is a scalar survival exp(-eps*f_leak*K*t) and a
# shot that leaves the one-magnon sector is gone. That is a
# single-error picture. The deep end has eps_step*STEPS = 0.08*60 =
# 4.85 EXPECTED errors, so multi-error is the regime, and a magnon
# number that random-walks under a depolarizing channel can walk BACK
# to 1. Post-selection then keeps that shot, because post-selection
# reads the final magnon number and cannot read the history.
#
# Measured as an exact Markov chain on the magnon number, built from
# the same uniform depolarizing model the f_leak count comes from:
#     depth  3: kept 0.841, never left 0.839, re-entrant  0.3%
#     depth 24: kept 0.325, never left 0.245, re-entrant 24.5%
#     depth 36: kept 0.228, never left 0.121, re-entrant 46.6%
#     depth 60: kept 0.147, never left 0.030, re-entrant 79.8%
# Two consequences and they point opposite ways. The kept COUNT is
# UNDER-stated (0.147 against the scalar model's 0.077 at f_leak =
# 8/15), which is margin the deep-end floor did not know it had. And
# four fifths of the deepest kept shots carry no beat at all, a
# pedestal that GROWS with depth, which dilutes every arm's split and
# therefore d.
#
# THE ASSUMPTION, stated because it is the weak point: a re-entrant
# has suffered at least two sector-changing errors, so it is modelled
# as UNIFORM over the six sites and phase-less. That is the
# maximum-entropy end and also the maximal-dilution end, so it reads
# as an upper bound on the systematic rather than a best estimate.
#
# WHAT IT ACTUALLY DOES TO THE ESTIMATOR, and it is not what "pedestal"
# suggests. A uniform site vector projects to EXACTLY ZERO on all three
# registered dyads, measured 1e-15: arm_rates reads yd = 6*Yp @ A with
# A[l, (i,j)] = W[l,i]*W[l,j], and sum_l W[l,i]W[l,j] is the Floquet
# modes' orthonormality, delta_ij, while every registered dyad (1,2),
# (2,3), (3,5) is OFF-diagonal. So the re-entrants are invisible in the
# numerator and present in the denominator: the dyad amplitude PER KEPT
# SHOT falls by (1 - share(t)). The two modes below therefore produce
# the SAME traces and differ only in the fit weights, which is why
# their measured effects agree to a couple of points.
# And the loss is not a constant added rate. Its effective rate,
# -ln(1-share)/t, climbs 0.0065 -> 0.0370 -> 0.0781 -> 0.1163 ->
# 0.1777 across the grid. A constant added rate would leave the
# VARIANCE of the three fitted rates untouched; a depth-growing one
# distorts the shape, which is why s2 moves at all and why rbar rises
# while s2 falls.
# Nowhere in the pre-registration does re-entry appear; the words
# "re-enter", "multi-error" and "Markov" are absent from all 3760
# lines, and section 10 uses the 4.85 exponent only inside
# exp(-4.85*f_leak).
REENTRY = False              # exploration switch; the freeze runs False
# "add" keeps the beat mass and adds the pedestal beside it (total
# counts rise, deep points gain fit weight); "split" keeps the total
# and changes only the composition. The truth is between them: the
# chain says kept = 0.147 and beat-carrying = 0.030 at the deep end,
# and the harness's own survival, 0.077, is neither. Running both
# brackets the systematic instead of picking the flattering end.
REENTRY_MODE = "add"         # "add" | "split"
# THE CORRECTION, separate from the injection on purpose. The
# re-entrants reach the estimator through exactly one channel, the
# normalisation, so folding the KNOWN curve (1 - f(t)) into the fit
# model beside the fitted exponential removes them without a free
# parameter. f(t) comes from the calibrated p2 through the chain.
# Kept switchable so the three-way read stays available: data clean,
# data contaminated, data contaminated AND corrected.
REENTRY_FIT = False          # fold (1 - f(t)) into the fit model


_CHAIN_CACHE = {}


def reentrant_share(eps_step=None, steps=None, k_grid=None):
    """Per grid depth, the fraction of the KEPT shots that left the
    one-magnon sector at some point and came back.

    1 - never_left/kept from the chain below. This is the only number
    taken from the chain, deliberately: see raw_counts_for_cell."""
    key = (8 * GPB * P2 if eps_step is None else eps_step,
           STEPS if steps is None else steps,
           K if k_grid is None else k_grid)
    # MEMOISED. The chain is a function of the configuration only, and
    # one_rep asked for it per rep: 33600 python-level iterations x 400
    # reps x 2 campaigns tripled the A/B's wall clock for nothing.
    if key not in _CHAIN_CACHE:
        kept, never = magnon_chain(*key)
        _CHAIN_CACHE[key] = 1.0 - never / np.maximum(kept, 1e-300)
    return _CHAIN_CACHE[key]


def magnon_chain(eps_step, steps, k_grid):
    """(kept, never_left) per grid depth under one depolarizing error
    per step at weight eps_step, as an exact chain on the magnon number.

    Returns two arrays over the grid depths, both starting at 1.0."""
    bonds = [(b, b + 1) for b in range(N - 1)]
    T = np.zeros((N + 1, N + 1))
    for n in range(N + 1):
        for s in range(2 ** N):
            if bin(s).count("1") != n:
                continue
            bits = [(s >> (N - 1 - l)) & 1 for l in range(N)]
            for (a, b) in bonds:
                for p0 in "IXYZ":
                    for p1 in "IXYZ":
                        if p0 == p1 == "I":
                            continue
                        nb = bits.copy()
                        if p0 in "XY":
                            nb[a] ^= 1
                        if p1 in "XY":
                            nb[b] ^= 1
                        T[n, sum(nb)] += 1.0
        T[n] /= T[n].sum()
    p = np.zeros(N + 1)
    p[1] = 1.0
    nl = 1.0
    kept, never = [1.0], [1.0]
    for t in range(1, steps + 1):
        p = (1 - eps_step) * p + eps_step * (p @ T)
        nl *= (1 - eps_step * (1 - T[1, 1]))
        if t % k_grid == 0:
            kept.append(float(p[1]))
            never.append(nl)
    return np.array(kept), np.array(never)


def raw_counts_for_cell(occ_t, surv, shots_bind, rng, Mfwd,
                        out_of_sector=None, reentry=None):
    """Per-binding RAW 64-dim counts for one (arm, prep) cell.

    occ_t: [T, M, N] per-binding one-magnon occupations from the gate's
    run_bindings. Returns [T, M, 64] integer counts, which is exactly
    what the runner's cell_raw_matrix produces from hardware.

    reentry, when given, is the (kept, never_left) pair from
    magnon_chain: the in-sector mass becomes kept[t] instead of surv[t],
    split into never_left[t] carrying the beat and the remainder
    carrying a flat, phase-less pedestal.
    """
    oos = OUT_OF_SECTOR if out_of_sector is None else out_of_sector
    T, M, _ = occ_t.shape
    out = np.zeros((T, M, D6), dtype=np.int64)
    off_p = np.zeros(D6)
    off_p[_OFF_SECTOR] = 1.0 / len(_OFF_SECTOR)
    flat = np.zeros(D6)
    flat[_ONE_MAGNON] = 1.0 / len(_ONE_MAGNON)
    for t in range(T):
        p = np.clip(occ_t[t], 0, None)
        p = p / np.maximum(p.sum(axis=1, keepdims=True), 1e-300)
        s = float(min(surv[t], 1.0))
        true = np.zeros((M, D6))
        if reentry is None:
            true[:, _ONE_MAGNON] = s * p
        else:
            # ISOLATE THE RE-ENTRY, CHANGE NOTHING ELSE. The chain
            # derives its own leak from the depolarizing model while
            # the harness treats f_leak as a free knob, so using the
            # chain's kept curve directly would move two things at
            # once and the A/B would measure their sum. Only the
            # chain's SHARE is taken, and it is applied to the
            # harness's own survival: the beat-carrying mass stays
            # exactly s*p, and a pedestal is added beside it in the
            # ratio the chain says.
            share = float(reentry[t])
            if REENTRY_MODE == "add":
                # keep the beat-carrying mass, ADD the pedestal beside
                # it. Total kept rises, so deep points gain fit weight:
                # anti-conservative on the weighting, an upper edge.
                true[:, _ONE_MAGNON] = s * p
                ped = s * share / max(1.0 - share, 1e-12)
                true += ped * flat[None, :]
                s = s + ped
            elif REENTRY_MODE == "split":
                # keep the TOTAL kept mass and split it. The weighting
                # is untouched and only the composition changes:
                # conservative on the weighting, a lower edge.
                true[:, _ONE_MAGNON] = s * (1.0 - share) * p
                true += s * share * flat[None, :]
            else:
                raise RuntimeError(
                    f"unknown re-entry mode {REENTRY_MODE!r}")
        if oos == "measured":
            true += (1.0 - s) * off_p[None, :]
        elif oos != "dropped":
            raise RuntimeError(f"unknown out-of-sector mode {oos!r}")
        meas = true @ Mfwd.T
        # an explicit "not recorded" column keeps the multinomial exact
        # when leaked shots are dropped (rows then sum to s, not 1)
        lost = np.clip(1.0 - meas.sum(axis=1), 0.0, None)
        pv = np.concatenate([meas, lost[:, None]], axis=1)
        pv = pv / pv.sum(axis=1, keepdims=True)
        draw = rng.multinomial(shots_bind, pv)
        out[t] = draw[:, :D6]
    return out


def arm_traces(profile, Wm, U, err_deph_dt, p_hop, surv, rng, Mfwd, Minv,
               floor, out_of_sector=None, counts_chain="analyze",
               reentry=None):
    """One arm, all three preparations.

    counts_chain switches ONLY the layer the refreeze exists to change,
    with identical physics and identical draws upstream of it:
      "analyze" -> per-binding RAW 64-dim counts, pooled, ONE confusion
                   inversion, ONE clip, kept-count floor (the flown
                   chain, the runner's own functions);
      "gate"    -> the committed gate's cell_counts, which applies
                   readout as a per-shot keep/reroute and never inverts
                   (the chain the frozen numbers were measured on).
    Running both is how a difference gets attributed to the CHAIN
    rather than to this file having a bug in it.
    """
    shots_bind = max(1, SHOTS // M_BIND)
    Yp, Wp, clipped = [], [], 0.0
    for (i, j) in DY:
        prep = (Wm[:, i] + Wm[:, j]) / np.sqrt(2)
        ph = gate.phase_tables(profile * JDT, STEPS, M_BIND, rng)
        occ = gate.run_bindings(U, JDT, err_deph_dt, p_hop, STEPS, K,
                                prep, ph, rng)
        if counts_chain == "gate":
            Y, W = gate.cell_counts(occ, surv, SHOTS, M_BIND, rng)
            Yp.append(Y); Wp.append(W)
            continue
        if counts_chain != "analyze":
            raise RuntimeError(f"unknown counts chain {counts_chain!r}")
        raw = raw_counts_for_cell(occ, surv, shots_bind, rng, Mfwd,
                                  out_of_sector, reentry)
        Y, W = [], []
        for t in range(raw.shape[0]):
            occ_t, w_t, cl = rcb.pool_invert_postselect(
                raw[t].astype(float), Minv, floor=floor)
            Y.append(occ_t); W.append(w_t); clipped += cl
        Yp.append(np.array(Y)); Wp.append(np.array(W))
    return Yp, Wp, clipped


def one_rep(rng, f_leak, hop_frac, null=False, out_of_sector=None,
            arms=("U", "C", "Cp"), counts_chain="analyze"):
    """One realization of the whole flight through the analyze chain.

    null selects WHICH null, and the two are not the same H0:
      False    -> H1, every arm on its own profile.
      "C"      -> the C arm flies the UNIFORM profile and is ANALYZED
                  with the C estimator. This is the H0 the committed
                  gate implements and it is sufficient for theta_D,
                  whose statistic reads {C, U} only.
      "both"   -> C AND C-prime both fly uniform. Section 8 registers
                  this as the H0 that the W false rate and
                  P(false CONFIRMED) require, "the BOTH-corner-arms-null
                  H0", and it is the only correct null for theta_W:
                  under the C-only null the C-prime arm still splits, so
                  d_W would be frozen against a live alternative rather
                  than against no class physics at all.
    """
    gbar = 1.0 / Q
    U = gate.strang(JDT)
    em, Wm = gate.floquet_modes(U, JDT)
    dyf = [em[i] - em[j] for (i, j) in DY]
    A = np.array([[Wm[l, i] * Wm[l, j] for (i, j) in DY] for l in range(N)])

    Tn = STEPS // K + 1
    eps_step = 8 * GPB * P2
    surv = np.exp(-eps_step * f_leak * K * np.arange(Tn))
    eps_in = eps_step * (1 - f_leak)
    # the within-sector residue, split by hop_frac between the Z-like
    # face and the hop face (the gate pins 0.5; the bracket is registered)
    err_deph_dt = gate.GCOUNT / gate.GCOUNT.mean() * (eps_in * (1 - hop_frac)) / 4
    p_hop = eps_in * hop_frac

    # The device's TRUE confusion, and the estimate the analysis has of
    # it. They are the same object only when CAL_SHOTS = 0.
    eps0 = EPS0 * EPS_SCALE
    eps1 = EPS1 * EPS_SCALE
    Mfwd = confusion_forward(eps0, eps1)
    if CAL_SHOTS:
        # one CAL estimate per rep: a flight buys one calibration, and
        # its error is a systematic across the whole grid, not noise
        # that averages away over depths
        e0 = eps0 + rng.normal(
            0.0, np.sqrt(eps0 * (1 - eps0) / CAL_SHOTS), size=N)
        e1 = eps1 + rng.normal(
            0.0, np.sqrt(eps1 * (1 - eps1) / CAL_SHOTS), size=N)
        Minv = np.linalg.inv(confusion_forward(np.clip(e0, 1e-6, 0.5),
                                               np.clip(e1, 1e-6, 0.5)))
    else:
        Minv = np.linalg.inv(Mfwd)
    floor = rcb.KEPT_COUNT_FLOOR
    reentry = reentrant_share() if REENTRY else None
    env_known = (1.0 - reentrant_share()) if REENTRY_FIT else None
    ts = rcb.fit_time_axis()
    r_max = rcb.R_MAX_FIT

    prof = {"U": gate.UNIF * gbar, "C": gate.CORN * gbar,
            "Cp": gate.CORNP * gbar, "N0": np.zeros(N)}
    if null == "C":
        prof["C"] = gate.UNIF * gbar
    elif null == "both":
        prof["C"] = gate.UNIF * gbar
        prof["Cp"] = gate.UNIF * gbar
    elif null == "offlocus":
        prof["C"] = OFF_LOCUS * gbar
    elif null:
        raise RuntimeError(
            f"unknown null {null!r}: use False, 'C', 'both', 'offlocus'")

    res = {}
    for a in arms:
        Yp, Wp, cl = arm_traces(prof[a], Wm, U, err_deph_dt, p_hop, surv,
                                rng, Mfwd, Minv, floor, out_of_sector,
                                counts_chain, reentry)
        ts_phase = (None if PHASE_AXIS == "dose" else nominal_time_axis())
        rates, dws, sses = rcb.arm_rates(Yp, Wp, A, dyf, ts, a, r_max,
                                         ts_phase, env_known)
        res[a] = {"rates": rates, "dws": dws, "sse": sses,
                  "s2": rcb.s2_of(rates), "rbar": float(np.mean(rates)),
                  "clipped": cl,
                  "floored": int(sum(int(w == 0.0) for W in Wp for w in W))}
    if "C" in res and "U" in res:
        res["d"] = res["C"]["s2"] - res["U"]["s2"]
    if "C" in res and "Cp" in res:
        res["d_W"] = res["C"]["s2"] - res["Cp"]["s2"]
    return res


def campaign(reps, f_leak, hop_frac, null, seed, out_of_sector=None,
             arms=("U", "C", "Cp"), progress=None,
             counts_chain="analyze"):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(reps):
        rows.append(one_rep(rng, f_leak, hop_frac, null=null,
                            out_of_sector=out_of_sector, arms=arms,
                            counts_chain=counts_chain))
        if progress and (i + 1) % progress == 0:
            print(f"      {i + 1}/{reps}", flush=True)
    return rows


def _stats(vals):
    v = np.array([x for x in vals if x is not None and np.isfinite(x)])
    if v.size == 0:
        return {"n": 0, "mean": float("nan"), "sd": float("nan")}
    return {"n": int(v.size), "mean": float(v.mean()),
            "sd": float(v.std(ddof=1)) if v.size > 1 else float("nan")}


def three_sigma(rows, key):
    """The registered shape: E[x|H0] + 3*sigma_hat, ONE formula."""
    s = _stats([r.get(key) for r in rows])
    return s["mean"] + 3 * s["sd"], s


def _col(rows, key, sub=None):
    """The per-rep column, persisted so a threshold can be re-read
    without re-running the campaign (the freeze re-evaluates every rate
    at the REGISTERED maximum, not at the corner's own value)."""
    out = []
    for r in rows:
        v = r.get(key)
        if sub is not None and isinstance(v, dict):
            v = v.get(sub)
        out.append(None if v is None or not np.isfinite(v) else float(v))
    return out


def _live(rows, key, sub=None):
    """The finite draws of a column: the rate's real denominator."""
    return [v for v in _col(rows, key, sub) if v is not None]


def _nan_count(rows, key, sub=None):
    return int(sum(1 for v in _col(rows, key, sub) if v is None))


# How many TRACES each arm's three channels are fitted on. U and N0 fit
# three single-trace channels; C and C-prime fit one single-trace
# channel (the decoupled dyad) and two channels jointly over two
# preparations. Read off section 6.3's channel structure and confirmed
# against run_corner_beat.arm_rates.
CHANNEL_TRACES = {"U": (1, 1, 1), "N0": (1, 1, 1),
                  "C": (1, 2, 2), "Cp": (1, 2, 2)}


def _fit_col(rows, what, traces=None):
    """Per-fit health readouts pooled over every arm and channel of
    every rep: the |dw| excursions, the weighted residuals, or the
    fitted rates themselves.

    traces, when given, keeps only the channels fitted on that many
    traces. An SSE is a SUM over its traces, so pooling the two
    populations makes a bimodal distribution no single bound fits."""
    out = []
    for r in rows:
        for a in ("U", "C", "Cp"):
            v = r.get(a)
            if not isinstance(v, dict):
                continue
            arity = CHANNEL_TRACES[a]
            for ch in range(3):
                if traces is not None and arity[ch] != traces:
                    continue
                if what == "dws":
                    out.append(abs(v["dws"][ch]))
                elif what == "sse":
                    out.append(v["sse"][ch])
                elif what == "rates":
                    out.append(v["rates"][ch])
    return out


def _quant(vals):
    """The distribution a margin gets frozen from, not just its top:
    a margin set at the maximum of a finite sample is a sample size in
    disguise."""
    v = np.array([x for x in np.asarray(vals, dtype=float)
                  if np.isfinite(x)])
    if v.size == 0:
        return None
    return {"n": int(v.size), "mean": float(v.mean()),
            "sd": float(v.std(ddof=1)) if v.size > 1 else None,
            "p50": float(np.percentile(v, 50)),
            "p99": float(np.percentile(v, 99)),
            "p999": float(np.percentile(v, 99.9)),
            "max": float(v.max())}


def _rate(vals, thr, above=True):
    """An exceedance rate over the LIVE draws only.

    _stats drops non-finite values from a threshold while np.mean over a
    boolean list silently counts a NaN rep as "did not exceed" and keeps
    it in the denominator, so threshold and rate ran over different
    populations. A dead rep is not a non-detection; it is not a draw.
    The count of what was dropped is reported beside the rate."""
    v = np.array([x for x in np.asarray(vals, dtype=float)
                  if np.isfinite(x)])
    if v.size == 0:
        return None
    return float(np.mean(v > thr if above else v < thr))


# =====================================================================
# THE s2(C) FLOOR, REGISTERED AS A FUNCTION OF THE ARM'S OWN s2(U)
# =====================================================================
# WHY IT IS NOT A SCALAR (Tom's decision, 2026-08-18). Section 7
# registers the floor as the MAXIMUM over the f_leak and hop brackets
# AND subject to P(s2C < floor | H1) <= 0.01. Measured on the
# 2026-08-18 four-corner run, the maximum (0.00389) breaches its own
# reachability at three of the four corners, 0.064 / 0.046 / 0.074
# against the registered 0.01. The cause is structural, not a tuning
# accident: the SIGNAL is corner-robust (s2C under H1 spans a factor
# 1.23 across the bracket) while the estimator NOISE spans 9.59, so a
# floor calibrated on the loudest corner sits inside the signal at the
# quiet ones. Taking the MINIMUM instead passes reachability and leaves
# the floor below the H0 level exactly where it is needed. Both routes
# fail, so the floor becomes a function of a quantity the flight
# measures.
#
# WHY s2(U) AND NOT THE BRACKET PARAMETERS THEMSELVES. The honest
# regressor would be (f_leak, hop), since those ARE the unknowns the
# bracket brackets. Measured and rejected, both from the code:
#   * they enter the model only through the PRODUCT
#     p_hop = eps_step*(1 - f_leak)*hop_frac (one_rep above), so no
#     observable separates them;
#   * the site resolution that would expose a hop is destroyed at
#     analyze step 2, `yd = [6 * Yp[p] @ A ...]` in arm_rates, which
#     compresses six site occupations into three dyads before the first
#     fit runs. Every verdict statistic is a function of 12 fitted
#     rates and none of them carries a site index.
# And partial identification does not help: with f_leak known and hop
# bracketed the worst reachability is 0.064, with hop known and f_leak
# bracketed 0.046. Only knowing BOTH reaches 0.004, which is what the
# per-corner floor already is.
#
# WHY THE COEFFICIENT IS NOT 1.67. The affine form fitted to the four
# corners' H0 EXPECTATIONS gives floor ~ 0.00046 + 1.67*s2U. Applied
# per rep it is worse than the scalar it replaces (reachability 0.070
# at the loud corner), because b > 1 AMPLIFIES s2(U)'s own draw noise,
# and a single s2(U) has two degrees of freedom and a CV of 0.43-0.99.
# The two things s2(U) carries pull in opposite directions: BETWEEN
# corners it tracks the noise level, which is what the floor must
# follow; WITHIN a corner it fluctuates, and a floor that follows THAT
# inverts the guard section 7 built it for (a low s2(U) raises d and
# would lower the floor at the same time). Fitting the coefficients to
# the registered CRITERIA instead of to the expectations is what keeps
# b small enough that the between-corner signal survives and the
# within-corner noise does not take over.
#
# THE SELECTION IS ITSELF PRE-REGISTERED, and this block is its text:
# a fixed grid, a fixed objective, a fixed admissibility criterion, no
# judgement at freeze time. The coefficients come out of the FIT
# campaigns and are read on HELD-OUT ones, because they are estimated
# quantities and section 8a's held-out rule applies to them exactly as
# it does to theta_D.
FLOOR_B_GRID = tuple(i * 0.25 for i in range(7))       # 0.00 .. 1.50
FLOOR_A_GRID = tuple(i * 2.5e-5 for i in range(161))   # 0 .. 0.004
FLOOR_REACH_CEILING = 0.01     # section 7's registered rate
FLOOR_REACH_ALPHA = 0.05       # the bound's confidence level


def floor_value(s2U, coef):
    """The registered floor at one rep's own measured s2(U).

    coef = (a, b); b = 0 reproduces a plain scalar floor exactly, so
    the old registration is the degenerate member of the new one."""
    a, b = coef
    return a + b * s2U


def binom_upper_bound(k, n, alpha=FLOOR_REACH_ALPHA):
    """One-sided Clopper-Pearson upper bound on a binomial rate.

    The root of P(X <= k; n, p) = alpha, found by bisection on the exact
    tail rather than through a special-function library, so the value
    can be checked against the tail it inverts (and, at k = 0, against
    the closed form 1 - alpha^(1/n))."""
    if n <= 0:
        raise ValueError("n must be positive")
    if k >= n:
        return 1.0

    def tail(p):
        return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i)
                   for i in range(k + 1))

    lo, hi = 0.0, 1.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if tail(mid) > alpha:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def reach_admissible(k, n, ceiling=FLOOR_REACH_CEILING,
                     alpha=FLOOR_REACH_ALPHA):
    """Is the reachability rate BELOW the ceiling, as a bound?

    Not as a point estimate: at n = 500 a true rate of 0.01 has a
    standard error of 0.0044, so an observed 0.010 is equally
    consistent with 0.02. The criterion is the confidence bound, which
    at n = 500 and these constants admits at most ONE breach in 500.
    The point-estimate short circuit is exact, never a speedup with a
    different answer: the bound is always at or above the estimate."""
    if n <= 0:
        return False
    if k / n > ceiling:
        return False
    return binom_upper_bound(k, n, alpha) <= ceiling


def min_n_for_ceiling(ceiling=FLOOR_REACH_CEILING,
                      alpha=FLOOR_REACH_ALPHA):
    """The smallest sample at which the ceiling is DEMONSTRABLE at all.

    Below it, even a perfect campaign (zero breaches in every corner)
    cannot put the upper bound under the ceiling, so `reach_admissible`
    refuses everything and the refusal says nothing about the floor.
    At the registered 0.01 / 0.05 this is 299, from
    1 - alpha^(1/n) <= ceiling, so the freeze's >= 500 held-out reps
    clear it with room for exactly one breach. A smoke run at 6 reps
    does not, and must say which of the two it hit."""
    n = 1
    while binom_upper_bound(0, n, alpha) > ceiling:
        n += 1
    return n


def floor_trip_counts(s2C, s2U, coef):
    """(trips, live draws) of the floor over one corner's per-rep pairs.

    A dead rep is dropped from BOTH, the same rule _rate applies: a NaN
    rep is a VOID, not a rep that failed to trip."""
    if len(s2C) != len(s2U):
        raise ValueError(
            f"paired columns must have equal length, got "
            f"{len(s2C)} and {len(s2U)}")
    k = n = 0
    for c, u in zip(s2C, s2U):
        if c is None or u is None:
            continue
        c, u = float(c), float(u)
        if not (math.isfinite(c) and math.isfinite(u)):
            continue
        n += 1
        if c < floor_value(u, coef):
            k += 1
    return k, n


def conjunction_rate(d_col, s2C_col, s2U_col, theta_D, coef):
    """P(d > theta_D AND s2C >= floor) over reps live in all THREE.

    The registered detection rule since round 9 is this conjunction, not
    the d exceedance alone, and it can only be evaluated rep by rep now
    that the floor reads the same rep's s2(U). A rep dead in ANY of the
    three columns is dropped from both numerator and denominator: the
    verdict for such a rep is VOID, and a VOID is not a non-detection."""
    if not (len(d_col) == len(s2C_col) == len(s2U_col)):
        raise ValueError(
            f"paired columns must have equal length, got "
            f"{len(d_col)}, {len(s2C_col)} and {len(s2U_col)}")
    k = n = 0
    for dd, sc, su in zip(d_col, s2C_col, s2U_col):
        if dd is None or sc is None or su is None:
            continue
        dd, sc, su = float(dd), float(sc), float(su)
        if not (math.isfinite(dd) and math.isfinite(sc)
                and math.isfinite(su)):
            continue
        n += 1
        if dd > theta_D and sc >= floor_value(su, coef):
            k += 1
    return (k / n) if n else None


def is_provably_inert(coef, theta_D):
    """Can this floor change ANY D-sign verdict, even in principle?

    THE LEMMA (2026-08-18, from the statistics review, re-derived here).
    A rep passes conjunct 1 and fails conjunct 2 exactly when
        s2U + theta_D < s2C < a + b*s2U,
    which needs a - theta_D > (1 - b)*s2U. On this design s2U is a
    VARIANCE, so s2U >= 0, so for b <= 1 the right side is >= 0 and the
    inequality forces a > theta_D. Below that, the floor is not merely
    weak, it is unsatisfiable: it cannot remove one detection, true or
    false, at any corner, on any sample.

    This matters because P(trip | H0) does NOT see it. That rate counts
    reps the floor stops, and under H0 nearly all of them had already
    failed conjunct 1, so a floor with a trip rate of 0.98 can be
    changing nothing at all. Measured on the 2026-08-18 run: over the
    whole registered grid, EVERY member admissible under reachability
    is inert, on 2000 H1 and 2000 held-out H0 draws, and the one member
    that does bite (the superseded scalar 0.00389 > theta_D) removes 64
    true detections to remove 1 false one."""
    a, b = coef
    return a <= theta_D and b <= 1.0


def incremental_guard_counts(d_col, s2C_col, s2U_col, theta_D, coef):
    """(removed, live): reps that pass conjunct 1 and die at the floor.

    THE floor's actual contribution, as opposed to its trip rate. Read
    on H0 draws it is the false detections the floor removes; on H1
    draws it is the true detections it costs."""
    if not (len(d_col) == len(s2C_col) == len(s2U_col)):
        raise ValueError("paired columns must have equal length")
    k = n = 0
    for dd, sc, su in zip(d_col, s2C_col, s2U_col):
        if dd is None or sc is None or su is None:
            continue
        dd, sc, su = float(dd), float(sc), float(su)
        if not (math.isfinite(dd) and math.isfinite(sc)
                and math.isfinite(su)):
            continue
        n += 1
        if dd > theta_D and sc < floor_value(su, coef):
            k += 1
    return k, n


def _weakest_h0_trip(h0_by_corner, coef):
    weakest = 1.0
    for c in h0_by_corner.values():
        k, n = floor_trip_counts(c["s2C"], c["s2U"], coef)
        weakest = min(weakest, (k / n) if n else 0.0)
    return weakest


def select_floor_coefficients(h1_by_corner, h0_by_corner,
                              b_grid=FLOOR_B_GRID, a_grid=FLOOR_A_GRID,
                              ceiling=FLOOR_REACH_CEILING,
                              alpha=FLOOR_REACH_ALPHA, theta_D=None):
    """Choose (a, b) on the FIT campaigns. Deterministic, no judgement.

    ADMISSIBLE: reachability holds at EVERY corner, as a bound. Note
    what this replaces and improves: the max-over-bracket rule needed a
    reading of section 7's "both bracket ends", where the producer
    sentence names two brackets and the reachability sentence one. A
    single function required to be admissible at every corner satisfies
    both readings and needs neither.
    OBJECTIVE: the strongest guard, i.e. the largest WEAKEST-corner
    P(trip | H0), which is section 7 (i)'s registered purpose ("under a
    TRUE NULL the floor trips by design"). Ties break to the smaller b
    and then the smaller a, the more conservative direction on
    reachability.

    Returns None only when no grid member is admissible. On real data
    that does not happen, because s2(C) is a variance and therefore
    (a, b) = (0, 0) trips nothing and is always admissible. So the
    operative check is NOT this refusal, it is the reported guard: a
    selection that lands on the degenerate member has produced a
    decorative floor, which is route (b)'s failure under another name,
    and `degenerate` in the returned dict says so out loud.

    WHY THE OBJECTIVE IS STILL THE MARGINAL TRIP RATE, against the
    statistics review's suggestion to make it the INCREMENTAL guard
    (the false detections the floor actually removes). Tried and
    measured before rejecting: on this design the incremental guard is
    identically 0 for EVERY admissible member (is_provably_inert above
    says why, and the sweep confirms it over the whole grid). An
    objective that is constant does not select, it ties, and the
    tie-break then walks to the smallest b and the smallest a, i.e. it
    registers the do-nothing floor. That is strictly worse than what
    this objective gives, which is a floor POSITIONED at the H0 level,
    the thing section 7 (i) actually asks for ("under a TRUE NULL the
    floor trips by design"). The review's point stands entirely as a
    reading, and it is answered by REPORTING: `provably_inert` and the
    incremental counts travel with every selection so the trip rate can
    never be read as protection it does not give."""
    best = None
    n_adm = 0
    for b in b_grid:
        for a in a_grid:
            coef = (a, b)
            if not all(reach_admissible(
                    *floor_trip_counts(c["s2C"], c["s2U"], coef),
                    ceiling=ceiling, alpha=alpha)
                    for c in h1_by_corner.values()):
                continue
            n_adm += 1
            weakest = _weakest_h0_trip(h0_by_corner, coef)
            key = (weakest, -b, -a)
            if best is None or key > best[0]:
                best = (key, coef, weakest)
    if best is None:
        return None
    _, coef, weakest = best
    worst_h1 = 0.0
    for c in h1_by_corner.values():
        k, n = floor_trip_counts(c["s2C"], c["s2U"], coef)
        worst_h1 = max(worst_h1, (k / n) if n else 1.0)
    return {"coefficients": coef,
            "weakest_h0_trip_fit": weakest,
            "worst_h1_trip_fit": worst_h1,
            "degenerate": coef[0] == 0.0 and weakest == 0.0,
            # NOT a property of the coefficients alone: a floor can be
            # non-degenerate (it trips constantly) and still change no
            # verdict. Both flags are reported because they fail
            # differently and the second is the one that matters.
            "provably_inert": (None if theta_D is None
                               else is_provably_inert(coef, theta_D)),
            "theta_D_used": theta_D,
            "n_admissible_members": n_adm,
            "reach_ceiling": ceiling, "reach_alpha": alpha,
            "b_grid": list(b_grid),
            "a_grid_min": min(a_grid), "a_grid_max": max(a_grid),
            "a_grid_n": len(a_grid)}


def validate_floor_coefficients(coef, h1_by_corner,
                                ceiling=FLOOR_REACH_CEILING,
                                alpha=FLOOR_REACH_ALPHA):
    """Read the chosen coefficients on a HELD-OUT H1 sample.

    The number that matters is `passes`. It is measured on draws the
    selection never saw, because the greedy optimum on this design
    reads 0.008 in-sample and 0.012 held out: the coefficients sit on
    the edge of the admissible region and a 500-draw sample places that
    edge with an error of its own. A freeze whose reachability was only
    ever read in-sample is the theta_D mistake one level up."""
    corners, ok = {}, True
    for tag, c in h1_by_corner.items():
        k, n = floor_trip_counts(c["s2C"], c["s2U"], coef)
        adm = reach_admissible(k, n, ceiling=ceiling, alpha=alpha)
        ok = ok and adm
        corners[tag] = {"trips": k, "n": n,
                        "rate": (k / n) if n else None,
                        "upper_bound": (binom_upper_bound(k, n, alpha)
                                        if n else None),
                        "admissible": adm}
    return {"passes": ok, "coefficients": list(coef),
            "reach_ceiling": ceiling, "reach_alpha": alpha,
            "corners": corners}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=500)
    ap.add_argument("--seed", type=int, default=20260817)
    ap.add_argument("--heldout", type=int, default=0,
                    help="reps in the HELD-OUT null campaigns that "
                         "measure the rates (0 = same as --reps). "
                         "Section 8a registers held-out splits; the "
                         "thresholds keep the full first sample")
    ap.add_argument("--quick", action="store_true",
                    help="one bracket corner only, for a smoke run")
    ap.add_argument("--oos", default=None, choices=["measured", "dropped"])
    ap.add_argument("--eps-scale", dest="eps_scale", type=float,
                    default=None,
                    help="scale the readout confusion; 1.0 is the "
                         "pinned ~1.5%%, 4/3 reaches the 2%% Class-1 "
                         "guard bound (section 8a's owed ladder)")
    ap.add_argument("--cal-shots", dest="cal_shots", type=int,
                    default=None,
                    help="shots per CAL state behind the inversion "
                         "(0 = invert with the truth, the old "
                         "behaviour; default 57344 = Amendment 1.5's "
                         "pooled head + per-job pairs)")
    ap.add_argument("--offlocus", action="store_true",
                    help="run ONLY the off-locus diagnostic: an "
                         "equal-dose inhomogeneous profile that is "
                         "neither a transversal nor on the R90 locus, "
                         "flown as C, to price how much of D-sign the "
                         "class physics actually carries")
    ap.add_argument("--theta-d", dest="theta_d", type=float,
                    default=0.00291,
                    help="the registered theta_D the diagnostic is "
                         "read against (default: the 2026-08-18 "
                         "refreeze value)")
    ap.add_argument("--phase-axis", dest="phase_axis", default=None,
                    choices=["dose", "nominal"],
                    help="which axis carries the beat phase; 'dose' is "
                         "the committed single-axis behaviour, 'nominal' "
                         "gives the cosine its own counter")
    # DEFAULTS TO THE REGISTERED AXIS (2026-08-18, review finding). It
    # used to default to None, which left rcb.TIME_AXIS at its
    # "nominal" module value: the axis section 5 RETIRED, on a harness
    # whose whole reason for existing is the realized-dose refit, with
    # the documented Run: line omitting the flag. Section 8a measured
    # what that costs: theta_D ~14% low and P(detect) ~0.03 high. The
    # reviewer who found this walked into it while smoke-testing, which
    # is the argument for closing it rather than documenting it.
    ap.add_argument("--axis", default="realized_dose",
                    choices=["nominal", "realized_dose"],
                    help="throw the fit-axis switch on BOTH sides at "
                         "once (section 5 registers the refit on "
                         "realized_dose; the gate half landed 2026-08-17)")
    ap.add_argument("--reentry", action="store_true",
                    help="add the RE-ENTRANT pedestal: shots that "
                         "left the one-magnon sector and walked "
                         "back, which post-selection keeps because "
                         "it reads the final magnon number and not "
                         "the history. Off by default because the "
                         "registered chain does not model it; the "
                         "switch exists to MEASURE what leaving it "
                         "out costs, which at the deep end is 80%% "
                         "of the kept shots")
    ap.add_argument("--curve", type=int, default=0, metavar="NPOINTS",
                    help="sweep NPOINTS on the COUNTED curve "
                         "(f_leak, hop) = leak_hop_at(q), q running 1 "
                         "down to 0, INSTEAD of the four corners of the "
                         "product box. The two parameters are not "
                         "independent: the Pauli count section 10 uses "
                         "for f_leak also fixes hop, and two of the "
                         "box's corners lie off the curve, including "
                         "the one where the floor's existence condition "
                         "fails. 5 is the sensible value and costs what "
                         "the four corners cost")
    ap.add_argument("--shots", type=int, default=None,
                    help="science shots per (arm, depth, prep), default "
                         "the runner's SHOTS_SCIENCE. EXPLORATION ONLY: "
                         "the freeze runs the flown value. It exists to "
                         "MEASURE the shot scaling of theta_D and of the "
                         "H1 s2C spread, which is what prices section "
                         "7's design-side remedy; extrapolating that "
                         "scaling instead of measuring it was how a "
                         "1.5x estimate got quoted for what the two "
                         "different exponents make ~1.2x")
    args = ap.parse_args()
    if args.reentry:
        global REENTRY
        REENTRY = True
    if args.shots is not None:
        global SHOTS
        if args.shots < M_BIND:
            raise SystemExit(
                f"--shots {args.shots} is below M_BIND {M_BIND}: the "
                "per-binding count would floor at 1 and the shot "
                "scaling being measured would be the flooring instead")
        SHOTS = args.shots
    if args.phase_axis:
        global PHASE_AXIS
        PHASE_AXIS = args.phase_axis
    if args.eps_scale is not None:
        global EPS_SCALE
        EPS_SCALE = args.eps_scale
    if args.cal_shots is not None:
        global CAL_SHOTS
        CAL_SHOTS = args.cal_shots
    if args.axis:
        rcb.TIME_AXIS = args.axis
        gate.TIME_AXIS = args.axis
        assert np.array_equal(
            rcb.fit_time_axis(),
            gate.fit_time_axis(STEPS, K, JDT)), "axis switch out of sync"

    ends = F_LEAK_ENDS[:1] if args.quick else F_LEAK_ENDS
    hops = HOP_FRACTIONS[:1] if args.quick else HOP_FRACTIONS

    # THE SWEEP SET: the flown product BOX, or the counted CURVE.
    # Same machinery downstream either way; the registered value stays
    # the envelope over whatever set is swept, so narrowing the set is
    # the only change and every threshold moves with it.
    if args.curve:
        qs = ([1.0] if args.quick
              else [i / (args.curve - 1) for i in range(args.curve)][::-1])
        POINTS = []
        for q in qs:
            flq, hq = leak_hop_at(q)
            POINTS.append((f"q={q:.3f},f_leak={flq:.4f},hop={hq:.4f}",
                           flq, hq, q))
    else:
        POINTS = [(f"f_leak={fl:.4f},hop={hop:g}", fl, hop, None)
                  for fl in ends for hop in hops]
    t0 = time.time()

    if args.offlocus:
        # A DIAGNOSTIC, never a registered threshold: how much of
        # D-sign survives when the flown profile is merely inhomogeneous
        # at equal dose, off the locus and not a transversal. Only the
        # {U, C} arms are needed, d being their difference.
        global OFF_LOCUS
        thr = args.theta_d
        res = {"mode": "offlocus_diagnostic",
               "timestamp": datetime.now().isoformat(),
               "profiles": {k: v.tolist()
                            for k, v in OFF_LOCUS_PROFILES.items()},
               "theta_D_tested": thr, "seed": args.seed,
               "reps": args.reps, "time_axis": rcb.TIME_AXIS,
               "phase_axis": PHASE_AXIS, "eps_scale": EPS_SCALE,
               "cal_shots": CAL_SHOTS,
               "out_of_sector": args.oos or OUT_OF_SECTOR, "variants": {}}
        odd_c = (gate.CORN - gate.CORN[::-1]) / 2
        nrm = float(odd_c @ odd_c)
        print(f"MIRROR-ODD PROJECTION LADDER, {args.reps} reps, tested "
              f"against theta_D = {thr}, axis = {rcb.TIME_AXIS}")
        print("  d as a function of the projection onto C's mirror-odd "
              "part, which G4 registers as the first-order coordinate\n")
        for name, prof in OFF_LOCUS_PROFILES.items():
            OFF_LOCUS = prof
            odd_p = (prof - prof[::-1]) / 2
            proj = float(odd_p @ odd_c / nrm)
            pairs = [float(prof[l] + prof[5 - l]) for l in range(3)]
            on_locus = all(abs(x - 2.0) < 1e-9 for x in pairs)
            res["variants"][name] = {
                "profile": prof.tolist(),
                "profile_variance": float(prof.var()),
                "projection_onto_C_odd": proj,
                "pair_sums": pairs, "on_locus": on_locus,
                "corners": {}}
            print(f"  --- {name}: {np.round(prof, 3)} "
                  f"proj {proj:.3f}, var {prof.var():.3f}, "
                  f"{'ON' if on_locus else 'OFF'} locus ---")
            for fl in ends:
                for hop in hops:
                    tag = f"f_leak={fl:.4f},hop={hop:g}"
                    print(f"  {tag} ...", flush=True)
                    rows = campaign(args.reps, fl, hop, "offlocus",
                                    args.seed + 300000, args.oos,
                                    arms=("U", "C"),
                                    progress=max(1, args.reps // 5))
                    st = _stats([r["d"] for r in rows])
                    rate = _rate([r["d"] for r in rows], thr, above=True)
                    res["variants"][name]["corners"][tag] = {
                        "d": st, "P_exceeds_theta_D": rate,
                        "draws_d": _col(rows, "d")}
                    print(f"    d = {st['mean']:+.5f} +- {st['sd']:.5f}"
                          f"   P(d > theta_D) = {rate:.3f}\n", flush=True)
        res["wall_clock_s"] = time.time() - t0
        p = HERE / "results" / (
            "corner_beat_offlocus_"
            + datetime.now().strftime("%Y%m%d_%H%M%S") + ".json")
        p.parent.mkdir(exist_ok=True)
        with open(p, "w") as f:
            json.dump(res, f, indent=2)
        print(f"  {res['wall_clock_s']:.1f} s -> {p}")
        return
    out = {"mode": "refreeze", "timestamp": datetime.now().isoformat(),
           "config": {"Q": Q, "JDT": JDT, "STEPS": STEPS, "K": K,
                      "SHOTS": SHOTS, "M_BIND": M_BIND, "GPB": GPB,
                      "p2": P2, "reps": args.reps, "seed": args.seed,
                      "time_axis": rcb.TIME_AXIS,
                      "phase_axis": PHASE_AXIS,
                      "R_MAX_FIT": rcb.R_MAX_FIT,
                      "kept_count_floor": rcb.KEPT_COUNT_FLOOR,
                      "eps_scale": EPS_SCALE, "cal_shots": CAL_SHOTS,
                      "out_of_sector": args.oos or OUT_OF_SECTOR,
                      "chain": "analyze-side (pool -> invert -> clip)"},
           "corners": {}}

    print(f"REFREEZE: {args.reps} reps per corner, "
          f"{len(ends)} f_leak x {len(hops)} hop corners, "
          f"axis = {rcb.TIME_AXIS}, phase axis = {PHASE_AXIS}, "
          f"chain = analyze-side\n")

    for tag, fl, hop, qval in POINTS:
        print(f"  {tag}")
        print("    H0-C (C flown as U, analyzed as C) ...", flush=True)
        h0 = campaign(args.reps, fl, hop, "C", args.seed,
                      args.oos, progress=max(1, args.reps // 5))
        print("    H0-both (C and C-prime flown as U) ...", flush=True)
        h0b = campaign(args.reps, fl, hop, "both", args.seed + 2,
                       args.oos, progress=max(1, args.reps // 5))
        print("    H1 ...", flush=True)
        h1 = campaign(args.reps, fl, hop, False, args.seed + 1,
                      args.oos, progress=max(1, args.reps // 5))
        # THE HELD-OUT NULLS (2026-08-18). Section 8a registers the
        # package as ">= 500 H0 reps with HELD-OUT SPLITS", and the
        # frozen 200-rep run honoured it (theta_D from the first 50,
        # the false rate on the second 50). This harness did not:
        # it built theta_D from h0 and then counted exceedances on
        # h0. The fraction of a sample above its own mean + 3 sd is
        # near zero by construction and says nothing about the
        # flight's type-I rate. Thresholds keep the full first
        # sample; every RATE is measured on fresh draws. Seeds are
        # far off the 0/1/2 stride so today's three campaigns stay
        # reproducible bit for bit.
        n_held = args.heldout if args.heldout else args.reps
        print(f"    H0-C held out ({n_held}) ...", flush=True)
        h0h = campaign(n_held, fl, hop, "C", args.seed + 100000,
                       args.oos, progress=max(1, n_held // 5))
        print(f"    H0-both held out ({n_held}) ...", flush=True)
        h0bh = campaign(n_held, fl, hop, "both", args.seed + 200000,
                        args.oos, progress=max(1, n_held // 5))
        # THE HELD-OUT H1 (2026-08-18, the floor decision). Until
        # the floor became a FUNCTION, nothing was estimated on the
        # H1 side: the thresholds all came from H0 chains and H1
        # only ever supplied a rate. The floor's coefficients are
        # chosen against H1 draws (reachability is an H1 quantity),
        # so H1 now carries an estimated object too and section
        # 8a's held-out rule reaches it. Without this campaign the
        # reachability that gates the freeze would be read on the
        # very sample it was fitted to, which is the theta_D
        # mistake of 2026-08-17 one level up. Measured cost of
        # NOT doing it, on the four-corner run: the greedy optimum
        # reads 0.008 in-sample and 0.012 held out.
        print(f"    H1 held out ({n_held}) ...", flush=True)
        h1h = campaign(n_held, fl, hop, False, args.seed + 400000,
                       args.oos, progress=max(1, n_held // 5))
        thD, sD = three_sigma(h0, "d")
        # theta_W freezes against the BOTH-arms null, section 8:
        # under the C-only null C-prime still splits, so d_W would
        # be frozen against a live alternative
        thW, sW = three_sigma(h0b, "d_W")
        fl_floor, sF = three_sigma(
            [{"x": r["C"]["s2"]} for r in h0], "x")
        d1 = _stats([r["d"] for r in h1])
        w1 = _stats([r["d_W"] for r in h1])
        s2C1 = _stats([r["C"]["s2"] for r in h1])
        det = _rate([r["d"] for r in h1], thD, above=True)
        false_rate = _rate([r["d"] for r in h0h], thD, above=True)
        floor_trip = _rate([r["C"]["s2"] for r in h1], fl_floor,
                           above=False)
        floor_false = _rate([r["C"]["s2"] for r in h0h], fl_floor,
                            above=False)
        w_false = _rate([r["d_W"] for r in h0bh], thW, above=True)
        out["corners"][tag] = {
            # the sweep point's own coordinates, so the artifact says
            # WHERE it was measured without a reader parsing the tag,
            # and so a curve run is distinguishable from a box run cell
            # by cell and not only by its header
            "f_leak": fl, "hop": hop, "q": qval,
            "theta_D": thD, "H0_d": sD,
            "theta_W": thW, "H0_dW_both_arms_null": sW,
            "W_false_rate_both_arms_null_own_threshold": w_false,
            "P_W_exceeds_given_H1_own_threshold": _rate(
                [r["d_W"] for r in h1], thW, above=True),
            "s2C_floor": fl_floor, "H0_s2C": sF,
            "H1_d": d1, "H1_dW": w1, "H1_s2C": s2C1,
            # NAMED _own_threshold, because these are the
            # SUPERSEDED copies: section 7 registers the operating
            # characteristics AT the registered threshold, and the
            # governing values live in at_registered below. The
            # first version of this repair computed the right
            # numbers and then left both copies sharing four key
            # names, with the superseded ones at the outer level,
            # which is how the far end's (1.000 / 0.020) pair got
            # into the document in the first place.
            "P_detect_own_threshold": det,
            "H0_false_rate_own_threshold": false_rate,
            "P_s2C_below_floor_given_H0_own_threshold": floor_false,
            "P_s2C_below_floor_given_H1_own_threshold": floor_trip,
            "power": (d1["mean"] / d1["sd"]) if d1["sd"] else None,
            "held_out_reps": n_held,
            # every rate above is measured on HELD-OUT draws; the
            # per-rep values are persisted so a re-evaluation at a
            # different threshold never costs another campaign
            # s2(U) IS PERSISTED IN ITS OWN COLUMN (2026-08-18), not
            # left to be recovered as s2C - d. The arithmetic is
            # exact, but _col maps each column's dead reps to None
            # independently, so two columns are rep-aligned only
            # while nothing died; a floor evaluated per rep must not
            # depend on that holding.
            "draws": {
                "H1_d": _col(h1, "d"), "H1_dW": _col(h1, "d_W"),
                "H1_s2C": _col(h1, "C", "s2"),
                "H1_s2U": _col(h1, "U", "s2"),
                # the FIT-side H0 pairs, persisted because the floor
                # coefficients are chosen against them: without
                # these the selection could not be reproduced from
                # the artifact, only re-run
                "H0_s2C_fit": _col(h0, "C", "s2"),
                "H0_s2U_fit": _col(h0, "U", "s2"),
                "H0_d_heldout": _col(h0h, "d"),
                "H0_s2C_heldout": _col(h0h, "C", "s2"),
                "H0_s2U_heldout": _col(h0h, "U", "s2"),
                "H0both_dW_heldout": _col(h0bh, "d_W"),
                # the held-out H1 sample the floor coefficients are
                # VALIDATED on, never fitted to
                "H1_d_heldout": _col(h1h, "d"),
                "H1_s2C_heldout": _col(h1h, "C", "s2"),
                "H1_s2U_heldout": _col(h1h, "U", "s2"),
            },
            # FIT HEALTH, the same registered UNIT (section 8a): the
            # dw excursion margin must be frozen from the IDEAL
            # construction's own excursion PLUS H0, never at the
            # search grid's edge, and the saturation fraction needs
            # its base. These fall out of the campaigns already run,
            # so they are measured here rather than in a second pass
            # over the same hours.
            "fit_health": {
                "dw_abs_max_H0": _quant(_fit_col(h0, "dws")),
                "dw_abs_max_H1": _quant(_fit_col(h1, "dws")),
                # RESIDUALS SPLIT BY CHANNEL ARITY, because
                # fit_residual_bound is applied PER CHANNEL
                # (run_corner_beat.arm_unhealthy) while the C and
                # C-prime +- channels sum their SSE over TWO traces
                # and every other channel over one. The pooled
                # distribution is bimodal with a ~7-10x separation,
                # so a bound frozen from its upper tail could never
                # fire on the five single-trace channels.
                "sse_H0_1trace": _quant(_fit_col(h0, "sse", 1)),
                "sse_H1_1trace": _quant(_fit_col(h1, "sse", 1)),
                "sse_H0_2trace": _quant(_fit_col(h0, "sse", 2)),
                "sse_H1_2trace": _quant(_fit_col(h1, "sse", 2)),
                # SATURATION, with its BASE named as one of the two
                # registered strings (section 8a: the freeze states
                # `r_max` or `1.1*r_max`, and the runner hard-aborts
                # on any third spelling). The reachable rate ceiling
                # is r_max*(1 + 1/12 + 1/60) = 1.1*r_max exactly,
                # the same refinement arithmetic that puts the |dw|
                # ceiling at 0.10; the first version of this block
                # applied that insight to dw only and measured the
                # rate against 0.9*r_max, which is neither base.
                "rate_dist_H1": _quant(_fit_col(h1, "rates")),
                "rate_dist_H0": _quant(_fit_col(h0, "rates")),
                "r_reachable_ceiling": 1.1 * rcb.R_MAX_FIT,
                "r_rail_fraction_H1": _rate(
                    _fit_col(h1, "rates"),
                    1.1 * rcb.R_MAX_FIT - 1e-9, above=True),
                "r_rail_fraction_H0": _rate(
                    _fit_col(h0, "rates"),
                    1.1 * rcb.R_MAX_FIT - 1e-9, above=True),
                "saturation_frac_base_r_max": _rate(
                    _fit_col(h1, "rates"), rcb.R_MAX_FIT - 1e-9,
                    above=True),
                "saturation_frac_base_1p1_r_max": _rate(
                    _fit_col(h1, "rates"),
                    1.1 * rcb.R_MAX_FIT - 1e-9, above=True),
                # The coarse grid ends at 0.06, but fit_rate then
                # refines twice by +-0.02, so the REACHABLE ceiling
                # is 0.10. Measured here because the H1 fits sit ON
                # it: a margin frozen from an excursion that is
                # railing measures the search, not the physics, and
                # section 8a already asks for a wider dw grid or a
                # two-frequency channel model so gate and analyze
                # share an INTERIOR optimum. The rail fraction is
                # what says whether that work is still owed.
                "dw_coarse_grid_edge": 0.06,
                "dw_reachable_ceiling": 0.10,
                "dw_rail_fraction_H1": _rate(
                    _fit_col(h1, "dws"), 0.0999, above=True),
                "dw_rail_fraction_H0": _rate(
                    _fit_col(h0, "dws"), 0.0999, above=True),
                "r_max": rcb.R_MAX_FIT,
            },
            # THE DENOMINATORS. _rate drops non-finite reps, which
            # is right (a dead rep is a VOID, not a non-detection),
            # but the first version of that repair reported the
            # dropped count only for the campaigns that build the
            # THRESHOLDS, while it had just moved every RATE onto
            # the held-out campaigns. Each population reports its
            # own live count now, named for the sample it is.
            "health": {
                "live_H1": len(_live(h1, "d")),
                "live_H1_heldout": len(_live(h1h, "d")),
                "live_H0_threshold": len(_live(h0, "d")),
                "live_H0_heldout": len(_live(h0h, "d")),
                "live_H0both_threshold": len(_live(h0b, "d_W")),
                "live_H0both_heldout": len(_live(h0bh, "d_W")),
                "nan_reps_H1": _nan_count(h1, "d"),
                "nan_reps_H1_heldout": _nan_count(h1h, "d"),
                "nan_reps_H0_threshold": _nan_count(h0, "d"),
                "nan_reps_H0_heldout": _nan_count(h0h, "d"),
                "nan_reps_H0both_heldout": _nan_count(h0bh, "d_W"),
                "floored_H1": int(sum(r[a]["floored"] for r in h1
                                      for a in ("U", "C", "Cp"))),
                "clipped_H1": float(sum(r[a]["clipped"] for r in h1
                                        for a in ("U", "C", "Cp"))),
            },
        }
        print(f"    theta_D = {thD:.5f}  (E[d|H0] = {sD['mean']:+.5f}, "
              f"sd {sD['sd']:.5f}, n {sD['n']})")
        print(f"    theta_W = {thW:.5f}   s2C floor = {fl_floor:.5f}")
        print(f"    H1: d = {d1['mean']:+.5f} +- {d1['sd']:.5f} "
              f"(power {d1['mean'] / d1['sd']:.2f}), "
              f"P(detect) = {det:.3f}, "
              f"H0 false (held out, n {n_held}) = {false_rate:.3f}")
        print(f"    P(s2C < floor | H1) = {floor_trip:.3f} "
              f"(section 7 asks <= 0.01)\n", flush=True)

    # THE SUMMARY IS GUARDED, THE ARTIFACT IS NOT OPTIONAL.
    # Everything above cost the campaigns; everything below only
    # reads them. A raise down here used to take the whole run
    # with it, so the dump moved into a finally.
    try:
        if out["corners"]:
            out["registered"] = {
                "theta_D": max(c["theta_D"] for c in out["corners"].values()),
                "theta_W": max(c["theta_W"] for c in out["corners"].values()),
                # THE FLOOR IS A LOWER CUT, SO IT FREEZES AT THE LOWER
                # ENVELOPE (2026-08-18). Not a new rule: the house already
                # owns it, decided after measuring what the other choice
                # costs. RECORD_PARITY_HARDWARE_PREDICTION.md, round 28,
                # verbatim: "the fewer-false-void rule thereby extends to
                # the LEVEL axis: a valid device lives anywhere in
                # [eta_min, 1], so the guard's envelope must too; the FLOORS
                # are unaffected, since for a LOWER cut the worst basis IS
                # the lower envelope, verified: the post-fix corner
                # breakdown shows no floor trips". That document reached the
                # rule the expensive way: a cut frozen at the single worst
                # admitted basis "terminally VOIDed the BEST admitted
                # device", VOID in 269 of 300 runs, and it counts the
                # episode as the seventh instance of an enumerated
                # "protection-interaction" class.
                # Our measured 0.046-0.074 reachability breach at the
                # MAXIMUM is the same disease, milder, and the same cure
                # applies. Note what the precedent also says: under the
                # correct freeze the floor does NOT trip. Inertness is the
                # expected state of a correctly frozen floor here, not a
                # defect to engineer away, which is what the function form
                # explored earlier this session was trying to do.
                "s2C_floor": min(
                    c["s2C_floor"] for c in out["corners"].values()),
                # the max is kept beside it because the amendment has to
                # quote what it supersedes, and a reader has to be able to
                # recompute the breach that forced the change
                "s2C_floor_max_superseded": max(
                    c["s2C_floor"] for c in out["corners"].values()),
            }
            print("REGISTERED over the swept set "
                  f"({len(out['corners'])} points, "
                  f"{'counted curve' if args.curve else 'product box'}):")
            for k, v in out["registered"].items():
                print(f"    {k} = {v:.5f}")
            print("    (theta_D and theta_W are UPPER cuts and freeze at the "
                  "upper envelope; s2C_floor is a LOWER cut and freezes at "
                  "the LOWER one, RECORD_PARITY round 28)")

            # ------------------------------------------------------------
            # WHAT THE REGISTERED FLOOR ACTUALLY DOES. It is frozen above,
            # by the lower-envelope rule; nothing is chosen here. This block
            # only MEASURES it, because the one number the pre-registration
            # never carried is the one that says whether the floor changes
            # any verdict at all.
            #
            # WHY THE MEASUREMENT AND NOT A FITTED FUNCTION (2026-08-18).
            # Earlier this session the floor was going to be registered as
            # an affine function of the arm's own s2(U), with coefficients
            # chosen by a grid search under the reachability criterion. The
            # machinery is built and tested (floor_value, binom_upper_bound,
            # select_floor_coefficients, validate_floor_coefficients) and it
            # is deliberately left in place as exploration, but it is NOT on
            # the freeze path, for two reasons found after building it.
            # First, it could not work: passing conjunct 1 and failing
            # conjunct 2 needs a > theta_D for b <= 1 (is_provably_inert),
            # and reachability forbids a anywhere near theta_D, so over the
            # whole registered grid every admissible member is inert, on
            # 2000 H1 and 2000 held-out H0 draws. Second, and this is the
            # part that mattered, the house already had the answer: a lower
            # cut freezes at the lower envelope, and the precedent that
            # settled it also verified that the floor then does not trip.
            # Choosing coefficients after seeing the data would have been
            # correction-shape 58 besides.
            # ------------------------------------------------------------
            tags = list(out["corners"])
            reg_thD = out["registered"]["theta_D"]
            reg_floor = out["registered"]["s2C_floor"]
            coef = (reg_floor, 0.0)          # a scalar IS the b = 0 member
            g_h0 = g_h1 = n_h0 = n_h1 = 0
            reach = {}
            for tag in tags:
                dr = out["corners"][tag]["draws"]
                k, n = incremental_guard_counts(
                    dr["H0_d_heldout"], dr["H0_s2C_heldout"],
                    dr["H0_s2U_heldout"], reg_thD, coef)
                g_h0 += k; n_h0 += n
                k, n = incremental_guard_counts(
                    dr["H1_d_heldout"], dr["H1_s2C_heldout"],
                    dr["H1_s2U_heldout"], reg_thD, coef)
                g_h1 += k; n_h1 += n
                kt, nt = floor_trip_counts(dr["H1_s2C_heldout"],
                                           dr["H1_s2U_heldout"], coef)
                reach[tag] = {"trips": kt, "n": nt,
                              "rate": (kt / nt) if nt else None,
                              "upper_bound": (binom_upper_bound(kt, nt,
                                                                FLOOR_REACH_ALPHA)
                                              if nt else None)}
            inert = is_provably_inert(coef, reg_thD)
            out["floor_report"] = {
                "registered_floor": reg_floor, "theta_D": reg_thD,
                "rule": "lower envelope over the swept set "
                        "(RECORD_PARITY round 28: for a LOWER cut the worst "
                        "basis IS the lower envelope)",
                "provably_inert": inert,
                "false_detections_removed": g_h0, "H0_live": n_h0,
                "true_detections_removed": g_h1, "H1_live": n_h1,
                "reachability_heldout": reach}
            print("")
            print(f"  s2C FLOOR = {reg_floor:.5f} (lower envelope)")
            print(f"     WHAT IT REMOVES, held out: {g_h0} false detections "
                  f"of {n_h0} H0 draws, {g_h1} true detections of {n_h1} "
                  f"H1 draws")
            print("     reachability P(s2C < floor | H1), held out:")
            for tag, cc in reach.items():
                print(f"       {tag}  {cc['trips']}/{cc['n']} = "
                      f"{cc['rate']:.4f}, 95% upper bound "
                      f"{cc['upper_bound']:.4f}")
            if inert:
                print(f"     PROVABLY INERT: floor {reg_floor:.5f} <= theta_D "
                      f"{reg_thD:.5f}, so no rep can pass d > theta_D and "
                      "fail the floor. This is the EXPECTED state of a "
                      "correctly frozen lower cut here (the precedent "
                      "verified 'no floor trips'); the amendment states it "
                      "rather than engineering it away.")

            # SECOND PASS, at the REGISTERED thresholds (2026-08-18).
            # Section 7: "Both ends' operating characteristics must be read
            # AT the registered threshold" -- written because the far end's
            # published (1.000 / 0.020) pair had been read at its own 3-sigma
            # value and was superseded. This harness IS that re-evaluation
            # and it repeated the mistake: every rate above is the corner's
            # own. The per-corner numbers stay as diagnostics; these are the
            # ones the flight is bought on.
            reg = out["registered"]
            print("\nAT THE REGISTERED THRESHOLDS (what the flight actually "
                  "reads; the per-corner rates above are diagnostics):")
            # the registered floor as the (a, b) pair the readers take: a
            # scalar floor is the b = 0 member, so one code path serves both
            coef = (reg["s2C_floor"], 0.0)
            for tag, c in out["corners"].items():
                dr = c["draws"]
                c["at_registered"] = {
                    # the FIRST conjunct alone, kept because every earlier
                    # record quotes it and the comparison with the
                    # conjunction is what shows what the floor costs
                    "P_detect_d_only": _rate(dr["H1_d"], reg["theta_D"], True),
                    "P_detect_d_only_heldout": _rate(
                        dr["H1_d_heldout"], reg["theta_D"], True),
                    "H0_false_rate_d_only": _rate(dr["H0_d_heldout"],
                                                  reg["theta_D"], True),
                    "P_W_exceeds_given_H1": _rate(dr["H1_dW"],
                                                  reg["theta_W"], True),
                    "W_false_rate": _rate(dr["H0both_dW_heldout"],
                                          reg["theta_W"], True),
                    # the reachability of the SUPERSEDED max-over-the-set
                    # floor, kept because it is the breach that forced the
                    # rule change and the amendment has to quote it. Read
                    # on the HELD-OUT H1 draws like everything else here;
                    # the first version of this line read the fit column
                    # while every neighbour read held-out.
                    "P_s2C_below_max_floor_given_H1_heldout": _rate(
                        dr["H1_s2C_heldout"], reg["s2C_floor_max_superseded"],
                        False),
                }
                a = c["at_registered"]
                if coef is not None:
                    # THE OPERATIVE RULE (round 9): the CONJUNCTION. Read on
                    # HELD-OUT draws for both hypotheses, because both the
                    # threshold and the floor coefficients are estimated.
                    a["P_detect_conjunction_heldout"] = conjunction_rate(
                        dr["H1_d_heldout"], dr["H1_s2C_heldout"],
                        dr["H1_s2U_heldout"], reg["theta_D"], coef)
                    a["H0_false_rate_conjunction_heldout"] = conjunction_rate(
                        dr["H0_d_heldout"], dr["H0_s2C_heldout"],
                        dr["H0_s2U_heldout"], reg["theta_D"], coef)
                    a["P_floor_trips_given_H1_heldout"] = (
                        floor_trip_counts(dr["H1_s2C_heldout"],
                                          dr["H1_s2U_heldout"], coef))
                    a["P_floor_trips_given_H0_heldout"] = (
                        floor_trip_counts(dr["H0_s2C_heldout"],
                                          dr["H0_s2U_heldout"], coef))
                print(f"  {tag}")
                print(f"     P(detect), d alone   = "
                      f"{a['P_detect_d_only_heldout']:.3f} held out   "
                      f"H0 false = {a['H0_false_rate_d_only']:.3f}   "
                      f"P(W|H1) = {a['P_W_exceeds_given_H1']:.3f}")
                if coef is not None:
                    kt, nt = a["P_floor_trips_given_H0_heldout"]
                    print(f"     P(detect), CONJUNCTION = "
                          f"{a['P_detect_conjunction_heldout']:.3f} held out   "
                          f"H0 false = "
                          f"{a['H0_false_rate_conjunction_heldout']:.3f}   "
                          f"floor trips on H0 {kt}/{nt}")
    finally:
        out["wall_clock_s"] = time.time() - t0
        _dump(out)


def _dump(out):
    """Persist the run. Called from a finally, so a broken print or a
    None rate cannot destroy two hours of campaigns.

    THE FAILURE THIS PREVENTS HAPPENED, on a four-rep smoke run minutes
    after a reviewer named it: the summary block raised a KeyError on a
    renamed key AFTER every campaign had finished and BEFORE anything
    was written, and the whole run was lost. On the freeze that is ~2 h."""
    path = HERE / "results" / (
        "corner_beat_refreeze_"
        + datetime.now().strftime("%Y%m%d_%H%M%S") + ".json")
    path.parent.mkdir(exist_ok=True)
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\n  -> {path}")


if __name__ == "__main__":
    main()
