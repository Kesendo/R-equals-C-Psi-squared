"""The corner-beat flight runner (G7 of CORNER_BEAT_HARDWARE_PREDICTION.md).

Pre-registration: experiments/CORNER_BEAT_HARDWARE_PREDICTION.md (v7.3,
R=CPsi-squared commits 05a3129 + 6d3b75a + 9e480a5). Committed gate:
simulations/corner_beat_gate.py (v2.1). The physics/estimator functions
below marked VERBATIM are copied from that gate unchanged; the G7 review
diffs them against the committed source.

FROZEN CONFIGURATION (Tom's two budget decisions, 2026-08-16):
Q = 10, J*dt = 0.15, depth grid 21 points x 3 Trotter steps (0..60),
16384 total shots per (arm, depth, prep) split over M = 1024 bindings
(16 shots/binding), FRACTIONAL-RZZ (2 two-qubit gates per XY block;
CZ = NO-FLIGHT), arms per depth in the pinned order N0 -> U -> C' -> C,
preparations = the three room dyads (1,2), (2,3), (3,5).

THE FLOWN CIRCUIT (verified from below before this file existed; the
--certify mode re-verifies all of it at machine precision):
  - One Strang step = XY blocks only (odd half / even full / odd half,
    each block rxx(2 theta) . ryy(2 theta) = exactly 2 rzz after basis
    change), PLUS the ZZ part restricted to the one-magnon sector,
    which is a single-particle potential: CONSTANT for the odd bonds
    (a global phase) and 2 n_0 + 2 n_5 for the even bonds (an rz layer
    on qubits 0 and 5, free). On the one-magnon sector the step equals
    the committed gate's strang(dt) to machine precision; off-sector
    components (from noise) see a different unitary and are removed by
    the pinned post-selection on total excitation = 1.
  - Preparation: 5 Givens rotations from |100000> building the orbital
    (psi_i + psi_j)/sqrt(2) of the GAUGE-PINNED Floquet modes.
  - Injection: per step, per site, one parameterized rz;
    phi ~ N(0, 4 gamma_l dt); tables frozen per (arm, depth, prep) from
    named SeedSequence keys (persisted as seeds + sha256, not values).

MODES
  --certify           local, no network: sector parity, prep parity,
                      trajectory parity vs the committed-gate simulator,
                      transpile-shape assertions on a local rzz basis,
                      dose certificates on the frozen tables.
  --calibrate         line selection on the pinned backend by the
                      concentrator uniform-line rule + the BOTH-TWINS
                      fractional-rzz check; writes corner_beat_chain_*.json.
  --aer [--full]      hardware-shaped SYNTHETIC artifact through a
                      hand-built noise model (never from_backend), same
                      JSON schema as --hardware; feeds the --analyze
                      end-to-end proof (G7). Default reduced M = 64.
  --hardware [--yes]  the flight: Class-1 guards (hard aborts, no
                      override), day-of snapshot, billing projection cap,
                      typed FLY confirmation, stub-before-wait persistence,
                      per-binding raw counts, job.usage() persisted.
  --analyze FILE      committed estimator on a persisted artifact:
                      per-binding raw counts pooled (over the resampled
                      bindings), ONE confusion inversion of the pooled
                      vector (exactly equal to per-binding inversion by
                      linearity), ONE clip, post-selection, eigenchannel
                      damped-cosine fits, s2 statistics, binding
                      bootstrap, verdict-rule evaluation against
                      corner_beat_constants.json. The pool-then-invert
                      order and the clip policy amend the section 5/6.1
                      letter (a numbered pre-data Amendment carries
                      them: per-binding invert-then-clip is nonlinear at
                      16 shots/binding and depth-biased, measured).
                      THE COMMITTED DOCUMENT GOVERNS, NOT THIS PRINTOUT.

PENDING (remaining committed-gate work, tracked in the pre-registration
section 8a "Still outstanding", SHARPENED by the G7 review round): the
N0 full 3x3 generator fit B-hat and var(B-hat) (G2), the dressed centre
functions f_C / f_Cp as code, the difference-null generator fits and
margins, family rates at >= 500 H0 reps, the per-M retention criterion,
asymmetric confusion, the hop-split bracket, AND (G7 round findings):
a theta_D REFREEZE through the analyze-side chain (fixed R_MAX_FIT,
confusion inversion in the loop, >= 500 H0 reps; the frozen 0.00253
came from 50 reps on the gate's counts path), G3(b)/G5 re-runs at the
frozen 21 x 3 grid (they ran at the retired 13 x 5), the A-arm and Pi3
centres as DRESSED functions (the noiseless flown construction misses
the ideal r-bar by -1.4%/-3.4% and Pi3(C) is not 0 under the committed
estimator: dressing, not error; both centres are corridor-dressing per
section 6, to be frozen as functions like f_C), the committed G1-stage
record (the s2(Cp)-s2(U) line exists only in uncommitted output), and
the per-binding billing overhead measured at M = 1024 (the concentrator
law was measured at M = 256). --analyze computes and prints every raw
statistic now and WITHHOLDS any verdict whose constant is not yet
frozen in corner_beat_constants.json; --hardware --yes REFUSES to
submit while any constant is unfrozen or the commit gate fails.
R_MAX_FIT = 0.9 is pinned here as the analyze-side fit bound and joins
the theta_D refreeze above.
"""

import argparse
import hashlib
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
RESULTS_DIR = SCRIPT_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

try:
    from dotenv import load_dotenv
    load_dotenv(SCRIPT_DIR / ".env")
except ImportError:
    pass
IBM_TOKEN = os.environ.get("IBM_QUANTUM_TOKEN")

# ----------------------------------------------------------------- pins
N = 6
JDT = 0.15                      # J*dt per Trotter step (J = 1 units)
Q = 10
GBAR = 1.0 / Q                  # engineered gamma-bar in J units
GBAR_DT = JDT / Q               # gamma-bar * dt
K_GRID = 3                      # Trotter steps per grid point
STEPS_MAX = 60                  # deep end
GRID_STEPS = list(range(0, STEPS_MAX + 1, K_GRID))   # 21 points
M_BIND = 1024
SHOTS_SCIENCE = 16384
SHOTS_BIND = SHOTS_SCIENCE // M_BIND                 # 16
R_BOOT = 500
R_MAX_FIT = 0.9                 # pinned analyze-side fit bound (see docstring)

DY = [(1, 2), (2, 3), (3, 5)]
ODD_BONDS = [(0, 1), (2, 3), (4, 5)]
EVEN_BONDS = [(1, 2), (3, 4)]
MIRROR_PAIRS = [(0, 5), (1, 4), (2, 3)]

T1_CLEAN_SCALE_BOUND = 0.25     # G4 v2, in gbar; section 9's scale face


def t1_clean_profile(tau_step_us, t1_us):
    """The T1-CLEAN us->gbar bridge of section 9, as code.

    gamma_l^hw = tau_step / (4 * T1_l * dt): the Gamma_l/4-equivalent
    dephasing in the simulation's J units. Gamma/4 is the LINDBLAD-book
    T1 contribution to Z-dephasing (a T1 process at rate Gamma costs
    Gamma/2 of coherence rate, and this repository's gamma is half a
    coherence rate); writing Gamma/2 would be the coherence book and
    would double every station.

    Returns the per-site profile in J units and in gbar, plus the SCALE
    the G4 bound is stated against: the profile MEAN in gbar units, the
    definition the committed gate draws from ("'scale' = Gam_mean/4 in
    units gbar", simulations/corner_beat_gate.py, uniform(0, 2*scale*gbar)).
    """
    if len(t1_us) != N:
        raise ValueError(f"expected {N} stations, got {len(t1_us)}")
    if any(t <= 0 for t in t1_us):
        raise ValueError(f"T1 must be positive on every station: {t1_us}")
    if tau_step_us <= 0:
        raise ValueError(f"tau_step_us must be positive: {tau_step_us}")
    g_j = [tau_step_us / (4.0 * t * JDT) for t in t1_us]
    g_gbar = [g / GBAR for g in g_j]
    return {"gamma_hw_j": g_j,
            "gamma_hw_gbar": g_gbar,
            "scale_gbar": sum(g_gbar) / len(g_gbar)}


def t1_list_from_snapshot(snapshot, chain):
    """Per-station T1 in LINE order, out of the same structure the day-of
    re-gate writes: snapshot["qubits"][str(q)]["T1_us"]. The order is the
    chain's, never the dict's; a missing station raises rather than
    defaulting, since a silent nominal T1 is how a bad line passes."""
    qubits = snapshot["qubits"]
    return [float(qubits[str(q)]["T1_us"]) for q in chain]


def t1_clean_scale_face(tau_step_us, t1_us):
    """The scale face of T1-CLEAN: the profile scale must not exceed the
    registered bound. Returns (passes, profile). The comparison is <=,
    per section 9's own wording; a boundary that excludes its own edge is
    an unregistered margin."""
    prof = t1_clean_profile(tau_step_us, t1_us)
    return prof["scale_gbar"] <= T1_CLEAN_SCALE_BOUND, prof


UNIF = np.ones(N)
CORN = np.zeros(N); CORN[[0, 3, 4]] = 2.0    # C  (maximizing)
CORNP = np.zeros(N); CORNP[[0, 1, 2]] = 2.0  # C' (non-maximizing)
ARM_NAMES = ["N0", "U", "Cp", "C"]           # pinned within-depth order
ARM_PROFILE = {"N0": None, "U": UNIF, "Cp": CORNP, "C": CORN}

BACKEND_NAME = "ibm_kingston"
TWIN_NAMES = ["ibm_kingston", "ibm_marrakesh"]
RULE_MIN_T2_US = 150.0
RULE_MAX_T2_RATIO = 2.0
RULE_MAX_READOUT = 0.02
# Amendment 1.8 (Tom, 2026-08-17): 0.006 -> 0.005. The >= 3x deep-end
# kept-count margin failed at 0.6% (1.74x) and holds at 0.5% (4.17x),
# and the bound now sits ON the measured p2 ladder's first rung instead
# of between rungs. Registered as a manifest key and sync-checked, the
# ruling round 12 made for kept_count_floor: a threshold that
# hard-aborts a paid submission is not a code constant. It stays an
# ISOLATED-gate number; the layered error is the day-of layer-fidelity
# gate's business (doc section 9), and >= 3x is equivalent to an
# EFFECTIVE p2 <= 0.538%.
RULE_MAX_P2_RZZ = 0.005          # Class-1 guard on the used edges
# Round 12: a REGISTERED threshold, not a code constant. A cell whose
# post-selected kept counts fall under it is FLOORED, which makes the
# grid incomplete, which VOIDS D-sign and W. It is sync-checked
# against the constants file like the billing cap and the dose
# criterion; the >= 3x deep-end margin over it is a machinery item.
KEPT_COUNT_FLOOR = 50
BILL_ANCHOR_MS_PER_SHOT = 0.327  # delay-bearing anchor (staircase flights)
BILL_ABORT_MIN = 25.0            # accepted band 21-25 QPU min (v7.2)
DOSE_CRITERION = 0.02            # provisional (per-M criterion pending, G3)

# named seeds; NEVER hash() (process-salted)
SEED_TABLES = 20260816
SEED_TRANSPILER = 42
SEED_AER = 7
SEED_BOOT = 1234

CONSTANTS_PATH = SCRIPT_DIR / "corner_beat_constants.json"
# In-repo copy: the repository root is two levels up from
# simulations/flight/, so this resolves for any clone. The live pipeline
# copy hardcodes an absolute path instead, because it sits outside the
# repository; see README.md in this directory.
RCPSI_REPO = Path(__file__).resolve().parents[2]
GATE_PATH = RCPSI_REPO / "simulations/corner_beat_gate.py"
RCPSI_COMMIT_PATHS = [
    "experiments/CORNER_BEAT_HARDWARE_PREDICTION.md",
    "simulations/corner_beat_gate.py",
    "simulations/results/corner_beat/gate_frozen.txt",
    "simulations/results/corner_beat/gate_g2345.txt",
]


# --------------------------------------------- the line error budget
# The Class-1 guard is a MAXIMUM over the used edges. The quantity it
# protects, the deep-end post-selected kept count, is a WEIGHTED SUM:
# a Strang step gives each ODD bond two XY blocks and each EVEN bond
# one, at 2 rzz per block, and the preparation puts 2 more rzz on every
# bond (five Givens on the overlapping pairs). So an odd edge costs
# exactly twice an even one, and the guard cannot see the difference.
# These weights are DERIVED, not chosen: they sum to the 970 two-qubit
# gates the pre-registration flies and asserts post-transpile.
CHAIN_BONDS = [(l, l + 1) for l in range(N - 1)]
DEEP_END_MARGIN = 3.0    # section 8a, over the kept-count floor
F_LEAK_WORST = 0.9       # the registered bracket's bad end


def edge_gate_weights(n_steps=None):
    """Two-qubit gates on each chain bond at the deepest grid point."""
    n_steps = GRID_STEPS[-1] if n_steps is None else n_steps
    w = {}
    for b in CHAIN_BONDS:
        blocks = 2 if list(b) in [list(x) for x in ODD_BONDS] else 1
        w[b] = blocks * 2 * n_steps + 2      # + 2 rzz of preparation
    return w


def line_error_budget(p2_by_bond, n_steps=None):
    """The deep-end error exponent of a line, per unit f_leak."""
    w = edge_gate_weights(n_steps)
    return float(sum(w[b] * p2_by_bond[b] for b in CHAIN_BONDS))


def projected_deep_end_kept(p2_by_bond, shots=None, n_steps=None):
    """Post-selected kept counts in the deepest cell at the worst
    f_leak end. This is what section 8a's >= 3x margin is ABOUT, so a
    line's margin can be projected from its own edges instead of
    assumed at the guard's worst case."""
    shots = SHOTS_SCIENCE if shots is None else shots
    return float(shots * np.exp(-line_error_budget(p2_by_bond, n_steps)
                                * F_LEAK_WORST))


# section 8a's >= 3x requirement, written as the budget it is. NOT a
# guard, and deliberately so: any line passing the per-edge ceiling has
# budget <= 970 * RULE_MAX_P2_RZZ = 4.85, which is under this bound, so
# a budget GUARD could never fire and would be machinery that only
# looks like a gate. It is used as the line-selection SCORE and as the
# projected margin the calibrate record carries. A test pins the
# redundancy, and fails if an amendment ever raises the ceiling past it.
P2_BUDGET_MAX = float(-np.log(
    DEEP_END_MARGIN * KEPT_COUNT_FLOOR / SHOTS_SCIENCE) / F_LEAK_WORST)


def _day_of_rule_text():
    """The day-of line rule, stated ONCE and built from the constants
    that enforce it. It used to be spelled a second time as a literal
    inside the hard-abort message, so tightening the p2 guard would
    have printed a refusal naming the rule the code no longer applies
    (Amendment 1.8)."""
    return (f"T2 >= {RULE_MIN_T2_US:g} us, "
            f"ratio <= {RULE_MAX_T2_RATIO:g}, "
            f"readout <= {RULE_MAX_READOUT * 100:g}%, "
            f"rzz p2 <= {RULE_MAX_P2_RZZ * 100:g}%")


def _require(cond, msg):
    """Class-1 structural guard that SURVIVES python -O (a bare assert
    is an override flag the pre-registration says does not exist)."""
    if not cond:
        raise RuntimeError(f"GUARD FAILED: {msg}")
DEFAULT_CONSTANTS = {
    "_source": "CORNER_BEAT_HARDWARE_PREDICTION.md v7.3 + "
               "simulations/results/corner_beat/gate_frozen.txt; entries "
               "with frozen=false await the remaining committed-gate work",
    "theta_D": {"value": 0.00253, "frozen": True,
                "refreeze_required": True,
                "note": "held-out-half H0, worst f_leak end, frozen grid; "
                        "exact gate_frozen.txt value (no rounding). "
                        "refreeze_required BLOCKS submission: the "
                        "freezing run used 50 H0 replicates, the gate's "
                        "counts path (no confusion inversion), a "
                        "config-dependent r_max, and non-nested tables; "
                        "the re-freeze runs the analyze-side chain, "
                        "fixed R_MAX_FIT, nested tables, delta-omega "
                        "bound pinned, >= 500 reps"},
    "dw_excursion_margin": {
        "value": None, "frozen": False,
        "note": "fit-health VOID threshold on |delta-omega| (G2). NOT "
                "the coarse-grid edge 0.06: the refinement reaches "
                "+-0.10, the +- channels are NOT eigenchannels of the "
                "Trotter-detuned block (two damped cosines per trace; "
                "|dw| up to 0.03 noiseless, measured), so the margin "
                "must be frozen from the IDEAL construction's own "
                "excursion + H0, inside the theta_D refreeze; until "
                "then the dw flags are informational only"},
    "kept_count_floor": {
        "value": 50, "frozen": False, "refreeze_required": True,
        "note": "the post-selected kept-count FLOOR; a cell under it "
                "is floored, which makes the grid incomplete, which "
                "VOIDS D-sign and W. Registered threshold, not a code "
                "constant (round 12); freezes with the >= 3x deep-end "
                "margin, which reads 4.17x at the p2 guard bound since "
                "Amendment 1.8 tightened it to 0.5% (it read 1.74x at "
                "0.6%, failing the requirement). The margin holds on "
                "the ISOLATED axis only; layered 1.3-2x takes it to "
                "1.12x and then to a floored deep end"},
    "p2_guard_bound": {
        "value": RULE_MAX_P2_RZZ, "frozen": True,
        "note": "Amendment 1.8 (Tom, 2026-08-17): the Class-1 guard on "
                "the calibrated 2q error of the used edges, 0.6% -> "
                "0.5%, taken as the third of the three levers priced "
                "for the failing >= 3x deep-end margin (more shots is "
                "~40 QPU min against a 25-min cap; a lower floor needs "
                "a G3 justification that does not exist). Sync-checked "
                "against RULE_MAX_P2_RZZ. ISOLATED-gate number: the "
                "layered error is the day-of layer-fidelity gate's, "
                "whose own bound is still owed and whose candidate "
                "value, an effective p2 <= 0.538%, is what >= 3x is "
                "equivalent to"},
    "r_saturation_frac": {
        "value": None, "base": "r_max", "frozen": False,
        "note": "fit-health VOID threshold on rate ceiling proximity "
                "(G2); the BASE is a FIELD of this entry, not prose "
                "(round 11): the search's true ceiling after "
                "refinement is 1.1*r_max, not r_max, and the runner "
                "multiplies the fraction by whichever of 'r_max' / "
                "'1.1*r_max' this field names. G2 sets both fields "
                "together; informational until frozen"},
    "dose_criterion": {"value": 0.02, "frozen": False,
                       "note": "transplanted from the concentrator at a "
                               "different dose; the per-M criterion is "
                               "outstanding gate work (G3); used "
                               "provisionally by certify/dry runs"},
    "billing_cap_min": {"value": 25.0, "frozen": True,
                        "note": "the accepted 21-25 band's top (Tom, "
                                "2026-08-16, second budget decision)"},
    "R_MAX_FIT": {"value": 0.9, "frozen": False,
                  "note": "analyze-side fit bound AND rate lattice "
                          "(terminal grid r_max/240); the gate froze "
                          "theta_D on a config-dependent r_max, so "
                          "this freezes WITH the theta_D refreeze"},
    "R_BOOT": {"value": 2000, "frozen": False,
               "note": "bootstrap replicates; the percentile endpoints "
                       "feed the theta_F containment tests, so 500 is "
                       "too jittery for the FALSIFIED rule; freezes "
                       "with theta_F"},
    "theta_W": {"value": None, "frozen": False},
    "theta_F": {"value": None, "frozen": False},
    "dmag_band_C": {"value": None, "frozen": False},
    "dmag_band_Cp": {"value": None, "frozen": False},
    # round 18: one key per A LINE (the partition reads the C line
    # alone; the (2,3) null voids only the C-prime line)
    "A_margin_C": {"value": None, "frozen": False},
    "A_margin_Cp": {"value": None, "frozen": False},
    "kappa": {"value": None, "frozen": False},
    "null_margin_12": {"value": None, "frozen": False},
    "null_margin_23": {"value": None, "frozen": False},
    "fit_residual_bound": {"value": None, "frozen": False},
    "time_axis": {"value": None, "frozen": False},
}

# ------------------------------------------------ physics (VERBATIM gate)
D6 = 1 << N


def bond_H_1m(bonds):
    H = np.zeros((D6, D6))
    for a in range(D6):
        for l in bonds:
            za = 1 - 2 * ((a >> l) & 1)
            zb = 1 - 2 * ((a >> (l + 1)) & 1)
            H[a, a] += za * zb
            if ((a >> l) & 1) != ((a >> (l + 1)) & 1):
                H[a, a ^ (1 << l) ^ (1 << (l + 1))] += 2
    cfg = [1 << l for l in range(N)]
    return np.array([[H[a, b] for b in cfg] for a in cfg])


Ho = bond_H_1m([0, 2, 4])
He = bond_H_1m([1, 3])
E0, V0 = np.linalg.eigh(Ho + He)


def strang(dt):
    Eo, Uo = np.linalg.eigh(Ho)
    Ee, Ue = np.linalg.eigh(He)
    Ah = Uo @ np.diag(np.exp(-1j * Eo * dt / 2)) @ Uo.T
    B = Ue @ np.diag(np.exp(-1j * Ee * dt)) @ Ue.T
    return Ah @ B @ Ah


def floquet_modes(U, dt):
    lam, W = np.linalg.eig(U)
    eps = -np.angle(lam) / dt
    order = []
    for kk in range(N):
        ov = np.abs(V0[:, kk] @ W.conj())
        for cand in np.argsort(-ov):
            if cand not in order:
                order.append(cand); break
    Wm = W[:, order]
    for kk in range(N):
        ph = Wm[np.argmax(np.abs(Wm[:, kk])), kk]
        Wm[:, kk] = Wm[:, kk] * np.exp(-1j * np.angle(ph))
    Wr = np.real(Wm)
    # GAUGE PIN (v2.1): mode sign fixed by overlap with the continuous-H
    # mode, NOT by the largest component.
    for kk in range(N):
        if V0[:, kk] @ Wr[:, kk] < 0:
            Wr[:, kk] = -Wr[:, kk]
    return eps[order], Wr


# ------------------------------------------------------------ estimator
# fit_rate's SEARCH is the committed gate's, unchanged; the return value
# is extended to (r, dw) because the pre-registration consumes the fitted
# detuning (the G5 in-analysis J estimator, section 6.4) and the fit-health
# VOID trigger needs a saturation readout. --certify asserts equality of
# the rate against the committed gate's fit_rate on shared traces.
def fit_rate(traces, ts, om, r_max, ts_phase=None, env_known=None):
    # ts_phase (2026-08-18): A MEASURED AND REJECTED ALTERNATIVE, kept
    # because the flight's estimator uses ONE axis for two counters and
    # a reader is entitled to know that was tested rather than assumed.
    # The observation is real: the envelope counts the dephasing DOSE,
    # n-1 injection layers at grid depth n (section 5), while the beat
    # phase counts n Trotter steps. With one axis for both, the
    # section-5 axis switch shifts the cosine too; every depth n >= 1
    # shifts by exactly one step and a uniform shift is absorbed by the
    # free amplitude and phase, so the ENTIRE effect sits on depth 0,
    # where the switch repairs a ~5% dose error and introduces a
    # 1 - cos(om*dt) = 4.5% phase error in its place. That is why
    # section 8a measures the registered axis correction as COSTING
    # detection.
    # Passing ts_phase = the NOMINAL axis fits
    # exp(-r*t_dose)*cos((om+dw)*t_nominal + phi), the circuit's own two
    # counters. Measured on the gate's deterministic NOISELESS
    # construction, where the rate triples are known in closed form: the
    # two-axis fit has the LARGER weighted residual on all three arms
    # (C 5.12e-2 vs 4.68e-2, Cp 1.68e-1 vs 1.55e-1, U 7.97e-3 vs
    # 6.43e-3) and inflates s2(C) further from the bare ideal (1.400x
    # vs 1.343x). It buys 15% of d and 8% of power in a 200-rep paired
    # H0/H1 run, which is exactly the section-6.3 inflation trap, not
    # accuracy. The mask and the step unitary do not commute, so the
    # true solution is a non-normal 3x3 generator rather than a pure
    # envelope times a pure cosine; both counters are individually
    # right and their naive separation is still the worse approximation.
    # ts_phase = None is the committed single-axis behaviour, bit for
    # bit, so --certify's parity against the committed gate is
    # untouched (verified exactly, both defaults).
    tp = ts if ts_phase is None else ts_phase
    # env_known (2026-08-18): a KNOWN, arm-independent multiplicative
    # envelope folded into the model beside the fitted exponential.
    # It exists for the RE-ENTRANT population: shots that left the
    # one-magnon sector and walked back, which post-selection keeps
    # because it reads the final magnon number and not the history. At
    # the deep end they are ~78% of the kept shots (4.85 expected
    # errors over the grid; exact 64-dimensional check agrees with the
    # magnon chain to 0.5%).
    # Such a shot is provably beat-less: a uniform release is rho =
    # I/6, which commutes with U and is untouched by the dephasing
    # mask, so it is a permanent fixed point; and a uniform site vector
    # projects to EXACTLY ZERO on the three registered dyads (1e-15, by
    # orthonormality of the Floquet columns). Its only channel into the
    # estimator is therefore the NORMALISATION, which multiplies every
    # dyad trace by (1 - f(t)).
    # Leaving it out is not neutral, and the reason is the shape:
    # (1 - f(t)) is strongly NON-exponential, its effective rate
    # climbing 0.010 -> 0.181 across the grid, so the free amplitude
    # cannot absorb it. A convex common envelope shifts SLOW channels
    # more than fast ones, the rate SPREAD compresses, and s2 is a
    # variance of rates. Measured on an exact NOISELESS construction:
    # d falls to 0.776 of baseline with it left out, and recovers to
    # 99.5% with it folded in.
    #
    # AND THAT RECOVERY DOES NOT SURVIVE THE FULL CHAIN, which is why
    # this stays OFF. Measured 2026-08-18, 400 reps, paired, at the
    # theta_D-setting corner, clean / contaminated / corrected:
    #   d      0.00512 / 0.00434 / 0.00502     the SIGNAL comes back
    #   rbar   0.5010  / 0.5587  / 0.4999      the bias goes entirely
    #   theta_D 0.00287 / 0.00287 / 0.00327    and the THRESHOLD rises
    #   power  2.93    / 2.55    / 2.53
    #   P(detect) at each variant's own threshold
    #          0.907   / 0.807   / 0.810
    # The correction acts under H0 too, where the signal is ~zero and
    # dividing by an envelope that falls to 0.202 at the deep end
    # amplifies the null's own noise; sd(d|H0) grows, theta_D walks up
    # 14%, and the recovered signal is eaten. It buys 0.003.
    # THE LESSON, worth more than the switch: a recovery of the MEAN,
    # measured without noise, says nothing about a decision that is a
    # signal-to-noise ratio.
    # Untried if this is ever revisited: down-weighting by the envelope
    # rather than only dividing by it, and re-opening the depth-grid
    # endpoint, whose NEGATIVE result was measured with this OFF.
    # NO free parameter is added: f(t) comes from the calibrated p2
    # through the chain. It does inherit p2's calibration error, so it
    # would register as a BRACKET over p2, never as a point.
    # None is the committed behaviour, bit for bit.

    def sse(r, dw):
        tot = 0.0
        for (y, w) in traces:
            ok = ~np.isnan(y) & (w > 0)
            if ok.sum() < 5:
                continue
            t = ts[ok]; yy = y[ok]; ww = w[ok] / max(w[ok].max(), 1.0)
            tph = tp[ok]
            sw = np.sqrt(ww)
            env = np.exp(-r * t)
            if env_known is not None:
                env = env * np.asarray(env_known, dtype=float)[ok]
            X = np.column_stack([np.ones_like(t),
                                 env * np.cos((om + dw) * tph),
                                 env * np.sin((om + dw) * tph)])
            try:
                coef, *_ = np.linalg.lstsq(X * sw[:, None], yy * sw,
                                           rcond=None)
            except np.linalg.LinAlgError:
                return np.inf
            tot += float(((yy - X @ coef) ** 2 * ww).sum())
        return tot
    if all((~np.isnan(y) & (w > 0)).sum() < 5 for (y, w) in traces):
        return np.nan, np.nan, np.nan  # dead cell: NaN, never a silent 0
    rs = np.linspace(0.0, r_max, 25)
    dws = np.linspace(-0.06, 0.06, 7)
    _, r0, dw0 = min(((sse(r, dw), r, dw) for r in rs for dw in dws))
    for span in (r_max / 12, r_max / 60):
        rs = np.linspace(max(0, r0 - span), r0 + span, 9)
        dws = np.linspace(dw0 - 0.02, dw0 + 0.02, 5)
        s0, r0, dw0 = min(((sse(r, dw), r, dw) for r in rs for dw in dws))
    return r0, dw0, s0


def s2_of(rates):
    r = np.asarray(rates, dtype=float)
    return float(((r[0] - r[1]) ** 2 + (r[0] - r[2]) ** 2
                  + (r[1] - r[2]) ** 2) / 6)


# the N0 arm is read through the U-arm channel shape (three bare dyads);
# the full 3x3 B-hat generator fit is G2-pending. The C-arm +- channels
# demodulate at om = mean of ALL THREE dyad frequencies, exactly as the
# committed gate does (a deliberate carry-over: the 0.0075 offset to the
# mixing pair's own mean sits inside the +-0.06 detuning search).
def arm_rates(Yp, Wp, A, dyf, ts, arm, r_max, ts_phase=None,
              env_known=None):
    """The three fitted channel (rate, detuning) pairs of an arm; the
    channel structure is section 6.3's, identical to the committed
    gate's arm_s2 branches. Returns (rates[3], dws[3]).

    ts_phase is fit_rate's second axis and defaults to None, which is
    the single-axis behaviour."""
    yd = [6 * Yp[p] @ A for p in range(3)]
    om = float(np.mean(dyf))
    tp = ts_phase
    if arm in ("U", "N0"):
        out = [fit_rate([(yd[p][:, p], Wp[p])], ts, dyf[p], r_max, tp, env_known)
               for p in range(3)]
    elif arm == "C":
        out = [fit_rate([(yd[0][:, 0], Wp[0])], ts, dyf[0], r_max, tp, env_known),
               fit_rate([(yd[1][:, 1] + yd[1][:, 2], Wp[1]),
                         (yd[2][:, 1] + yd[2][:, 2], Wp[2])],
                        ts, om, r_max, tp, env_known),
               fit_rate([(yd[1][:, 1] - yd[1][:, 2], Wp[1]),
                         (yd[2][:, 1] - yd[2][:, 2], Wp[2])],
                        ts, om, r_max, tp, env_known)]
    elif arm == "Cp":
        omp = float(np.mean([dyf[0], dyf[2]]))
        out = [fit_rate([(yd[1][:, 1], Wp[1])], ts, dyf[1], r_max, tp, env_known),
               fit_rate([(yd[0][:, 0] + yd[0][:, 2], Wp[0]),
                         (yd[2][:, 0] + yd[2][:, 2], Wp[2])],
                        ts, omp, r_max, tp, env_known),
               fit_rate([(yd[0][:, 0] - yd[0][:, 2], Wp[0]),
                         (yd[2][:, 0] - yd[2][:, 2], Wp[2])],
                        ts, omp, r_max, tp, env_known)]
    else:
        raise ValueError(f"unknown arm {arm!r}")
    return ([r for r, _, _ in out], [d for _, d, _ in out],
            [s for _, _, s in out])


# ------------------------------------------------------ circuit builders
from qiskit import QuantumCircuit, transpile
from qiskit.circuit import Parameter, ParameterExpression

LOCAL_BASIS = ["rzz", "rz", "sx", "x", "measure", "delay"]


def _rzz_pos(qc, a, b, phi):
    """rzz with the angle folded into (0, pi/2] for fractional hardware:
    rzz(-phi) = X_a rzz(phi) X_a. All pinned angles are <= pi/2 by
    construction (max 2*JDT = 0.3 for the steps, |Givens| <= pi/2)."""
    _require(abs(phi) <= np.pi / 2 + 1e-12, f"rzz angle out of range: {phi}")
    _require(abs(phi) > 1e-15,
             "zero rzz angle would silently break the pinned 2q budget")
    if phi < 0:
        qc.x(a); qc.rzz(-phi, a, b); qc.x(a)
    else:
        qc.rzz(phi, a, b)


def _rxx(qc, a, b, phi):
    qc.h(a); qc.h(b)
    _rzz_pos(qc, a, b, phi)
    qc.h(a); qc.h(b)


def _ryy(qc, a, b, phi):
    qc.s(a); qc.s(b); qc.h(a); qc.h(b)
    _rzz_pos(qc, a, b, phi)
    qc.h(a); qc.h(b); qc.sdg(a); qc.sdg(b)


def xy_block(qc, a, b, theta):
    """exp(-i theta (XX + YY)) on bond (a,b): exactly 2 rzz."""
    _rxx(qc, a, b, 2 * theta)
    _ryy(qc, a, b, 2 * theta)


def strang_step_circuit(qc, dt):
    """One flown Strang step (see module docstring): XY blocks + the
    even-bond ZZ sector potential exp(-i dt (2 n_0 + 2 n_5)) as
    rz(-2 dt) on qubits 0 and 5 (the odd-bond ZZ part is constant on
    the sector). 8 XY blocks x 2 rzz = 16 two-qubit gates."""
    for (a, b) in ODD_BONDS:
        xy_block(qc, a, b, dt / 2)
    for (a, b) in EVEN_BONDS:
        xy_block(qc, a, b, dt)
    qc.rz(-2 * dt, 0)
    qc.rz(-2 * dt, 5)
    for (a, b) in ODD_BONDS:
        xy_block(qc, a, b, dt / 2)


def givens_angles(c):
    """Descending Givens chain preparing sum_l c_l |l> from |100000>,
    all angles within [-pi/2, pi/2] (the remainder sign is chosen by
    lookahead so cos never has to go negative)."""
    c = np.asarray(c, float)
    angs, rem = [], 1.0
    for l in range(N - 1):
        cl = np.clip(c[l] / rem, -1.0, 1.0) if abs(rem) > 1e-14 else 0.0
        mag = float(np.sqrt(max(0.0, 1.0 - cl * cl)))
        # lookahead: give the remainder the sign of the next target
        nxt = c[l + 1]
        s = 1.0 if (nxt >= 0) == (rem >= 0) else -1.0
        th = float(np.arctan2(s * mag, cl))
        angs.append(th)
        rem = rem * float(np.sin(th))
        if abs(rem) < 1e-14:
            rem = 1e-14 * (1 if rem >= 0 else -1)
    return angs


def givens_gate(qc, a, b, th):
    """|a> -> cos th |a> + sin th |b> on the one-magnon pair (measured
    convention: the rz-sandwiched rxx.ryy pair at +th gives cos - sin,
    so the rotation angles are negated)."""
    qc.rz(-np.pi / 2, a)
    _rxx(qc, a, b, -th)
    _ryy(qc, a, b, -th)
    qc.rz(np.pi / 2, a)


def prep_circuit_ops(qc, coeffs):
    # the orbital is defined up to a global sign; pick the sign with
    # c_0 >= 0 so every Givens angle stays within [-pi/2, pi/2]
    c = np.asarray(coeffs, float)
    nz = np.flatnonzero(np.abs(c) > 1e-12)
    if c[nz[0]] < 0:
        c = -c
    qc.x(0)
    for l, th in enumerate(givens_angles(c)):
        givens_gate(qc, l, l + 1, th)


def science_skeleton(n_steps, prep_coeffs):
    """The shared skeleton for one (depth, prep) cell: prep + n_steps
    Strang steps, each followed by one parameterized injection rz per
    site (the last layer is a Z-basis no-op, flown for skeleton
    uniformity), then measure_all. All four arms bind THIS circuit."""
    qc = QuantumCircuit(N)
    prep_circuit_ops(qc, prep_coeffs)
    qc.barrier()
    params = {site: [Parameter(f"phi_q{site}_s{s:02d}")
                     for s in range(n_steps)] for site in range(N)}
    for s in range(n_steps):
        strang_step_circuit(qc, JDT)
        for site in range(N):
            qc.rz(params[site][s], site)
        qc.barrier()
    qc.measure_all()
    return qc, params


def build_preps():
    em, Wm = floquet_modes(strang(JDT), JDT)
    dyf = [em[i] - em[j] for (i, j) in DY]
    A = np.array([[Wm[l, i] * Wm[l, j] for (i, j) in DY] for l in range(N)])
    preps = [(Wm[:, i] + Wm[:, j]) / np.sqrt(2) for (i, j) in DY]
    return em, Wm, dyf, A, preps


def demod_basis_record():
    """The demodulation basis AS ARRAYS for the artifact: on paid data
    'the flown basis governs' is enforceable only if the basis is in
    the record (a 16-hex digest cannot be fallen back on)."""
    em, Wm, dyf, A, preps = build_preps()
    return {"em": em.tolist(), "Wm": Wm.tolist(),
            "dyf": list(dyf), "A": A.tolist(),
            "preps": [p.tolist() for p in preps]}


def demod_basis_sha():
    """sha of the demodulation basis (Floquet modes, dyad frequencies,
    site->dyad map, preps): the analyze side recomputes these via a
    greedy eigensolver ordering, so an artifact records what was flown
    and the analysis audits against it."""
    em, Wm, dyf, A, preps = build_preps()
    h = hashlib.sha256()
    for arr in (em, Wm, np.array(dyf), A, np.array(preps)):
        h.update(np.round(np.ascontiguousarray(arr), 12).tobytes())
    return h.hexdigest()[:16]


# ------------------------------------------------------ frozen bindings
def phase_table_full(arm, prep_idx, m_bind):
    """ONE frozen channel realization per (arm, prep, binding), drawn at
    the full depth [m_bind, STEPS_MAX, N]; every grid point slices its
    prefix. NESTED across depth is load-bearing: binding m then means
    THE SAME realization at every depth, which is (a) what the
    committed gate's run_bindings models (one trajectory recorded at
    grid points) and (b) what makes 'the binding, resampled jointly
    across depths' (section 6.6) a meaningful resampling unit at all.
    sigma_l = sqrt(4 gamma_l dt); N0 gets the zero table."""
    prof = ARM_PROFILE[arm]
    if prof is None:
        return np.zeros((m_bind, STEPS_MAX, N))
    rng = np.random.default_rng(
        [SEED_TABLES, ARM_NAMES.index(arm), prep_idx])
    sig = np.sqrt(4 * prof * GBAR_DT)
    return (rng.normal(0.0, 1.0, size=(m_bind, STEPS_MAX, N))
            * sig[None, None, :])


_TABLE_CACHE = {}


def phase_table(arm, tp_idx, prep_idx, m_bind):
    """The [m_bind, n_steps, N] prefix slice for one grid cell."""
    key = (arm, prep_idx, m_bind)
    if key not in _TABLE_CACHE:
        _TABLE_CACHE[key] = phase_table_full(arm, prep_idx, m_bind)
    return _TABLE_CACHE[key][:, :GRID_STEPS[tp_idx], :]


def table_sha(tab):
    return hashlib.sha256(np.ascontiguousarray(tab).tobytes()).hexdigest()[:16]


def dose_certificates(tab, arm):
    """Per-site and two-site realized retention on a flown table vs the
    exact e^{-2 gamma dt} / e^{-2(gamma_l+gamma_m) dt}; lit/lit,
    lit/dark, dark/dark pairs for the corner arms."""
    prof = ARM_PROFILE[arm]
    if prof is None or tab.shape[1] == 0:
        return None
    g = prof * GBAR_DT
    out = {"per_site_max_dev": 0.0, "pairs": {}}
    for l in range(N):
        ret = np.abs(np.mean(np.exp(1j * tab[:, :, l]), axis=0))
        dev = float(np.max(np.abs(ret - np.exp(-2 * g[l]))))
        out["per_site_max_dev"] = max(out["per_site_max_dev"], dev)
    lit = [l for l in range(N) if prof[l] > 0]
    dark = [l for l in range(N) if prof[l] == 0]
    pairs = {}
    if len(lit) >= 2:
        pairs["lit/lit"] = (lit[0], lit[1])
    if lit and dark:
        pairs["lit/dark"] = (lit[0], dark[0])
    if len(dark) >= 2:
        pairs["dark/dark"] = (dark[0], dark[1])
    if arm == "U":
        pairs = {"lit/lit": (0, 3)}
    for name, (l, m) in pairs.items():
        ret = np.abs(np.mean(np.exp(1j * (tab[:, :, l] - tab[:, :, m])),
                             axis=0))
        dev = float(np.max(np.abs(ret - np.exp(-2 * (g[l] + g[m])))))
        # CUMULATIVE realized retention over the trajectory (what the
        # estimator actually reads): the M-binding sample mean floors
        # at ~1/sqrt(M) while the target keeps decaying, an arm- and
        # depth-dependent SYSTEMATIC the per-step number cannot see
        # (round-3 finding: +3100% at the deep end; the round-4
        # recompute showed the weighted FIT does NOT import the
        # floor - rates unmoved - so this is a certificate gap, not
        # a verdict bias); REPORTED here, criterion ownership stays
        # with G3's per-M work
        cum = np.abs(np.mean(np.exp(
            1j * np.cumsum(tab[:, :, l] - tab[:, :, m], axis=1)),
            axis=0))
        tgt = np.exp(-2 * (g[l] + g[m]) * np.arange(1, tab.shape[1] + 1))
        out["pairs"][name] = {
            "sites": [l, m], "max_dev": dev,
            "cum_realized_deep": float(cum[-1]),
            "cum_target_deep": float(tgt[-1]),
            "cum_floor_ratio_deep": float(cum[-1] / max(tgt[-1], 1e-300)),
        }
    devs = [out["per_site_max_dev"]] + [p["max_dev"]
                                        for p in out["pairs"].values()]
    out["max_dev"] = max(devs)
    out["pass_0.02"] = bool(out["max_dev"] < DOSE_CRITERION)
    return out


# -------------------------------------------------- transversal classes
def profile_class(profile):
    """Mirror-transversal check + class on LOGICAL sites: exactly one
    lit site per mirror pair; maximizing iff sigma_outer != sigma_middle
    (pairs {0,5} and {2,3})."""
    lit = {l for l in range(N) if profile[l] > 0}
    for (a, b) in MIRROR_PAIRS:
        if (a in lit) == (b in lit):
            return "NOT_TRANSVERSAL"
    s_out = +1 if 0 in lit else -1
    s_mid = +1 if 2 in lit else -1
    return "maximizing" if s_out != s_mid else "non-maximizing"


# ------------------------------------------------ transpile + assertions
def transpile_local(circs):
    return transpile(circs, basis_gates=LOCAL_BASIS, optimization_level=1,
                     seed_transpiler=SEED_TRANSPILER)


def assert_skeleton_invariants(tc, n_steps, chain=None, label=""):
    """Class-1 structural assertions on ONE transpiled skeleton:
    (a) every 2q op is rzz (CZ = NO-FLIGHT), angles in (0, pi/2];
    (b) 2q budget = 10 (prep) + 16 per step, exactly;
    (c) exactly one PARAMETERIZED rz per (site, step), matched by name;
    (d) routing added nothing (active qubits inside the chain)."""
    phys = list(range(N)) if chain is None else list(chain)
    site_of = {p: l for l, p in enumerate(phys)}
    ops2 = []
    par_count = {}
    active = set()
    for inst in tc.data:
        if inst.operation.name in ("barrier",):
            continue
        qs = [tc.find_bit(q).index for q in inst.qubits]
        active |= set(qs)
        if inst.operation.num_qubits == 2:
            ops2.append((inst.operation.name, tuple(qs),
                         inst.operation.params))
        if inst.operation.name == "rz":
            pars = [p for p in inst.operation.params
                    if isinstance(p, ParameterExpression) and p.parameters]
            for p in pars:
                for par in p.parameters:
                    par_count[par.name] = par_count.get(par.name, 0) + 1
                    # the phase must enter with coefficient EXACTLY 1
                    # (an offset from a merged fixed rz is fine)
                    _require(len(p.parameters) == 1,
                             f"{label}: rz expression carries "
                             f"{len(p.parameters)} parameters (merged "
                             f"injections?)")
                    slope = (float(p.assign(par, 1.0))
                             - float(p.assign(par, 0.0)))
                    _require(abs(slope - 1.0) < 1e-12,
                             f"{label}: parameter {par.name} carries "
                             f"coefficient {slope}, not 1")
    bad2 = [o for o in ops2 if o[0] != "rzz"]
    _require(not bad2,
             f"{label}: non-rzz 2q gate (CZ = NO-FLIGHT): {bad2[:4]}")
    n_expect = 10 + 16 * n_steps
    _require(len(ops2) == n_expect,
             f"{label}: 2q budget broken: {len(ops2)} != {n_expect}")
    for name, qs, params in ops2:
        ang = float(params[0])
        _require(0 < ang <= np.pi / 2 + 1e-9,
                 f"{label}: rzz angle out of fractional range: {ang}")
    expected = {f"phi_q{site}_s{s:02d}" for site in range(N)
                for s in range(n_steps)}
    got = {p.name for p in tc.parameters}
    _require(got == expected,
             f"{label}: parameter set mismatch "
             f"({len(got)} vs {len(expected)})")
    for nm in expected:
        _require(par_count.get(nm, 0) == 1,
                 f"{label}: parameter {nm} appears {par_count.get(nm, 0)}x"
                 f" (exactly-one-injected-RZ assertion)")
    meas_target = set(phys)
    _require(active <= meas_target,
             f"{label}: routing left the line: "
             f"{sorted(active - meas_target)}")
    # DD-off face: the transpiled science circuit must contain NO delay
    # (an inserted DD sequence carries delays + extra 1q pulses; the 2q
    # budget above catches 2q insertions, this catches the schedule)
    for inst in tc.data:
        _require(inst.operation.name != "delay",
                 f"{label}: delay in science circuit (DD/scheduling ran?)")
    # layout identity: logical site l must sit on chain[l]; a
    # permutation WITHIN the line would relabel sites and swap the
    # transversal classes silently
    if chain is not None:
        fin = None
        try:
            fin = list(tc.layout.final_index_layout())
        except Exception:
            pass
        if fin is not None:
            _require(fin == list(chain),
                     f"{label}: final layout {fin} != pinned chain "
                     f"{list(chain)}")
        else:
            raise RuntimeError(
                f"{label}: cannot read final_index_layout; the identity-"
                f"layout assertion is mandatory on hardware transpiles")
    return {"n_2q": len(ops2), "params": len(got)}


def assert_shared_skeleton(pubs, meta):
    """The four arms of every (depth, prep) cell must fly the IDENTICAL
    transpiled circuit object, differing only in bound phase tables
    (section 3's shared-skeleton invariant, asserted in code)."""
    by_cell = {}
    for info in meta:
        key = (info["tp_idx"], info["prep_idx"])
        by_cell.setdefault(key, []).append(id(pubs[info["pub_index"]][0]))
    for key, ids in by_cell.items():
        _require(len(ids) == 4 and len(set(ids)) == 1,
                 f"shared-skeleton assertion FAILED at cell {key}")


# --------------------------------------------------- counts -> analysis
def counts_to_vec(counts):
    v = np.zeros(D6)
    for bs, c in counts.items():
        b = bs.replace(" ", "")
        v[int(b, 2)] += c        # bit l of int(b,2) IS qubit l (measure_all)
    return v


def confusion_matrices(cal0, cal1):
    """Per-qubit asymmetric confusion from the CAL PUBs; C_q maps
    true -> measured."""
    v0, v1 = counts_to_vec(cal0), counts_to_vec(cal1)
    n0, n1 = v0.sum(), v1.sum()
    Cs = []
    for q in range(N):
        bit = (np.arange(D6) >> q) & 1
        eps0 = float(v0[bit == 1].sum() / n0)   # P(read 1 | true 0)
        eps1 = float(v1[bit == 0].sum() / n1)   # P(read 0 | true 1)
        Cs.append(np.array([[1 - eps0, eps1], [eps0, 1 - eps1]]))
    return Cs


def inverse_confusion_full(Cs):
    Minv = np.array([[1.0]])
    for q in range(N - 1, -1, -1):   # qubit 5 = most significant bit
        Minv = np.kron(Minv, np.linalg.inv(Cs[q]))
    return Minv


CFG_1M = [1 << l for l in range(N)]


def cell_raw_matrix(counts_list):
    """Per-binding RAW count vectors [M, 64]; no reduction here. The
    binding stays the resampling unit as raw data; the confusion
    inversion is applied ONCE to the (resampled) pooled vector, so the
    chain stays linear until the single post-pooling clip (the
    per-binding invert-then-clip variant is nonlinear at 16
    shots/binding and biases depth-dependently; measured in the G7
    review round and verified: Minv(sum v) = sum(Minv v) to fp
    precision, ~1e-12)."""
    return np.stack([counts_to_vec(cd) for cd in counts_list])


def pool_invert_postselect(raw, Minv, idx=None, floor=None):
    """(Resampled) pooled raw vector -> one confusion inversion -> one
    clip -> post-selected one-magnon occupations. Returns (occ[N],
    kept_weight, clipped_mass) where clipped_mass counts the negative
    quasi-counts over the FULL 64-dim inverted vector (the inversion's
    strain lives mostly outside the one-magnon sector).

    `floor` is the kept-count floor IN FORCE (round 13: the module
    constant was still the one the estimator read, so a re-freeze of
    kept_count_floor would have been honoured by the submit gate and
    ignored by every analysis, including a re-analysis of an already
    flown artifact). Analyze passes the frozen value."""
    v = raw.sum(axis=0) if idx is None or raw.shape[0] == 1 \
        else raw[idx].sum(axis=0)
    q = Minv @ v
    clipped = float(-np.clip(q, None, 0).sum())
    w = np.clip(q[CFG_1M], 0, None)
    tot = float(w.sum())
    if tot < (KEPT_COUNT_FLOOR if floor is None else floor):
        return np.full(N, np.nan), 0.0, clipped
    return w / tot, tot, clipped


# ------------------------------------------------------------ PUB plan
def build_pub_plan():
    """The deterministic PUB order: depth ascending; within each depth
    arm N0 -> U -> Cp -> C; preps (1,2),(2,3),(3,5) within each arm."""
    plan = []
    for tp_idx in range(len(GRID_STEPS)):
        for arm in ARM_NAMES:
            for prep_idx in range(3):
                plan.append({"arm": arm, "tp_idx": tp_idx,
                             "prep_idx": prep_idx})
    return plan


def make_pubs(transpiled, plan, m_bind, shots_bind, shots_science):
    """PUB per plan entry; every arm at depth > 0 is SWEPT over the
    m_bind bindings of ONE shared transpiled skeleton, bound by
    parameter NAME (post-transpile assertion). N0's table is the zero
    table (the pre-registration's "phi = 0 bindings", plural: N0 keeps
    the binding axis, so var(B-hat) and the bootstrap reach it; note
    N0's bindings all bind the identical zero table, so its binding
    axis carries SHOT noise only, no channel-realization variance -
    its bootstrap spread is not on the other arms' footing). Depth 0
    has no parameters; all four arms fly the identical circuit there."""
    pubs, meta = [], []
    for entry in plan:
        arm, tp_idx, prep_idx = entry["arm"], entry["tp_idx"], entry["prep_idx"]
        tc = transpiled[(tp_idx, prep_idx)]
        n_steps = GRID_STEPS[tp_idx]
        tab = phase_table(arm, tp_idx, prep_idx, m_bind)
        info = dict(entry)
        info["table_sha"] = table_sha(tab)
        info["pub_index"] = len(pubs)
        if n_steps == 0:
            pubs.append((tc, None, shots_science))
            info["swept"] = False
        else:
            tpars = sorted(tc.parameters, key=lambda p: p.name)
            # SamplerV2 binds a bare ndarray in circuit.parameters
            # order; the sweep columns are placed by sorted name. These
            # coincide on qiskit 2.x (ParameterView is name-sorted and
            # s%02d zero-padding makes the order well-defined), and a
            # divergence would make the flown dataset silently
            # worthless, so it is a hard guard, not a convention:
            _require([p.name for p in tc.parameters]
                     == [p.name for p in tpars],
                     "circuit.parameters is not name-sorted: sweep "
                     "columns would bind to the wrong parameters")
            name_to_col = {p.name: c for c, p in enumerate(tpars)}
            sweep = np.zeros((m_bind, len(tpars)))
            for site in range(N):
                for s in range(n_steps):
                    col = name_to_col.get(f"phi_q{site}_s{s:02d}")
                    _require(col is not None,
                             f"post-transpile parameter lost: q{site} s{s}")
                    sweep[:, col] = tab[:, s, site]
            pubs.append((tc, sweep, shots_bind))
            info["swept"] = True
        meta.append(info)
    return pubs, meta


def rebuild_table_from_pub(pub, tc, n_steps):
    """Reconstruct the per-binding phase table [M, n_steps, N] from the
    sweep array actually bound into a PUB (certifying WHAT FLIES, not
    the pre-binding table; catches binding/transposition bugs)."""
    sweep = pub[1]
    if sweep is None:
        return None
    tpars = sorted(tc.parameters, key=lambda p: p.name)
    name_to_col = {p.name: c for c, p in enumerate(tpars)}
    tab = np.zeros((sweep.shape[0], n_steps, N))
    for site in range(N):
        for s in range(n_steps):
            tab[:, s, site] = sweep[:, name_to_col[f"phi_q{site}_s{s:02d}"]]
    return tab


def cal_circuits():
    c0 = QuantumCircuit(N)
    c0.measure_all()
    c1 = QuantumCircuit(N)
    c1.x(range(N))
    c1.measure_all()
    return [c0, c1]


SHOTS_AUX = 4096      # aux T1/T2 PUBs need survival fractions, not fits


def t1t2_station_circuits(d_t1, d_echo):
    """ONE station's aux set, flown identically at EVERY station so
    drift and delay are deconfounded (the G7 round-2 blocker: one
    delay per station made mid-batch T1 telegraphing, the recorded
    staircase killer, algebraically invisible):
    - T1 at two delays (short, long); CAL1 is the shared zero point;
    - a fixed-delay Ramsey (T2*-class): its absolute value is
      detuning-scrambled, but the SAME circuit repeated across
      stations is a clean drift canary;
    - a Hahn echo at the same total delay: detuning-immune, fittable
      dephasing readout."""
    circs, meta = [], []
    for d in (d_t1, 3 * d_t1):
        qt = QuantumCircuit(N)
        qt.x(range(N))
        for q in range(N):
            qt.delay(d, q, unit="us")
        qt.measure_all()
        circs.append(qt); meta.append({"kind": "T1", "delay_us": d})
    qr = QuantumCircuit(N)
    qr.sx(range(N))
    for q in range(N):
        qr.delay(d_echo, q, unit="us")
    qr.sx(range(N))
    qr.measure_all()
    circs.append(qr); meta.append({"kind": "Ramsey_canary",
                                   "delay_us": d_echo})
    qe = QuantumCircuit(N)
    qe.sx(range(N))
    for q in range(N):
        qe.delay(d_echo / 2, q, unit="us")
    qe.x(range(N))
    for q in range(N):
        qe.delay(d_echo / 2, q, unit="us")
    qe.sx(range(N))
    qe.measure_all()
    circs.append(qe); meta.append({"kind": "T2echo", "delay_us": d_echo})
    # a CAL pair per station: readout confusion is the one calibration
    # input every cell consumes, and the head-of-batch pair alone
    # cannot see mid-batch readout drift (the station pairs are drift
    # DIAGNOSTICS; the head pair stays the section-5 inversion source)
    for kind, cal in zip(("CAL0_station", "CAL1_station"),
                         cal_circuits()):
        circs.append(cal); meta.append({"kind": kind, "delay_us": 0.0})
    return circs, meta


# ------------------------------------------------------------ save/load
def _save_json(payload, mode):
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = RESULTS_DIR / f"corner_beat_{mode}_{ts}.json"
    with open(path, "w") as f:
        json.dump(payload, f, indent=2, default=str)
    print(f"\n  saved -> {path}")
    return path


def load_constants():
    if CONSTANTS_PATH.exists():
        return json.load(open(CONSTANTS_PATH))
    return DEFAULT_CONSTANTS


# ------------------------------------------------------------- certify
def cmd_certify():
    print("CERTIFY (local, no network)\n")
    if not __debug__:
        raise SystemExit("ABORT: never run --certify under python -O "
                         "(its record feeds the Class-1 gate)")
    em, Wm, dyf, A, preps = build_preps()

    from qiskit.quantum_info import Operator, Statevector
    cfg = [1 << l for l in range(N)]

    # (1) Strang-step sector parity
    qc = QuantumCircuit(N)
    strang_step_circuit(qc, JDT)
    Uc = Operator(qc).data[np.ix_(cfg, cfg)]
    Ug = strang(JDT)
    kk = np.unravel_index(np.argmax(np.abs(Ug)), Ug.shape)
    ph = Uc[kk] / Ug[kk]
    dev = np.max(np.abs(Uc - ph * Ug))
    print(f"  (1) Strang-step sector parity: max|dev| = {dev:.3e}")
    _require(dev < 1e-10, 'sector parity')
    Ufull = Operator(qc).data
    pop = np.array([bin(a).count('1') for a in range(D6)])
    leak = np.max(np.abs(Ufull)[pop[:, None] != pop[None, :]])
    print(f"      number conservation of the step: {leak:.3e}")
    _require(leak < 1e-12, 'number conservation')

    # (2) prep parity, all three dyads
    prep_devs = []
    for p, (i, j) in enumerate(DY):
        qp = QuantumCircuit(N)
        prep_circuit_ops(qp, preps[p])
        sv = Statevector.from_instruction(qp).data[cfg]
        ref = preps[p]
        php = sv[np.argmax(np.abs(ref))] / ref[np.argmax(np.abs(ref))]
        d = np.max(np.abs(sv - php * ref))
        print(f"  (2) prep dyad {DY[p]}: max|dev| = {d:.3e}")
        _require(d < 1e-10, 'prep parity')
        prep_devs.append(float(d))

    # (3) trajectory parity vs the committed-gate simulator: one frozen
    # binding, no noise, occupations at every grid point to depth 12
    steps_chk = 12
    tab = phase_table("C", GRID_STEPS.index(12), 0, 4)[:1]  # 1 binding
    Ug1 = strang(JDT)
    Psi = preps[0].astype(complex).copy()
    occ_ref = [np.abs(Psi) ** 2]
    for s in range(steps_chk):
        Psi = Ug1 @ Psi
        Psi = Psi * np.exp(1j * tab[0, s, :])
        occ_ref.append(np.abs(Psi) ** 2)
    qt = QuantumCircuit(N)
    prep_circuit_ops(qt, preps[0])
    worst = 0.0
    for s in range(steps_chk):
        strang_step_circuit(qt, JDT)
        for site in range(N):
            qt.rz(float(tab[0, s, site]), site)
        sv = Statevector.from_instruction(qt).data[cfg]
        worst = max(worst, float(np.max(np.abs(np.abs(sv) ** 2
                                               - occ_ref[s + 1]))))
    print(f"  (3) trajectory parity vs gate simulator (12 steps, one "
          f"frozen binding): max|dev| = {worst:.3e}")
    _require(worst < 1e-10, 'trajectory parity')

    # (4) transpile-shape assertions on the local rzz basis
    stats = {}
    for tp_idx in [0, 1, len(GRID_STEPS) - 1]:
        n_steps = GRID_STEPS[tp_idx]
        for prep_idx in range(3):
            qs, _ = science_skeleton(n_steps, preps[prep_idx])
            tc = transpile_local(qs)
            st = assert_skeleton_invariants(
                tc, n_steps, label=f"t{tp_idx}/p{prep_idx}")
            stats[(tp_idx, prep_idx)] = st
    deep = stats[(len(GRID_STEPS) - 1, 0)]
    print(f"  (4) transpile shape: all-rzz, angles in (0, pi/2], one "
          f"parameterized rz per site/step; deep-end 2q = {deep['n_2q']} "
          f"(expected {10 + 16 * STEPS_MAX})")

    # (5) profile classes
    for arm in ("C", "Cp"):
        cls = profile_class(ARM_PROFILE[arm])
        want = "maximizing" if arm == "C" else "non-maximizing"
        print(f"  (5) profile {arm}: {cls}")
        _require(cls == want, f'profile class {cls} != {want}')

    # (6) dose certificates on ALL 180 swept cells: 9 nested tables
    # (arm x prep) read at their 20 depth prefixes; the "worst of
    # 180" is a running maximum along 9 realizations
    worst_dev, worst_cell, table_report = 0.0, None, {}
    for arm in ("U", "Cp", "C"):
        for tp_idx in range(1, len(GRID_STEPS)):
            for prep_idx in range(3):
                tab = phase_table(arm, tp_idx, prep_idx, M_BIND)
                cert = dose_certificates(tab, arm)
                if cert["max_dev"] > worst_dev:
                    worst_dev = cert["max_dev"]
                    worst_cell = (arm, GRID_STEPS[tp_idx], prep_idx)
                    table_report["worst"] = cert
    print(f"  (6) dose certificates, ALL 180 swept tables, M = {M_BIND}:"
          f" worst max|dev| = {worst_dev:.4f} at {worst_cell} "
          f"({'PASS' if worst_dev < DOSE_CRITERION else 'FAIL'} "
          f"@ {DOSE_CRITERION})")
    _require(worst_dev < DOSE_CRITERION, 'dose certificate')

    # (7) fit parity against the COMMITTED gate: same traces, same
    # r_max, the rate must agree exactly (the runner's fit_rate extends
    # the return value, never the search)
    import importlib.util as _ilu
    _require(GATE_PATH.exists(), f"committed gate not found: {GATE_PATH}")
    spec = _ilu.spec_from_file_location("cbgate", GATE_PATH)
    gate = _ilu.module_from_spec(spec)
    spec.loader.exec_module(gate)
    rng7 = np.random.default_rng(7)
    ts7 = fit_time_axis()         # the parity test runs the FLOWN axis
    worst7 = 0.0
    for trial in range(6):
        r_true = 0.2 + 0.5 * rng7.random()
        om7 = -2.0 + 0.1 * rng7.standard_normal()
        y = (0.05 + 0.2 * np.exp(-r_true * ts7) * np.cos(om7 * ts7 + 0.3)
             + 0.01 * rng7.standard_normal(len(ts7)))
        w = 1000 * np.exp(-ts7) + 50
        # Round 11 (MAJOR): every parity case was a SINGLE trace, and
        # with one trace the per-trace peak normalization section 6.3
        # pins is a no-op. The joint branch is the only place the
        # weighting does anything (the C and C' +- channels are fitted
        # jointly), so the certified 0.0 rested on evidence that never
        # touched the claim. Second trace: a different rate, phase and
        # a peak weight 40x smaller, so unequal-scale pooling is what
        # is compared.
        y2 = (0.04 + 0.15 * np.exp(-(r_true + 0.3) * ts7)
              * np.cos(om7 * ts7 - 0.7) + 0.01 * rng7.standard_normal(len(ts7)))
        w2 = 25 * np.exp(-ts7) + 2
        for rm in (0.8, R_MAX_FIT):
            for traces in ([(y, w)], [(y, w), (y2, w2)]):
                rg = gate.fit_rate(traces, ts7, -2.0, rm)
                rr, _, _ = fit_rate(traces, ts7, -2.0, rm)
                worst7 = max(worst7, abs(rg - rr))
    print(f"  (7) fit parity vs committed gate (6 x 2 r_max x "
          f"[single, joint two-trace at 40x unequal peak weight]): "
          f"max|dev| = {worst7:.2e}")
    _require(worst7 == 0.0, "fit_rate diverged from the committed gate")

    # (8) THE AXIS ITSELF, not the fits on a given axis. Check (7) hands
    # ONE array to both implementations, so it can certify that the two
    # FITS agree and never that the gate CONSTRUCTS the same axis; the
    # pre-registration said so in as many words and left the gate half
    # to be verified by reading the file. Both sides now derive the axis
    # from a switch, and this compares them on BOTH settings, so a
    # refreeze that throws the switch on one side only is caught here
    # instead of in the fits. Exact equality: both compute the same
    # product of the same integers, so any residual is a defect in the
    # construction rather than rounding (house rule).
    axis_devs = {}
    saved_axis = TIME_AXIS
    try:
        for axis in ("nominal", "realized_dose"):
            globals()["TIME_AXIS"] = axis
            mine = fit_time_axis()
            theirs = gate.fit_time_axis(STEPS_MAX, K_GRID, JDT, axis=axis)
            same = (len(mine) == len(theirs)
                    and bool(np.array_equal(mine, theirs)))
            axis_devs[axis] = same
            _require(same, f"gate and runner disagree on the {axis} fit "
                           f"axis: gate {theirs} vs runner {mine}")
    finally:
        globals()["TIME_AXIS"] = saved_axis
    print(f"  (8) fit-axis parity vs committed gate, both settings "
          f"(nominal, realized_dose): EXACT; axis in force = {TIME_AXIS}")
    _require(gate.TIME_AXIS == TIME_AXIS,
             f"gate TIME_AXIS {gate.TIME_AXIS!r} != runner {TIME_AXIS!r}: "
             f"the refreeze must throw the switch on BOTH sides")

    print("\nCERTIFY: ALL CHECKS PASS")
    _save_json({"mode": "certify",
                "debug_mode": bool(__debug__),
                "basis_sha": demod_basis_sha(),
                "sector_parity_max_dev": float(dev),
                "number_conservation_leak": float(leak),
                "prep_parity_max_devs": prep_devs,
                "trajectory_parity_max_dev": float(worst),
                "worst_dose_dev": worst_dev,
                "worst_dose_cell": worst_cell,
                "fit_parity_max_dev": worst7,
                "fit_axis_parity": axis_devs,
                "time_axis_runner": TIME_AXIS,
                "time_axis_gate": gate.TIME_AXIS,
                "transpile_stats": {f"{k}": v for k, v in stats.items()},
                "dose": table_report}, "certify")


# ---------------------------------------------------------------- aer
# ONE declaration of the synthetic noise point, because it was two: the
# function defaults and the artifact's "noise" record were independent
# literals, so a change to the simulation would have been recorded
# wrongly in the artifact that proves the chain. Note what Amendment 1.8
# did to its MEANING: p2 = 0.005 used to sit one guard-step BELOW the
# 0.6% bound, and now it IS the bound, so the synthetic proof runs at
# the worst device the Class-1 guard admits.
AER_NOISE = {"p2": 0.005, "readout": 0.015, "t1_us": 250.0,
             "t2_us": 200.0}


def build_noise_model(p2=AER_NOISE["p2"], ro=AER_NOISE["readout"],
                      t1_us=AER_NOISE["t1_us"], t2_us=AER_NOISE["t2_us"],
                      dur_2q_ns=100.0):
    """Hand-built (never from_backend): depolarizing + thermal on every
    rzz orientation of the line, asymmetric-capable readout."""
    from qiskit_aer.noise import (NoiseModel, depolarizing_error,
                                  thermal_relaxation_error, ReadoutError)
    nm = NoiseModel(basis_gates=["rzz", "rz", "sx", "x"])
    dep = depolarizing_error(4.0 / 3.0 * p2, 2)
    for a in range(N - 1):
        t = thermal_relaxation_error(t1_us * 1e-6,
                                     min(t2_us, 2 * t1_us) * 1e-6,
                                     dur_2q_ns * 1e-9)
        te = t.expand(thermal_relaxation_error(
            t1_us * 1e-6, min(t2_us, 2 * t1_us) * 1e-6, dur_2q_ns * 1e-9))
        nm.add_quantum_error(dep.compose(te), "rzz", [a, a + 1])
        nm.add_quantum_error(dep.compose(te), "rzz", [a + 1, a])
    for q in range(N):
        nm.add_readout_error(ReadoutError([[1 - ro, ro], [ro, 1 - ro]]), [q])
    return nm


def cmd_aer(full=False):
    m_bind = M_BIND if full else 64
    shots_bind = SHOTS_SCIENCE // m_bind
    print(f"AER synthetic artifact: M = {m_bind}, "
          f"{shots_bind} shots/binding (hardware-shaped schema"
          f"{'' if full else ', reduced M'})\n")
    em, Wm, dyf, A, preps = build_preps()
    plan = build_pub_plan()
    transpiled = {}
    for tp_idx in range(len(GRID_STEPS)):
        for prep_idx in range(3):
            qs, _ = science_skeleton(GRID_STEPS[tp_idx], preps[prep_idx])
            transpiled[(tp_idx, prep_idx)] = transpile_local(qs)
    # invariants on every skeleton (cheap, and the flight runs the same)
    for (tp_idx, prep_idx), tc in transpiled.items():
        assert_skeleton_invariants(tc, GRID_STEPS[tp_idx],
                                   label=f"t{tp_idx}/p{prep_idx}")
    pubs, meta = make_pubs(transpiled, plan, m_bind, shots_bind,
                           SHOTS_SCIENCE)
    assert_shared_skeleton(pubs, meta)
    cal_tc = transpile_local(cal_circuits())
    for c in cal_tc:
        pubs.append((c, None, SHOTS_SCIENCE))
    print(f"  {len(pubs)} PUBs (252 science + 2 CAL; T1/T2* skipped on "
          f"the simulator)")

    from qiskit_aer.primitives import SamplerV2 as AerSampler
    nm = build_noise_model()
    sampler = AerSampler(options={
        "backend_options": {"noise_model": nm,
                            "seed_simulator": SEED_AER}})
    print("  running Aer (grouped submission: aer's batch postprocessing "
          "chokes on the full 254-PUB list at M = 1024) ...")
    GROUP = 12
    results = []
    for lo in range(0, len(pubs), GROUP):
        chunk = pubs[lo:lo + GROUP]
        res = sampler.run(chunk).result()
        results.extend(_pub_counts_list(pr) for pr in res)
        print(f"    PUBs {lo}..{lo + len(chunk) - 1} done")

    raw = {"science": {}, "cal": {}, "t1t2": []}
    for info in meta:
        arm, tp, pi = info["arm"], info["tp_idx"], info["prep_idx"]
        counts = results[info["pub_index"]]
        raw["science"].setdefault(arm, {}).setdefault(
            str(GRID_STEPS[tp]), {})[str(pi)] = counts
    raw["cal"]["cal0"] = results[len(meta)][0]
    raw["cal"]["cal1"] = results[len(meta) + 1][0]

    payload = {
        "mode": "aer", "synthetic": True,
        "timestamp": datetime.now().isoformat(),
        "basis_sha": demod_basis_sha(),
        "basis": demod_basis_record(),
        "config": {"N": N, "JDT": JDT, "Q": Q, "K_GRID": K_GRID,
                   "GRID_STEPS": GRID_STEPS, "M_BIND": m_bind,
                   "SHOTS_BIND": shots_bind,
                   "SHOTS_SCIENCE": SHOTS_SCIENCE,
                   "seed_tables": SEED_TABLES, "seed_aer": SEED_AER,
                   "numpy": np.__version__},
        "noise": dict(AER_NOISE),
        "pub_meta": meta,
        "raw_counts": raw,
    }
    path = _save_json(payload, "aer")
    print("  analyze with: python run_corner_beat.py --analyze "
          f"\"{path}\"")


def _pub_counts_list(pub_result):
    """1 counts dict for an unswept PUB, M dicts for a swept one
    (the concentrator idiom, reading the 'meas' register)."""
    databin = pub_result.data
    bitarr = getattr(databin, "meas", None)
    if bitarr is None:
        for nm in list(getattr(databin, "keys", lambda: [])()):
            cand = getattr(databin, nm, None)
            if hasattr(cand, "get_counts"):
                bitarr = cand
                break
    if bitarr is None:
        raise RuntimeError("PUB result has no readable classical register")
    # PER BINDING via get_counts(k): a bare get_counts() POOLS across the
    # sweep axis (the concentrator's recorded instrument deviation, which
    # the pre-registration pins this runner to bypass)
    shape = tuple(getattr(bitarr, "shape", ()) or ())
    if not shape:
        c = bitarr.get_counts()
        return [{str(k): int(v) for k, v in c.items()}]
    m = int(np.prod(shape))
    return [{str(k): int(v) for k, v in bitarr.get_counts(i).items()}
            for i in range(m)]


# ------------------------------------------------------------- analyze
def cmd_analyze(path):
    print(f"ANALYZE {path}")
    print("  (the committed pre-registration governs; this printout is "
          "not the verdict)\n")
    art = json.load(open(path))
    cfgd = art["config"]
    for key, want in (("GRID_STEPS", GRID_STEPS), ("JDT", JDT), ("Q", Q),
                      ("M_BIND", None), ("SHOTS_SCIENCE", SHOTS_SCIENCE),
                      ("seed_tables", SEED_TABLES)):
        if want is not None:
            _require(cfgd[key] == want,
                     f"artifact config mismatch: {key} = {cfgd[key]} "
                     f"!= {want}")
    m_bind = cfgd["M_BIND"]
    if art.get("mode") == "hardware":
        _require(m_bind == M_BIND,
                 f"hardware artifact at M_BIND = {m_bind}, frozen "
                 f"config is {M_BIND}")
    if "basis" in art:
        # THE FLOWN BASIS GOVERNS, from the record itself
        b = art["basis"]
        em = np.array(b["em"]); Wm = np.array(b["Wm"])
        dyf = list(b["dyf"]); A = np.array(b["A"])
        preps = [np.array(p) for p in b["preps"]]
        print("  demodulation basis: loaded from the artifact "
              "(the flown basis governs)")
    else:
        em, Wm, dyf, A, preps = build_preps()
        print("  demodulation basis: RECOMPUTED (artifact predates "
              "basis persistence)")
    ts = fit_time_axis()          # governed by the frozen time_axis
    consts = load_constants()

    # frozen-table audit: recompute every cell's table and check the
    # persisted sha. A WARNING, never an abort: the tables are archival
    # (the estimator consumes counts only), a hardware flight persists
    # them as .npz, and an analysis of paid data must not die on a
    # numpy-stream bookkeeping check.
    sha_bad = []
    for info in art.get("pub_meta", []):
        tab = phase_table(info["arm"], info["tp_idx"], info["prep_idx"],
                          m_bind)
        if table_sha(tab) != info["table_sha"]:
            sha_bad.append((info["arm"], info["tp_idx"], info["prep_idx"]))
    print(f"  frozen-table sha audit: "
          + ("OK (all cells)" if not sha_bad else
             f"WARNING: {len(sha_bad)} cells mismatch (numpy stream "
             f"drift?); the flown .npz record governs: "
             f"{art.get('tables_npz', 'ABSENT')}"))
    bsha = demod_basis_sha()
    art_bsha = art.get("basis_sha")
    print(f"  demodulation-basis sha audit: "
          + ("OK" if art_bsha == bsha else
             f"WARNING: artifact {art_bsha} != recomputed {bsha} "
             f"(eigensolver/library drift; the flown basis governs)"))

    cal_d = art["raw_counts"].get("cal", {})
    cal0_keys = sorted(k for k in cal_d if k.startswith("cal0"))
    cal1_keys = sorted(k for k in cal_d if k.startswith("cal1"))
    if not cal0_keys or not cal1_keys:
        raise SystemExit("ABORT: CAL PUBs missing from artifact; "
                         "the confusion inversion cannot be built")

    def _sum_counts(keys):
        tot = {}
        for k in keys:
            for bs, c in cal_d[k].items():
                tot[bs] = tot.get(bs, 0) + c
        return tot
    # ALL CAL pairs pool into the inversion (head + one pair per job,
    # Amendment 1); the per-pair readings are the drift record
    cal0, cal1 = _sum_counts(cal0_keys), _sum_counts(cal1_keys)
    if len(cal0_keys) > 1:
        print(f"  CAL pairs pooled into the inversion: {cal0_keys}")
    # per-pair drift record (Amendment 1.5's registered promise;
    # round-8 repair: it was recoverable from the raw artifact but
    # never computed): per-qubit (eps0, eps1) for each pair
    cal_pairs_drift = {}
    for k0 in cal0_keys:
        k1 = k0.replace("cal0", "cal1")
        if k1 in cal_d:
            try:
                Cp = confusion_matrices(cal_d[k0], cal_d[k1])
                cal_pairs_drift[k0] = [
                    (float(C[1, 0]), float(C[0, 1])) for C in Cp]
            except Exception as e:
                cal_pairs_drift[k0] = f"unreadable: {e!r}"
    # registered estimator constants: the frozen file value governs
    # once frozen; the module default is the pre-freeze provisional
    r_max_fit = (consts["R_MAX_FIT"]["value"]
                 if consts.get("R_MAX_FIT", {}).get("frozen")
                 else R_MAX_FIT)
    r_boot = int(consts["R_BOOT"]["value"]
                 if consts.get("R_BOOT", {}).get("frozen") else R_BOOT)
    kept_floor = float(consts["kept_count_floor"]["value"]
                       if consts.get("kept_count_floor", {}).get("frozen")
                       else KEPT_COUNT_FLOOR)
    Cs = confusion_matrices(cal0, cal1)
    Minv = inverse_confusion_full(Cs)
    eps_ro = [(float(C[1, 0]), float(C[0, 1])) for C in Cs]
    print("  CAL confusion (eps0, eps1) per qubit: "
          + ", ".join(f"({a:.3f},{b:.3f})" for a, b in eps_ro))

    # per-cell per-binding RAW vectors; binding-count hard assertion
    raws = {}      # (arm, tp_idx, prep_idx) -> [M_or_1, 64]
    missing = []
    for arm in ARM_NAMES:
        for tp_idx, n_steps in enumerate(GRID_STEPS):
            for pi in range(3):
                try:
                    cl = art["raw_counts"]["science"][arm][str(n_steps)][
                        str(pi)]
                except KeyError:
                    missing.append((arm, n_steps, pi))
                    raws[(arm, tp_idx, pi)] = None
                    continue
                raw = cell_raw_matrix(cl)
                if n_steps > 0 and raw.shape[0] != m_bind:
                    # per-binding structure lost (get_counts pooling?):
                    # EXCLUDE the cell loudly rather than kill the
                    # analysis of an un-repeatable dataset
                    missing.append((arm, n_steps, pi,
                                    f"bindings {raw.shape[0]}"))
                    raws[(arm, tp_idx, pi)] = None
                    continue
                raws[(arm, tp_idx, pi)] = raw
    if missing:
        print(f"  WARNING: {len(missing)} cells excluded (missing or "
              f"binding structure lost; traces degrade to NaN there): "
              f"{missing[:6]} ...")

    clip_log = []
    boot_floor_trips = []   # replicate-level floor trips (round 13)
    floored = []  # cells that FLEW but fell under the kept-count
                  # floor on the full (un-resampled) data: they are
                  # grid-incompleteness like a missing cell (T1 /
                  # leakage / readout pathologies kill post-selection
                  # preferentially at the deep end, exactly the rate
                  # lever arm, so a silent drop would bias the fit
                  # with no trace in the record)

    def arm_traces(arm, idx=None):
        Yp, Wp = [], []
        for pi in range(3):
            Y, W = [], []
            for tp_idx in range(len(GRID_STEPS)):
                raw = raws[(arm, tp_idx, pi)]
                if raw is None:
                    Y.append(np.full(N, np.nan)); W.append(0.0)
                    continue
                y, w, clip = pool_invert_postselect(raw, Minv, idx=idx,
                                                    floor=kept_floor)
                if idx is None:
                    clip_log.append(clip)
                    if w == 0.0:
                        floored.append((arm, GRID_STEPS[tp_idx], pi))
                elif w == 0.0:
                    # round 13: a RESAMPLED pool can dip under the floor
                    # at the deep end (~208 kept expected at the guard
                    # bound since Amendment 1.8; ~87 before it),
                    # dropping that grid point from that
                    # replicate's fit silently. Same mechanism as the
                    # NaN-replicate exclusion section 7 registers, one
                    # level down, and it narrows the CI the same way.
                    boot_floor_trips.append((arm, GRID_STEPS[tp_idx], pi))
                Y.append(y); W.append(w)
            Yp.append(np.array(Y)); Wp.append(np.array(W))
        return Yp, Wp

    stats = {}
    for arm in ARM_NAMES:
        Yp, Wp = arm_traces(arm)
        rates, dws, sses = arm_rates(Yp, Wp, A, dyf, ts, arm, r_max_fit)
        # INFORMATIONAL saturation flag only (round 15: 0.95 is a
        # display heuristic and is registered nowhere; the REGISTERED
        # trigger is r_saturation_frac with its declared base, applied
        # in arm_unhealthy once G2 freezes it. Never gate on this.)
        sat = [bool(np.isfinite(r) and r > 0.95 * r_max_fit)
               for r in rates]
        dead = [bool(not np.isfinite(r)) for r in rates]
        # frequency-offset excursion, INFORMATIONAL until G2 freezes
        # the margin (the +- channels are not eigenchannels of the
        # detuned block, so |dw| up to ~0.03 is DESIGNED-IN dressing;
        # 0.06 is only the coarse-grid edge, the true reachable rail
        # is +-0.10)
        dw_edge = [bool(np.isfinite(d) and abs(d) >= 0.06 - 1e-9)
                   for d in dws]
        dw_rail = [bool(np.isfinite(d) and abs(d) >= 0.10 - 1e-9)
                   for d in dws]
        pi3 = float(np.prod(np.asarray(rates) - np.mean(rates))) \
            if all(np.isfinite(rates)) else float("nan")
        stats[arm] = {"rates": rates, "dws": dws, "fit_sse": sses,
                      "s2": s2_of(rates),
                      "rbar": float(np.mean(rates)),
                      "pi3_centred": pi3,
                      "fit_saturated": sat, "fit_dead": dead,
                      "dw_edge": dw_edge, "dw_rail": dw_rail}
        flag = ""
        if any(sat):
            flag += "  [rate near r_max ceiling]"
        if any(dead):
            flag += "  [FIT-HEALTH: dead channel]"
        if any(dw_rail):
            flag += "  [delta-omega at TRUE rail +-0.10]"
        elif any(dw_edge):
            flag += "  [delta-omega past coarse grid edge]"
        print(f"  {arm:2}: rates = "
              + ", ".join(f"{r:.4f}" for r in rates)
              + "  dw = " + ", ".join(f"{d:+.4f}" for d in dws)
              + f"  s2 = {stats[arm]['s2']:.5f}"
                f"  rbar = {stats[arm]['rbar']:.4f}"
                f"  Pi3 = {pi3:+.2e}{flag}")
    print(f"  inversion negative quasi-count mass, full 64-dim vector "
          f"(only the 6 sector entries are clipped; pooled, all cells): "
          f"{sum(clip_log):.1f}")

    d_sign = stats["C"]["s2"] - stats["U"]["s2"]
    d_W = stats["C"]["s2"] - stats["Cp"]["s2"]
    d_rbar_C = stats["C"]["rbar"] - stats["U"]["rbar"]
    d_rbar_Cp = stats["Cp"]["rbar"] - stats["U"]["rbar"]
    print(f"\n  d      = s2(C) - s2(U)  = {d_sign:+.5f}")
    print(f"  d_W    = s2(C) - s2(Cp) = {d_W:+.5f}")
    print(f"  rbar(C) - rbar(U)  = {d_rbar_C:+.5f}")
    print(f"  rbar(Cp) - rbar(U) = {d_rbar_Cp:+.5f}")

    # bootstrap over bindings, jointly across depths and preps WITHIN
    # each arm (section 6.6): each arm draws its own index vector; N0
    # is swept and participates; depth-0 cells are single-PUB and
    # contribute no resampling variance (stated, not hidden)
    rngb = np.random.default_rng(SEED_BOOT)
    boot = {k: [] for k in ("d", "dW", "s2C", "s2Cp", "s2U", "s2N0",
                            "drbarC", "drbarCp")}
    print(f"\n  bootstrap over bindings (R = {r_boot}, per-arm draws; "
          f"depth-0 cells carry no binding axis) ...")
    for _ in range(r_boot):
        s2b, rbarb = {}, {}
        for arm in ("N0", "U", "Cp", "C"):
            idx = rngb.integers(0, m_bind, size=m_bind)
            Yp, Wp = arm_traces(arm, idx=idx)
            rr, _, _ = arm_rates(Yp, Wp, A, dyf, ts, arm, r_max_fit)
            s2b[arm] = s2_of(rr); rbarb[arm] = float(np.mean(rr))
        boot["d"].append(s2b["C"] - s2b["U"])
        boot["dW"].append(s2b["C"] - s2b["Cp"])
        boot["s2C"].append(s2b["C"]); boot["s2Cp"].append(s2b["Cp"])
        boot["s2U"].append(s2b["U"]); boot["s2N0"].append(s2b["N0"])
        boot["drbarC"].append(rbarb["C"] - rbarb["U"])
        boot["drbarCp"].append(rbarb["Cp"] - rbarb["U"])
    # NaN-safe: one dead channel in one replicate must not turn every
    # CI into NaN; the NaN replicate count is itself a health number
    # Round 11 (MAJOR): the dropped-replicate count was computed from
    # boot["d"] alone while nanpercentile silently drops NaN
    # replicates from EVERY key. W's own count and A's C'-line count
    # were therefore invisible, and section 7 registers this number
    # precisely because silent exclusion narrows a CI in the one
    # direction that turns INCONCLUSIVE into FALSIFIED. Count per key.
    n_nan_by_key = {k: int(sum(1 for v in vs if not np.isfinite(v)))
                    for k, vs in boot.items()}
    n_nan = n_nan_by_key["d"]          # D-sign's own count, as before
    if boot_floor_trips:
        print(f"    WARNING: {len(boot_floor_trips)} replicate-level "
              f"kept-count floor trips (resampled pools under the "
              f"floor {kept_floor:g}); the affected grid points drop "
              f"from those replicates' fits, which narrows the CI")
    for k, cnt in n_nan_by_key.items():
        if cnt:
            print(f"    WARNING: {cnt}/{r_boot} bootstrap replicates "
                  f"of {k} carried a dead fit (excluded from that CI; "
                  f"exclusion narrows the CI toward FALSIFIED)")
    ci = {k: (float(np.nanpercentile(v, 2.5)),
              float(np.nanpercentile(v, 97.5)),
              float(np.nanstd(v, ddof=1))) for k, v in boot.items()}
    for k in ("d", "dW", "s2C", "s2Cp", "s2U", "s2N0", "drbarC",
              "drbarCp"):
        lo, hi, se = ci[k]
        print(f"    {k:7}: 95% CI [{lo:+.5f}, {hi:+.5f}]  SE {se:.5f}")

    # in-situ T1/T2* PUBs: raw survival report (the G4 CLEAN thresholds
    # are pending; the RECORD carries the numbers either way)
    t1t2 = art["raw_counts"].get("t1t2", [])
    if t1t2:
        print("\n  in-situ T1/T2* stations (per-qubit one-fraction):")
        for rec in t1t2:
            v = counts_to_vec(rec["counts"])
            tot = v.sum()
            ones = [float(v[(np.arange(D6) >> q) & 1 == 1].sum() / tot)
                    for q in range(N)]
            print(f"    {rec['kind']:13} delay {rec['delay_us']:7.1f} us "
                  f"station {rec.get('station', '?')}: "
                  + " ".join(f"{o:.3f}" for o in ones))
    else:
        print("\n  (no in-situ T1/T2* PUBs in this artifact)")

    # verdict-rule evaluation against the frozen constants
    print("\n  VERDICT-RULE EVALUATION (frozen constants only; every "
          "pending constant is named):")
    # round 12: a bare index would KeyError on a constants file that
    # lost the entry, i.e. analyze dying on paid data over a missing
    # THRESHOLD. A missing entry reads as unfrozen, which is what the
    # manifest audit above already reports it as.
    thD = consts.get("theta_D", {"value": None, "frozen": False})

    def arm_unhealthy(arms):
        """Dead channels always disqualify; dw/saturation disqualify
        ONLY once their G2 margins are frozen (until then they are
        designed-in dressing readouts, not registered VOID triggers;
        the round-3 reviews measured that the raw 0.06 edge fires on
        every artifact including the ideal construction)."""
        bad = any(any(stats[a]["fit_dead"]) for a in arms)
        cdw = consts.get("dw_excursion_margin", {})
        if cdw.get("frozen") and cdw.get("value") is not None:
            bad |= any(any(np.isfinite(d) and abs(d) >= cdw["value"]
                           for d in stats[a]["dws"]) for a in arms)
        csat = consts.get("r_saturation_frac", {})
        if csat.get("frozen") and csat.get("value") is not None:
            # Round 11 (MAJOR): the entry's own note says the frozen
            # value states WHICH BASE the fraction multiplies (the
            # search's true ceiling after refinement is 1.1*r_max, not
            # r_max), but the schema had no field for it and the code
            # hardwired r_max. If G2 froze the fraction against
            # 1.1*r_max, the runner applied a threshold 10% tighter
            # than registered, on a VOID trigger. The base is now read.
            sat_base = str(csat.get("base", "r_max"))
            base_val = SAT_BASES.get(sat_base)
            if base_val is None:
                # analyze never dies on paid data (the registered rule
                # for its audits): warn loudly and drop the term. The
                # submit gate hard-aborts on the same string, which is
                # where an unreadable base belongs, before the spend.
                print(f"    ** CONSTANTS AUDIT: r_saturation_frac base "
                      f"{sat_base!r} unknown (expected one of "
                      f"{sorted(SAT_BASES)}); the saturation term is "
                      f"NOT evaluated in this run.")
            else:
                bad |= any(any(np.isfinite(r) and r >= csat["value"]
                               * base_val(r_max_fit)
                               for r in stats[a]["rates"])
                           for a in arms)
        cres = consts.get("fit_residual_bound", {})
        if cres.get("frozen") and cres.get("value") is not None:
            bad |= any(any(np.isfinite(s) and s >= cres["value"]
                           for s in stats[a]["fit_sse"])
                       for a in arms)
        return bad

    grid_ok = not missing and not floored
    # one verdict label for a brake-truncated artifact (round 8: the
    # hw status says BRAKE-TRUNCATED, the analyze verdict says grid
    # incomplete; the RECORD needs them as one line)
    brake_note = ("; brake-truncated per Amendment 1.6"
                  if str(art.get("status", "")).startswith("BRAKE")
                  else "")
    if floored:
        print(f"  WARNING: {len(floored)} flown cells under the "
              f"kept-count floor (grid-incomplete, Amendment 1.2): "
              f"{floored[:6]} ...")
    if not grid_ok:
        print(f"    D-sign: VOID (grid incomplete: {len(missing)} "
              f"science cells missing or binding-broken, "
              f"{len(floored)} flown but kept-count-floored"
              f"{brake_note}; theta_D "
              f"is frozen ON the 21x3 grid and no point is dropped "
              f"after data exists)")
    elif not np.isfinite(d_sign):
        # round 11, the W blocker's mirror image: a non-finite d would
        # have read as "NOT DETECTED", which is the falsification-
        # shaped half of the partition, produced by instrument failure.
        print("    D-sign: VOID (fit health: d is not finite, i.e. a "
              "dead channel in C or U)")
    elif arm_unhealthy(("C", "U")):
        # round 12: the "dead channel" half of this message is
        # unreachable, since a dead channel makes d non-finite and the
        # branch above takes it. What survives here is the frozen-
        # margin half.
        print("    D-sign: FIT-HEALTH VOID CANDIDATE (a frozen G2 "
              "margin exceeded in C or U; dead channels take the VOID "
              "branch above)")
    elif thD.get("frozen"):
        ok = d_sign > thD["value"]
        # round 14: the partition's labels split in round 13 and this
        # string still named one INCONCLUSIVE
        verdict_word = ("PASS" if ok else
                        "NOT DETECTED (the FALSIFIED / Anti-D / "
                        "INCONCLUSIVE (underpowered) / INCONCLUSIVE "
                        "(indeterminate) partition awaits theta_F "
                        "and kappa)")
        print(f"    D-sign: d = {d_sign:+.5f} (bootstrap SE "
              f"{ci['d'][2]:.5f}) vs theta_D = {thD['value']:.5f} -> "
              f"{verdict_word} "
              f"[kappa power condition PENDING"
              + ("; theta_D REFREEZE REQUIRED pre-data, this line is "
                 "provisional" if thD.get("refreeze_required") else "")
              + "]")
    else:
        print("    D-sign: theta_D PENDING")
    for name, val, key in (("W", d_W, "theta_W"),):
        c = consts.get(key, {"value": None, "frozen": False})
        if not grid_ok:
            # the same branch D-sign takes: an incomplete grid VOIDS
            # W outright (Amendment 1.2), never a PASS/FAIL with a
            # suffix
            print(f"    {name}: VOID (grid incomplete: {len(missing)} "
                  f"missing or binding-broken, {len(floored)} "
                  f"kept-count-floored{brake_note})")
            continue
        # round 11 (BLOCKER): W had the fit-health state as a SUFFIX
        # while D-sign had it as a BRANCH, so a dead channel in C or
        # C' left d_W = nan, and `nan > theta_W` is False: the line
        # printed FAIL. A falsification-shaped output produced by
        # instrument failure is exactly what section 3 forbids, and
        # section 7 registers the opposite (a NaN verdict statistic on
        # a complete grid IS fit-health VOID). W now takes D-sign's
        # branch structure. The non-finite guard is separate from
        # arm_unhealthy on purpose: arm_unhealthy gates the dw and
        # saturation margins behind their G2 freeze (round 3), while a
        # statistic that is not a number is never comparable to a
        # threshold, frozen or not.
        if not np.isfinite(val):
            print(f"    {name}: VOID (fit health: the {name} statistic "
                  f"is not finite, i.e. a dead channel in C or C')")
            continue
        if arm_unhealthy(("C", "Cp")):
            print(f"    {name}: FIT-HEALTH VOID CANDIDATE (a frozen "
                  f"G2 margin exceeded in C or C'; dead channels take "
                  f"the VOID branch above); measured {val:+.5f}")
        elif c["frozen"]:
            print(f"    {name}: {val:+.5f} vs {key} = {c['value']} -> "
                  f"{'PASS' if val > c['value'] else 'FAIL'}")
        else:
            print(f"    {name}: measured {val:+.5f}; {key} PENDING "
                  f"(gate work outstanding)")
    # round 11 (MAJOR): this was a hardcoded subset of ten keys, so
    # two thirds of the manifest (the s2(C) floor, the time
    # axis, R_MAX_FIT/R_BOOT, the three T1-CLEAN faces, the dressing
    # functions, the layer-fidelity bound, the |2>-leakage price, ...)
    # never appeared in the analysis printout at any freeze state. The
    # executor reading the paid-data output could not see which of the
    # constants governing the flight were frozen. The manifest is the
    # registered list; iterate IT, and report all three defect states
    # the submit gate knows (missing / unfrozen / frozen-but-null).
    # round 12: the submit gate knows FOUR defect states and this
    # report knew three, so a refreeze_required entry counted as
    # frozen here and as blocking there (theta_D today: the printed
    # pending count came out one short of the gate's).
    # round 19: the submit gate knows FIVE defect states since the
    # value-form check landed and this report knew four, the same
    # drift round 12 repaired when it knew three
    unfrozen, missing_keys, null_frozen, stale, badform = [], [], [], [], []
    for key in CONSTANTS_MANIFEST:
        if key not in consts:
            missing_keys.append(key)
        elif not consts[key].get("frozen"):
            unfrozen.append(key)
        elif consts[key].get("value") is None:
            null_frozen.append(key)
        elif not _value_admissible(key, consts[key]["value"]):
            badform.append(key)
        elif consts[key].get("refreeze_required"):
            stale.append(key)
    for label, keys in (("PENDING (not frozen)", unfrozen),
                        ("MISSING from the constants file", missing_keys),
                        ("frozen with a NULL value", null_frozen),
                        ("frozen with a value of the WRONG FORM",
                         badform),
                        ("frozen but REFREEZE REQUIRED (blocks "
                         "submission like an unfrozen one)", stale)):
        if keys:
            print(f"    {len(keys)}/{len(CONSTANTS_MANIFEST)} manifest "
                  f"constants {label}: {', '.join(keys)}")
    # round 12: printed conditionally, so in the DEFAULT state (code
    # "nominal", entry null) neither the axis nor its audit appeared
    # anywhere the reader could see, while section 5 claims analyze
    # reports it. The axis the fits ran on is never a detail.
    print(f"    time axis in force: {TIME_AXIS} "
          f"(constants entry: "
          f"{consts.get('time_axis', {}).get('value')})")
    # The submit gate HARD-ABORTS on this pair; analyze only warns and
    # records, because analyze runs on paid data and the registered
    # rule for its audits (basis sha, table shas) is warn-never-abort.
    axis_ok = consts.get("time_axis", {}).get("value") in (None, TIME_AXIS)
    if not axis_ok:
        print(f"    ** AXIS AUDIT: constants time_axis "
              f"{consts['time_axis']['value']!r} != code {TIME_AXIS!r}; "
              f"the fits above ran on {TIME_AXIS}. Every rate and every "
              f"s2 on this page is off-axis: re-run analyze with the "
              f"code axis set to the frozen entry before reading any "
              f"verdict line.")
    print("    B-hat (N0 3x3 generator fit) + difference-nulls + "
          "N0-CLEAN + T1-CLEAN thresholds + the s2(C) floor + the "
          "band-validity window + the G5 in-analysis J estimator "
          "(a VOID trigger, the G5 J pass band; the fitted dw above are its raw "
          "input): PENDING machinery (G1/G2/G4/G5); the raw inputs "
          "for all of them are in this artifact")

    out = {"mode": "analyze", "artifact": str(path),
           "timestamp": datetime.now().isoformat(),
           "artifact_synthetic": bool(art.get("synthetic", False)),
           "constants_sha": hashlib.sha256(
               json.dumps(consts, sort_keys=True).encode()
           ).hexdigest()[:16],
           "grid_complete": grid_ok,
           "time_axis": TIME_AXIS,
           "axis_audit": ("ok" if axis_ok else "MISMATCH"),
           "manifest_unfrozen": unfrozen,
           "manifest_missing": missing_keys,
           "manifest_frozen_null": null_frozen,
           "manifest_bad_form": badform,
           "seed_boot": SEED_BOOT,
           "table_sha_mismatches": len(sha_bad),
           "basis_sha_audit": ("ok" if art_bsha == bsha
                               else f"{art_bsha}!={bsha}"),
           "constants": consts, "stats": stats,
           "d_sign": d_sign, "d_W": d_W,
           "d_rbar_C": d_rbar_C, "d_rbar_Cp": d_rbar_Cp,
           "boot_ci": ci, "clipped_mass": float(sum(clip_log)),
           "missing_cells": missing,
           "floored_cells": floored,
           "cal_pairs_drift": cal_pairs_drift,
           "boot_nan_replicates": n_nan,
           "boot_nan_replicates_by_key": n_nan_by_key,
           "kept_count_floor_in_force": kept_floor,
           "boot_floor_trips": len(boot_floor_trips),
           # round 14: d is a BETWEEN-ARM difference, so asymmetric
           # trips between C and U bias it and a scalar total cannot
           # show that. Per-arm, as the NaN counts are per key.
           "boot_floor_trips_by_arm": {
               a: sum(1 for t in boot_floor_trips if t[0] == a)
               for a in ARM_NAMES},
           "boot_replicates": r_boot,
           "confusion": eps_ro, "R_MAX_FIT": r_max_fit}
    _save_json(out, "analysis")


# ----------------------------------------------------------- calibrate
def _connect(name=None, fractional=True):
    from qiskit_ibm_runtime import QiskitRuntimeService
    if not IBM_TOKEN:
        print("ERROR: set IBM_QUANTUM_TOKEN in .env")
        sys.exit(1)
    try:
        svc = QiskitRuntimeService(channel="ibm_quantum_platform",
                                   token=IBM_TOKEN)
    except Exception:
        QiskitRuntimeService.save_account(channel="ibm_quantum_platform",
                                          token=IBM_TOKEN, overwrite=True)
        svc = QiskitRuntimeService(channel="ibm_quantum_platform")
    if name is None:
        return svc
    return svc, svc.backend(name, use_fractional_gates=fractional)


def _edge_rzz_error(target, a, b):
    try:
        props = target["rzz"].get((a, b)) or target["rzz"].get((b, a))
        return getattr(props, "error", None)
    except Exception:
        return None


def cmd_calibrate():
    print("CALIBRATE: line selection + BOTH-TWINS fractional check\n")
    svc = _connect()
    twin_ok = {}
    for name in TWIN_NAMES:
        try:
            bk = svc.backend(name, use_fractional_gates=True)
            has = "rzz" in bk.target.operation_names
        except Exception as e:
            has = False
            print(f"  {name}: backend fetch failed ({e})")
        twin_ok[name] = has
        print(f"  {name}: fractional rzz exposed = {has}")
    if not twin_ok.get(BACKEND_NAME):
        print(f"\n  ** {BACKEND_NAME} does not expose fractional rzz: "
              f"the candidate is PARKED by the Class-1 rule.")
        sys.exit(1)
    backend = svc.backend(BACKEND_NAME, use_fractional_gates=True)
    props = backend.properties()
    target = backend.target
    edges = set()
    for (a, b) in backend.coupling_map.get_edges():
        edges.add((a, b)); edges.add((b, a))
    adj = {}
    for (a, b) in edges:
        adj.setdefault(a, set()).add(b)

    def qubit_ok(q):
        try:
            t2 = props.t2(q) * 1e6
            ro = props.readout_error(q)
        except Exception:
            return None
        if t2 is None or ro is None:
            return None
        return (t2, ro)

    lines = []
    def dfs(path):
        if len(path) == N:
            lines.append(list(path)); return
        for nxt in sorted(adj.get(path[-1], ())):
            if nxt not in path:
                path.append(nxt); dfs(path); path.pop()
    for q0 in sorted(adj):
        dfs([q0])
    print(f"  {len(lines)} candidate 6-qubit lines")
    scored = []
    for line in lines:
        if line[0] > line[-1]:
            continue                     # dedupe reversals; orientation rule
        vals = [qubit_ok(q) for q in line]
        if any(v is None for v in vals):
            continue
        t2s = [v[0] for v in vals]; ros = [v[1] for v in vals]
        if min(t2s) < RULE_MIN_T2_US or max(t2s) / min(t2s) > RULE_MAX_T2_RATIO:
            continue
        if max(ros) > RULE_MAX_READOUT:
            continue
        p2s = []
        for a, b in zip(line, line[1:]):
            e = _edge_rzz_error(target, a, b)
            if e is None:
                break
            p2s.append(e)
        if len(p2s) < N - 1 or max(p2s) > RULE_MAX_P2_RZZ:
            continue
        # SCORE = the deep-end error budget, not the unweighted mean it
        # used to be. The mean treats all five edges alike while an odd
        # bond carries twice the gates of an even one, so it can rank a
        # worse line first: two admissible lines exist whose mean order
        # and whose deep-end kept counts disagree (pinned in the test
        # suite). The budget IS the exponent, so ranking by it ranks by
        # the quantity section 8a's margin is about. Readout keeps its
        # place as the secondary term, on the same footing as before.
        budget = line_error_budget(dict(zip(CHAIN_BONDS, p2s)))
        score = float(budget / 970.0 + np.mean(ros))
        scored.append((score, line, p2s, t2s, ros))
    if not scored:
        print("\n  ** NO line passes the rule today: do not fly.")
        sys.exit(1)
    scored.sort()
    score, chain, p2s, t2s, ros = scored[0]
    print(f"\n  chosen line (deterministic low-index orientation): {chain}")
    print(f"    rzz p2 per edge: {['%.4f' % p for p in p2s]}")
    print(f"    T2 (us): {['%.0f' % t for t in t2s]}, readout: "
          f"{['%.3f' % r for r in ros]}")
    # The margin section 8a requires, PROJECTED FROM THIS LINE'S OWN
    # EDGES instead of assumed at the guard's worst case. Reported, not
    # gated: the per-edge ceiling already implies a passing budget (see
    # P2_BUDGET_MAX), and the >= 3x evaluator itself is G3's, still on
    # the machinery ledger. What this gives the executor is the number
    # the requirement is stated in, for the line actually chosen.
    edges = dict(zip(CHAIN_BONDS, p2s))
    budget = line_error_budget(edges)
    kept = projected_deep_end_kept(edges)
    margin = kept / KEPT_COUNT_FLOOR
    print(f"    deep-end error budget {budget:.3f} "
          f"(gate-weighted mean p2 {budget / 970:.5f}; "
          f"a uniform line at the guard would read {970 * RULE_MAX_P2_RZZ:.3f})")
    print(f"    projected deepest-cell kept counts {kept:.0f} of "
          f"{SHOTS_SCIENCE} at f_leak {F_LEAK_WORST:g} = {margin:.2f}x the "
          f"floor {KEPT_COUNT_FLOOR:g} (section 8a asks "
          f">= {DEEP_END_MARGIN:g}x; ISOLATED-gate p2, the layered "
          f"question is the day-of layer-fidelity gate's)")
    if margin < DEEP_END_MARGIN:
        print(f"    ** NOTE: below the registered {DEEP_END_MARGIN:g}x on "
              f"this line's own edges. Not an abort (the >= 3x evaluator "
              f"is G3's, machinery ledger), but the day-of addendum must "
              f"carry it and the executor should expect a deep-end VOID "
              f"risk from grid-incompleteness.")
    info = {"backend": BACKEND_NAME, "chain": chain,
            "timestamp": datetime.now().isoformat(),
            "twin_fractional": twin_ok, "score": score,
            "rzz_p2": p2s, "t2_us": t2s, "readout": ros,
            "deep_end": {"error_budget": budget,
                         "edge_gate_weights":
                             [edge_gate_weights()[b] for b in CHAIN_BONDS],
                         "f_leak_worst": F_LEAK_WORST,
                         "projected_kept": kept,
                         "margin_over_floor": margin,
                         "required_margin": DEEP_END_MARGIN},
            "rule": {"min_t2_us": RULE_MIN_T2_US,
                     "max_ratio": RULE_MAX_T2_RATIO,
                     "max_readout": RULE_MAX_READOUT,
                     "max_p2_rzz": RULE_MAX_P2_RZZ}}
    path = SCRIPT_DIR / ("corner_beat_chain_"
                         + datetime.now().strftime("%Y%m%d") + ".json")
    with open(path, "w") as f:
        json.dump(info, f, indent=2)
    print(f"  saved -> {path}")


def load_chain_info():
    cands = sorted(SCRIPT_DIR.glob("corner_beat_chain_*.json"))
    if not cands:
        print("ERROR: no corner_beat_chain_*.json; run --calibrate first")
        sys.exit(1)
    return json.load(open(cands[-1]))


# ------------------------------------------------------------ hardware
BOUND_CIRC_MAX = 80_000  # bound circuits per job (the concentrator
                         # precedent flew 6144 in one job: 24 sink
                         # PUBs x 256 bindings, no-sink PUBs unbound)
PAYLOAD_MAX_MB = 48.0   # per-job parameter payload cap; ~11 jobs at the
                        # frozen grid (372 MB in one job is untested
                        # scale; measured plan: 11 jobs)
STATION_AFTER_TP = [3, 9, 15, 20]    # aux stations after these GRID
                                     # INDICES = steps 9/27/45/60
QUEUE_WARN = 50


# the COMPLETE section-8 freeze manifest: a constants file that lacks a
# key is not "all frozen", it is incomplete (the gate must not read
# green just because the unfinished quantities have no entry yet)
CONSTANTS_MANIFEST = [
    "theta_D", "theta_W", "theta_F", "dmag_band_C", "dmag_band_Cp",
    # round 18: A has TWO independently decisive lines (the
    # partition reads its C line alone; the (2,3) null voids only its
    # C-prime line), so one scalar could not carry both, the same
    # ruling round 6 made for T1-CLEAN's three faces
    "A_margin_C", "A_margin_Cp", "kappa",
    "null_margin_12", "null_margin_23",
    "dw_excursion_margin", "r_saturation_frac", "dose_criterion",
    "billing_cap_min", "R_MAX_FIT", "R_BOOT",
    # section-8 quantities still without machinery (entries must be
    # created by the gate work before they can freeze):
    "s2C_floor", "band_validity_window", "kappa_projected_SE",
    "family_error_rates", "N0_clean_thresholds",
    # T1-CLEAN's three registered faces (doc section 9, round 6:
    # one scalar cannot carry three bands)
    "T1_clean_projection_band", "T1_clean_scale",
    "T1_clean_station_band",
    "f_C_dressing", "f_Cp_dressing",
    "A_centres_dressed", "pi3_centre_dressed",
    "pi3_report_band", "J_pass_band",
    "dayof_width_scaling", "fit_residual_bound", "time_axis",
    # round 12: two thresholds that VOID a verdict and were code
    # constants. kept_count_floor is the sole gate on the
    # grid-incompleteness trigger (a cell under it is floored, which
    # voids D-sign and W), so Amendment 1.5's "no verdict consumes
    # them" was false for it; boot_nan_replicate_bound is the bound
    # section 7 says freezes with the fit-health margins, which none
    # of the three fit-health keys carries.
    "kept_count_floor", "deep_end_kept_profile",
    "boot_nan_replicate_bound",
    # round 10: the day-of layer-fidelity gate's frozen BOUND (the
    # machinery had a ledger entry, the threshold had no key) and
    # the |2>-leakage pricing (section 11 calls it the largest
    # unpriced systematic; it could never block without a key)
    "layer_fidelity_bound", "leakage_ket2_price",
    # Amendment 1.8: the Class-1 p2 guard bound itself. It hard-aborts
    # a paid submission on both the calibrate and the day-of path and
    # was a bare module constant, outside the manifest and outside the
    # sync-check loop, so an edit to it was invisible to every gate.
    "p2_guard_bound",
]

# the fit time axis the analysis runs on. The pre-registration
# (section 5, round-3 repair) registers the REALIZED-DOSE axis
# t_eff = max(n-1, 0)*dt for the refreeze; this constant is
# sync-checked against the frozen "time_axis" manifest entry at
# submission, so freezing that entry to "realized_dose" FORCES the
# code change HERE before any flight. Not in the committed gate,
# though (round-12 correction): certify hands ONE axis array to both
# fit_rate implementations, so the gate's own axis is never compared.
# The gate side rides the machinery entry realized_dose_time_axis_refit.
TIME_AXIS = "nominal"


# the bases r_saturation_frac's fraction may multiply (round 11: the
# entry's note named them in prose and the schema had no field). Each
# is a function of the analyze-side r_max in force.
SAT_BASES = {
    "r_max": lambda r_max: r_max,
    "1.1*r_max": lambda r_max: 1.1 * r_max,
}


def fit_time_axis():
    """The fit time axis, DERIVED from TIME_AXIS.

    Round 11 (MAJOR): TIME_AXIS was sync-checked against the frozen
    constants entry but never read by the estimator, which built its
    axis as n*dt inline. Freezing the entry to "realized_dose" would
    have forced an edit to this one string and nothing else, leaving
    the fit on the nominal axis while the manifest recorded the
    repair as landed. The section-5 off-by-one it hides is worth
    ~16.6% of s2 and biases the rates ~-5%, so the freeze must move
    the axis itself.

    "nominal"       t = n*dt              (today's committed behaviour)
    "realized_dose" t = max(n-1, 0)*dt    (section 5, round-3 repair;
                                           freezes with the theta_D
                                           refreeze, gate and runner
                                           together)
    """
    if TIME_AXIS == "nominal":
        return np.array([s * JDT for s in GRID_STEPS])
    if TIME_AXIS == "realized_dose":
        return np.array([max(s - 1, 0) * JDT for s in GRID_STEPS])
    raise RuntimeError(f"unknown TIME_AXIS: {TIME_AXIS!r}")

# The manifest keys whose CONSUMER ALREADY EXISTS in this file, so that
# what they still owe is a NUMBER and not code. Declared explicitly
# because section 8a's promise ("every machinery entry names the keys
# it discharges, so that 'empty the ledger by its names' is an
# operation an executor can actually perform") only holds if every key
# is reachable one way or the other. A key named by no entry and with
# no live consumer is a number nobody has to build anything for, which
# is exactly how a freeze releases a flight with a dead arm. Each is
# consumed where named:
#   billing_cap_min      -> the pre-submit projection and the brake
#   p2_guard_bound       -> the Class-1 line rule, calibrate + day-of
#   theta_D              -> the D-sign comparison in cmd_analyze
#   theta_W              -> the W comparison in cmd_analyze
#   R_MAX_FIT, R_BOOT    -> read by cmd_analyze before the fits
#   dose_criterion       -> the dose certificates in certify / dry run
#   dw_excursion_margin  -> fit health, arm_unhealthy
#   r_saturation_frac    -> fit health, arm_unhealthy (+ its base gate)
#   fit_residual_bound   -> fit health, arm_unhealthy
#   kept_count_floor     -> the floor in pool_invert_postselect
#   time_axis            -> fit_time_axis(), sync-checked at submit
# Several ALSO appear in the ledger, legitimately: kept_count_floor has
# a live floor comparison AND an owed >= 3x margin evaluator, and the
# time axis is read here while the gate-side refit is still owed.
KEYS_WITH_LIVE_CONSUMERS = frozenset({
    "billing_cap_min", "p2_guard_bound", "theta_D", "theta_W",
    "R_MAX_FIT", "R_BOOT", "dose_criterion", "dw_excursion_margin",
    "r_saturation_frac", "fit_residual_bound", "kept_count_floor",
    "time_axis",
})

# verdict machinery still to be BUILT (distinct from unfrozen
# numbers): submission hard-aborts while this list is non-empty;
# the gate work removes entries as the consumers land.
MACHINERY_PENDING = [
    "B_hat_3x3_generator_fit (G2)",
    "difference_null_evaluators_12_23 (G2) [null_margin_12, null_margin_23]",
    "N0_clean_evaluator (G2) [N0_clean_thresholds]",
    "T1_clean_three_faces_evaluator (G4) [T1_clean_projection_band, T1_clean_scale, T1_clean_station_band]",
    "s2C_floor_adjudication (G1) [s2C_floor]",
    "band_validity_window_check (G1) [band_validity_window]",
    "f_C_f_Cp_dressing_functions (G1/G2) [f_C_dressing, f_Cp_dressing]",
    "A_margin_C_and_Cp_tost + pi3_report_band (G1) [A_margin_C, A_margin_Cp, pi3_report_band]",
    "J_pass_band_evaluator (G5) [J_pass_band]",
    "dayof_width_scaling_consumer (G1) [dayof_width_scaling]",
    "realized_dose_time_axis_refit: SWITCH + certify axis parity BUILT "
    "2026-08-17 (gate and runner derive the axis from the same "
    "definition; certify check 8 compares both settings exactly). "
    "WHAT REMAINS is EXECUTION: re-fit at >= 500 reps with the "
    "switch thrown (theta_D refreeze) [time_axis]",
    "layer_fidelity_dayof_gate (section 9, round 3) [layer_fidelity_bound]",
    "depth0_shot_resample_or_CI_inflation (theta_D refreeze, round 6)",
    "dose_certificates_on_flown_tables_in_analyze (G3, round 8) [dose_criterion]",
    "s2C_floor_forces_notW_and_notDsign_evaluator (G1, rounds 8-9) [s2C_floor]",
    # round 11: the ledger held every owed NUMBER's consumer but not
    # the two pieces of machinery that turn those consumers into the
    # RECORD's one line. Emptying the list by name would have released
    # a flight whose falsification arm has no evaluator and whose
    # verdict has no composer.
    # The composer itself is BUILT and tested: corner_beat_verdict.py
    # beside this file, with test_corner_beat_verdict.py (21 tests,
    # including an exhaustive walk of all 972 arm/line states across
    # both D-mag regimes, asserting that every state yields exactly
    # one line and that no line leaves section 7's five-form
    # vocabulary). What remains is the WIRING, which cannot land
    # before the arm evaluators below produce TRUE/FALSE/VOID per
    # line: a composer with no inputs would be a frozen number
    # without its evaluator, one storey up.
    "verdict_line_composer_WIRED_into_analyze (needs the arm evaluators)",
    "notDsign_partition_evaluator_thetaF_antiD_kappa (section 7/G1) [theta_F, kappa, kappa_projected_SE]",
    # the entry ABOVE also owns kappa_projected_SE, the manifest key
    # whose only consumer is the power condition inside the partition
    # (round 16; comment moved in round 17, where it sat above the
    # composer entry and read as if the composer discharged it)
    # round 11 (writer, while confirming Open question 2): the parking
    # rule's third condition is the only one with no machine hook.
    # Conditions 1-2 (fractional exposure on the flown backend AND the
    # twin, day-of p2 <= RULE_MAX_P2_RZZ on the used edges) hard-abort
    # in cmd_hardware; "the committed gate fails to reproduce section
    # 8a" is today a human act, since --certify checks PARITY against
    # the gate, never that the gate reproduces the governing numbers.
    "committed_gate_reproduces_section_8a_dayof (G1 table mode)",
    # round 12: the deep-end margin is REGISTERED (>= 3x the kept-count
    # floor at the worst modeled corner). It FAILED at the old 0.6%
    # guard bound (1.74x) and Amendment 1.8 met it at 0.5% (4.17x), but
    # the EVALUATOR is still missing, and it is what would catch the
    # layered case, where the same arithmetic gives 1.12x and then a
    # floored deep end.
    "deep_end_kept_margin_3x_evaluator (G3, round 5) [deep_end_kept_profile, kept_count_floor]",
    "boot_nan_replicate_bound_evaluator (G2 fit health) [boot_nan_replicate_bound]",
    # round 14: three manifest keys had NO reader in code and no entry
    # here, so freezing them would have released submission silently.
    # The D-mag bands are a named conjunct of CONFIRMED and the
    # dressing entry covers the CENTRES, not the containment test;
    # the |2>-leakage price is what section 11 calls the largest
    # unpriced systematic, and round 10 added its key so it "could
    # now block", which it could not without a consumer.
    "dmag_band_containment_evaluator_C_and_Cp (G1) [dmag_band_C, dmag_band_Cp]",
    "leakage_ket2_price_consumer (G1 counts model) [leakage_ket2_price]",
    # round 15: the same sweep one key further. family_error_rates had
    # no reader and no entry either, so freezing it would have
    # released submission silently; A_centres_dressed and
    # pi3_centre_dressed were covered only implicitly by the
    # A_margin_tost entry, and are named here so the map is total.
    "family_error_rates_reporter (G1 joint false-VOID accounting) [family_error_rates]",
    "A_centres_and_pi3_centres_as_dressed_functions (G1/G2) [A_centres_dressed, pi3_centre_dressed]",
]


def _require_frozen_constants():
    """A constant blocks submission if ABSENT from the file, frozen=
    false, refreeze_required=true (a stale freeze is not a freeze), or
    frozen=true with a null value (a freeze without a number).

    ONE registered non-number is admissible (round 12): the sentinel
    "REPORTED", meaning the quantity is deliberately NOT registered by
    pre-data amendment and its arm prints as a reported reading. It
    exists because the document tells the executor to EXPECT that
    regime for dmag_band_Cp (today's gate numbers admit only the
    reverted band), and there was no way to record it: a null blocked,
    a missing entry blocked, and inventing a number for a band the
    pre-registration says will not exist is the worse option. Any
    evaluator consuming such a key must print the reading, never a
    verdict."""
    consts = load_constants()
    pending = []
    for k in CONSTANTS_MANIFEST:
        v = consts.get(k)
        if not isinstance(v, dict) or "frozen" not in v:
            pending.append(f"{k} (MISSING ENTRY)")
        elif not v["frozen"] or v.get("refreeze_required"):
            pending.append(k)
        elif v.get("value") is None:
            pending.append(f"{k} (frozen but value null)")
        elif not _value_admissible(k, v["value"]):
            pending.append(f"{k} (frozen with an inadmissible value "
                           f"{v['value']!r})")
    return consts, pending


REPORTED_SENTINEL = "REPORTED"
# Round 18 (BLOCKER): the admissibility guard demanded a NUMBER on
# every manifest key, while the pre-registration registers nine of
# them as objects that cannot be scalars. Section 7 says of the
# dressed centres, in as many words, "a frozen NUMBER cannot exist
# pre-flight"; the day-of width SCALINGS are functions too; the
# family error rates are a PAIR; N0-CLEAN is THREE thresholds; the
# band-validity window has dimensions; the deep-end kept profile is
# one number per deep cell. The guard would therefore have rejected
# the correctly frozen objects at the last gate before a paid
# submission, leaving no route forward except inventing scalars,
# which is post-data freezing under another name. Each key now
# declares its FORM and the guard checks that form.
#   scalar   a real number
#   interval [lo, hi], lo <= hi
#   vector   a non-empty list of reals (a profile)
#   map      a non-empty dict of named reals
#   code     {"source": <path>, "symbol": <name>, "sha256": <hex>},
#            the committed function IS the frozen object
VALUE_FORMS = {
    "f_C_dressing": "code",
    "f_Cp_dressing": "code",
    "A_centres_dressed": "code",
    "pi3_centre_dressed": "code",
    "dayof_width_scaling": "code",
    "family_error_rates": "map",
    "N0_clean_thresholds": "map",
    # (band_validity_window is declared once, below, as
    # map_of_intervals; round 19 left a "map" entry here as well and
    # Python kept the later one, which is this arc's own signature
    # failure, a rule stated twice and repaired once, sitting inside
    # the guard that exists to prevent it)
    "deep_end_kept_profile": "vector",
    # Round 19: "scalar_or_interval" accepted two INCOMPATIBLE
    # readings of a CONFIRMED conjunct (a half-width around a frozen
    # centre versus two absolute endpoints), which is not a guard.
    # Section 7 says "inside its frozen band AROUND f_C(B-hat)", so
    # the half-width is the registered reading and these are scalars.
    "dmag_band_C": "scalar",
    "dmag_band_Cp": "scalar",
    "pi3_report_band": "scalar",
    # a pass band around the frozen J: half-width, same reading, and
    # NOT interval-only, which would have made it a tenth key unable
    # to hold a scalar and hard-aborted a perfectly good freeze
    "J_pass_band": "scalar",
    # "band" in a T1-CLEAN name is likewise a half-width
    "T1_clean_projection_band": "scalar",
    "T1_clean_station_band": "scalar",
    # the window has DIMENSIONS: named axes, each an interval
    "band_validity_window": "map_of_intervals",
    "R_BOOT": "positive_integer",
}


def _form_ok(form, v):
    num = (int, float)
    if isinstance(v, bool):
        return False
    if form == "scalar":
        return isinstance(v, num)
    if form == "interval":
        return (isinstance(v, (list, tuple)) and len(v) == 2
                and all(isinstance(x, num) and not isinstance(x, bool)
                        for x in v) and v[0] <= v[1])
    if form == "positive_integer":
        return isinstance(v, int) and v > 0
    if form == "map_of_intervals":
        return (isinstance(v, dict) and len(v) > 0
                and all(_form_ok("interval", x) for x in v.values()))
    if form == "vector":
        return (isinstance(v, (list, tuple)) and len(v) > 0
                and all(isinstance(x, num) and not isinstance(x, bool)
                        for x in v))
    if form == "map":
        return (isinstance(v, dict) and len(v) > 0
                and all(isinstance(x, num) and not isinstance(x, bool)
                        for x in v.values()))
    if form == "code":
        # Round 19: this accepted any non-empty strings, so
        # {"source":"todo","symbol":"todo","sha256":"todo"} passed the
        # last gate before a paid submission. Same hole round 14
        # closed for time_axis, one key over. The frozen object IS the
        # committed function, so the file must exist and hash to the
        # recorded digest.
        import hashlib as _h
        import re as _re
        if not (isinstance(v, dict)
                and all(isinstance(v.get(k), str) and v.get(k)
                        for k in ("source", "symbol", "sha256"))):
            return False
        if not _re.fullmatch(r"[0-9a-f]{64}", v["sha256"]):
            return False
        src = (SCRIPT_DIR / v["source"]) if not Path(v["source"]).is_absolute() \
            else Path(v["source"])
        if not src.exists():
            return False
        return _h.sha256(src.read_bytes()).hexdigest() == v["sha256"]
    raise RuntimeError(f"unknown value form: {form!r}")
# the ONLY keys whose freeze may be the sentinel instead of a number,
# and the reason each is eligible: the pre-registration tells the
# executor to EXPECT the reverted-band regime for the C' magnitude, in
# which the band does not exist and its centre function is owed only
# as the REPORTED centre.
# Round 19 (BLOCKER): round 18 added A_margin_Cp here and said so only
# in a JSON note. Section 7 registers TWO sentinel keys "and no
# others", and the note asserted something the pre-registration never
# does: that the D-mag C-prime BAND reverting to REPORTED also makes
# the A arm's C-prime line REPORTED. Those are different statistics
# (s2(C') at SNR ~ 1 versus rbar(C') - rbar(U)) and A's C-prime
# discriminability is assessed nowhere. A new pre-data registration
# cannot be made in a code comment, so the key comes back out; if the
# coupling is wanted it is a documented amendment with its own test.
SENTINEL_ELIGIBLE = ("dmag_band_Cp", "f_Cp_dressing")
# registered STRING keys, each with its closed enumeration. Round 14:
# this was a type test that accepted any non-empty string, so freezing
# time_axis to "realised_dose" or "realized-dose" passed the last
# guard before a paid submission (the sync check only compares it to
# the module constant, and the same typo would sit in both) and then
# raised inside fit_time_axis on the paid data. That is the exact
# failure mode the sentinel guard was built to close, one key over.
_STRING_KEYS = {"time_axis": ("nominal", "realized_dose")}


def _value_admissible(key, value):
    """Round 13: round 12 DOCUMENTED a single admissible non-number
    and never implemented it, so the gate accepted ANY non-null value
    on ANY key. A string in theta_D would have passed the last guard
    before a paid submission and then raised TypeError on the paid
    data at `d > thD["value"]`. Round 18: but "everything is a number"
    was equally wrong in the other direction, and blocked nine keys
    the pre-registration registers as functions, pairs, profiles or
    threshold sets. A key is admissible iff it matches ITS OWN
    registered form."""
    if key in _STRING_KEYS:
        return value in _STRING_KEYS[key]
    if value == REPORTED_SENTINEL:
        return key in SENTINEL_ELIGIBLE
    return _form_ok(VALUE_FORMS.get(key, "scalar"), value)


def _rcpsi_git_gate():
    """(ok, head, detail): the pre-registration, gate, and frozen
    records must be tracked and clean in the R=CPsi-squared repo."""
    import subprocess
    detail = []
    ok = True
    for rel in RCPSI_COMMIT_PATHS:
        tracked = subprocess.run(
            ["git", "-C", str(RCPSI_REPO), "ls-files",
             "--error-unmatch", rel],
            capture_output=True).returncode == 0
        dirty = subprocess.run(
            ["git", "-C", str(RCPSI_REPO), "status", "--porcelain",
             "--", rel], capture_output=True, text=True).stdout.strip()
        detail.append((rel, tracked, dirty == ""))
        ok &= tracked and dirty == ""
    r = subprocess.run(["git", "-C", str(RCPSI_REPO), "rev-parse",
                        "HEAD"], capture_output=True, text=True)
    head = r.stdout.strip() if r.returncode == 0 else None
    return ok, head, detail


def _git_state(paths):
    """(relpath, tracked, clean) for each path, in the pipeline repo."""
    import subprocess
    out = []
    for p in paths:
        rel = os.path.relpath(str(p), str(SCRIPT_DIR)).replace("\\", "/")
        tracked = subprocess.run(
            ["git", "ls-files", "--error-unmatch", rel], cwd=SCRIPT_DIR,
            capture_output=True).returncode == 0
        dirty = subprocess.run(
            ["git", "status", "--porcelain", "--", rel], cwd=SCRIPT_DIR,
            capture_output=True, text=True).stdout.strip()
        out.append((rel, tracked, dirty == ""))
    return out


def _git_head():
    import subprocess
    r = subprocess.run(["git", "rev-parse", "HEAD"], cwd=SCRIPT_DIR,
                       capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else None


def _align_delay_us(delay_us, backend):
    """Round a delay to the backend's timing granularity (an unaligned
    delay in an auxiliary PUB must never be able to reject the job)."""
    try:
        tgt = backend.target
        dt = float(tgt.dt)
        gran = int(getattr(tgt, "granularity", 16) or 16)
    except Exception:
        return round(delay_us, 1)
    samples = max(gran, int(round(delay_us * 1e-6 / dt / gran)) * gran)
    return samples * dt * 1e6


def cmd_hardware(yes=False):
    print("HARDWARE (Class-1 guards are hard aborts; no override flags)\n")
    if not __debug__:
        raise SystemExit("ABORT: never run --hardware under python -O "
                         "(guards must not be strippable)")
    today = datetime.now().strftime("%Y%m%d")

    # Class 1: a passing --certify record from TODAY (sector parity,
    # prep parity, trajectory parity, transpile shape, dose, fit parity)
    certs = sorted(RESULTS_DIR.glob(f"corner_beat_certify_{today}_*.json"))
    if not certs:
        print(f"  no --certify record from today ({today})")
        if yes:
            print("\n  ** HARD ABORT: run --certify first (same day).")
            sys.exit(1)
    else:
        crec = json.load(open(certs[-1]))
        if yes and (not crec.get("debug_mode", False)
                    or crec.get("basis_sha") != demod_basis_sha()):
            print("\n  ** HARD ABORT: certify record was produced under "
                  "python -O or against a different demodulation basis.")
            sys.exit(1)
        print(f"  certify record: {certs[-1].name} (debug_mode ok)")

    # Class 1: constants all frozen (submission only; dry run reports)
    consts, pending = _require_frozen_constants()
    # the code-enforced module constants must MATCH the constants file
    # (either alone can be edited; the pair may not diverge)
    # round-11 nit: this was the one Class-1 guard that surfaced as a
    # RuntimeError traceback in a DRY run instead of the printed
    # hard-abort every other guard uses. It still aborts; it now
    # reads like its siblings.
    csat_gate = consts.get("r_saturation_frac", {})
    if str(csat_gate.get("base", "r_max")) not in SAT_BASES:
        print(f"\n  ** HARD ABORT: r_saturation_frac base "
              f"{csat_gate.get('base')!r} is not one of "
              f"{sorted(SAT_BASES)}; the fit-health VOID trigger has "
              f"no defined threshold.")
        sys.exit(1)
    for key, mod_val in (("billing_cap_min", BILL_ABORT_MIN),
                         ("dose_criterion", DOSE_CRITERION),
                         ("kept_count_floor", KEPT_COUNT_FLOOR),
                         ("p2_guard_bound", RULE_MAX_P2_RZZ),
                         ("time_axis", TIME_AXIS)):
        if key in consts and consts[key].get("value") is not None:
            same = (str(consts[key]["value"]) == str(mod_val)
                    if key == "time_axis"
                    else float(consts[key]["value"]) == float(mod_val))
            if not same:
                print(f"\n  ** HARD ABORT: {key}: constants file "
                      f"{consts[key]['value']!r} != code {mod_val!r}; "
                      f"the frozen entry and the code that reads it "
                      f"have diverged.")
                sys.exit(1)
    # Class 1: the manifest freezes NUMBERS; this list freezes
    # MACHINERY (round-5 gap: writing a number into the JSON must
    # not release submission while the evaluating code is absent).
    # The gate work DELETES entries here as it lands the consumers.
    if MACHINERY_PENDING:
        print(f"  verdict machinery NOT all implemented "
              f"({len(MACHINERY_PENDING)}): {MACHINERY_PENDING}")
        if yes:
            print("\n  ** HARD ABORT: frozen numbers without their "
                  "evaluating machinery cannot govern a verdict.")
            sys.exit(1)
    if pending:
        print(f"  constants NOT all frozen ({len(pending)} pending): "
              f"{pending}")
        if yes:
            print("\n  ** HARD ABORT: submission with unfrozen verdict "
                  "constants is post-data freezing (section 13). "
                  "Finish the gate work first.")
            sys.exit(1)

    # Class 1: pre-registration commit gate (runner, constants, chain,
    # day-of addendum committed at a real hash BEFORE the Batch opens)
    chain_info = load_chain_info()
    chain = chain_info["chain"]
    if chain_info.get("backend") != BACKEND_NAME:
        print(f"\n  ** HARD ABORT: chain file is for "
              f"{chain_info.get('backend')}, flying {BACKEND_NAME}.")
        sys.exit(1)
    chain_date = chain_info.get("timestamp", "")[:10].replace("-", "")
    if chain_date != today:
        print(f"\n  ** {'HARD ABORT' if yes else 'WARNING'}: chain file "
              f"is from {chain_date}, not today ({today}); re-run "
              f"--calibrate.")
        if yes:
            sys.exit(1)
    addendum = SCRIPT_DIR / f"corner_beat_dayof_addendum_{today}.json"
    commit_paths = [Path(__file__), CONSTANTS_PATH,
                    SCRIPT_DIR / f"corner_beat_chain_{chain_date}.json",
                    addendum]
    dayof = None
    if yes:
        if not addendum.exists():
            print(f"\n  ** HARD ABORT: day-of addendum {addendum.name} "
                  f"missing (the section 9 re-gate output, committed "
                  f"BEFORE the Batch opens).")
            sys.exit(1)
        dayof = json.load(open(addendum))
        # tau_step_us: the transpiled per-Strang-step wall duration,
        # the T1-CLEAN us->gbar bridge's shared conversion (doc
        # section 9; round 6: registered but previously unenforced)
        # round 14: job_plan joins the schema, because Amendment 1.5
        # tells the executor to carry the dry run's realized plan into
        # this addendum and check it against the pre-registered range,
        # and the addendum had no key for it.
        for key in ("calibration_snapshot", "band_width_factors",
                    "tau_step_us", "job_plan"):
            if key not in dayof:
                print(f"\n  ** HARD ABORT: day-of addendum lacks "
                      f"'{key}' (an empty file is not a re-gate).")
                sys.exit(1)
        # The T1-CLEAN us->gbar bridge, COMPUTED (2026-08-18). Until now
        # tau_step_us was checked for PRESENCE and never read, so section
        # 9's conversion would have been done by hand, from prose, across
        # two repositories, on flight day. It is printed, not enforced:
        # T1-CLEAN is a registered GLOBAL VOID trigger (section 7, rank 3
        # of the void order), NOT a Class-1 abort, and a guard must not
        # convert a registered void into instrument failure.
        try:
            _t1 = t1_list_from_snapshot(dayof["calibration_snapshot"], chain)
            _ok, _prof = t1_clean_scale_face(float(dayof["tau_step_us"]), _t1)
            print(f"\n  T1-CLEAN scale face (computed, not enforced here): "
                  f"scale = {_prof['scale_gbar']:.4f} gbar against the "
                  f"registered <= {T1_CLEAN_SCALE_BOUND} "
                  f"-> {'within' if _ok else 'OUTSIDE'}")
            print("    per-station gamma_hw/gbar: "
                  + "  ".join(f"Q{q}:{g:.4f}"
                              for q, g in zip(chain, _prof["gamma_hw_gbar"])))
            if not _ok:
                print("    ** the scale face is outside its registered "
                      "bound: this is a VOID condition to carry into the "
                      "verdict, per section 7, not a reason to abort here.")
        except (KeyError, ValueError, TypeError) as exc:
            print(f"\n  ** WARNING: the T1-CLEAN bridge could not be "
                  f"computed from the addendum ({exc}); section 9's face "
                  f"must then be evaluated by hand and recorded.")
        bad = [(rel, tr, cl) for (rel, tr, cl) in _git_state(commit_paths)
               if not (tr and cl)]
        if bad:
            print("\n  ** HARD ABORT: pre-registration commit gate: the "
                  "following must be committed and unmodified at HEAD:")
            for rel, tr, cl in bad:
                print(f"       {rel}: tracked={tr} clean={cl}")
            sys.exit(1)
        rc_ok, rc_head, rc_detail = _rcpsi_git_gate()
        if not rc_ok:
            print("\n  ** HARD ABORT: R=CPsi-squared commit gate (the "
                  "pre-registration, gate, and frozen records must be "
                  "committed and clean):")
            for rel, tr, cl in rc_detail:
                print(f"       {rel}: tracked={tr} clean={cl}")
            sys.exit(1)
        # the runner's analysis chain amends section 5/6.1 (pool-then-
        # invert, clip policy, R_MAX_FIT, kept-count floor): the
        # committed document must CARRY that amendment (section 13);
        # a docstring claim is not a pre-registration
        doc_text = (RCPSI_REPO / RCPSI_COMMIT_PATHS[0]).read_text(
            encoding="utf-8")
        if not any(ln.strip().startswith(("## Amendment 1",
                                          "### Amendment 1"))
                   for ln in doc_text.splitlines()):
            print("\n  ** HARD ABORT: the pre-registration carries no "
                  "'Amendment 1' (the pool-then-invert/clip/R_MAX_FIT "
                  "amendment must be appended and committed pre-data).")
            sys.exit(1)
        print(f"  commit gate: PASS (pipeline HEAD {_git_head()}, "
              f"R=CPsi-squared HEAD {rc_head}, Amendment 1 present)")

    svc, backend = _connect(BACKEND_NAME)
    target = backend.target

    # Class 1: fractional exposed and taken, on the flown backend AND
    # the twin (a twin without it drops out of the pinned pair: record)
    _require("rzz" in target.operation_names,
             "fractional rzz not exposed: PARKED (CZ = NO-FLIGHT)")
    twin_status = {}
    for nm in TWIN_NAMES:
        try:
            bk = backend if nm == BACKEND_NAME else \
                svc.backend(nm, use_fractional_gates=True)
            twin_status[nm] = "rzz" in bk.target.operation_names
        except Exception as e:
            twin_status[nm] = f"unavailable: {e!r}"
    print(f"  twin fractional status: {twin_status}")
    if twin_status.get(BACKEND_NAME) is not True:
        print("\n  ** HARD ABORT: flown backend lost fractional rzz.")
        sys.exit(1)

    # Class 1: day-of properties snapshot + line rule + p2 guard
    props = backend.properties(refresh=True) if _has_refresh(backend) \
        else backend.properties()
    snap = {"timestamp": datetime.now().isoformat(), "qubits": {},
            "edges": {}}
    rule_ok, t2s = True, []
    for q in chain:
        try:
            # props.t2() is the properties-file T2 (echo-calibrated on
            # IBM's pipeline); the pre-registration's T2echo rule reads
            # this as its proxy, and the first aux station's Hahn echo
            # is the in-data cross-check of the same quantity
            t2 = props.t2(q) * 1e6
            ro = props.readout_error(q)
        except Exception:
            print(f"    Q{q}: properties missing"); rule_ok = False
            continue
        snap["qubits"][str(q)] = {"T1_us": props.t1(q) * 1e6,
                                  "T2_us": t2, "readout_error": ro}
        t2s.append(t2)
        if t2 < RULE_MIN_T2_US or ro > RULE_MAX_READOUT:
            rule_ok = False
    if len(t2s) == N and max(t2s) / min(t2s) > RULE_MAX_T2_RATIO:
        rule_ok = False
    p2_today = []
    for a, b in zip(chain, chain[1:]):
        e = _edge_rzz_error(target, a, b)
        snap["edges"][f"{a}-{b}"] = {"rzz_error": e}
        if e is None or e > RULE_MAX_P2_RZZ:
            rule_ok = False
        if e is not None:
            p2_today.append(e)
    print(f"  day-of line {chain}: T2 {['%.0f' % t for t in t2s]}, "
          f"rzz p2 {['%.4f' % p for p in p2_today]}")
    if not rule_ok:
        print(f"\n  ** HARD ABORT: the line violates the day-of rule "
              f"({_day_of_rule_text()}). No override flag by design.")
        sys.exit(1)
    print("  day-of line rule + p2 guard: PASS")
    if len(p2_today) == N - 1:
        edges_today = dict(zip(CHAIN_BONDS, p2_today))
        kept_today = projected_deep_end_kept(edges_today)
        print(f"  day-of deep-end projection: budget "
              f"{line_error_budget(edges_today):.3f}, kept "
              f"{kept_today:.0f} of {SHOTS_SCIENCE} = "
              f"{kept_today / KEPT_COUNT_FLOOR:.2f}x the floor "
              f"(section 8a asks >= {DEEP_END_MARGIN:g}x; reported, and "
              f"isolated-gate: the layered face is the layer-fidelity "
              f"gate's)")
        snap["deep_end_projection"] = {
            "error_budget": line_error_budget(edges_today),
            "projected_kept": kept_today,
            "margin_over_floor": kept_today / KEPT_COUNT_FLOOR}

    # build + transpile
    em, Wm, dyf, A, preps = build_preps()
    plan = build_pub_plan()
    print("  transpiling 63 skeletons ...")
    skeletons = {}
    for tp_idx in range(len(GRID_STEPS)):
        for prep_idx in range(3):
            qs, _ = science_skeleton(GRID_STEPS[tp_idx], preps[prep_idx])
            skeletons[(tp_idx, prep_idx)] = qs
    keys = sorted(skeletons)
    tcs = transpile([skeletons[k] for k in keys], backend=backend,
                    initial_layout=chain, optimization_level=1,
                    seed_transpiler=SEED_TRANSPILER)
    transpiled = {k: tc for k, tc in zip(keys, tcs)}
    assert_report = {}
    for (tp_idx, prep_idx), tc in transpiled.items():
        st = assert_skeleton_invariants(tc, GRID_STEPS[tp_idx], chain=chain,
                                        label=f"t{tp_idx}/p{prep_idx}")
        assert_report[f"t{tp_idx}/p{prep_idx}"] = st
    print("  post-transpile assertions (all-rzz, budget, angles, "
          "one-injected-RZ, no routing, no delay, layout identity): PASS")

    pubs, meta = make_pubs(transpiled, plan, M_BIND, SHOTS_BIND,
                           SHOTS_SCIENCE)
    assert_shared_skeleton(pubs, meta)

    # realized-profile class assertion, on WHAT FLIES: rebuild each
    # swept cell's table from the bound sweep columns; the lit LOGICAL
    # sites must match the pinned profile, whose class is asserted, and
    # the layout-identity assertion above ties logical site l to
    # physical qubit chain[l]
    for info in meta:
        if not info["swept"]:
            continue
        arm, tp_idx, pi = info["arm"], info["tp_idx"], info["prep_idx"]
        tc = transpiled[(tp_idx, pi)]
        tab = rebuild_table_from_pub(pubs[info["pub_index"]], tc,
                                     GRID_STEPS[tp_idx])
        lit = {l for l in range(N) if np.std(tab[:, :, l]) > 0}
        prof = ARM_PROFILE[arm]
        want_lit = (set() if prof is None
                    else {l for l in range(N) if prof[l] > 0})
        _require(lit == want_lit,
                 f"realized-profile assertion FAILED at "
                 f"{arm}/t{tp_idx}/p{pi}: lit sites {sorted(lit)} != "
                 f"pinned {sorted(want_lit)}")
    _require(profile_class(ARM_PROFILE["C"]) == "maximizing", "C class")
    _require(profile_class(ARM_PROFILE["Cp"]) == "non-maximizing",
             "Cp class")
    print("  realized-profile assertion (bound sweeps -> lit sites -> "
          "class, on the identity layout): PASS")

    # dose certificates on ALL swept cells, computed FROM THE BOUND
    # SWEEPS (certifies what flies; catches binding/transposition bugs)
    dose, worst_dose, worst_cell = {}, 0.0, None
    for info in meta:
        if not info["swept"] or info["arm"] == "N0":
            continue
        arm, tp_idx, pi = info["arm"], info["tp_idx"], info["prep_idx"]
        tc = transpiled[(tp_idx, pi)]
        tab = rebuild_table_from_pub(pubs[info["pub_index"]], tc,
                                     GRID_STEPS[tp_idx])
        cert = dose_certificates(tab, arm)
        dose[f"{arm}/t{tp_idx}/p{pi}"] = cert
        if cert["max_dev"] > worst_dose:
            worst_dose, worst_cell = cert["max_dev"], (arm, tp_idx, pi)
    _require(worst_dose < DOSE_CRITERION,
             f"dose certificate FAILED pre-submit: {worst_dose:.4f} at "
             f"{worst_cell}")
    print(f"  dose certificates (ALL {sum(1 for i in meta if i['swept'] and i['arm'] != 'N0')} swept cells, from bound "
          f"sweeps): worst {worst_dose:.4f} at {worst_cell} PASS")

    cal_tc = transpile(cal_circuits(), backend=backend,
                       initial_layout=chain, optimization_level=1,
                       seed_transpiler=SEED_TRANSPILER)
    # aux stations: the IDENTICAL 6-PUB set (T1 short/long, Ramsey
    # canary, Hahn echo, CAL0/CAL1) flies after EVERY station depth,
    # so drift and
    # delay are deconfounded (round-2 blocker: one delay per station
    # made mid-batch T1 telegraphing invisible). Delays are scaled to
    # the deep circuit's estimated duration and ALIGNED to the backend
    # timing granularity (an unaligned delay in an auxiliary PUB must
    # not be able to reject a science job).
    dur_rzz = _typical_rzz_duration(target, chain)
    # LAYERED duration model (round-8 repair: the previous serial
    # gate-count model gave ~145 us, ~4x the physical wall time, and
    # its 1.5x T1 point at ~218 us sat past Heron T1, exactly what
    # the comment below says to avoid). The three disjoint rzz of a
    # sub-layer schedule in PARALLEL: 6 rzz sub-layers per Strang
    # step + the rz layer; prep ~5 sub-layers; 1.5x buffer for 1q
    # overhead and scheduling gaps. The doc (section 10) carries the
    # same model; the transpiled schedule at submit is the authority
    # and the day-of addendum records the measured tau_step_us.
    # round 12: the prep books TEN rzz sub-layers, not five. The five
    # Givens act on the OVERLAPPING pairs (l, l+1), so they cannot
    # schedule in parallel, and each is two rzz. 1.4% of t_deep, but
    # the stations derive from this number. The injection rz layer
    # stays at zero duration on purpose: RZ is a virtual frame change
    # on Heron. Section 10 dropped the same allowance in round 12, so
    # doc and runner are ONE model at two buffers; this comment had
    # section 10 still charging ~100 ns per step and running above
    # this one, superseded text the executor reads at submit (round
    # 15). The retired gap figure is not restated here: section 10
    # carries it, and a number in two places is how it goes stale.
    n_sublayers = 10 + 6 * STEPS_MAX
    t_deep_us = n_sublayers * dur_rzz * 1e6 * 1.5
    # T1 points at 0.5x and 1.5x the deep-circuit duration (3x would
    # sit past T1 on Heron and read noise); the echo delay is built as
    # 2 x an ALIGNED half so both halves stay granularity-aligned
    d_t1 = _align_delay_us(max(10.0, 0.5 * t_deep_us), backend)
    d_echo = 2 * _align_delay_us(max(10.0, 0.25 * t_deep_us), backend)
    st_c, st_meta = t1t2_station_circuits(d_t1, d_echo)
    st_tc = transpile(st_c, backend=backend, initial_layout=chain,
                      optimization_level=1,
                      seed_transpiler=SEED_TRANSPILER,
                      scheduling_method=None)
    aux_by_station = {}
    t1t2_meta = []
    for k, tp_after in enumerate(STATION_AFTER_TP):
        recs = [dict(st_meta[i], station=k, after_tp=tp_after)
                for i in range(len(st_tc))]
        t1t2_meta.extend(recs)
        aux_by_station[tp_after] = list(zip(st_tc, recs))
    n_aux = sum(len(v) for v in aux_by_station.values())

    # final PUB order: CAL0, CAL1, then depth blocks ascending with the
    # aux station after its depth block; job chunking respects depth-
    # block boundaries (never splits inside a depth, never by arm)
    final_pubs, final_kind, blocks = [], [], []
    blk = []
    for c in cal_tc:
        final_kind.append(("cal", "cal0" if not final_pubs else "cal1"))
        final_pubs.append((c, None, SHOTS_SCIENCE))
        blk.append(len(final_pubs) - 1)
    blocks.append(("cal", blk))
    by_tp = {}
    for info in meta:
        by_tp.setdefault(info["tp_idx"], []).append(info)
    for tp_idx in range(len(GRID_STEPS)):
        blk = []
        for info in by_tp[tp_idx]:
            info["pub_index_final"] = len(final_pubs)
            final_kind.append(("science", info))
            final_pubs.append(pubs[info["pub_index"]])
            blk.append(len(final_pubs) - 1)
        for (c, mrec) in aux_by_station.get(tp_idx, []):
            final_kind.append(("aux", mrec))
            final_pubs.append((c, None, SHOTS_AUX))
            blk.append(len(final_pubs) - 1)
        blocks.append((f"tp{tp_idx}", blk))
    n_sci = len(meta)
    print(f"  PUBs: {n_sci} science + 2 CAL + {n_aux} aux "
          f"(identical station set after depths {STATION_AFTER_TP}; "
          f"T1 {d_t1:.1f}/{3 * d_t1:.1f} us, Ramsey+echo {d_echo:.1f} us,"
          f" {SHOTS_AUX} shots each)")

    # payload guard: chunk into jobs on the parameter-payload size,
    # respecting depth-block boundaries (max_circuits alone counts PUBs
    # and would submit a ~280 MB single job)
    def pub_mb(p):
        return (p[1].size * 8 / 1e6) if p[1] is not None else 0.0

    def pub_bc(p):
        return p[1].shape[0] if p[1] is not None else 1
    job_plan, cur = [], []
    cur_mb, cur_bc = 0.0, 0
    for desc, blk in blocks:
        mb = sum(pub_mb(final_pubs[i]) for i in blk)
        bc = sum(pub_bc(final_pubs[i]) for i in blk)
        if cur and (cur_mb + mb > PAYLOAD_MAX_MB
                    or cur_bc + bc > BOUND_CIRC_MAX):
            job_plan.append(cur); cur, cur_mb, cur_bc = [], 0.0, 0
        cur.extend(blk); cur_mb += mb; cur_bc += bc
    if cur:
        job_plan.append(cur)
    # a CAL pair flies IN EVERY JOB (the head pair alone would leave
    # the deep half of the grid inverted with a head-of-batch Minv);
    # the analysis pools every CAL pair into the inversion (Amendment 1)
    # per-job CALs fly at SHOTS_AUX: at SHOTS_SCIENCE the ten extra
    # pairs alone would push the projection to the 25-min cap; pooled
    # over all pairs the confusion precision stays ~0.05%/qubit
    for jn in range(1, len(job_plan)):
        for kind, cal in zip((f"cal0_job{jn}", f"cal1_job{jn}"), cal_tc):
            final_kind.append(("cal", kind))
            final_pubs.append((cal, None, SHOTS_AUX))
            job_plan[jn].append(len(final_pubs) - 1)
    total_mb = sum(pub_mb(p) for p in final_pubs)
    max_c = getattr(backend, "max_circuits", None)
    if max_c is not None:
        _require(all(len(j) <= max_c for j in job_plan),
                 f"a job exceeds max_circuits={max_c}; shrink "
                 f"PAYLOAD_MAX_MB")
    print(f"  parameter payload {total_mb:.0f} MB -> {len(job_plan)} "
          f"job(s) inside one Batch (cap {PAYLOAD_MAX_MB:.0f} MB/job, "
          f"depth blocks never split)")
    # round 14: Amendment 1.5 tells the executor to carry the REALIZED
    # job plan into the day-of addendum and check it against the
    # pre-registered 8-to-11 range, and the per-job numbers were
    # computed and never printed or persisted, so there was nothing to
    # carry.
    job_plan_record = [
        {"job": jn, "pubs": len(blk),
         "bound_circuits": int(sum(pub_bc(final_pubs[i]) for i in blk)),
         "payload_mb": round(sum(pub_mb(final_pubs[i]) for i in blk), 1)}
        for jn, blk in enumerate(job_plan)]
    for rec in job_plan_record:
        print(f"    job{rec['job']}: {rec['pubs']} PUBs, "
              f"{rec['bound_circuits']} bound circuits "
              f"({100.0 * rec['bound_circuits'] / BOUND_CIRC_MAX:.0f}% "
              f"of cap), {rec['payload_mb']:.1f} MB")

    # Class 1: billing projection under the cap. Per-binding overhead
    # RE-MEASURED (2026-08-17, read-only job.usage() query, twice):
    # concentrator job d99a970tcv6s73dn2atg (24 sink PUBs x 256
    # bindings = 6144 bound circuits, 12 no-sink science PUBs unbound,
    # 376,832 shots) billed 119 s =
    # 0.316 ms/shot, INSIDE the delay-bearing anchor band 0.309-0.327
    # despite its 6144 bound circuits (the job carried 8 delay-bearing
    # aux PUBs; an earlier same-day note quoted 69 s = 0.183, refuted
    # by the API re-query). No gross per-circuit surcharge is visible,
    # but the 0.327 anchor's margin over the measured 0.316 is only
    # ~3%: the SERIAL submit brake below is the binding backstop, not
    # this projection.
    n_jobcal = sum(1 for k, _ in final_kind
                   if k == "cal" and _ not in ("cal0", "cal1"))
    total_shots = ((n_sci + 2) * SHOTS_SCIENCE
                   + (n_aux + n_jobcal) * SHOTS_AUX)
    proj_min = total_shots * BILL_ANCHOR_MS_PER_SHOT / 1000 / 60
    print(f"  billing projection: {total_shots/1e6:.2f}M shots ~ "
          f"{proj_min:.1f} QPU min (cap {BILL_ABORT_MIN})")
    if proj_min > BILL_ABORT_MIN:
        raise SystemExit(f"ABORT: projection {proj_min:.1f} min > cap")

    # pending queue (wall-clock protection; validity is not affected)
    try:
        st = backend.status()
        print(f"  backend status: operational={st.operational}, "
              f"pending={st.pending_jobs}"
              + (f"  [WARNING: queue > {QUEUE_WARN}]"
                 if st.pending_jobs > QUEUE_WARN else ""))
        if not st.operational:
            raise SystemExit("ABORT: backend not operational")
    except SystemExit:
        raise
    except Exception as e:
        # round 14: this handler used to print and CONTINUE, so a
        # status() exception skipped the operational check entirely
        # while section 9 lists it as a Class-1 guard that prevents
        # the spend. A guard that cannot read its input fails CLOSED
        # on a paid submission; the dry run still reports and goes on.
        print(f"  backend status unavailable: {e!r}")
        if yes:
            print("\n  ** HARD ABORT: the operational check is a "
                  "Class-1 guard and it could not be read.")
            sys.exit(1)

    if not yes:
        print("\n  DRY STOP (no --yes): nothing submitted. Re-run with "
              "--hardware --yes after the day-of addendum is committed.")
        _save_json({"mode": "hardware_dry", "chain": chain,
                    "properties_before": snap, "pub_meta": meta,
                    "dose_worst": {"dev": worst_dose,
                                   "cell": worst_cell},
                    "assertions": assert_report, "twin": twin_status,
                    "constants_pending": pending,
                    "jobs_planned": len(job_plan),
                    "job_plan": job_plan_record,
                    "payload_mb": total_mb,
                    "projection_min": proj_min}, "hw_dry")
        return

    ans = input("Type FLY to submit (anything else aborts): ")
    if ans != "FLY":
        raise SystemExit("Aborted: FLY not confirmed.")

    from qiskit_ibm_runtime import SamplerV2, Batch
    ts_run = datetime.now().strftime("%Y%m%d_%H%M%S")
    stub_path = RESULTS_DIR / f"corner_beat_hw_stub_{ts_run}.json"
    # persist the flown phase tables themselves (seeds + sha alone
    # hold paid data hostage to numpy stream stability)
    tables_npz = RESULTS_DIR / f"corner_beat_hw_tables_{ts_run}.npz"
    np.savez_compressed(
        tables_npz,
        **{f"{arm}_p{pi}": phase_table_full(arm, pi, M_BIND)
           for arm in ("U", "Cp", "C") for pi in range(3)},
        # N0's swept ZERO tables too: Amendment 1.4 says "the flown
        # tables are persisted", and the record must not make the
        # reader reconstruct even a deterministic zero
        **{f"N0_p{pi}": np.zeros_like(phase_table_full("U", pi, M_BIND))
           for pi in range(3)})
    print(f"  flown phase tables persisted -> {tables_npz.name}")
    payload = {"mode": "hardware", "timestamp": ts_run,
               "backend": BACKEND_NAME, "chain": chain,
               "chain_info": chain_info, "git_head": _git_head(),
               "rcpsi_head": _rcpsi_git_gate()[1],
               "dayof_addendum": dayof,
               "tables_npz": tables_npz.name,
               "basis_sha": demod_basis_sha(),
               "basis": demod_basis_record(),
               "constants_sha": hashlib.sha256(json.dumps(
                   consts, sort_keys=True).encode()).hexdigest()[:16],
               "config": {"N": N, "JDT": JDT, "Q": Q, "K_GRID": K_GRID,
                          "GRID_STEPS": GRID_STEPS, "M_BIND": M_BIND,
                          "SHOTS_BIND": SHOTS_BIND,
                          "SHOTS_SCIENCE": SHOTS_SCIENCE,
                          "seed_tables": SEED_TABLES,
                          "seed_transpiler": SEED_TRANSPILER,
                          "seed_boot": SEED_BOOT,
                          "payload_max_mb": PAYLOAD_MAX_MB,
                          "bound_circ_max": BOUND_CIRC_MAX,
                          # round 13: section 10 says all four
                          # code-pinned flight constants are recorded
                          # in every hardware payload; the anchor and
                          # the station placement were not
                          "bill_anchor_ms_per_shot": BILL_ANCHOR_MS_PER_SHOT,
                          # round 15: the realized plan was persisted
                          # in the DRY artifact only and reached the
                          # flight record just via the copied
                          # addendum, i.e. on executor discipline
                          "job_plan": job_plan_record,
                          "station_after_tp": list(STATION_AFTER_TP),
                          "kept_count_floor": KEPT_COUNT_FLOOR,
                          "numpy": np.__version__},
               "pub_meta": meta, "dose": dose,
               "assertions": assert_report, "twin": twin_status,
               "properties_before": snap,
               "t1t2_meta": t1t2_meta,
               "delays_us": {"t1": [d_t1, 3 * d_t1], "echo": d_echo},
               "jobs": {},
               "status": "SUBMITTED (awaiting results)"}

    def write_stub(with_counts=False):
        """Persist the in-flight payload. with_counts=True on the two
        TERMINAL paths (brake, job failure) so the stub is a readable
        artifact and not just a list of job IDs; the routine
        after-every-submit calls stay counts-free, since raw is ~100 MB
        at flight scale and rewriting it 11 times would cost minutes of
        wall clock inside the Batch (round 11)."""
        try:
            if with_counts:
                payload["raw_counts"] = raw
                payload["job_usage"] = usage
                # analyze keys the brake state off startswith("BRAKE");
                # the stub said "BILLING BRAKE ...", so a stub read
                # after a crash lost the brake label (round 11 minor)
                st = str(payload.get("status", ""))
                if st.startswith("BILLING BRAKE"):
                    payload["brake"] = st
                    payload["status"] = ("BRAKE-TRUNCATED (flight VOID "
                                         "per Amendment 1.6; partial "
                                         "counts follow)")
            with open(stub_path, "w") as f:
                json.dump(payload, f, indent=2, default=str)
        except Exception:
            pass

    raw = {"science": {}, "cal": {}, "t1t2": []}
    usage = {}
    billed_s = 0.0
    with Batch(backend=backend) as batch:
        sampler = SamplerV2(mode=batch)
        _disable_mitigation(sampler)
        # SERIAL submit -> wait -> read usage per job: the billing cap
        # must be able to BRAKE mid-flight (submitting all jobs up
        # front makes the cap a projection with no hands)
        fallback_s = 0.0
        for jn, idxs in enumerate(job_plan):
            job_shots = sum(final_pubs[i][2] * pub_bc(final_pubs[i])
                            for i in idxs)
            proj_next = job_shots * BILL_ANCHOR_MS_PER_SHOT / 1000.0
            spent = max(billed_s, fallback_s)
            # PREDICTIVE brake on max(measured, shots-fallback): the
            # runtime's usage() can be None/unreadable, and a brake
            # with no hands is worse than none
            if (spent + proj_next) / 60.0 > BILL_ABORT_MIN:
                payload["status"] = (
                    f"BILLING BRAKE before job{jn}: spent "
                    f"{spent / 60.0:.1f} min (measured "
                    f"{billed_s / 60.0:.1f}, shots-fallback "
                    f"{fallback_s / 60.0:.1f}) + next-job projection "
                    f"{proj_next / 60.0:.1f} > cap {BILL_ABORT_MIN}; "
                    f"remaining jobs NOT submitted")
                print(f"\n  ** {payload['status']}")
                write_stub(with_counts=True)
                break
            desc = f"job{jn}"
            job = sampler.run([final_pubs[i] for i in idxs])
            payload["jobs"][desc] = job.job_id()
            write_stub()          # job IDs persisted after EVERY submit
            print(f"  {desc}: {len(idxs)} PUBs -> {job.job_id()}")
            try:
                res = job.result()
            except Exception as e:
                # Round 11 (MAJOR): raw_counts was attached to the
                # payload only AFTER the submit loop, so this path
                # saved a paid artifact with no counts key at all and
                # --analyze died on it at art["raw_counts"]. Every
                # count retrieved before the failure was on disk and
                # unreadable by the committed estimator, and the state
                # HAS a registered reading (VOID, grid incomplete,
                # with D-mag/A as REPORTED excluded-cell fits). Attach
                # what was retrieved, then re-raise.
                payload["status"] = f"JOB FAILED ({desc}): {e!r}"
                write_stub(with_counts=True)
                _save_json(payload, "hw")
                raise
            u_raw = _safe(lambda: job.usage())
            usage[desc] = {"usage": u_raw,
                           "metrics_usage": _safe(
                               lambda: job.metrics().get("usage"))}
            fallback_s += proj_next
            if isinstance(u_raw, (int, float)):
                billed_s += float(u_raw)
                print(f"    billed so far: {billed_s / 60.0:.2f} min "
                      f"(cap {BILL_ABORT_MIN})")
            else:
                print(f"    ** WARNING: job.usage() unreadable "
                      f"({u_raw!r}); the brake runs on the shots "
                      f"fallback ({fallback_s / 60.0:.2f} min so far)")
            for local_i, global_i in enumerate(idxs):
                kind, info = final_kind[global_i]
                try:
                    counts = _pub_counts_list(res[local_i])
                except Exception as e:
                    raw.setdefault("errors", []).append(
                        {"pub": global_i, "error": repr(e)})
                    continue
                if kind == "science":
                    if info["swept"] and len(counts) != M_BIND:
                        raw.setdefault("errors", []).append(
                            {"pub": global_i,
                             "error": f"per-binding structure lost: "
                                      f"{len(counts)} != {M_BIND}"})
                    raw["science"].setdefault(info["arm"], {}).setdefault(
                        str(GRID_STEPS[info["tp_idx"]]), {})[
                        str(info["prep_idx"])] = counts
                elif kind == "cal":
                    raw["cal"][info] = counts[0]
                else:
                    mrec = dict(info)
                    mrec["counts"] = counts[0]
                    raw["t1t2"].append(mrec)
            # partial dump after every job (a dump failure cannot kill
            # the flight). Round 11 (MAJOR): this dumped the bare
            # counts dict, which --analyze cannot read at all (it dies
            # on art["config"] before it reaches the counts). The
            # partial is the ONLY artifact if the process dies between
            # jobs, so it carries the full artifact schema.
            try:
                part = {k: v for k, v in payload.items()
                        if k != "raw_counts"}
                part["raw_counts"] = raw
                part["job_usage"] = usage
                part["partial_after"] = desc
                part["status"] = (f"PARTIAL (counts through {desc}; "
                                  f"grid incomplete by construction)")
                with open(RESULTS_DIR /
                          f"corner_beat_hw_partial_{ts_run}_{desc}.json",
                          "w") as f:
                    json.dump(part, f, default=str)
            except Exception:
                pass
    payload["raw_counts"] = raw
    payload["job_usage"] = usage
    # the brake's status string is the RECORD's evidence that the
    # flight was parked by Amendment 1.6's registered rule; it must
    # never be overwritten by the completion label (round-7 repair:
    # the fall-through stamped COUNTS RETRIEVED unconditionally)
    # round 11: write_stub(with_counts=True) may ALREADY have
    # normalized the label on the brake path, so the test must accept
    # both spellings; keying only on "BILLING BRAKE" would have let
    # the fall-through stamp COUNTS RETRIEVED over a brake-truncated
    # flight, which is the exact overwrite round 7 repaired.
    _st = str(payload.get("status", ""))
    if _st.startswith("BILLING BRAKE"):
        payload["brake"] = _st
        payload["status"] = ("BRAKE-TRUNCATED (flight VOID per "
                             "Amendment 1.6; partial counts follow)")
    elif not _st.startswith("BRAKE"):
        payload["status"] = "COUNTS RETRIEVED"
    payload["properties_after"] = snapshot_after(backend, chain)
    path = _save_json(payload, "hw")
    print(f"\n  counts persisted BEFORE any analysis -> {path}")
    print("  analyze with: python run_corner_beat.py --analyze "
          f"\"{path}\"")


def snapshot_after(backend, chain):
    snap = {"timestamp": datetime.now().isoformat(), "qubits": {}}
    try:
        props = backend.properties()
        for q in chain:
            snap["qubits"][str(q)] = {
                "T1_us": props.t1(q) * 1e6, "T2_us": props.t2(q) * 1e6,
                "readout_error": props.readout_error(q)}
    except Exception as e:
        snap["error"] = repr(e)
    return snap


def _has_refresh(backend):
    try:
        import inspect
        return "refresh" in inspect.signature(backend.properties).parameters
    except Exception:
        return False


def _typical_rzz_duration(target, chain):
    for a, b in zip(chain, chain[1:]):
        try:
            props = target["rzz"].get((a, b)) or target["rzz"].get((b, a))
            if props is not None and props.duration:
                return float(props.duration)
        except Exception:
            continue
    return 100e-9


def _disable_mitigation(sampler):
    """DD and twirling OFF, FAIL-CLOSED: DD would rewrite the very
    dephasing profile under test, so failing to set these flags is a
    hard abort, not a shrug (the transpiled-circuit no-delay assertion
    covers the compile side; this covers the runtime side)."""
    try:
        sampler.options.dynamical_decoupling.enable = False
        sampler.options.twirling.enable_gates = False
        sampler.options.twirling.enable_measure = False
    except Exception as e:
        raise SystemExit(
            f"ABORT: cannot disable runtime DD/twirling ({e!r}); "
            f"flying with mitigation state unknown is forbidden (the "
            f"section 5 no-runtime-mitigation requirement)")
    got = (sampler.options.dynamical_decoupling.enable,
           sampler.options.twirling.enable_gates,
           sampler.options.twirling.enable_measure)
    _require(got == (False, False, False),
             f"mitigation flags read {got}")
    print("  sampler mitigation: dd.enable=False, twirling.gates=False, "
          "twirling.measure=False (asserted)")


def _safe(fn):
    try:
        return fn()
    except Exception as e:
        return f"unavailable: {e!r}"


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--certify", action="store_true")
    g.add_argument("--calibrate", action="store_true")
    g.add_argument("--aer", action="store_true")
    g.add_argument("--hardware", action="store_true")
    g.add_argument("--analyze", metavar="FILE")
    ap.add_argument("--full", action="store_true",
                    help="--aer at the full M = 1024")
    ap.add_argument("--yes", action="store_true",
                    help="actually submit in --hardware")
    a = ap.parse_args()
    np.set_printoptions(precision=5, suppress=True)
    if not CONSTANTS_PATH.exists():
        with open(CONSTANTS_PATH, "w") as f:
            json.dump(DEFAULT_CONSTANTS, f, indent=2)
        print(f"(wrote default constants -> {CONSTANTS_PATH})\n")
    if a.certify:
        cmd_certify()
    elif a.calibrate:
        cmd_calibrate()
    elif a.aer:
        cmd_aer(full=a.full)
    elif a.hardware:
        cmd_hardware(yes=a.yes)
    elif a.analyze:
        cmd_analyze(a.analyze)


if __name__ == "__main__":
    main()
