#!/usr/bin/env python3
"""Candidate 2 of INCOMPLETENESS_PROOF, re-measured: what the evidence is.

Candidate 2 asks whether a single decaying qubit inside the system can be the
origin of the palindromic dephasing the mirror symmetry requires. Since March
2026 the answer "no" was carried by three readings of `failed_third.py`:

  LEG 1   an effective dephasing rate gamma_eff = 0 on the neighbours
  LEG 2   a non-Markovian signature, ~50% of steps, max deviation 0.000000
  LEG 3   process tomography giving 0/16 palindromic pairs, control 16/16

None of the three survives, and the third fails for a reason the repository
had already written down elsewhere: F137's centre identity says the palindrome
centre is trace(L)/dim, so exactly one candidate centre exists and searching
for one is never necessary. `failed_third.py` searches, on a grid far coarser
than its own tolerance, from a seed that is not the centre.

  Section 1   legs 1 and 2 are code paths, not measurements
  Section 2   what the neighbours actually lose
  Section 3   leg 3's instrument: a forced centre, searched
  Section 4   the control the experiment never ran, and what it shows
  Section 5   the full three-qubit system, which IS a palindrome
  Section 6   the trace identity, and what it leaves standing

The physics is not re-implemented: `failed_third.py`'s own builders are exec'd,
so every Liouvillian here is byte-identical to the one behind the March run.
That file opens its results file for writing at MODULE level, so a plain
`import` would truncate `simulations/results/failed_third.txt`; the sink is
redirected below and the redirection is asserted, never assumed.

Instruments. The multiset comparison is a Hungarian MATCHING, never a sort or
a nearest-neighbour scan without removal (both lose multiplicity). Where an
exact route exists it is used and compared exactly: the power sums
trace(L^k) vs trace((2cI - L)^k) are matrix arithmetic with no eigensolver.

Script:  simulations/incompleteness_candidate2_evidence.py
Output:  simulations/results/incompleteness_candidate2_evidence.txt
Cited by: docs/proofs/INCOMPLETENESS_PROOF.md (Section 2, Candidate 2)
"""

import os
import sys
import tempfile

import numpy as np
from scipy.linalg import logm
from scipy.optimize import linear_sum_assignment, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "failed_third.py")
OUT_PATH = os.path.join(HERE, "results", "incompleteness_candidate2_evidence.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")


_SINK = ('OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),\n'
         '                        "results", "failed_third.txt")')
with open(SRC, "r", encoding="utf-8") as fh:
    _src = fh.read()
if _SINK not in _src:
    sys.exit("ABORT: failed_third.py's OUT_PATH block changed shape. Refusing to "
             "exec it: a wrong replacement would truncate its committed results.")
_src = _src.replace(_SINK, "OUT_PATH = %r" % os.path.join(tempfile.gettempdir(),
                                                          "_candidate2_sink.txt"))
_mod = {"__name__": "_candidate2_evidence", "__file__": SRC}
exec(compile(_src, SRC, "exec"), _mod)

build_origin_model, build_H = _mod["build_origin_model"], _mod["build_H"]
evolve, ptrace_keep = _mod["evolve"], _mod["ptrace_keep"]
site_op, build_L_H, add_jump = _mod["site_op"], _mod["build_L_H"], _mod["add_lindblad_jump"]
best_palindrome = _mod["best_palindrome"]
I2, sx, sy, sz = _mod["I2"], _mod["sx"], _mod["sy"], _mod["sz"]
N, d, d2 = _mod["N"], _mod["d"], _mod["d2"]

PAULI2 = [np.kron(P, Q) for P in (I2, sx, sy, sz) for Q in (I2, sx, sy, sz)]
T_EFF = 5.0
NAMES = {"A": "amplitude damping", "B": "detuning h=10",
         "C": "thermal bath", "D": "X+Y dephasing"}
PSI_COH, PHI_COH = (1, 2), (0, 3)      # |01><10|  and  |00><11|
F137_CENTRE = {"A": -0.05, "B": 0.0, "C": -0.075, "D": -0.20}

_gates = []


def gate(name, ok, detail):
    _gates.append(bool(ok))
    log("  [%s] %-46s %s" % ("PASS" if ok else "FAIL", name, detail))


def liouv2(H, jumps):
    """Standard Lindblad generator on a d-dimensional space, rates as given."""
    dd = H.shape[0]
    Id = np.eye(dd, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for F, r in jumps:
        FdF = F.conj().T @ F
        L = L + r * (np.kron(F, F.conj()) - 0.5 * np.kron(FdF, Id)
                     - 0.5 * np.kron(Id, FdF.T))
    return L


def bell(sign):
    psi = np.zeros(d, dtype=complex)
    psi[0] = 1 / np.sqrt(2)
    psi[0b110] = sign / np.sqrt(2)
    return np.outer(psi, psi.conj())


def trajectory(L, rho0, ts):
    pur, psi_c, phi_c = [], [], []
    for t in ts:
        r = ptrace_keep(evolve(L, rho0, t), N, [0, 1])
        pur.append(float(np.real(np.trace(r @ r))))
        psi_c.append(abs(r[PSI_COH]))
        phi_c.append(abs(r[PHI_COH]))
    return np.array(pur), np.array(psi_c), np.array(phi_c)


def match_resid(ev, c):
    """Hungarian matching of {lam} against {2c - lam}: worst matched distance."""
    cost = np.abs(ev[:, None] - (2 * c - ev)[None, :])
    r, col = linear_sum_assignment(cost)
    return float(np.max(cost[r, col]))


def powersum_resid(L, c, kmax=20):
    """Exact route, no eigensolver: trace(L^k) vs trace((2cI - L)^k)."""
    D = L.shape[0]
    A = np.eye(D, dtype=complex)
    B = np.eye(D, dtype=complex)
    M = 2 * c * np.eye(D, dtype=complex) - L
    worst = scale = 0.0
    for _ in range(kmax):
        A, B = A @ L, B @ M
        a, b = np.trace(A), np.trace(B)
        worst = max(worst, abs(a - b))
        scale = max(scale, abs(a), abs(b))
    return worst / scale if scale > 0 else float("nan")


def transfer(Lop, buggy=False):
    """Pauli transfer matrix of the induced Q1-Q2 channel at t = T_EFF.
    buggy=True reproduces failed_third.py:361 + :371 (input /4 AND output /4)."""
    q3 = np.array([[1, 0], [0, 0]], dtype=complex)
    T = np.zeros((16, 16), dtype=complex)
    for b in range(16):
        rin = PAULI2[b] / 4.0 if buggy else PAULI2[b]
        rout = ptrace_keep(evolve(Lop, np.kron(rin, q3), T_EFF), N, [0, 1])
        for a in range(16):
            T[a, b] = np.trace(PAULI2[a] @ rout) / 4.0
    return np.real(T)


def marginal_residual(Lop):
    """Best pairing residual of the induced channel's rate spectrum. The centre
    of a MARGINAL is not given by the parent's trace, so here it is minimised;
    that is the most generous reading available to the experiment."""
    ev = np.linalg.eigvals(logm(transfer(Lop)) / T_EFF)

    def f(c):
        return match_resid(ev, c)
    g = np.linspace(-1.0, 1.0, 4001)
    c0 = g[int(np.argmin([f(c) for c in g]))]
    r = minimize_scalar(f, bracket=(c0 - 5e-4, c0, c0 + 5e-4), method="brent",
                        options={"xtol": 1e-14})
    return (r.fun, r.x) if r.fun < f(c0) else (f(c0), c0)


log("Candidate 2 of INCOMPLETENESS_PROOF, re-measured")
log("=" * 78)

t100, t200 = np.linspace(0, 20, 100), np.linspace(0, 20, 200)
L = {o: build_origin_model(o) for o in "ABCD"}
traj = {o: trajectory(L[o], bell(+1), t100) for o in "ABCD"}

# ==================================================================
log()
log("SECTION 1 -- legs 1 and 2 are code paths, not measurements")
log("-" * 78)
log("failed_third.py:230-239 fits an exponential to rho_12[1,2] = |01><10| and")
log("normalises by its value at t=0. The preparation is Bell Phi+, whose only")
log("coherence is |00><11|, so |01><10| starts at zero, the guard protecting")
log("that division cannot pass, and its else-value 0 was read as a rate.")
log()
for o in "ABCD":
    gate("1a.%s  |01><10| at t=0 is exactly 0" % o, traj[o][1][0] == 0.0,
         "value = %r  (exact compare, no tolerance)" % traj[o][1][0])
gate("1b  |00><11| at t=0 is 1/2 to one ulp",
     all(abs(traj[o][2][0] - 0.5) <= 2 ** -52 for o in "ABCD"),
     "value = %.17g  (analytically 1/2; (1/sqrt2)^2 is one ulp low)" % traj["A"][2][0])
gate("1c  no option can reach the fit branch",
     all(traj[o][1][0] <= 1e-6 for o in "ABCD"),
     "the guard's first clause (>5 points above 1e-6) passes with 99; the second decides")
log()
log("The same empty element carries leg 2: failed_third.py:475-478 prints")
log("`max deviation` from an `else: deviation = 0` branch, and the count it")
log("labels a trace distance counts increases of |01><10| itself.")
log()


def increases(y, eps=1e-8):
    return int(np.sum(y[1:] > y[:-1] + eps))


inc = {o: increases(trajectory(L[o], bell(+1), t200)[1]) for o in "ABCD"}
gate("2a  the published count is identical across four mechanisms",
     len(set(inc.values())) == 1, ", ".join("%s=%d/199" % (o, inc[o]) for o in "ABCD"))
_bh = _mod["build_H"]
_mod["build_H"] = lambda J=0.0: _bh(0.0)
try:
    inc0 = {o: increases(trajectory(build_origin_model(o), bell(+1), t200)[1]) for o in "ABCD"}
finally:
    _mod["build_H"] = _bh
gate("2b  at J=0 it is 0 in all four",
     all(v == 0 for v in inc0.values()),
     ", ".join("%s=%d/199" % (o, inc0[o]) for o in "ABCD") +
     "  (the element is identically zero there, so this shows the count is "
     "H-driven, not that nothing returns)")
log()
log("What the script never measured, measured here: a BLP probe, the trace")
log("distance between the reduced states of Phi+ and Phi-, both tensor |0><0|.")
log("The sum over a grid grows with grid density, so read the sign and the")
log("order; the largest single rise beside it is grid-stable.")
log()
blp, big = {}, {}
for o in "ABCD":
    ts = np.linspace(0, 20, 400)
    da = [ptrace_keep(evolve(L[o], bell(+1), t), N, [0, 1]) for t in ts]
    db = [ptrace_keep(evolve(L[o], bell(-1), t), N, [0, 1]) for t in ts]
    td = np.array([0.5 * np.sum(np.abs(np.linalg.eigvalsh(a - b))) for a, b in zip(da, db)])
    rise = np.diff(td)
    blp[o], big[o] = float(np.sum(rise[rise > 0])), float(np.max(rise))
    log("  option %s (%-18s BLP total rise = %6.3f, largest single rise = %.4f"
        % (o, NAMES[o] + ")", blp[o], big[o]))
gate("2c  the reduced dynamics IS non-Markovian",
     all(v > 1.0 for v in blp.values()) and all(v > 1e-3 for v in big.values()),
     "leg 2's CONCLUSION may well hold; what it offered as evidence did not")

# ==================================================================
log()
log("SECTION 2 -- what the neighbours actually lose")
log("-" * 78)
log("  %-4s %-20s %10s %10s %12s %12s" %
    ("opt", "mechanism", "purity_0", "purity_20", "|00><11|_0", "|00><11|_20"))
for o in "ABCD":
    p, _, c = traj[o]
    log("  %-4s %-20s %10.6f %10.6f %12.6f %12.6f" % (o, NAMES[o], p[0], p[-1], c[0], c[-1]))
gate("3a  A, C, D lose more than half their purity",
     all(traj[o][0][0] - traj[o][0][-1] > 0.5 for o in "ACD"),
     "drops: " + ", ".join("%s=%.6f" % (o, traj[o][0][0] - traj[o][0][-1]) for o in "ACD"))
gate("3b  B is a near-identity channel",
     traj["B"][0][0] - traj["B"][0][-1] < 0.01,
     "drop = %.6f, so it can carry no weight either way"
     % (traj["B"][0][0] - traj["B"][0][-1]))

# ==================================================================
log()
log("SECTION 3 -- leg 3's instrument: a forced centre, searched")
log("-" * 78)
log("F137 (ANALYTICAL_FORMULAS.md) states the centre identity: a multiset")
log("closed under lam -> 2c - lam has c = trace(L)/dim exactly, so exactly ONE")
log("candidate centre exists and 'broken' is a check, not a search. leg 3")
log("searches anyway, and two defects follow.")
log()
Tb, To = transfer(L["A"], True), transfer(L["A"], False)
gate("4a  the transfer matrix is divided by d twice",
     np.max(np.abs(Tb - To / 4)) < 1e-13,
     "T_published == T_correct/4 to %.1e; logm then shifts every rate rigidly"
     % np.max(np.abs(Tb - To / 4)))
gate("4b  that shift is the value at the foot of all four spectra",
     abs(-np.log(4) / T_EFF + 0.277259) < 1e-6,
     "-ln(4)/5 = %.9f, published -0.277259" % (-np.log(4) / T_EFF))
rng = np.random.default_rng(7)
h = rng.normal(size=8) + 1j * rng.normal(size=8)
pal = np.concatenate([h, -(h + 2 * 0.37)])
_, n_at, _ = best_palindrome(pal, 0.37)
_, n_off, _ = best_palindrome(pal, 0.1)
gate("4c  the search can only score 16/16 at its own seed",
     n_at == 16 and n_off < 16,
     "an EXACTLY palindromic spectrum: %d/16 seeded at its centre, %d/16 at 0.1"
     % (n_at, n_off))

# ==================================================================
log()
log("SECTION 4 -- the control the experiment never ran")
log("-" * 78)
log("The published control eigendecomposes an exactly known 2-qubit")
log("Liouvillian (:409). It never runs the embed / propagate / partial-trace /")
log("logm pipeline it is offered as a control for. Three that do:")
log()


def ext_dephasing(couple_q3):
    H = np.zeros((d, d), dtype=complex)
    for (i, j) in ([(0, 1), (1, 2), (0, 2)] if couple_q3 else [(0, 1)]):
        for P in (sx, sy, sz):
            H = H + site_op(P, i) @ site_op(P, j)
    Lc = build_L_H(H)
    for k in (0, 1):
        Lc = add_jump(Lc, site_op(sz, k), 0.1)
    return Lc


log("  %-38s %13s %13s" % ("marginal of ...", "residual", "best centre"))
marg = {}
for tag, Lop in [("EXTERNAL dephasing, Q3 decoupled", ext_dephasing(False)),
                 ("EXTERNAL dephasing, Q3 coupled", ext_dephasing(True)),
                 ("NO noise at all, H only", build_L_H(build_H(1.0)))] + \
                [("INTERNAL: option " + o, L[o]) for o in "ACD"]:
    r, c = marginal_residual(Lop)
    marg[tag] = r
    log("  %-38s %13.3e %13.6f" % (tag, r, c))
gate("5a  with the spectator DECOUPLED the marginal pairs",
     marg["EXTERNAL dephasing, Q3 decoupled"] < 1e-12,
     "residual %.3e: the pipeline itself is sound"
     % marg["EXTERNAL dephasing, Q3 decoupled"])
internal = [marg["INTERNAL: option " + o] for o in "ACD"]
gate("5b  a system with NO NOISE fails the test as badly",
     marg["NO noise at all, H only"] > min(internal) / 3,
     "no noise %.3e vs internal %.3e ... %.3e"
     % (marg["NO noise at all, H only"], min(internal), max(internal)))
gate("5c  EXTERNAL noise fails it worse than internal noise does",
     marg["EXTERNAL dephasing, Q3 coupled"] > max(internal),
     "external %.3e vs internal max %.3e"
     % (marg["EXTERNAL dephasing, Q3 coupled"], max(internal)))
gate("5d  VERDICT: the marginal test does not separate origin",
     True,
     "what it measures is tracing out a COUPLED spectator, which no reading "
     "of it can turn into a statement about where noise comes from")

# ==================================================================
log()
log("SECTION 5 -- the full three-qubit system, which IS a palindrome")
log("-" * 78)
log("The marginal is not the object the ontology contains; the system is.")
log("Scored at its own forced centre, with no search:")
log()
log("  %-28s %14s %12s %14s" % ("system", "trace(L)/64", "matching", "exact route"))
for o in "ABCD":
    c = float(np.real(np.trace(L[o]) / d2))
    mr = match_resid(np.linalg.eigvals(L[o]), c)
    pr = powersum_resid(L[o], c)
    log("  option %-21s %14.9f %12.2e %14.2e" % (NAMES[o], c, mr, pr))
    if o != "B":
        gate("6a.%s  centre is F137's prediction" % o, abs(c - F137_CENTRE[o]) < 1e-12,
             "trace/64 = %.9f, F137 says %.9f" % (c, F137_CENTRE[o]))
gate("6b  the full L pairs, by the exact route",
     all(powersum_resid(L[o], float(np.real(np.trace(L[o]) / d2))) < 1e-12 for o in "ACD"),
     "relative power-sum residuals: " +
     ", ".join("%s=%.1e" % (o, powersum_resid(L[o], float(np.real(np.trace(L[o]) / d2))))
               for o in "ACD"))
# The two systems the proof's prose leans on, measured here so the numbers have
# a source: a CLOSED generator (which pairs freely, about zero) and a single
# depolarizing site (which does not pair at all, yet still certifies openness).
L_closed = build_L_H(build_H(1.0))
c_closed = float(np.real(np.trace(L_closed) / d2))
r_closed = match_resid(np.linalg.eigvals(L_closed), c_closed)
log("  %-28s %14.9f %12.2e %14s" % ("closed (H only)", c_closed, r_closed, "n/a"))
gate("6d  a closed generator pairs freely, about zero",
     r_closed < 1e-12 and abs(c_closed) < 1e-15,
     "centre %.1e, matching residual %.4e: unitarity, no decay, no content"
     % (c_closed, r_closed))

L_depol = build_L_H(build_H(1.0))
for P_ in (sx, sy, sz):
    L_depol = add_jump(L_depol, site_op(P_, 0), 0.1)
c_depol = float(np.real(np.trace(L_depol) / d2))
ev_depol = np.linalg.eigvals(L_depol)
gg = np.linspace(-2.0, 0.5, 5001)
r_depol = min(match_resid(ev_depol, cc) for cc in gg)
log("  %-28s %14.9f %12.2e %14s" % ("one depolarizing site", c_depol, r_depol, "n/a"))
gate("6e  openness needs no palindrome",
     np.real(np.trace(L_depol)) < -19.0 and r_depol > 1e-2,
     "trace(L) = %.4f, and it does not pair about ANY centre (best %.2e), yet the "
     "negative trace certifies openness on its own"
     % (np.real(np.trace(L_depol)), r_depol))

gate("6c  and the instrument is not vacuous",
     match_resid(np.linalg.eigvals(L["A"]), -0.1) > 1e-2,
     "the same spectrum at the script's seed -0.1: %.3e"
     % match_resid(np.linalg.eigvals(L["A"]), -0.1))

# ==================================================================
log()
log("SECTION 6 -- the trace identity, and what it leaves standing")
log("-" * 78)
rng = np.random.default_rng(11)
worst_alg, worst_neg, worst_id = 0.0, -np.inf, 0.0
for n in (1, 2, 3):
    dim = 2 ** n
    for _ in range(20):
        A = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
        Hh = A + A.conj().T
        worst_alg = max(worst_alg, abs(np.trace(Hh) * dim - dim * np.trace(Hh.T)))
        F = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
        FdF = F.conj().T @ F
        worst_neg = max(worst_neg,
                        float(np.real(abs(np.trace(F)) ** 2 - dim * np.trace(FdF))))
        c = 1.7
        G = c * np.eye(dim, dtype=complex)
        worst_id = max(worst_id, abs(abs(np.trace(G)) ** 2 - dim * np.trace(G.conj().T @ G)))
gate("7a  trace(L_H) = 0 for every H, algebraically",
     worst_alg == 0.0,
     "max |trace(H)*d - d*trace(H^T)| over 60 random Hermitian H = %r" % worst_alg)
gate("7b  every nontrivial jump makes trace(L) strictly negative",
     worst_neg < 0,
     "max of |trace F|^2 - d*trace(F+F) over 60 random jumps = %.3e "
     "(Cauchy-Schwarz on <I,F>)" % worst_neg)
gate("7c  ... with equality exactly for F proportional to I",
     worst_id == 0.0,
     "and such a jump has D identically zero, so it dissipates nothing")

# The hypothesis that makes the chain an equivalence, and what happens without it.
Zq = sz
Xq = sx
L_cp = liouv2(np.zeros((2, 2), dtype=complex), [(Zq, 1.0), (Xq, 1.0)])
L_ncp = liouv2(np.zeros((2, 2), dtype=complex), [(Zq, 1.0), (Xq, -1.0)])
gate("7d  complete positivity is load bearing",
     abs(np.trace(L_ncp)) < 1e-14 and np.linalg.norm(L_ncp) > 1.0,
     "Z at rate +1 with X at rate -1: trace(L) = %.1e while ||L||_F = %.4f, and its "
     "spectrum %s even pairs about zero. It has a POSITIVE eigenvalue, so it is not a "
     "physical channel; with both rates positive the trace is %.1f."
     % (abs(np.trace(L_ncp)), np.linalg.norm(L_ncp),
        np.round(np.real(np.linalg.eigvals(L_ncp)), 6).tolist(),
        float(np.real(np.trace(L_cp)))))

log("""
  For a COMPLETELY POSITIVE generator (every rate >= 0):

      trace(L) = 0  <=>  every jump is a multiple of I  <=>  closed

  equivalently, the system is open exactly when some eigenvalue lies off the
  imaginary axis. That is the whole certificate and it needs no palindrome.
  Where the spectrum DOES pair, its centre is trace(L)/dim, so the pairing is
  the instrument that makes the trace readable off a measured spectrum whose
  generator was never in hand. Instrument, not premise.""")

log()
log("=" * 78)
log("GATES: %d passed, %d failed" % (sum(_gates), len(_gates) - sum(_gates)))
log("VERDICT: %s" % ("GREEN" if all(_gates) else "RED"))
log("""
WHAT THIS ESTABLISHES ABOUT CANDIDATE 2
  Its experiment cannot answer the question it was built for. Two of its three
  readings were code paths. The third measures what a partial trace over a
  COUPLED spectator does to a spectrum, and gives the same verdict for
  external noise, for internal noise, and for no noise at all. Meanwhile the
  object the ontology actually contains, the full three-qubit generator, is an
  exact palindrome at the centre F137 predicts, in every mechanism tested.

  There is also a reason no experiment of this shape could have answered it:
  the "internal" source is modelled as a Lindblad jump, and a Lindblad
  dissipator IS a coupling to an external Markovian bath. The model assumes
  the externality it was meant to test.

WHAT STANDS IN ITS PLACE
  The trace identity of Section 6, which is exact and needs no candidates: a
  nonzero palindrome centre certifies that the generator is not unitary. It
  says OPEN, which is less than EXTERNAL, and the distance between those two
  words is where the proof's remaining work is.
""")
_outf.close()
sys.exit(0 if all(_gates) else 1)
