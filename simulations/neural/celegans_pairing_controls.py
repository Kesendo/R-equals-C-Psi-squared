"""What the C. elegans cavity analysis measured, and what its own model does.

WHAT THE SWEEP RETURNED
-----------------------
The mechanism that breaks the pairing claim is the repository's own, and this
script is its fourth sighting rather than its discovery. Named stores:

  * `docs/ANALYTICAL_FORMULAS.md` (F137): "The centre is an identity and carries
    no evidence ... it is equally well-defined for a spectrum that does not
    pair. ... The pairing is the claim."
  * `experiments/CHAIN_SELECTION_TEST.md` (correction of 2026-08-05): the same
    defect found and repaired on the qubit side, an absolute tolerance against
    level gaps mostly smaller than it, verdict "a measure of how a greedy
    first-fit scrambles in a clustered spectrum". Also
    `experiments/CONCENTRATOR_MAPPING.md`, a saturating score: "these
    percentages carry no information".
  * `docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md`: the degree-matched null ALREADY
    RUN on this animal, on a different instrument, Erdos-Renyi 0.108 against
    degree-preserving rewiring 0.013 = C. elegans's own 0.013, concluding "the
    degree distribution fully explains the palindrome advantage".
  * `docs/CAUGHT_ERRORS.md` returns TWO entries, not nothing: A5 is about this
    very document (a 40 Hz claim asserted as fact against its own Result 1), and
    the greedy-matcher orbit entry is the failure shape STEP 3 measures here.
  * `docs/neural/proofs/PROOF_PALINDROME_NEURAL.md`: condition (a) is selective
    damping, tau_E != tau_I, which this Jacobian does not have.
  * `simulations/neural/` itself returns `hopf_threshold.py`, which integrates
    this same model with solve_ivp and brentq. The Hopf was already understood
    in this folder while the cavity script hand-rolled a Picard iteration; STEP
    7 is that sibling's method applied to the cavity script's own parameters.
  * The OpenArcs registry returns `substrate_q_provenance` (5)(b), which parks
    the neural row's provenance, and `benzene_center_tier_upgrade`, closing with
    "Reopen only if someone finds a route from a trace to a PAIRING".
  * `fw.Confirmations` returns nothing neural; that registry is hardware-only.
  * `docs/GLOSSARY.md` returns the word warnings honoured here: "concentration"
    belongs to the gamma profile (F9), so the quantity below is the Re-spread
    about the centre, and a mode on the centre is self-paired.

WHAT THIS SCRIPT ESTABLISHES
----------------------------
  1. the pairing score is a reading of the tolerance against the spectral scale
     (STEP 1), of the normalisation constant (STEP 2), and of nothing about Dale's law (2b);
  2. the strict matcher is a function of the eigenvalue LIST rather than of the
     spectrum (STEP 3), which is what the 18-unpaired mode list rests on;
  3. the wiring's degeneracy at zero is real and is mostly a matching fact
     rather than a spectral one (STEP 5);
  4. and the model reaches the gamma band at its own committed parameters
     (STEP 7). The cavity document's surviving null said it does not.

Output: simulations/results/celegans_pairing_controls.txt
"""

import json
import numpy as np
from pathlib import Path
import sys, os

if sys.platform == "win32":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    sys.stdout.reconfigure(encoding="utf-8")

NEURAL_DIR = Path(__file__).parent
RESULTS_DIR = NEURAL_DIR.parent / "results"

SEED = 20260825
TOL = 0.01          # the committed matching tolerance, both scripts
F_PRIME = 0.3       # the committed linearisation slope, neural_gamma_cavity.py:215
PRIMES = (2147483647, 2147483629, 104729)
R_ENS = 200         # ensemble size for every null below

out = []


def log(msg=""):
    print(msg)
    out.append(msg)


failures = []


def gate(name, ok, detail):
    if not ok:
        failures.append(name)
    log(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")


# -------------------------------------------------------------------
# The two matchers, reproduced verbatim (vectorised, same answers)
# -------------------------------------------------------------------

def loose_pairing(eig, tol=TOL):
    """neural_gamma_cavity.py:248-254. Self-pairing is NOT excluded.

    Reproduced including the original's quirk that the target's imaginary part
    is +ev.imag: the -ev.imag computed into partner_target is then discarded.
    """
    center = np.mean(eig.real)
    targets = (2 * center - eig.real) + 1j * eig.imag
    d = np.abs(eig[None, :] - targets[:, None])
    return int(np.sum(d.min(axis=1) < tol))


def exclusive_pairing(eig, tol=TOL):
    """neural_gamma_cavity_unpaired.py:83-109. Forbids j == i.

    Greedy and first-index-biased, so it is a function of the eigenvalue LIST,
    not of the spectrum. STEP 3 measures how much that costs.
    """
    n = len(eig)
    center = np.mean(eig.real)
    re, im = eig.real, np.abs(eig.imag)
    free = np.ones(n, dtype=bool)
    pairs, unpaired = [], []
    for i in range(n):
        if not free[i]:
            continue
        d = np.hypot(re - (2 * center - re[i]), im - im[i])
        d[i] = np.inf
        d[~free] = np.inf
        j = int(np.argmin(d))
        if d[j] < tol:
            pairs.append((i, j))
            free[i] = False
            free[j] = False
        else:
            unpaired.append(i)
            free[i] = False
    return pairs, unpaired


def self_paired_count(eig, tol=TOL):
    """Modes on the centre, which are their own reflection.

    |ev - ((2c - Re) + i Im)| = 2|Re - c|, so this is |Re - c| < tol/2 written
    the way the matcher writes it. A definition, not a measurement.
    """
    return int(np.sum(np.abs(eig.real - np.mean(eig.real)) < tol / 2))


def jacobian(W, f_prime=F_PRIME):
    """The committed C. elegans construction, neural_gamma_cavity.py:210-216.

    Note what is NOT here: no tau. The committed line is a bare -I. Reading its
    unit damping as 1/tau with tau = 1 ms is an interpretation; tau_E = tau_I =
    1.0 at :48-49 belongs to the Wilson-Cowan block in the same file.
    """
    m = np.max(np.abs(W))
    return -np.eye(W.shape[0]) + f_prime * (W / m if m > 0 else W)


# -------------------------------------------------------------------
# Controls
# -------------------------------------------------------------------

def degree_preserving_swap(W_raw, rng, passes=10):
    """Directed double-edge swap on the unsigned matrix.

    Preserves every in-degree and out-degree exactly, and destroys only which
    neuron connects to which. Dale's law is applied afterwards (a per-
    presynaptic-neuron sign, so it follows neuron identity rather than the
    edge) and is preserved exactly too. STEP 9 checks the mixing rather than
    asserting it.
    """
    W = W_raw.copy()
    edges = [(int(i), int(j)) for i, j in zip(*np.nonzero(W))]
    present = set(edges)
    n_edges = len(edges)
    done = 0
    for _ in range(passes * n_edges):
        p = int(rng.integers(n_edges))
        q = int(rng.integers(n_edges))
        a, b = edges[p]
        c, d = edges[q]
        if len({a, b, c, d}) < 4:
            continue
        if (a, d) in present or (c, b) in present:
            continue
        W[a, d] = W[a, b]
        W[c, b] = W[c, d]
        W[a, b] = 0.0
        W[c, d] = 0.0
        present.discard((a, b))
        present.discard((c, d))
        present.add((a, d))
        present.add((c, b))
        edges[p] = (a, d)
        edges[q] = (c, b)
        done += 1
    return W, done


def spread_matched(eig, target_sd):
    """Rescale a spectrum about its centre to a chosen Re-spread.

    The matcher can only see the tolerance against the spectral scale, and the
    committed construction sets that scale by dividing by max|W|, one synapse
    count. Matching the spread is the only honest way to place a control.
    """
    c = np.mean(eig.real)
    return c + (eig - c) * (target_sd / np.std(eig.real))


def gf_rank(M, p):
    """Exact rank of an integer matrix over GF(p). No tolerance, no eigensolver."""
    A = (np.rint(M).astype(np.int64) % p)
    rows, cols = A.shape
    r = 0
    for col in range(cols):
        nz = np.nonzero(A[r:, col])[0]
        if nz.size == 0:
            continue
        piv = r + int(nz[0])
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        A[r] = (A[r] * pow(int(A[r, col]), p - 2, p)) % p
        below = np.nonzero(A[r + 1:, col])[0]
        if below.size:
            idx = r + 1 + below
            A[idx] = (A[idx] - np.outer(A[idx, col], A[r])) % p
        r += 1
        if r == rows:
            break
    return r


def zero_multiplicity(M, p):
    """Algebraic multiplicity of the eigenvalue 0, exactly, by the rank chain.

    nullity(M^k) rises until the Jordan structure is exhausted, and where it
    stops is the algebraic multiplicity. That is the count of modes the wiring
    does nothing to, which is what the self-paired column counts with a
    tolerance. The nullity alone (k = 1) counts eigenVECTORS and is a smaller,
    different number.
    """
    n = M.shape[0]
    A = (np.rint(M).astype(np.int64) % p)
    prev, P = -1, A.copy()
    for _ in range(12):
        cur = n - gf_rank(P, p)
        if cur == prev:
            return cur
        prev = cur
        P = (P @ A) % p
    return prev


def structural_rank(M):
    """Term rank: the largest set of nonzeros with distinct rows and columns.

    A property of the zero PATTERN alone, decided by bipartite matching with no
    field and no weights, and an upper bound on the true rank. The part of a
    rank deficit it already explains is combinatorial, not numerical.
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import maximum_bipartite_matching
    m = maximum_bipartite_matching(csr_matrix((M != 0).astype(np.int8)),
                                   perm_type='column')
    return int(np.sum(m >= 0))


# -------------------------------------------------------------------

log("=" * 75)
log("WHAT THE C. ELEGANS CAVITY ANALYSIS MEASURED")
log("=" * 75)
log()
log(f"Seed {SEED}; tolerance {TOL} and f' = {F_PRIME} are the committed values.")
log()

cdata = json.loads((NEURAL_DIR / "celegans_connectome.json").read_text())
W_raw = np.array(cdata['chemical'], dtype=float)
sign = np.array(cdata['chemical_sign'], dtype=float)
N = W_raw.shape[0]


def apply_dale(W):
    W = W.copy()
    for i in range(N):
        if sign[i] < 0:
            W[i, :] *= -1
    return W


W_real = apply_dale(W_raw)
rng = np.random.default_rng(SEED)
eig_real = np.linalg.eigvals(jacobian(W_real))
sd_real = float(np.std(eig_real.real))

log("-" * 75)
log("STEP 0: the centre the test reflects across is forced, not measured")
log("-" * 75)
log()
log("F137 states this in general: the centre is trace/dim, an identity, equally")
log("well defined for a spectrum that does not pair. NEURAL_CLOCK_TWO_HANDS:42-52")
log("proves the neural case is wiring-INDEPENDENT (its own value is -0.150000 at")
log("its own time constants; it is -1 here only because this Jacobian's damping")
log("is 1). Here the wiring-independence is exact rather than approximate.")
log()
n_selfsyn = int(np.count_nonzero(np.diag(W_raw)))
log(f"  neurons with a synapse onto themselves: {n_selfsyn}   (exact count)")
gate("G0 the wiring never touches the diagonal", n_selfsyn == 0,
     "so trace(J)/N = -1 exactly, for this connectome and for every rewiring "
     "of it")
log(f"  eigensolver mean Re(lambda) = {np.mean(eig_real.real):.12f}   (read; no")
log("    exact route through an eigensolver, the exact statement is G0)")
log()

log("-" * 75)
log("STEP 1: the score is a reading of the tolerance, and it was never swept")
log("-" * 75)
log()
gaps = np.diff(np.sort(eig_real.real))
re_width = float(eig_real.real.max() - eig_real.real.min())
log(f"  Re window                            : [{eig_real.real.min():.4f}, "
    f"{eig_real.real.max():.4f}], width {re_width:.4f}")
log(f"  mean nearest-neighbour spacing in Re : {np.mean(gaps):.6f}")
log(f"  committed tolerance                  : {TOL}")
log(f"  so each eigenvalue has of order {TOL / np.mean(gaps):.0f} candidate "
    "partners inside")
log("  the window before any structure is consulted.")
log()
log(f"  {'tol':>9s} {'loose':>8s} {'strict':>8s} {'self-paired':>12s}")
log("  " + "-" * 42)
sweep = {}
for tol in (0.1, TOL, 0.003, 0.001, 3e-4, 1e-4):
    lo = loose_pairing(eig_real, tol) / N * 100
    st = len(exclusive_pairing(eig_real, tol)[0]) * 2 / N * 100
    sf = self_paired_count(eig_real, tol) / N * 100
    sweep[tol] = (lo, st, sf)
    log(f"  {tol:9.4f} {lo:7.1f}% {st:7.1f}% {sf:11.1f}%")
log()
lo_fine, st_fine, sf_fine = sweep[1e-4]
log(f"  At the finest tolerance {lo_fine:.1f}% of modes still find a partner and")
log(f"  {sf_fine:.1f}% do it by being their own. The difference is a LOWER bound")
log("  on genuine cross-pairing, not an estimate: a self-paired mode that also")
log("  has a real partner is charged to the self column. STEP 4 puts a null")
log("  beside it, which is the only way to read a residue.")
log()
gate("G1 the committed 97.3 % reproduces", abs(sweep[TOL][0] - 97.3) < 0.2,
     f"at the committed tolerance the loose matcher gives {sweep[TOL][0]:.1f}%")
gate("G2 the score collapses as the tolerance leaves the clustering scale",
     lo_fine < 40.0,
     f"loose falls {sweep[TOL][0]:.1f}% -> {lo_fine:.1f}% between tol {TOL} "
     f"and {1e-4}")
log()

log("-" * 75)
log("STEP 2: the only thing the matcher can see is tolerance over Re-spread")
log("-" * 75)
log()
log("The committed construction fixes the scale by dividing W by max|W|, one")
log("synapse count. A control at a different scale is therefore not a control.")
log()
eigW = np.linalg.eigvals(W_real)
nz = W_real[W_real != 0]
normalisations = [
    ("max|w| (committed)", float(np.max(np.abs(W_real)))),
    ("spectral radius", float(np.max(np.abs(eigW)))),
    ("max row sum", float(np.max(np.abs(W_real).sum(axis=1)))),
    ("95th-percentile weight", float(np.percentile(np.abs(nz), 95))),
]
log(f"  {'normalisation by':26s} {'value':>9s} {'loose':>8s} {'sd Re':>8s}")
log("  " + "-" * 54)
div_scores = {}
for label, dv in normalisations:
    e = np.linalg.eigvals(-np.eye(N) + F_PRIME * (W_real / dv))
    div_scores[label] = loose_pairing(e) / N * 100
    log(f"  {label:26s} {dv:9.1f} {div_scores[label]:7.1f}% {np.std(e.real):8.4f}")
e_bin = np.linalg.eigvals(-np.eye(N) + F_PRIME * np.sign(W_real))
div_scores["binary, unweighted"] = loose_pairing(e_bin) / N * 100
log(f"  {'binary, unweighted':26s} {'-':>9s} "
    f"{div_scores['binary, unweighted']:7.1f}% {np.std(e_bin.real):8.4f}")
log()
log(f"  {len(nz)} directed edges, weight median {np.median(np.abs(nz)):.0f} and "
    f"max {np.max(np.abs(nz)):.0f}:")
log("  the committed normalisation is by one outlier synapse.")
gate("G3 the score is set by the normalisation constant, which nothing justifies",
     max(div_scores.values()) - min(div_scores.values()) > 50.0,
     f"the same connectome and the same matcher give "
     f"{min(div_scores.values()):.1f}% to {max(div_scores.values()):.1f}% "
     "across five defensible normalisations")
log()
log("  Matched on Re-spread instead, which is the only quantity the matcher")
log("  responds to, noise with no structure at all is not beaten:")
log()
dens = np.count_nonzero(W_raw) / W_raw.size
matched = {}
for label, M in [("dense iid gaussian", rng.normal(size=(N, N))),
                 ("sparse iid (matched density)",
                  rng.normal(size=(N, N)) * (rng.random((N, N)) < dens))]:
    e = np.linalg.eigvals(jacobian(M))
    matched[label] = loose_pairing(spread_matched(e, sd_real)) / N * 100
    log(f"    {label:30s} {matched[label]:6.1f}%")
log(f"    {'real connectome':30s} {sweep[TOL][0]:6.1f}%")
gate("G4 spread-matched noise is not beaten by the connectome",
     min(matched.values()) >= sweep[TOL][0],
     f"the weakest matched control scores {min(matched.values()):.1f}% against "
     f"the connectome's {sweep[TOL][0]:.1f}%")
log()

log("-" * 75)
log("STEP 2b: Dale's law makes no difference to the score")
log("-" * 75)
log()
log("NEURAL_GAMMA_CAVITY's thesis was that the excitatory/inhibitory")
log("classification 'creates the same SWAP structure that the Pi operator")
log("creates in the qubit chain'. That is a claim about the signs, so it is")
log("tested by changing them.")
log()
n_inh = int(np.sum(sign < 0))
r1 = np.ones(N)
r1[rng.choice(N, n_inh, replace=False)] = -1
variants = [
    (f"committed Dale ({n_inh} inhibitory)", sign),
    ("no Dale at all, every neuron excitatory", np.ones(N)),
    (f"{n_inh} inhibitory, chosen at random", r1),
    ("50/50 random signs", rng.choice([-1.0, 1.0], N)),
]
log(f"  {'sign assignment':42s} {'loose':>8s}")
log("  " + "-" * 52)
dale_scores = {}
for label, sg_ in variants:
    Wv = W_raw.copy()
    for i in range(N):
        if sg_[i] < 0:
            Wv[i, :] *= -1
    dale_scores[label] = loose_pairing(np.linalg.eigvals(jacobian(Wv))) / N * 100
    log(f"  {label:42s} {dale_scores[label]:7.1f}%")
log()
spread_dale = max(dale_scores.values()) - min(dale_scores.values())
log(f"  The four differ by {spread_dale:.1f} points, "
    f"{spread_dale / 100 * N:.0f} modes out of {N}. The claim to make is that")
log("  the signs do NOTHING, which is what kills the thesis. Saying the score")
log("  RISES without them would be a directional claim on that same noise.")
gate("G5 the signs make no difference either way", spread_dale < 3.0,
     f"every sign assignment lands within {spread_dale:.1f} points of every "
     "other, deleting Dale's law included")
log()

log("-" * 75)
log("STEP 3: the strict matcher is a function of the LIST, not of the spectrum")
log("-" * 75)
log()
log("It scans indices in order and takes the first best unused partner, so")
log("permuting the same eigenvalues changes its answer. The committed 18-")
log("unpaired mode list and its three biological categories are built on it.")
log()
native_unp = len(exclusive_pairing(eig_real)[1])
perm_unp = [len(exclusive_pairing(eig_real[rng.permutation(N)])[1])
            for _ in range(50)]
log(f"  LAPACK's own ordering                  : {native_unp} unpaired")
log(f"  50 permutations of the same eigenvalues: mean {np.mean(perm_unp):.1f}, "
    f"sd {np.std(perm_unp, ddof=1):.2f}, range {min(perm_unp)}-{max(perm_unp)}")
log()
log("  For scale: the committed run reports 18")
log(f"  (results/neural_gamma_cavity_unpaired.txt:13) and the same script on the")
log(f"  same data reports {native_unp} today.")
log()
gate("G6 ordering moves the count further than the run-to-run drift does",
     np.std(perm_unp, ddof=1) > 1.0 and abs(np.mean(perm_unp) - native_unp) > 2,
     f"reordering shifts the mean to {np.mean(perm_unp):.1f} from {native_unp}, "
     "while the committed 18 and today's 20 differ by 2")
log()
log("  A second reason is often given for distrusting per-mode readings here,")
log("  and it has to be scoped rather than asserted globally. The Jacobian is")
log("  non-normal, but that is manufactured by the defective cluster ON the")
log("  centre; the per-eigenvalue condition numbers say the rest is fine.")
log()
try:
    import scipy.linalg as sla
    kappa_max = float('nan')
    eig_error_bound = float('nan')
    _ev, _VL, _VR = sla.eig(jacobian(W_real), left=True, right=True)
    kappa = 1.0 / np.abs(np.sum(np.conj(_VL) * _VR, axis=0))
    off = np.abs(_ev + 1) > 1e-4
    k_off = np.sort(kappa[off])
    normJ = float(np.linalg.norm(jacobian(W_real), 2))
    motion = float(np.median(k_off)) * float(np.finfo(float).eps) * normJ
    log(f"  eigenvalues on the centre to 1e-4 : {int(np.sum(~off))}")
    log(f"  condition numbers of the other {int(np.sum(off)):3d} : median "
        f"{np.median(k_off):.0f}, 90th pct {np.percentile(k_off, 90):.0f}")
    log(f"  median induced eigenvalue motion  : {motion:.2e}, against a mean "
        f"spacing of {np.mean(gaps):.2e}")
    kappa_max = float(np.max(kappa))
    eig_error_bound = kappa_max * float(np.finfo(float).eps) * normJ
    log(f"  worst-conditioned eigenvalue      : kappa = {kappa_max:.1e}, so no")
    log(f"    eigenvalue of this matrix is resolved better than "
        f"{eig_error_bound:.1e}")
    gate("G7 the off-centre eigenvalues ARE well resolved",
         motion < np.mean(gaps) / 1e6,
         "so it is the ordering result above, not a conditioning argument, "
         "that carries the withdrawal of the mode list")
except ImportError:
    log("  (scipy unavailable; per-eigenvalue conditioning not computed)")
log()

log("-" * 75)
log("STEP 4: against degree-matched rewiring, by rank rather than by z")
log("-" * 75)
log()
log("ALGEBRAIC_PALINDROME_NEURAL:246-252 already ran this null on the algebraic")
log("residual and found it identical. Here it is at the eigenvalue level. The")
log("statistics are bounded discrete counts, so the read is the rank inside the")
log("ensemble. Every test is one-sided and the statistics were chosen after")
log("looking, so no multiplicity correction is applied and none is claimed.")
log()
ens = {k: [] for k in ("loose", "strict", "self-paired", "cross-fine",
                       "nullity", "multiplicity", "structural nullity")}
for _ in range(R_ENS):
    Wr_raw, _ = degree_preserving_swap(W_raw, rng)
    Wr = apply_dale(Wr_raw)
    er = np.linalg.eigvals(jacobian(Wr))
    ens["loose"].append(loose_pairing(er) / N * 100)
    ens["strict"].append(len(exclusive_pairing(er)[0]) * 2 / N * 100)
    ens["self-paired"].append(self_paired_count(er) / N * 100)
    ens["cross-fine"].append((loose_pairing(er, 1e-4)
                              - self_paired_count(er, 1e-4)) / N * 100)
    ens["nullity"].append(N - gf_rank(Wr_raw, PRIMES[0]))
    ens["multiplicity"].append(zero_multiplicity(Wr_raw, PRIMES[0]))
    ens["structural nullity"].append(N - structural_rank(Wr_raw))

real_vals = {
    "loose": sweep[TOL][0],
    "strict": sweep[TOL][1],
    "self-paired": sweep[TOL][2],
    "cross-fine": lo_fine - sf_fine,
    "nullity": N - gf_rank(W_raw, PRIMES[0]),
    "multiplicity": zero_multiplicity(W_raw, PRIMES[0]),
    "structural nullity": N - structural_rank(W_raw),
}

log(f"  R = {R_ENS} degree-matched rewirings; p is the one-sided rank p-value")
log("  (R - #{null <= real} + 1)/(R + 1).")
log()
log(f"  {'statistic':20s} {'real':>8s} {'null mean':>10s} {'null range':>13s} "
    f"{'p':>8s}")
log("  " + "-" * 64)
pvals = {}
for k in ("loose", "strict", "cross-fine", "self-paired", "nullity",
          "multiplicity", "structural nullity"):
    a = np.array(ens[k], dtype=float)
    ge = int(np.sum(real_vals[k] >= a))
    pvals[k] = (R_ENS - ge + 1) / (R_ENS + 1)
    log(f"  {k:20s} {real_vals[k]:8.1f} {np.mean(a):10.1f} "
        f"{f'{np.min(a):.1f}-{np.max(a):.1f}':>13s} {pvals[k]:8.3f}")
log()
gate("G8 the strict pairing score is ordinary against degree-matched rewiring",
     pvals["strict"] > 0.05,
     f"p = {pvals['strict']:.3f}; the document reads this matcher's unpaired "
     "modes as connectome structure")
n_tests = 7
bonf = 0.05 / n_tests
gate("G9 the tolerance-free residue is small and not established",
     (real_vals["cross-fine"] - float(np.mean(ens["cross-fine"]))) < 3.0
     and pvals["cross-fine"] > bonf,
     f"cross-pairing at tol 1e-4 is {real_vals['cross-fine']:.1f}% against a "
     f"null mean of {np.mean(ens['cross-fine']):.1f}%, an excess of "
     f"{real_vals['cross-fine'] - float(np.mean(ens['cross-fine'])):.1f} points "
     f"({(real_vals['cross-fine'] - float(np.mean(ens['cross-fine']))) / 100 * N:.0f} "
     f"modes) at p = {pvals['cross-fine']:.3f}, which does not clear "
     f"{bonf:.4f}, the level {n_tests} one-sided tests need. It is the one "
     "quantity here that is above its null at all, and it is not enough to "
     "carry a claim")
log()

log("-" * 75)
log("STEP 5: what IS different about this connectome, and what kind of fact")
log("-" * 75)
log()
log("The self-paired modes are eigenvalues of f'*W_norm at zero: modes the")
log("wiring does nothing to. Counting them needs a tolerance and an eigensolver.")
log("The synapse matrix is an integer matrix, so the same content is exact. Two")
log("different exact objects live here and only the second is the mode count:")
log()
log(f"  nullity(W)        = {real_vals['nullity']}   "
    "(geometric: independent eigenvectors at 0)")
log(f"  multiplicity of 0 = {real_vals['multiplicity']}   "
    "(algebraic, by the rank chain: the MODES)")
n_float = int(np.sum(np.abs(eig_real + 1) < 1e-4))
log(f"  float check, |lambda + 1| < 1e-4 : {n_float}")
gate("G10 the exact multiplicity matches the float count",
     real_vals['multiplicity'] == n_float,
     f"{real_vals['multiplicity']} == {n_float}, so the rank chain counts what "
     "the tolerance was counting, exactly and without one")
ranks_agree = len({gf_rank(W_raw, p) for p in PRIMES})
gate("G11 the exact rank agrees at every prime tried", ranks_agree == 1,
     f"{len(PRIMES)} primes, {ranks_agree} distinct value; rank over GF(p) "
     "bounds the rational rank from below, so agreement is the check")
nullity_signed = N - gf_rank(W_real, PRIMES[0])
gate("G12 the signs are irrelevant here, as row scaling must be",
     nullity_signed == real_vals["nullity"],
     f"nullity is {real_vals['nullity']} unsigned and {nullity_signed} with "
     "Dale's law applied, so the word 'signed' does not belong in this claim")
log()
log("  And most of the deficit is not spectral. The structural rank is decided")
log("  by bipartite matching on the zero pattern alone, no field, no weights,")
log("  no signs, and it upper-bounds the true rank:")
log()
excess_total = real_vals["nullity"] - float(np.mean(ens["nullity"]))
excess_struct = (real_vals["structural nullity"]
                 - float(np.mean(ens["structural nullity"])))
log(f"  {'':22s} {'real':>8s} {'null mean':>10s} {'excess':>8s}")
log(f"  {'structural nullity':22s} {real_vals['structural nullity']:8.0f} "
    f"{np.mean(ens['structural nullity']):10.2f} {excess_struct:8.2f}")
log(f"  {'exact nullity':22s} {real_vals['nullity']:8.0f} "
    f"{np.mean(ens['nullity']):10.2f} {excess_total:8.2f}")
log()
log(f"  Of the {excess_total:.1f}-point excess, {excess_struct:.1f} is already in")
log("  the zero pattern. The honest name for most of this is a small maximum")
log("  matching, not a numerical degeneracy.")
gate("G13 the connectome is more degenerate at zero than degree-matched rewiring",
     pvals["multiplicity"] <= 0.05,
     f"multiplicity {real_vals['multiplicity']} against a null mean of "
     f"{np.mean(ens['multiplicity']):.1f} (range "
     f"{min(ens['multiplicity'])}-{max(ens['multiplicity'])}), p = "
     f"{pvals['multiplicity']:.3f}")
log()

log("-" * 75)
log("STEP 6: Q_max is arithmetic on the damping, not a measured loss")
log("-" * 75)
log()
log("NEURAL_GAMMA_CAVITY:157 reads Q_max = 0.1 against the qubit cavity's 68 to")
log("75 as 'the biological cavity is extremely lossy'. The eigenvalues of")
log("-I/tau + f'*W_n are EXACTLY -1/tau + f'*mu, so Q is a function of the graph")
log("spectrum and the two chosen constants and of nothing else. That identity is")
log("an exact route, so it is compared rather than tolerated.")
log()
Wn = W_real / np.max(np.abs(W_real))
mu = np.linalg.eigvals(Wn)
log(f"  {'tau (ms)':>9s} {'f_prime':>8s} {'centre':>9s} {'Q_max':>9s}")
log("  " + "-" * 40)
worst = 0.0
q_vals = {}
for tau in (1.0, 2.0, 5.0, 10.0):
    for fp in (0.3, 0.6):
        A = -np.eye(N) / tau + fp * Wn
        lam_direct = np.sort_complex(np.linalg.eigvals(A))
        lam_identity = np.sort_complex(-1.0 / tau + fp * mu)
        worst = max(worst, float(np.max(np.abs(lam_direct - lam_identity))))
        m = np.abs(lam_direct.imag) > 1e-6
        q = float(np.max(np.abs(lam_direct[m].imag)
                         / np.abs(lam_direct[m].real))) if m.any() else 0.0
        q_vals[(tau, fp)] = q
        log(f"  {tau:9.1f} {fp:8.2f} {-1.0 / tau:9.4f} {q:9.3f}")
log()
log("  THIS SECTION IS A READ, NOT A GATE, and the reason is worth stating.")
log("  Q_max = f'|Im mu| / |1/tau - f' Re mu| is an identity, exact in exact")
log("  arithmetic. Both sides of any numerical check of it are eigensolver")
log("  output, so there is no exact route to compare along; and a bound built")
log("  from this matrix's worst conditioning is vacuous, because the defective")
log("  cluster on the centre pushes it past 1e40. A gate here would pass on")
log("  anything. What the identity says needs no gate:")
log()
log(f"  worst deviation, assembled vs -1/tau + f'*mu : {worst:.2e} (read)")
log(f"  Q_max across the eight rows                  : factor "
    f"{max(q_vals.values()) / min(q_vals.values()):.0f}")
log()
log("  Q_max is therefore arithmetic on two chosen constants and the graph")
log("  spectrum. The provenance is the finding, and it is checkable rather than")
log("  statistical: tau = 1 ms is stipulated bare at neural_gamma_cavity.py:48-49")
log("  for the Wilson-Cowan block, and the C. elegans Jacobian at :216 contains")
log("  no tau at all, only a bare -I. G0 above already pins that diagonal.")
log()

log("-" * 75)
log("STEP 7: the model reaches the gamma band at its own committed parameters")
log("-" * 75)
log()
log("neural_gamma_cavity.py:59 finds fixed points by Picard iteration, which")
log("converges only where the point is stable AS A MAP. That is not stability of")
log("the ODE; here it lands on the quiescent branch and reports its ringing. The")
log("sibling hopf_threshold.py in this same folder already integrates this model")
log("properly. Doing that, with NOTHING changed:")
log()
wEE, wEI, wIE, wII, alpha_, theta_ = 16.0, 12.0, 15.0, 3.0, 1.3, 4.0
sg = lambda x: 1.0 / (1.0 + np.exp(-alpha_ * (x - theta_)))


def picard(I_ext, n=1000):
    E, Iv = 0.5, 0.5
    for _ in range(n):
        E, Iv = sg(wEE * E - wEI * Iv + I_ext), sg(wIE * E - wII * Iv)
    return np.array([E, Iv])


def fp_residual(v, I_ext):
    E, Iv = v
    return np.array([sg(wEE * E - wEI * Iv + I_ext) - E,
                     sg(wIE * E - wII * Iv) - Iv])


try:
    from scipy.integrate import solve_ivp

    def rhs(t, y, I_ext):
        E, Iv = y
        return [-E + sg(wEE * E - wEI * Iv + I_ext),
                -Iv + sg(wIE * E - wII * Iv)]

    log(f"  {'I_ext':>6s} {'amplitude':>10s} {'freq (Hz)':>10s} {'band':>10s}")
    log("  " + "-" * 40)
    bands = {}
    for I_ext in (0.5, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0):
        s = solve_ivp(rhs, [0, 2000], [0.3, 0.3], args=(I_ext,), max_step=0.05,
                      dense_output=True, rtol=1e-9, atol=1e-11)
        t = np.linspace(1500, 2000, 200000)
        E = s.sol(t)[0]
        amp = float(E.max() - E.min())
        if amp < 1e-4:
            log(f"  {I_ext:6.1f} {amp:10.4f} {'-':>10s} {'fixed pt':>10s}")
            continue
        m = E - E.mean()
        zc = np.where((m[:-1] < 0) & (m[1:] >= 0))[0]
        f_hz = 1000.0 / float(np.mean(np.diff(t[zc]))) if len(zc) > 2 else np.nan
        band = ("gamma" if 30 <= f_hz <= 100
                else "alpha" if 8 <= f_hz <= 13
                else "above" if f_hz > 100 else "other")
        bands[I_ext] = (f_hz, band)
        log(f"  {I_ext:6.1f} {amp:10.4f} {f_hz:10.1f} {band:>10s}")
    log()
    in_gamma = [I for I, (f_hz, b) in bands.items() if b == "gamma"]
    gate("G15 the gamma band is reached with no parameter change",
         len(in_gamma) > 0,
         f"limit cycles inside 30-100 Hz at I_ext = "
         f"{', '.join(f'{I:g}' for I in in_gamma)}, at the same w, alpha, theta "
         "and tau the document calls insufficient")
    log()
    log("  What the committed script reports instead is a linearisation at the")
    log("  quiescent branch, which is where its iteration converges:")
    v = picard(2.0)
    res2 = float(np.linalg.norm(fp_residual(v, 2.0)))
    log(f"    at I_ext = 2.0 the iteration returns ({v[0]:.3f}, {v[1]:.3f}), "
        f"residual {res2:.4f}")
    gate("G16 the reported fixed points are not fixed points", res2 > 0.01,
         "a residual that size against a sigmoid whose whole range is 1 means "
         "the point was never found, so every stability reading taken there is "
         "undetermined")
except ImportError:
    log("  (scipy unavailable; the limit cycle is not integrated here)")
log()

log("-" * 75)
log("STEP 8: the isolated second cavity is an artifact of a discarded matrix")
log("-" * 75)
log()
PHARYNX = {'I1L', 'I1R', 'I2L', 'I2R', 'I3', 'I4', 'I5', 'I6', 'M1', 'M2L',
           'M2R', 'M3L', 'M3R', 'M4', 'M5', 'MCL', 'MCR', 'NSML', 'NSMR', 'MI'}
ids = []
with open(NEURAL_DIR / "celegans_neuron_ids.txt", encoding="utf-8") as fh:
    for line in fh:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            ids.append(parts[1])
Elec = np.array(cdata['electrical'], dtype=float)
pi = [i for i, n in enumerate(ids) if n in PHARYNX]
si = [i for i, n in enumerate(ids) if n not in PHARYNX]
chem_cross = int(np.count_nonzero(W_raw[np.ix_(pi, si)])
                 + np.count_nonzero(W_raw[np.ix_(si, pi)]))
elec_cross = int(np.count_nonzero(Elec[np.ix_(pi, si)])
                 + np.count_nonzero(Elec[np.ix_(si, pi)]))
log(f"  pharyngeal neurons identified  : {len(pi)} of {N}")
log(f"  electrical entries in the file : {int(np.count_nonzero(Elec))}, "
    "none of them used by the Jacobian")
log(f"  pharynx <-> soma, chemical     : {chem_cross}")
log(f"  pharynx <-> soma, electrical   : {elec_cross}")
for a in pi:
    for b in si:
        if Elec[a, b]:
            log(f"      {ids[a]} -- {ids[b]}")
gate("G17 the boundary reported as uncoupled is coupled in the discarded matrix",
     chem_cross == 0 and elec_cross > 0,
     f"{chem_cross} chemical and {elec_cross} electrical connections across it")
log()

log("-" * 75)
log("STEP 9: the null model's mixing, checked rather than asserted")
log("-" * 75)
log()
log(f"  {'passes':>8s} {'nullity mean':>14s} {'sd':>7s}")
log("  " + "-" * 32)
mixing = {}
for passes in (10, 100, 500):
    vals = [N - gf_rank(degree_preserving_swap(W_raw, rng, passes=passes)[0],
                        PRIMES[0]) for _ in range(12)]
    mixing[passes] = (float(np.mean(vals)), float(np.std(vals, ddof=1)))
    log(f"  {passes:8d} {mixing[passes][0]:14.2f} {mixing[passes][1]:7.2f}")
means = [m for m, _ in mixing.values()]
sds = [s for _, s in mixing.values()]
gate("G18 passes = 10 is already mixed",
     max(means) - min(means) < 2 * max(sds),
     f"the three mixing levels differ by {max(means) - min(means):.2f} against "
     f"a within-level sd of up to {max(sds):.2f}")
log()

log("=" * 75)
n_gates = 18
log(f"{len(failures)} of {n_gates} gates FAILED" if failures
    else f"ALL {n_gates} GATES PASS")
for f in failures:
    log(f"  FAILED: {f}")
log("=" * 75)
log()
log("WHAT THIS DOES NOT SHOW")
log("-----------------------")
log("It does not test condition (b) of PROOF_PALINDROME_NEURAL, where that")
log("proof's content sits. It does not identify WHICH modes go unpaired, and")
log("STEP 3 argues that question is not well posed as the committed script asks")
log("it. STEP 5's excess is measured against a null that does not hold the")
log("maximum matching fixed, so the residual numerical part is an upper bound; a")
log("matching-conditioned null would decide it. STEP 7 integrates one trajectory")
log("per input from one initial condition, so it finds the attractor that basin")
log("reaches and would miss a coexisting one.")

RESULTS_DIR.mkdir(exist_ok=True)
(RESULTS_DIR / "celegans_pairing_controls.txt").write_text('\n'.join(out),
                                                           encoding='utf-8')
print(f"\n>>> Results saved to: {RESULTS_DIR / 'celegans_pairing_controls.txt'}")
sys.exit(1 if failures else 0)
