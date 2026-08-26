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
  * `docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md`: the degree-null ALREADY RUN on this
    animal, on a different instrument. Read with care: as of 2026-08-26 BOTH
    rows of that table are withdrawn. Its residual has the closed form
    sqrt(2)*||W_eff||/||J|| on blocks this sparse, so it reads coupling
    magnitude and not wiring; the Erdos-Renyi arm compared two normalisation
    constants, and the degree-preserving arm cannot move under such a metric and
    mostly does not even run (mean 2 edges per block). The VERDICT there, no
    palindromic advantage, stands; its evidence does not.
  * `docs/CAUGHT_ERRORS.md` returns FOUR dated entries, not nothing: A5 is about this
    very document (a 40 Hz claim asserted as fact against its own Result 1), and
    the greedy-matcher orbit entry is the failure shape STEP 3 measures here.
  * `docs/neural/proofs/PROOF_PALINDROME_NEURAL.md`: condition (a) asks that Q
    send every seat to one of the OPPOSITE type, which holds at any tau and is
    automatic at uniform tau, so it imposes nothing here and (b) is the whole
    hypothesis. G0b below decides (b) on a count. Selective damping tau_E !=
    tau_I does not enable the pairing; it is what makes (a) bite on Q.
  * `simulations/neural/` itself returns `hopf_threshold.py`, which bisects on
    max Re lambda rather than reading stability off a point an iteration
    returned, and which carries tau_E != tau_I. Checked, not assumed: it is NOT
    the same model (random balanced N = 100..5000, a_E != a_I, theta_E != theta_I)
    and it hand-rolls a 3000-step Picard iteration of its own, so Picard is not
    the difference; its brentq is imported and never called. STEP 7 takes the
    habit of root-finding on the eigenvalue, not the script.
  * The OpenArcs registry returns `substrate_q_provenance` (5)(b), which parks
    the neural row's provenance, and `benzene_center_tier_upgrade`, closing with
    "Reopen only if someone finds a route from a trace to a PAIRING".
  * `fw.Confirmations` returns nothing neural; that registry is hardware-only.
  * `docs/GLOSSARY.md` returns NOTHING that fences the words used here. Its only
    bare "concentration" is inside the sum-gamma row and is used positively; the
    rest are the Concentrator, a different object. The Re-spread naming below is
    this file's own choice and not a convention inherited from there.

WHAT THIS SCRIPT ESTABLISHES
----------------------------
  1. the pairing score is a reading of the tolerance against the spectral scale
     (STEP 1), of the normalisation constant (STEP 2), and of nothing about Dale's law (2b);
  2. the strict matcher is a function of the eigenvalue LIST rather than of the
     spectrum (STEP 3), which is what the 18-unpaired mode list rests on;
  3. the wiring's degeneracy at zero is real and is mostly a matching fact
     rather than a spectral one (STEP 5);
  4. and the model runs on a limit cycle at its own committed parameters
     (STEP 7), where the cavity document's surviving null said it does not.
     NOT that it reaches the gamma band: the integrated equations carry no time
     constant, so the Hz are a stipulation and G15c gates that.

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
F_PRIME = 0.3       # the committed linearisation slope, neural_gamma_cavity.py:221
PRIMES = (2147483647, 2147483629, 104729)
R_ENS = 200         # ensemble size for every null below

out = []


def log(msg=""):
    print(msg)
    out.append(msg)


failures = []


gate_names = []


def gate(name, ok, detail):
    gate_names.append(name)
    if not ok:
        failures.append(name)
    log(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")


# -------------------------------------------------------------------
# The two matchers, reproduced verbatim (vectorised, same answers)
# -------------------------------------------------------------------

def loose_pairing(eig, tol=TOL):
    """neural_gamma_cavity.py:254-260. Self-pairing is NOT excluded.

    Reproduced including the original's quirk that the target's imaginary part
    is +ev.imag: the -ev.imag computed into partner_target is then discarded.
    """
    center = np.mean(eig.real)
    targets = (2 * center - eig.real) + 1j * eig.imag
    d = np.abs(eig[None, :] - targets[:, None])
    return int(np.sum(d.min(axis=1) < tol))


def exclusive_pairing(eig, tol=TOL):
    """neural_gamma_cavity_unpaired.py:88-114. Forbids j == i.

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
    """The committed C. elegans construction, neural_gamma_cavity.py:216-222.

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


def zero_multiplicity(M, p, return_chain=False):
    """Algebraic multiplicity of the eigenvalue 0, exactly, by the rank chain.

    nullity(M^k) rises until the Jordan structure is exhausted, and where it
    stops is the algebraic multiplicity. Only the k = 1 term, the nullity, counts
    modes the wiring does nothing to; on the rest the wiring acts nilpotently
    (Wv != 0 while W^k v = 0), so they are feedforward depth, not idleness.

    Arithmetic note: the product is taken in Python ints via dtype=object, not
    int64. It is the sum of PRODUCTS that overflows, not the sum of residues:
    two products of residues below p ~ 2**31 already reach 2 * 2**62 = 9.22e18,
    the int64 ceiling, so a 300-term row of P @ A overflows silently and returns
    a nullity chain that DECREASES, which is impossible (see
    reference_modp_int_overflow_near_2_30). gf_rank's own int64 path is safe: it
    scales rows elementwise and never sums a row of products.
    """
    n = M.shape[0]
    A = np.rint(M).astype(np.int64).astype(object) % p
    prev, P, chain = -1, A.copy(), []
    for _ in range(n + 1):
        cur = n - gf_rank(np.asarray(P % p, dtype=np.int64), p)
        if cur == prev:
            return (cur, chain) if return_chain else cur
        prev = cur
        chain.append(cur)
        P = (P @ A) % p
    raise RuntimeError("rank chain did not stabilise")


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
# Condition (b) of PROOF_PALINDROME_NEURAL, at uniform tau, asks for an
# involution Q with Q W Q = -W. Under Dale a neuron's whole outgoing row carries
# one sign, so Q would have to send every non-empty EXCITATORY row to a
# non-empty INHIBITORY one (an empty image would force the source row empty
# too). That is a count, and it is decidable here rather than left open.
_ne_rows = [i for i in range(N) if W_real[i, :].any()]
_ne_exc = [i for i in _ne_rows if sign[i] >= 0]
_ne_inh = [i for i in _ne_rows if sign[i] < 0]
rng = np.random.default_rng(SEED)
eig_real = np.linalg.eigvals(jacobian(W_real))
sd_real = float(np.std(eig_real.real))

log("-" * 75)
log("STEP 0: the centre the test reflects across is forced, not measured")
log("-" * 75)
log()
log(f"  Condition (b) of PROOF_PALINDROME_NEURAL asks for an involution Q with")
log(f"  Q W Q = -W. At uniform tau the proof's Step 3 imposes nothing (S is")
log(f"  (1/tau)*I for EVERY permutation), so (b) is the whole hypothesis, and")
log(f"  under Dale it is a COUNT: a neuron's outgoing row carries one sign, so Q")
log(f"  must send each non-empty excitatory row to a non-empty inhibitory one.")
log(f"    non-empty rows              : {len(_ne_rows)}")
log(f"    of them excitatory          : {len(_ne_exc)}")
log(f"    of them inhibitory          : {len(_ne_inh)}")
gate("G0b condition (b) FAILS on this connectome, and on a count rather than "
     "an untested hypothesis",
     len(_ne_exc) != len(_ne_inh),
     f"an involution with Q W Q = -W would have to send each of the "
     f"{len(_ne_exc)} non-empty excitatory rows to a non-empty row of the "
     f"opposite sign, and there are {len(_ne_inh)} of those. An empty image "
     f"will not do either: it would force the source row empty too. So no such "
     f"Q exists and the palindrome theorem does not apply here, for a reason "
     f"that is decided rather than left open. What (b) needs is a BIJECTION "
     f"between the two sets, so any inequality kills it and the predicate is "
     f"'!=', not '>'. The theorem's PREMISE fails too, on a different count "
     f"(274 excitatory against 26 over ALL neurons, not just non-empty rows); "
     f"the two conditions are independent and this gate is about (b). "
     f"Falsifier: an EQUAL count would make this gate fail and leave (b) open "
     f"to be tested directly. It is unreachable here only because "
     f"{len(_ne_rows)} is odd")
log()
log("-" * 75)
log("STEP 0c: the March 8x, reproduced from the connectome file and taken apart")
log("-" * 75)
log()
log("  The 8x enrichment of docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md was")
log("  withdrawn 2026-08-26 as a difference of normalisation constants. The")
log("  numbers that replaced it lived only in prose, so they are computed here")
log("  from the same JSON, in that page's protocol rather than this page's:")
log("  W[i,j] = sign[j]*chemical[j,i] divided by the GLOBAL max, tau = 10/20,")
log("  alpha = 0.3, RandomState(trial+100), 200 balanced blocks per size.")
log()

_TE, _TI, _AL = 10.0, 20.0, 0.3
_Wsig = np.zeros((N, N))
for _i in range(N):
    for _j in range(N):
        _Wsig[_i, _j] = sign[_j] * W_raw[_j, _i]
_Wglob = _Wsig / np.max(np.abs(_Wsig))
_exc_all = np.where(sign > 0)[0]
_inh_all = np.where(sign < 0)[0]


def _ap_jac(W, s):
    """algebraic_palindrome.py's build_jacobian, verbatim in effect."""
    n = len(s)
    J = np.zeros((n, n))
    for i in range(n):
        t = _TE if s[i] > 0 else _TI
        J[i, i] = -1.0 / t
        for j in range(n):
            if i != j:
                J[i, j] = _AL * W[i, j] / t
    return J


def _ap_swap(s):
    n = len(s)
    e = np.where(s > 0)[0]
    k = np.where(s < 0)[0]
    perm = np.arange(n)
    for a in range(min(len(e), len(k))):
        perm[e[a]] = k[a]
        perm[k[a]] = e[a]
    Q = np.zeros((n, n))
    for i in range(n):
        Q[i, perm[i]] = 1.0
    return perm, Q


def _ap_roff(J, Q):
    """S is fitted per seat, so the DIAGONAL residual is zero by construction
    and only this off-diagonal number carries anything."""
    QJQ = Q @ J @ Q.T
    S = -(np.diag(QJQ) + np.diag(J)) / 2.0
    R = QJQ + J + 2 * np.diag(S)
    R = R - np.diag(np.diag(R))
    return float(np.linalg.norm(R) / np.linalg.norm(J))


def _ap_weff(W, s):
    J = _ap_jac(W, s)
    return float(np.linalg.norm(J - np.diag(np.diag(J))))


_rows = {}
_maxrel = 0.0
for _nh in (5, 10, 13):
    _n = 2 * _nh
    _un, _ma, _ct, _we_c, _we_r, _cov = [], [], [], [], [], []
    _coll = _exact = 0
    for _t in range(200):
        _rng = np.random.RandomState(_t + 100)
        _e = _rng.choice(_exc_all, _nh, replace=False)
        _k = _rng.choice(_inh_all, _nh, replace=False)
        _idx = np.concatenate([_e, _k])
        _s = sign[_idx]
        _W = _Wglob[np.ix_(_idx, _idx)]
        _perm, _Q = _ap_swap(_s)
        _J = _ap_jac(_W, _s)
        _un.append(_ap_roff(_J, _Q))
        _we_c.append(_ap_weff(_W, _s))
        _mx = np.max(np.abs(_W))
        _ma.append(_ap_roff(_ap_jac(_W / _mx if _mx > 0 else _W, _s), _Q))
        _dens = np.count_nonzero(_W) / (_n * (_n - 1))
        _Wr = np.zeros((_n, _n))
        for _a in range(_n):
            for _b in range(_n):
                if _a != _b and _rng.random() < max(_dens, 0.01):
                    _Wr[_a, _b] = _s[_b] * _rng.exponential(0.3)
        _m2 = np.max(np.abs(_Wr))
        if _m2 > 0:
            _Wr = _Wr / _m2
        _ct.append(_ap_roff(_ap_jac(_Wr, _s), _Q))
        _we_r.append(_ap_weff(_Wr, _s))
        _nzs = set(zip(*np.nonzero(_W)))
        _is_cov = not any((_perm[_i], _perm[_j]) in _nzs for (_i, _j) in _nzs)
        _cov.append(_is_cov)
        if _is_cov:
            _coll += 1
            _pred = np.sqrt(2) * _we_c[-1] / np.linalg.norm(_J)
            _rel = abs(_pred - _un[-1]) / _un[-1] if _un[-1] else 0.0
            _Wm = _W / _mx if _mx > 0 else _W
            _Jm = _ap_jac(_Wm, _s)
            _predm = np.sqrt(2) * _ap_weff(_Wm, _s) / np.linalg.norm(_Jm)
            _relm = abs(_predm - _ma[-1]) / _ma[-1] if _ma[-1] else 0.0
            _maxrel = max(_maxrel, _rel, _relm)
            if _rel < 1e-12 and _relm < 1e-12:
                _exact += 1
    _d = np.array(_ma) - np.array(_ct)
    _cm = np.array(_cov)
    _MA, _CT = np.array(_ma), np.array(_ct)
    _rows[_n] = dict(
        un=float(np.mean(_un)), ma=float(np.mean(_ma)), ct=float(np.mean(_ct)),
        wr=float(np.mean(_we_r) / np.mean(_we_c)), coll=_coll, exact=_exact,
        t=float(_d.mean() / (_d.std(ddof=1) / np.sqrt(len(_d)))),
        r_cov=float(_MA[_cm].mean() / _CT[_cm].mean()) if _cm.any() else float("nan"),
        r_unc=float(_MA[~_cm].mean() / _CT[~_cm].mean()) if (~_cm).any() else float("nan"))

log("     N   unmatched   weff-ratio   matched   paired t   collapse blocks")
for _n, _v in _rows.items():
    log(f"  {_n:4d}   {_v['ct'] / _v['un']:9.3f}   {_v['wr']:10.3f}   "
        f"{_v['ma'] / _v['ct']:7.4f}   {_v['t']:8.2f}   {_v['coll']:3d}/200")
log()

gate("G0c the unmatched enrichment IS the arms' coupling magnitude, as an "
     "arithmetic consequence rather than as an independent test",
     all(abs((_v['ct'] / _v['un']) / _v['wr'] - 1.0) < 0.01
         for _v in _rows.values()),
     "at every block size the residual ratio and the ratio of the two arms' "
     "mean ||W_eff|| agree to better than 1 percent: "
     + ", ".join(f"N={_n} {_v['ct'] / _v['un']:.3f} against {_v['wr']:.3f}"
                 for _n, _v in _rows.items())
     + ". The 8.46 of the withdrawn table is the first of these. Read this as a "
       "consistency check on the arithmetic and NOT as evidence: given G0e's "
       "closed form the residual is sqrt(2)*x/sqrt(d*d + x*x) with x = ||W_eff|| "
       "and d = ||diag J|| identical in both arms by construction, so with "
       "x << d the two ratios are FORCED to agree to under a percent and no "
       "wiring effect could separate them. Falsifier: a block set where x is "
       "comparable to d, where the agreement is no longer forced and can be "
       "tested")

gate("G0d matched, the two arms are at parity ONLY at the smallest block size",
     (abs(_rows[10]['t']) < 2.0 and _rows[20]['t'] < -5.0
      and _rows[26]['t'] < -5.0
      and all(_v['ma'] / _v['ct'] < 1.0 for _v in _rows.values())),
     "give both arms the control's own rule and the ratio runs "
     + ", ".join(f"{_v['ma'] / _v['ct']:.4f} at N={_n} (paired t {_v['t']:.2f})"
                 for _n, _v in _rows.items())
     + ". So 0.960 is the N=10 answer and NOT the verdict: at N=20 and N=26 the "
       "worm's matched residual stays below the control's at more than 5 sigma. "
       "What that residue IS, this gate does not decide, and the split by G0e's "
       "condition is the reason it does not: on the blocks the closed form "
       "covers, the residual is a function of the weight multiset and cannot "
       "register wiring, but on the ones it does not cover it can, and the gap "
       "is on both sides ("
     + ", ".join(f"N={_n} covered {_v['r_cov']:.4f} against uncovered "
                 f"{_v['r_unc']:.4f}" for _n, _v in _rows.items())
     + "). Falsifier: parity at every size, or a ratio above 1")

gate("G0e where no Q-partner pair of edges is present the residual IS "
     "sqrt(2)*||W_eff||/||J||, with nothing left over for wiring",
     (all(_v['exact'] == _v['coll'] for _v in _rows.values())
      and _maxrel < 1e-14),
     "the closed form reproduces the measured residual on every block that "
     "meets the condition, on BOTH arms of the comparison, the globally "
     "normalised one and the matched one, "
     + ", ".join(f"{_v['exact']}/{_v['coll']} at N={_n}"
                 for _n, _v in _rows.items())
     + f", and the largest relative deviation over all of them is {_maxrel:.2e}, "
       f"or {_maxrel / 2.220446049250313e-16:.1f} ulp. The bound asserted here is "
       "1e-14, about 45 ulp, and it is a ceiling on float rounding for two "
       "routes to one number rather than a value read off the data; what would "
       "make it a law is a deviation that grows with N, and it does not. The condition holds in "
     + ", ".join(f"{_v['coll']}/200 blocks at N={_n}" for _n, _v in _rows.items())
     + ", so the instrument rarely gets to see wiring at all. Falsifier: a "
       "block meeting the condition whose residual departs from the closed form")
log()
log("-" * 75)
log()
log("F137 states this in general: the centre is trace/dim, an identity, equally")
log("well defined for a spectrum that does not pair. NEURAL_CLOCK_TWO_HANDS:42-59")
log("proves the neural case is wiring-INDEPENDENT (its own value is -0.150000 at")
log("its own time constants; it is -1 here only because this Jacobian's damping")
log("is 1). Here the wiring-independence is exact rather than approximate.")
log()
n_selfsyn = int(np.count_nonzero(np.diag(W_raw)))
log(f"  neurons with a synapse onto themselves: {n_selfsyn}   (exact count)")
_diag_ok = all(np.count_nonzero(np.diag(degree_preserving_swap(
    W_raw, np.random.default_rng(SEED + 100 + k), passes=1)[0])) == 0
    for k in range(20))
gate("G0 the wiring never touches the diagonal", n_selfsyn == 0 and _diag_ok,
     f"so trace(J)/N = -1 exactly, for this connectome and for every rewiring of "
     f"it. Be exact about what each half is worth. The zero self-synapse count "
     f"is DATA and could have come out otherwise. The 20 rewirings are a "
     f"REGRESSION GUARD, not a measurement: degree_preserving_swap rejects any "
     f"swap with a repeated endpoint, so the two cells it writes are "
     f"off-diagonal by construction and no draw could have failed. It guards "
     f"the swap rule against a future edit, which is worth having and is not "
     "evidence")
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
binary_score = loose_pairing(e_bin) / N * 100
log(f"  {'binary, unweighted':26s} {'-':>9s} "
    f"{binary_score:7.1f}% {np.std(e_bin.real):8.4f}")
log("    (listed for reference, NOT counted in G3: sign(W) discards every")
log("     weight, so it is a different matrix, not a rescaled one.)")
log()
log(f"  {len(nz)} directed edges, weight median {np.median(np.abs(nz)):.0f} and "
    f"max {np.max(np.abs(nz)):.0f}:")
log("  the committed normalisation is by one outlier synapse.")
log()
def _cand_complex(e, tol=TOL):
    """How many eigenvalues lie within tol of a given one's MIRROR.

    Not |e_i - e_j|: the matcher never compares an eigenvalue to another
    eigenvalue. It reflects one about the centre and asks what the spectrum has
    there, so the candidate set is |e_j - mirror(e_i)| < tol.
    """
    ctr = float(e.real.mean())
    mirror = (2.0 * ctr - e.real) + 1j * e.imag
    return (np.abs(e[None, :] - mirror[:, None]) < tol).sum(axis=1)


_cand_cplx = _cand_complex(eig_real)   # mirror-based, as the matcher works
log(f"  partners within TOL in the COMPLEX plane: median "
    f"{np.median(_cand_cplx):.0f}, mean {_cand_cplx.mean():.1f}")
gate("G2b every eigenvalue has many candidate partners before any structure",
     np.median(_cand_cplx) > 5,
     f"a median of {np.median(_cand_cplx):.0f} eigenvalues lie within the "
     f"committed tolerance of a given one's MIRROR, which is what the matcher "
     f"compares against, in the plane it works in. The distribution is skewed, "
     f"mean {_cand_cplx.mean():.1f}, so the median is the number to read. "
     "TOL/mean-gap is a looser and different quantity")
_span = max(div_scores.values()) - min(div_scores.values())
gate("G3 the score is set by the normalisation constant, which nothing justifies",
     _span > 25.0,
     f"the same connectome and the same matcher give "
     f"{min(div_scores.values()):.1f}% to {max(div_scores.values()):.1f}%, a "
     f"span of {_span:.1f} points, across the FOUR genuine normalisation "
     f"constants (the binary row is excluded; including it would inflate this "
     f"to {max(max(div_scores.values()), binary_score) - min(min(div_scores.values()), binary_score):.1f})")
binary_matched = loose_pairing(spread_matched(e_bin, sd_real)) / N * 100
log(f"  {'binary, AT the connectome sd':26s} {'-':>9s} "
    f"{binary_matched:7.1f}% {sd_real:8.4f}")
gate("G3b the binary row is a SPREAD effect, not a weight effect",
     binary_matched >= sweep[TOL][0] - 1.0,
     f"raw, the binary matrix scores {binary_score:.1f}% with a Re-spread of "
     f"{np.std(e_bin.real):.4f} against the connectome's {sd_real:.4f}. Rescaled "
     f"to the connectome's own spread it scores {binary_matched:.1f}% against "
     f"{sweep[TOL][0]:.1f}%. So discarding EVERY weight costs nothing once scale "
     "is controlled, and the 27% is the same measurement as G2 and G3, not a "
     "separate finding about which neuron connects to which")
log()
log("  Matched on Re-spread, noise with no structure at all is not beaten.")
log("  Read the fence with the result: sd is NOT the quantity the matcher sees.")
log("  What it sees is the tolerance against the LOCAL spacing, and the")
log("  connectome's Re-distribution is heavy-tailed, so matching the sd MIGHT")
log("  be thought to hand the control a denser spectrum. Measured, it runs the")
log("  other way: the counts below, and G4b, gate that.")
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
def _mean_gap(e):
    """CAUTION: mean(diff(sort(x))) is identically (max - min)/(N - 1).

    It is the range divided by a count. No clustering can move it, so it is a
    fine scale to quote and a useless one to compare densities with. The
    quantity the matcher actually works with is _candidates below.
    """
    r = np.sort(e.real)
    return float(np.mean(np.diff(r)))


def _candidates(e, tol=TOL):
    """How many other eigenvalues lie within tol of a given one, in Re."""
    r = np.sort(e.real)
    return np.array([int(np.sum(np.abs(r - x) < tol)) - 1 for x in r])

_gap_real = _mean_gap(np.linalg.eigvals(jacobian(W_real)))
_rng_aux = np.random.default_rng(SEED + 1)   # OWN stream: drawing from `rng`
# here would shift every downstream ensemble and silently change committed numbers.
# ONE draw, reused everywhere below: a score from one matrix and candidate counts
# from another are not statements about the same control.
_e_ctrl = spread_matched(np.linalg.eigvals(
    jacobian(_rng_aux.normal(size=(N, N)))), sd_real)
_gap_ctrl = _mean_gap(_e_ctrl)
_cand_real = _candidates(np.linalg.eigvals(jacobian(W_real)))
_cand_ctrl = _candidates(_e_ctrl)
# The fence below reads candidate counts and Re-width off _e_ctrl, while the
# table above scored a DIFFERENT gaussian draw (from `rng`). Joining the two
# would be exactly the seam this file's own comment warns about, so the fence's
# own matrix is scored here and gated to agree. The table's draws are left
# untouched: removing them would shift `rng` for every ensemble downstream.
_ctrl_score = loose_pairing(_e_ctrl) / N * 100
log(f"  the fence's OWN gaussian draw scores   : {_ctrl_score:.1f}% "
    f"(the table above scored a different draw at "
    f"{matched['dense iid gaussian']:.1f}%; every count and width in the fence "
    f"belongs to this one)")
gate("G4c the matrix the fence measures is the matrix the fence scores",
     _ctrl_score >= sweep[TOL][0],
     f"the gaussian draw whose candidate counts and Re-width the fence reports "
     f"scores {_ctrl_score:.1f}% itself, at or above the connectome's "
     f"{sweep[TOL][0]:.1f}%, so the 'a structureless matrix scores higher with "
     f"fewer partners' reading is one statement about one matrix rather than a "
     f"score borrowed from a sibling draw. The table's own draw is kept "
     f"separately at {matched['dense iid gaussian']:.1f}% and is not what the "
     "fence describes")
log()
log(f"    Re-width, connectome {np.ptp(eig_real.real):.4f}   control "
    f"{np.ptp(_e_ctrl.real):.4f}   (at EQUAL sd {sd_real:.4f})")
log(f"    Re-range/(N-1), connectome : {_gap_real:.3e}   control {_gap_ctrl:.3e}")
log("      (an identity for the range; it cannot see clustering, so it is NOT")
log("       the fence. The fence is the next two lines.)")
_cx_real = _cand_complex(eig_real)
_cx_ctrl = _cand_complex(_e_ctrl)
log(f"    partners within TOL, in Re     : connectome "
    f"median {np.median(_cand_real):.0f}, control {np.median(_cand_ctrl):.0f}")
log(f"    partners within TOL, in the plane the matcher uses : connectome "
    f"median {np.median(_cx_real):.0f}, control {np.median(_cx_ctrl):.0f}")
gate("G4 spread-matched noise scores at least as high as the connectome",
     min(matched.values()) >= sweep[TOL][0],
     f"the weakest matched control scores {min(matched.values()):.1f}% against "
     f"the connectome's {sweep[TOL][0]:.1f}%. NOTE this gate saturates: the "
     "control's score sits at the ceiling, so it cannot fail by much and is "
     "reported as an observation, not as the load-bearing control")
gate("G4b the control's 100% is not bought with more partners to choose from",
     np.median(_cx_ctrl) <= np.median(_cx_real)
     and np.median(_cand_ctrl) <= np.median(_cand_real),
     f"in the complex plane the connectome has a median of "
     f"{np.median(_cx_real):.0f} eigenvalues within TOL of a given one against "
     f"the control's {np.median(_cx_ctrl):.0f} (in Re alone, "
     f"{np.median(_cand_real):.0f} against {np.median(_cand_ctrl):.0f}; the axis "
     f"has to be named, the two counts differ). The control reaches 100% with FEWER "
     "candidates, so sd-matching did not hand it an easier spectrum. The "
     "Re-range cannot decide this either way: it is an identity that cannot see "
     "clustering. The degree-matched rewiring below is what carries Result 2")
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
log(f"  (results/neural_gamma_cavity_unpaired.txt:20) and the same script on the")
log(f"  same data reports {native_unp} today.")
_rng_perm = np.random.default_rng(SEED + 3)   # OWN stream
_perm200 = [len(exclusive_pairing(eig_real[_rng_perm.permutation(N)])[1])
            for _ in range(200)]
_sorted_re = len(exclusive_pairing(eig_real[np.argsort(eig_real.real)])[1])
log(f"  200 permutations                       : mean {np.mean(_perm200):.1f}, "
    f"range {min(_perm200)}-{max(_perm200)}")
log(f"  eigenvalues sorted by real part        : {_sorted_re} unpaired")
_rng_str = np.random.default_rng(SEED + 4)   # OWN stream
_strict_perm = [(N - len(exclusive_pairing(eig_real[_rng_str.permutation(N)])[1]))
                / N * 100 for _ in range(60)]
_strict_sorted = (N - len(exclusive_pairing(eig_real[np.argsort(eig_real.real)])[1])) / N * 100
_strict_native = sweep[TOL][1]
log()
log("  The same order-dependence reaches STEP 2's strict column and STEP 4's")
log("  strict p-value, which is not a separate statistic: it is this matcher.")
log(f"    LAPACK order      : {_strict_native:.1f}%")
log(f"    60 permutations   : {min(_strict_perm):.1f}-{max(_strict_perm):.1f}%")
log(f"    sorted by real    : {_strict_sorted:.1f}%")
_orbit = max(_strict_perm + [_strict_sorted, _strict_native]) - min(
    _strict_perm + [_strict_sorted, _strict_native])
gate("G6c the strict SCORE carries the same defect, so only the loose row is "
     "order-free",
     _orbit > 3.0,
     f"reordering moves the strict score over {_orbit:.1f} points "
     f"({min(_strict_perm + [_strict_sorted, _strict_native]):.1f} to "
     f"{max(_strict_perm + [_strict_sorted, _strict_native]):.1f}), which G8b "
     "below compares against the degree-matched null spread for the same "
     "statistic. "
     "Each null draw is scored in its own arbitrary order too, so the strict "
     "p-value compares one ordering against R others. The loose matcher has no "
     "such orbit and is what Result 2 should be read on")
gate("G6b the orbit is wider than any one sample of it shows",
     min(_perm200) < min(perm_unp) or max(_perm200) > max(perm_unp)
     or _sorted_re < native_unp,
     f"50 draws span {min(perm_unp)}-{max(perm_unp)}, 200 draws span "
     f"{min(_perm200)}-{max(_perm200)}, and a deterministic sort by real part "
     f"gives {_sorted_re}, below LAPACK's {native_unp}. The count is a property "
     "of the LIST, and its orbit has no natural width")
log()
gate("G6 ordering moves the count further than the run-to-run drift does",
     np.std(perm_unp, ddof=1) > 1.0
     and abs(np.mean(perm_unp) - native_unp) > max(abs(native_unp - 18), 1),
     f"reordering shifts the mean to {np.mean(perm_unp):.1f} from {native_unp}, "
     f"while the committed 18 and today's {native_unp} differ by "
     f"{abs(native_unp - 18)}")
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
    # Error model, stated because a threshold without one is a number: LAPACK's
    # backward error for the nonsymmetric eigenproblem is O(n)*eps*||A||, and the
    # induced eigenvalue motion is that times the eigenvalue's condition number.
    # The n is the part an earlier reading dropped, and it costs a factor of 300.
    motion = float(np.max(k_off)) * N * float(np.finfo(float).eps) * normJ
    motion_1 = float(np.max(k_off)) * float(np.finfo(float).eps) * normJ
    log(f"  eigenvalues on the centre to 1e-4 : {int(np.sum(~off))}")
    log(f"  condition numbers of the other {int(np.sum(off)):3d} : median "
        f"{np.median(k_off):.0f}, 90th pct {np.percentile(k_off, 90):.0f}, "
        f"max {np.max(k_off):.3g}")
    _e_off = _ev[off]
    _D = np.abs(_e_off[:, None] - _e_off[None, :])
    np.fill_diagonal(_D, np.inf)
    _nn = _D.min(axis=1)
    _min_sep = float(_nn.min())
    log(f"  WORST induced eigenvalue motion   : {motion:.2e}, against a mean "
        f"Re-spacing of {np.mean(gaps):.2e}")
    log(f"  closest two off-centre eigenvalues: {_min_sep:.2e} apart in the "
        f"complex plane (NOT in Re alone: conjugate pairs share a real part "
        f"exactly, so the smallest Re-gap is 0.0 by construction and measures "
        f"nothing)")
    kappa_max = float(np.max(kappa))
    eig_error_bound = kappa_max * float(np.finfo(float).eps) * normJ
    log(f"  worst-conditioned eigenvalue      : kappa = {kappa_max:.1e}, so no")
    log(f"    eigenvalue of this matrix is resolved better than "
        f"{eig_error_bound:.1e}")
    log(f"  READ, not a gate: motion {motion:.2e} against the mean Re-spacing "
        f"{np.mean(gaps):.2e}. That spacing is the range-over-count identity")
    log("  denounced 150 lines above and no eigenvalue in a clustered spectrum")
    log(f"  has anything like it; the closest actual separation is {_min_sep:.2e},")
    log("  the smaller and the honest number. G7 gates that one.")
    log()
    log("  THE MARGIN IS A READING OF THE CUTOFF, so it is not the statement.")
    log("  Which modes count as 'on the centre' is decided by a literal, and the")
    log("  worst conditioning outside it moves by orders when that literal does:")
    log(f"    {'cutoff':>8s} {'n_off':>6s} {'max kappa':>11s} {'motion':>10s} "
        f"{'vs min sep':>11s}")
    for _cut in (1e-4, 3e-4, 1e-3):
        _o = np.abs(_ev + 1) > _cut
        _e2 = _ev[_o]
        _D2 = np.abs(_e2[:, None] - _e2[None, :]); np.fill_diagonal(_D2, np.inf)
        _m2 = float(kappa[_o].max()) * N * float(np.finfo(float).eps) * normJ
        log(f"    {_cut:8.0e} {int(_o.sum()):6d} {kappa[_o].max():11.2e} "
            f"{_m2:10.2e} {_D2.min() / _m2:11.2f}")
    _straggler = float(np.abs(_ev[np.argmax(kappa * off)] + 1))
    log(f"  the binding eigenvalue sits at |lambda+1| = {_straggler:.2e}, i.e. "
        f"{_straggler / 1e-4:.2f}x the cutoff")
    log()
    log("  The cutoff-free statement pairs each eigenvalue's OWN motion against")
    log("  its OWN nearest neighbour, which is the distance its own pairing has")
    log("  to survive:")
    _mot_i = kappa[off] * N * float(np.finfo(float).eps) * normJ
    _nn_i = _nn
    _ratio = float((_mot_i / _nn_i).max())
    _j = int(np.argmax(_mot_i / _nn_i))
    log(f"    worst motion/gap ratio {_ratio:.4f}  -> margin {1 / _ratio:.2f}")
    log(f"    bound by kappa {kappa[off][_j]:.2e}, own gap {_nn_i[_j]:.2e}")
    log(f"    the CLOSEST pair, {_nn_i.min():.2e} apart, has kappa "
        f"{kappa[off][int(np.argmin(_nn_i))]:.0f} and moves "
        f"{_mot_i[int(np.argmin(_nn_i))]:.1e}")
    # And this margin is NOT constant in the cutoff either. What IS true is a
    # monotonicity: drop a mode from the set and every survivor's nearest
    # neighbour can only move away, so each ratio can only fall. The margin is
    # therefore non-decreasing in the cutoff, and the value at the FINEST cut is
    # a lower bound for every coarser one. That is the cutoff-robust statement.
    _marg = {}
    for _cut in (1e-8, 1e-6, 1e-5, 1e-4, 3e-4, 1e-3):
        _o = np.abs(_ev + 1) > _cut
        _e2 = _ev[_o]
        _D2 = np.abs(_e2[:, None] - _e2[None, :]); np.fill_diagonal(_D2, np.inf)
        _marg[_cut] = 1.0 / float(
            (kappa[_o] * N * float(np.finfo(float).eps) * normJ
             / _D2.min(axis=1)).max())
    log("    the same margin across cutoffs: "
        + ", ".join(f"{c:.0e}->{m:.2f}" for c, m in _marg.items()))
    _mvals = list(_marg.values())
    _monotone = all(a <= b + 1e-9 for a, b in zip(_mvals, _mvals[1:]))
    # The nearest neighbour above is taken WITHIN the surviving set, which is
    # itself a reading of the cutoff. Against the FULL spectrum, one mode loses.
    _Dfull = np.abs(_e_off[:, None] - _ev[None, :])
    for _k, _i in enumerate(np.where(off)[0]):
        _Dfull[_k, _i] = np.inf
    _nn_full = _Dfull.min(axis=1)
    _loses = int(np.sum(_mot_i > _nn_full))
    _ok = _mot_i < _nn_full
    _worst_ok = float((_nn_full[_ok] / _mot_i[_ok]).min())
    log(f"    against the FULL spectrum, not just the survivors: {_loses} of "
        f"{int(off.sum())} loses to its own nearest neighbour; the other "
        f"{int(_ok.sum())} clear theirs by at least {_worst_ok:.2f}x")
    gate("G7 all but ONE off-centre eigenvalue outruns its true nearest "
         "neighbour, and the exception is the straggler already convicted above",
         _loses == 1 and float(np.abs(_ev[np.argmax(kappa * off)] + 1)) < 2e-4
         and _monotone and abs(_marg[1e-4] - 1 / _ratio) < 1e-9,
         f"pairing each of the {int(np.sum(off))} off-centre eigenvalues with its "
         f"nearest neighbour AMONG THE SURVIVORS, the worst ratio of "
         f"rounding-level motion to that gap is {_ratio:.4f}, a margin of "
         f"{1 / _ratio:.2f}. That set is itself a reading of the cutoff, so the "
         f"honest version counts neighbours in the FULL spectrum, and there "
         f"exactly {_loses} mode loses: the kappa = {kappa[off][_j]:.1e} "
         f"straggler, whose nearest eigenvalue is a member of the excluded "
         f"cluster {_nn_full[_j]:.2e} away, below its own motion of "
         f"{_mot_i[_j]:.2e}. It is the same straggler the section above "
         f"convicts for sitting at 1.12x the cutoff. The other {int(_ok.sum())} "
         f"clear their true nearest neighbour by at least {_worst_ok:.2f}x. This margin is NOT independent of the cutoff "
         f"either, and saying so would repeat the defect one paragraph above: "
         f"it runs {', '.join(f'{c:.0e}->{m:.2f}' for c, m in _marg.items())}. "
         f"What is cutoff-robust is the DIRECTION, and it is a theorem rather "
         f"than a reading: drop a mode and every survivor's nearest neighbour "
         f"can only move away, so each ratio can only fall and the margin is "
         f"non-decreasing in the cutoff (checked over six cuts here). The "
         f"{1 / _ratio:.2f} is the value at the finest cut and therefore a "
         f"LOWER BOUND on the margin at every coarser one, which the "
         f"spectrum-wide comparison cannot offer: it moves from "
         f"{_min_sep / motion:.2f} to hundreds with no direction at all. "
         f"The error model is O(n) backward error, n*eps*||J||, and is named "
         f"because it is a choice: dropping the n would inflate this by about "
         f"{motion / motion_1:.0f}. So it is the ordering result above, not a "
         "conditioning argument, that carries the withdrawal of the mode list")
except ImportError:
    log("  (scipy unavailable; per-eigenvalue conditioning not computed)")
log()

log("-" * 75)
log("STEP 4: against degree-matched rewiring, by rank rather than by z")
log("-" * 75)
log()
log("ALGEBRAIC_PALINDROME_NEURAL:277-315 already ran this null on the algebraic")
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
log("  (R - #{null <= real} + 1)/(R + 1) = (#{null > real} + 1)/(R + 1).")
log("  Ties count in the connectome's favour, which is anti-conservative against")
log("  the usual #{null >= real} convention. Direction checked: G8 and G9 use p")
log("  as an UPPER bound, so the reported value is the harsher one there; G13")
log("  uses it as a lower bound and its statistic, multiplicity, has no ties.")
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
     0.05 < pvals["strict"] < 0.95,
     f"p = {pvals['strict']:.3f}; the document reads this matcher's unpaired "
     "modes as connectome structure")
n_tests = len(ens)   # the family itself, not a literal
bonf = 0.05 / n_tests
gate("G9 the tolerance-free residue is small and not established",
     (real_vals["cross-fine"] - float(np.mean(ens["cross-fine"]))) < 3.0
     and pvals["cross-fine"] > bonf,
     f"cross-pairing at tol 1e-4 is {real_vals['cross-fine']:.1f}% against a "
     f"null mean of {np.mean(ens['cross-fine']):.1f}%, an excess of "
     f"{real_vals['cross-fine'] - float(np.mean(ens['cross-fine'])):.1f} points "
     f"({(real_vals['cross-fine'] - float(np.mean(ens['cross-fine']))) / 100 * N:.0f} "
     f"modes) at p = {pvals['cross-fine']:.3f}, which does not clear "
     f"{bonf:.4f}, the level {n_tests} one-sided tests need. It is not the only "
     "statistic above its null: the loose score, self-paired, nullity, "
     "multiplicity and structural nullity are too. The last three are Result "
     "5's one finding measured three ways; self-paired is a partial proxy for "
     "it, since about half of what it counts is not at zero")
_null_strict_width = float(np.max(ens["strict"]) - np.min(ens["strict"]))
gate("G8b the ordering orbit is comparable to the null spread it is scored "
     "against",
     _orbit > 0.5 * _null_strict_width,
     f"reordering one spectrum moves the strict score by {_orbit:.1f} points, "
     f"against a degree-matched null range of {_null_strict_width:.1f} points "
     f"({np.min(ens['strict']):.1f} to {np.max(ens['strict']):.1f}) for the "
     f"same statistic, i.e. "
     f"{_orbit / _null_strict_width * 100:.0f}% of the spread the p-value is "
     "read against. Each null draw is scored in its own arbitrary order too, so "
     "the strict p-value compares one ordering against R others. A matcher "
     "whose orbit were small against this spread would fail this gate, and its "
     "strict row would be worth reading")
log()

log("-" * 75)
log("STEP 5: what IS different about this connectome, and what kind of fact")
log("-" * 75)
log()
log("At tolerances fine enough to resolve levels, the self-paired modes are")
log("eigenvalues of f'*W_norm at zero. Only the NULLITY counts modes the wiring")
log("does nothing to; on the rest it acts nilpotently, which is depth, not")
log("idleness. Counting them needs a tolerance and an eigensolver.")
log("The synapse matrix is an integer matrix, so the same content is exact. Two")
log("different exact objects live here and only the second is the mode count:")
log()
log(f"  nullity(W)        = {real_vals['nullity']}   "
    "(geometric: independent eigenvectors at 0)")
log(f"  multiplicity of 0 = {real_vals['multiplicity']}   "
    "(algebraic, by the rank chain: the MODES)")
n_float = int(np.sum(np.abs(eig_real + 1) < 1e-4))
log(f"  float check, |lambda + 1| < 1e-4 : {n_float}")
_decades = {c: int(np.sum(np.abs(eig_real + 1) < c))
            for c in (1e-4, 1e-5, 1e-6, 1e-7, 1e-8)}
log("  the same count over four decades of cut  : "
    + ", ".join(f"{c:.0e}->{v}" for c, v in _decades.items()))
# The two modes the committed tolerance adds over the exact count. They are NOT
# on the centre line: Re(mu) is nonzero, they are merely within tol/2 of it, and
# they leave at the next cut. That is the tolerance argument, not an exception.
# self_paired_count uses tol/2 on |Re - c|, so the line test must match it.
_ctr_line = np.abs(eig_real.real - float(eig_real.real.mean())) < 1e-4 / 2
_surplus = eig_real[_ctr_line & (np.abs(eig_real + 1) >= 1e-4)]
_sp_dec = {c: self_paired_count(eig_real, c) for c in (1e-4, 1e-5, 1e-6, 1e-8)}
log(f"  within tol/2 of the centre line, not at zero, at 1e-4: {len(_surplus)}"
    + (f" (|Im| = {', '.join(f'{abs(z.imag):.4f}' for z in _surplus)}, "
       f"|Re - c| = {', '.join(f'{abs(z.real - float(eig_real.real.mean())):.2e}' for z in _surplus)})"
       if len(_surplus) else ""))
log("  self-paired as the cut tightens        : "
    + ", ".join(f"{c:.0e}->{v}" for c, v in _sp_dec.items()))
# The other half of the identification, which the rank chain does NOT cover: a
# NONZERO purely imaginary mu would be self-paired at every tolerance and would
# break "self-paired == at zero" in the limit. Checked here rather than asserted.
_mu = (eig_real + 1.0) / F_PRIME
_osc = _mu[np.abs(_mu.imag) > 1e-6]
_min_re_osc = float(np.abs(_osc.real).min())
log(f"  smallest |Re mu| among modes with |Im mu| > 1e-6 : {_min_re_osc:.2e}")
gate("G10c no NONZERO purely imaginary mode exists, which is the other half of "
     "the identification",
     _min_re_osc > 1e-8,
     f"the smallest |Re mu| among the {len(_osc)} modes with |Im mu| > 1e-6 is "
     f"{_min_re_osc:.2e}, four orders above the 1e-8 a purely imaginary mode "
     f"would have to beat. Such a mode would sit on the centre line at EVERY "
     f"tolerance and would make the limit self-paired count exceed the count at "
     f"zero, which is what Result 5 rests on. The GF(p) chain cannot see this: "
     f"it counts multiplicity AT ZERO and is silent about the rest of the "
     f"imaginary axis. An exact route exists over the integers, through the gcd "
     f"of the characteristic polynomial with its reflection, and is not run "
     "here")

gate("G10b the float count is a plateau, not a reading of one cut",
     len(set(_decades.values())) == 1,
     f"the count of eigenvalues at zero is {n_float} at every cut from 1e-4 "
     f"down to 1e-8 ({', '.join(f'{c:.0e}:{v}' for c, v in _decades.items())}), "
     "so it is not tracking the literal the way the pairing score does. A "
     "spectrum whose cluster merely APPROACHED zero would shed modes as the cut "
     f"tightened and this gate would fail. The self-paired count does shed: "
     + ", ".join(f"{c:.0e}->{v}" for c, v in _sp_dec.items())
     + f". The {len(_surplus)} mode(s) it adds at 1e-4 are not ON the centre "
     "line, which would make them self-paired at every cut; their Re(lambda) - c "
     f"is {', '.join(f'{z.real + 1:.2e}' for z in _surplus)} (so Re(mu), which "
     f"is that over f_prime, is "
     f"{', '.join(f'{(z.real + 1) / F_PRIME:.2e}' for z in _surplus)}), "
     f"nonzero, so they are "
     "merely WITHIN tol/2 of it and are gone by 1e-5. Below about 1e-8 the "
     "eigensolver stops resolving the cluster at all, so the limit statement "
     "belongs to the exact chain and not to these counts")
gate("G10 the exact multiplicity matches the float count",
     real_vals['multiplicity'] == n_float,
     f"{real_vals['multiplicity']} == {n_float}, so the rank chain counts what "
     "the tolerance was counting, exactly and without one")
mults = {zero_multiplicity(W_raw, p) for p in PRIMES}
mults_signed = {zero_multiplicity(W_real, p) for p in PRIMES}
gate("G11b the MULTIPLICITY, not only the rank, agrees at every prime tried",
     len(mults) == 1 and len(mults_signed) == 1,
     f"{len(PRIMES)} primes give {len(mults)} distinct value "
     f"({sorted(mults)[0]}) for the unsigned matrix and {len(mults_signed)} "
     f"({sorted(mults_signed)[0]}) for the Dale-signed one. Note the direction: "
     "rank over GF(p) is at most the rational rank, so the GF(p) nullity and "
     "multiplicity are UPPER bounds on the rational ones, and agreement across "
     "primes is the check. G11 below covers the rank alone, which is the weaker "
     "statement")
ranks_agree = len({gf_rank(W_raw, p) for p in PRIMES})
gate("G11 the exact rank agrees at every prime tried", ranks_agree == 1,
     f"{len(PRIMES)} primes ({', '.join(str(q) for q in PRIMES)}), "
     f"{ranks_agree} distinct value; rank over GF(p) "
     "bounds the rational rank from below, so agreement is the check")
nullity_signed = N - gf_rank(W_real, PRIMES[0])
mult_signed = zero_multiplicity(W_real, PRIMES[0])
log(f"  nullity unsigned {real_vals['nullity']}, with Dale's law applied "
    f"{nullity_signed}")
log(f"  multiplicity unsigned {real_vals['multiplicity']}, with Dale's law "
    f"applied {mult_signed}")
log("    A left diagonal sign matrix preserves RANK, so the nullity row is")
log("    indeed a theorem and no data could move it. The MULTIPLICITY row is")
log("    not: it depends on rank((DW)^k) and (DW)^k != D W^k. Two-by-two")
log("    counterexample, checked below, so the invariance is measured here")
log("    rather than assumed. Either way it decides a WORD: 'signed' does not")
log("    belong in the claim.")
_cW = np.array([[1.0, 1.0], [1.0, 1.0]])
_cD = np.diag([1.0, -1.0])
_cm_plain = zero_multiplicity(_cW, PRIMES[0])
_cm_scaled = zero_multiplicity(_cD @ _cW, PRIMES[0])
gate("G22 row scaling can change the multiplicity, so its invariance here is "
     "a measurement",
     _cm_plain != _cm_scaled and mult_signed == real_vals['multiplicity'],
     f"on [[1,1],[1,1]] a diag(1,-1) moves the multiplicity of 0 from "
     f"{_cm_plain} to {_cm_scaled} while the nullity is 1 either way. That half "
     f"runs on a fixed literal and cannot come out otherwise; it is here to "
     f"show the phenomenon is real, not to measure it. The measured half is the "
     f"second: on this connectome the multiplicity happens NOT to move "
     f"({real_vals['multiplicity']} signed and unsigned), and that is a fact "
     "about this matrix which a different one could break")
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
log()
log("  WHERE THE ZERO MODES SIT, AND IN WHICH MATRIX. Rows are PREsynaptic and")
log("  columns POSTsynaptic (ASIL, PVDR, M4: non-empty row, exactly empty")
log("  column; DA7 the reverse), so the drive reaching neuron j is")
log("  sum_i W[i,j] x_i = (W^T x)_j and the model's Jacobian is -I + f'W^T,")
log("  while jacobian() above builds -I + f'W. Every SPECTRAL quantity in this")
log("  file is untouched: a matrix and its transpose share eigenvalues, nullity,")
log("  rank chain and per-eigenvalue condition numbers. EIGENVECTORS are not,")
log("  and zero modes are an eigenvector question, so the decomposition below")
log("  is computed on W^T. A kernel vector of W^T is a perturbation that DRIVES")
log("  NOTHING; the cheapest way to be one is to have no target inside the")
log("  model, which is what an empty row means. For the MOTOR neurons among")
log("  them that target is muscle, which the 300-neuron matrix does not hold;")
log("  that is not the only way a row can be empty, and STEP 8 counts how many")
log("  of them are nevertheless wired inside the 300 by gap junctions.")
log()
_names = [l.split('	')[1].strip()
          for l in (NEURAL_DIR / "celegans_neuron_ids.txt").read_text().splitlines()
          if l.strip()]
_WT = W_raw.T


def _explicit_kernel(M):
    """Kernel vectors forced by an empty column or by two identical columns."""
    zero_cols = [j for j in range(N) if not M[:, j].any()]
    seen, twins = {}, []
    for j in range(N):
        col = M[:, j]
        if not col.any():
            continue
        k = col.tobytes()
        if k in seen:
            twins.append((seen[k], j))
        else:
            seen[k] = j
    V = [np.eye(N)[:, j] for j in zero_cols]
    for a, b in twins:
        v = np.zeros(N); v[a] = 1.0; v[b] = -1.0
        V.append(v)
    V = np.array(V).T
    return zero_cols, twins, V, float(np.abs(M @ V).max()), int(
        np.linalg.matrix_rank(V))


_zero_cols, _twin_cols, _V, _resid, _indep = _explicit_kernel(_WT)
_w_cols, _w_twins, _, _, _w_indep = _explicit_kernel(W_raw)
_zero_rows = [i for i in range(N) if not W_raw[i, :].any()]
log(f"    empty ROWS of W = empty columns of W^T (no target in the model, so")
log(f"      the perturbation drives nothing): {len(_zero_cols)}  "
    f"{', '.join(_names[j] for j in _zero_cols)}")
log(f"    duplicate columns of W^T          : "
    f"{', '.join(f'({_names[a]},{_names[b]})' for a, b in _twin_cols)}")
log(f"    ON W INSTEAD (the transpose, i.e. NOT the model's Jacobian) the same")
log(f"      construction gives {len(_w_cols)} + {len(_w_twins)} = {_w_indep} and "
    f"the mirror-image story, neurons with")
log(f"      no chemical INPUT inside the 300 (not the same set as 'sensory':")
log(f"      the 15 include VC6, M4, MCL, MCR, DVB and AINL). That reading")
log(f"      belongs to W^T's")
log(f"      transpose and does not describe the model's zero modes.")
log(f"    distinct boundary neurons         : "
    f"{len(set(_zero_rows) | set(_w_cols))}  (VC6 is both an empty row and an "
    f"empty column, so {len(_zero_rows)} + {len(_w_cols)} double-counts it)")
log()
log(f"    explicit kernel vectors           : {_indep} of nullity "
    f"{real_vals['nullity']}  ({100 * _indep / real_vals['nullity']:.0f}%)")
log(f"    (the {_indep} vectors are exactly in the kernel, worst |W v| = "
    f"{_resid}, and independent by construction: zero columns and a star of "
    f"twin differences on disjoint indices. Neither is a measurement. What IS "
    f"a measurement, and what G23 gates, is how far short of the nullity they "
    f"fall.)")
gate("G23 the boundary and the twins do NOT explain the whole kernel",
     _indep < real_vals['nullity'],
     f"{len(_zero_cols)} zero columns and {len(_twin_cols)} duplicate-column "
     f"pairs give {_indep} independent kernel vectors, EXACTLY in the kernel "
     f"(worst |W v| = {_resid}), which is {_indep} of {real_vals['nullity']}: "
     f"the boundary and the twins explain "
     f"{100 * _indep / real_vals['nullity']:.0f}% of the nullity and the "
     f"remaining {real_vals['nullity'] - _indep} are other column dependencies")
log(f"    READ, not a gate: the {len(_w_cols)} neurons with no incoming chemical")
log("    edge inside the matrix sit in the LEFT kernel of the model's Jacobian,")
log("    so they do not enter this count. Both sides give the same NULLITY (a")
log("    matrix and its transpose have equal rank), which is why the statistic")
log("    Result 5 reports is orientation-proof while its localisation is not.")
_mult, _chain = zero_multiplicity(W_raw, PRIMES[0], return_chain=True)
_inc = [_chain[0]] + [_chain[i] - _chain[i - 1] for i in range(1, len(_chain))]
log()
log(f"    nullity(W^k) chain : {_chain}")
log(f"    increments         : {_inc}  -> nilpotency index {len(_chain)}")
gate("G24 the multiplicity is two mechanisms, not one",
     len(_chain) > 1 and _mult > _chain[0],
     f"{_chain[0]} dimensions of kernel plus "
     f"{_mult - _chain[0]} of Jordan structure to depth {len(_chain)}, on which "
     "the wiring acts before annihilating, so only that second part is what "
     "'defective' means. NO graph reading is claimed: the longest-path theorem "
     "needs a NILPOTENT matrix, W has 236 nonzero eigenvalues, and "
     "[[0,1,0,1],[0,0,1,0],[0,1,0,1],[1,0,1,0]] is strongly connected with a "
     "size-2 Jordan block at 0. Depth is a fact about the operator, not the "
     "wiring")
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
log("NEURAL_GAMMA_CAVITY, Result 2b, reads Q_max = 0.1 against the qubit")
log("cavity's 68 to 75 as 'extremely lossy'. The eigenvalues of")
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
    f"{max(q_vals.values()) / min(q_vals.values()):.1f}")
log()
log("  AND THE tau AXIS OF THAT FACTOR IS AN ARTIFACT OF THE EXPRESSION.")
log("  -I/tau + f'W puts tau on the leak term only. That is the Jacobian of no")
log("  Wilson-Cowan model: for tau*xdot = -x + S(Wx) the whole right-hand side")
log("  carries it, J = (-I + f'W)/tau, so lambda = (-1 + f'mu)/tau and")
log("      Q = |Im lambda| / |Re lambda| = f'|Im mu| / |1 - f' Re mu|,")
log("  in which tau has CANCELLED. There is an exact route here, so it is")
log("  compared and not tolerated: dividing by a power of two is exact in IEEE,")
log("  so for dyadic tau the two Q values must be bit-identical, and any")
log("  deviation would be a fact about the construction, not a rounding.")
log()
log(f"  {'tau (ms)':>9s} {'f_prime':>8s} {'Q_max (Wilson-Cowan)':>22s} {'vs tau=1':>12s}")
log("  " + "-" * 54)
q_wc, dyadic_dev, all_dev = {}, 0.0, 0.0
for tau in (1.0, 2.0, 4.0, 8.0, 10.0):
    for fp in (0.3, 0.6):
        lam = (-1.0 + fp * mu) / tau
        m = np.abs(lam.imag) > 1e-6 / tau
        q = float(np.max(np.abs(lam[m].imag) / np.abs(lam[m].real)))
        q_wc[(tau, fp)] = q
        d = abs(q - q_wc[(1.0, fp)])
        all_dev = max(all_dev, d)
        if tau in (1.0, 2.0, 4.0, 8.0):
            dyadic_dev = max(dyadic_dev, d)
        log(f"  {tau:9.1f} {fp:8.2f} {q:22.6f} {d:12.1e}")
log()
log("  READ, not a gate: tau's cancellation from Q is ALGEBRA, so no spectrum")
log("  could make it fail and gating it would be a theorem about the operation.")
log("  STEP 6 above refuses to gate Q_max for the neighbouring reason; the same")
log("  refusal applies here. What IS measured: for dyadic tau the Q values are")
log(f"  bit-identical to tau = 1 (deviation exactly {dyadic_dev}), and including")
log(f"  the non-dyadic tau = 10 the worst deviation is {all_dev:.1e}, at or below")
log("  one ulp. Reported, not thresholded: any threshold here would be a number")
log("  with no derivation sitting orders of magnitude above the only value the")
log("  arithmetic can produce.")
log()
log("  What the DATA decide, and what G21 therefore gates, is the SIZE of the")
log("  artifact on this connectome: how much of the factor above the malformed")
log("  expression owes to tau, which the correct form says should be none.")
_tau_bad = max(q_vals[(10.0, f)] / q_vals[(1.0, f)] for f in (0.3, 0.6))
_tau_wc = max(q_wc[(10.0, f)] / q_wc[(1.0, f)] for f in (0.3, 0.6))
gate("G21 the tau axis of the reported factor is an artifact of the expression",
     _tau_bad > 5.0,
     f"on this spectrum the malformed -I/tau + f'W moves Q by up to "
     f"{_tau_bad:.1f}x along tau alone. The Wilson-Cowan form moves it by "
     f"exactly {_tau_wc:.1f}x, which is REPORTED and not gated, because the "
     f"cancellation is algebra and no spectrum could make it fail; gating it "
     f"would be a theorem about the operation, as the lines above say. A connectome whose spectrum made the leak term "
     f"negligible would show no artifact and this gate would fail. The "
     f"factor of {max(q_vals.values()) / min(q_vals.values()):.1f} above is "
     f"an artifact: across the same grid the Wilson-Cowan form moves by "
     f"{max(q_wc.values()) / min(q_wc.values()):.2f}, all of it along f_prime. "
     f"In the malformed expression the tau axis alone moves Q by "
     f"{q_vals[(10.0, 0.3)] / q_vals[(1.0, 0.3)]:.1f} at f_prime = 0.3 and "
     f"{q_vals[(10.0, 0.6)] / q_vals[(1.0, 0.6)]:.1f} at 0.6, against "
     f"{q_vals[(1.0, 0.6)] / q_vals[(1.0, 0.3)]:.1f} along f_prime, so almost "
     "the whole of that factor was the artifact and not a sensitivity of Q")
log()
log("  AND THE TWO AXES ARE ONE. In the malformed form Q depends on tau and")
log("  f' only through their PRODUCT: Q_bad(tau, f') = f'|Im mu| / |1/tau -")
log("  f' Re mu| is a function of tau*f' alone, which the table proves exactly,")
log(f"  Q(1, 0.6) = {q_vals[(1.0, 0.6)]:.6f} = Q(2, 0.3) = {q_vals[(2.0, 0.3)]:.6f}.")
log("  So the factor above is one range, not two sensitivities multiplied, and")
log("  in the correct form the product collapses to f' alone.")
log()
log("  Q_max is therefore arithmetic on ONE chosen constant and the graph")
log("  spectrum. The provenance is the finding, and it is checkable rather than")
log("  statistical: tau = 1 ms is stipulated bare at neural_gamma_cavity.py:54-55")
log("  for the Wilson-Cowan block, and the C. elegans Jacobian at :222 contains")
log("  no tau at all, only a bare -I. G0 above already pins that diagonal.")
log()

log("-" * 75)
log("STEP 7: the model runs on a limit cycle at its committed parameters, and")
log("        the conversion of its period to Hz is a stipulation, not a result")
log("-" * 75)
log()
log("neural_gamma_cavity.py:65 finds fixed points by Picard iteration, which")
log("converges only where the point is stable AS A MAP. That is not stability of")
log("the ODE; here it lands on the quiescent branch and reports its ringing. The")
log("sibling hopf_threshold.py in this same folder does integrate properly, but")
log("NOT this model: it runs a random balanced network of N = 100 to 5000 with")
log("tau_E = 5 and tau_I = 10, so it is a sibling METHOD, not a second reading")
log("of the two-population block. The integration below is this file's own.")
log("Doing it, with NOTHING changed:")
log()
wEE, wEI, wIE, wII, alpha_, theta_ = 16.0, 12.0, 15.0, 3.0, 1.3, 4.0
log(f"  the two-population parameters, printed so the page need not declare "
    f"them: w_EE/w_EI/w_IE/w_II = {wEE:g}/{wEI:g}/{wIE:g}/{wII:g}, "
    f"alpha = {alpha_:g}, theta = {theta_:g}")
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

    log(f"  {'I_ext':>6s} {'amplitude':>10s} {'period T':>9s} "
        f"{'tau=1ms':>9s} {'tau=10ms':>9s} {'band(1ms)':>10s}")
    log("  " + "-" * 58)
    bands = {}
    for I_ext in (0.5, 1.0, 1.12, 1.15, 1.2, 1.5, 2.0, 2.5, 3.0, 3.05, 3.07,
                  3.5, 4.0):
        s = solve_ivp(rhs, [0, 2000], [0.3, 0.3], args=(I_ext,), max_step=0.05,
                      dense_output=True, rtol=1e-9, atol=1e-11)
        t = np.linspace(1500, 2000, 200000)
        E = s.sol(t)[0]
        amp = float(E.max() - E.min())
        if amp < 1e-4:
            log(f"  {I_ext:6.2f} {amp:10.4f} {'n/a':>9s} {'n/a':>9s} "
                f"{'n/a':>9s} {'fixed pt':>10s}")
            continue
        m = E - E.mean()
        zc = np.where((m[:-1] < 0) & (m[1:] >= 0))[0]
        f_hz = 1000.0 / float(np.mean(np.diff(t[zc]))) if len(zc) > 2 else np.nan
        band = ("gamma" if 30 <= f_hz <= 100
                else "alpha" if 8 <= f_hz <= 13
                else "above" if f_hz > 100 else "other")
        bands[I_ext] = (f_hz, band)
        log(f"  {I_ext:6.2f} {amp:10.4f} {1000.0 / f_hz:9.1f} {f_hz:9.1f} "
            f"{f_hz / 10.0:9.1f} {band:>10s}")
    # The refined scans the page's Result 1 rests on, run here rather than
    # declared as figures outside the run. A helper so the coarse table above
    # and these use one measurement, not two.
    def _period(I_ext):
        _s = solve_ivp(rhs, [0, 2000], [0.3, 0.3], args=(I_ext,), max_step=0.05,
                       dense_output=True, rtol=1e-9, atol=1e-11)
        _t = np.linspace(1500, 2000, 200000)
        _E = _s.sol(_t)[0]
        if float(_E.max() - _E.min()) < 1e-4:
            return None
        _m = _E - _E.mean()
        _zc = np.where((_m[:-1] < 0) & (_m[1:] >= 0))[0]
        if len(_zc) <= 2:
            return None
        return float(np.mean(np.diff(_t[_zc])))

    _fold_lo, _fold_hi = _period(1.126), _period(3.062)
    _fine = {}
    _I = 1.6
    while _I <= 2.7 + 1e-9:
        _p = _period(round(_I, 4))
        if _p is not None:
            _fine[round(_I, 4)] = _p
        _I += 0.025
    _Imin = min(_fine, key=_fine.get)
    _Tmin = _fine[_Imin]
    log(f"  refined scans: T = {_fold_lo:.1f} at I_ext = 1.126 and "
        f"T = {_fold_hi:.1f} at 3.062, against a minimum of {_Tmin:.2f} at "
        f"I_ext = {_Imin:g} over I in [1.6, 2.7] step 0.025 "
        f"(T at I_ext = 2.00 is {_fine[2.0]:.4f}, which the table above "
        f"displays rounded)")
    gate("G15e the period grows by orders towards both folds and bottoms out "
         "in between, so the coarse table's span is a floor and not the range",
         _fold_lo > 20 * _Tmin and _fold_hi > 10 * _Tmin
         and _Tmin < min(1000.0 / f for f, _ in bands.values()),
         f"at the lower fold T = {_fold_lo:.1f} and at the upper T = "
         f"{_fold_hi:.1f}, against a scanned minimum of {_Tmin:.4f} at "
         f"I_ext = {_Imin:g}: ratios of {_fold_lo / _Tmin:.0f} and "
         f"{_fold_hi / _Tmin:.0f}. The minimum is SAMPLED, so it bounds the "
         f"true infimum from above, and it already sits below the coarse "
         f"table's smallest entry of "
         f"{min(1000.0 / f for f, _ in bands.values()):.4f}. What would fail "
         "this gate is a period that varied by less than 20x toward the lower "
         "fold or 10x toward the upper one, measured against the scanned "
         "minimum; BOUNDEDNESS is not what is tested, and a bounded period "
         "with a large enough ratio passes")
    log()
    in_gamma = [I for I, (f_hz, b) in bands.items() if b == "gamma"]
    f_all = [f_hz for f_hz, _ in bands.values()]
    log()
    log("  READ THE GRID, NOT ITS HULL. The frequency is not monotone in I_ext")
    log("  and the period diverges toward the upper edge, so the sampled")
    log("  frequencies are points on a curve, not an interval the model")
    log(f"  occupies. Sampled here: {min(f_all):.1f} to {max(f_all):.1f} Hz, and a")
    log("  finer grid finds lower values at both ends. What is robust is that")
    log("  the model OSCILLATES at all. The period diverges at both folds, so")
    log("  the attained set is unbounded above, and whether it is an INTERVAL needs a continuity the sampling cannot give and this grid understates it; and")
    log("  whether any of it is 30-100 Hz is decided by tau, not here (G15c).")
    log()
    log("  THE Hz AXIS IS A STIPULATION, NOT A MEASUREMENT. The right-hand side")
    log("  integrated above is [-E + S(...), -I + S(...)]: there is no tau in it.")
    log("  f_hz = 1000/T asserts that one time unit is one millisecond. Insert a")
    log("  real tau, tau*xdot = -x + S(...), and the trajectory is unchanged in")
    log("  shape while every frequency scales as 1/tau:")
    log()

    def rhs_tau(t, y, I_ext, tau):
        E, Iv = y
        return [(-E + sg(wEE * E - wEI * Iv + I_ext)) / tau,
                (-Iv + sg(wIE * E - wII * Iv)) / tau]

    def measure_tau(I_ext, tau):
        s_ = solve_ivp(rhs_tau, [0, 2000 * tau], [0.3, 0.3], args=(I_ext, tau),
                       max_step=0.05 * tau, dense_output=True, rtol=1e-9,
                       atol=1e-11)
        t_ = np.linspace(1500 * tau, 2000 * tau, 200000)
        E_ = s_.sol(t_)[0]
        amp = float(E_.max() - E_.min())
        if amp < 1e-4:
            return amp, float('nan')
        m_ = E_ - E_.mean()
        zc_ = np.where((m_[:-1] < 0) & (m_[1:] >= 0))[0]
        if len(zc_) <= 2:
            return amp, float('nan')
        return amp, 1000.0 / float(np.mean(np.diff(t_[zc_])))

    log(f"  {'I_ext':>6s} {'amp tau=1':>11s} {'amp tau=10':>11s} "
        f"{'f tau=1':>9s} {'f tau=10':>9s} {'ratio':>7s}")
    log("  " + "-" * 58)
    tau_rows = []
    for I_ext in (1.15, 1.20, 3.00):
        a1, f1 = measure_tau(I_ext, 1.0)
        a10, f10 = measure_tau(I_ext, 10.0)
        tau_rows.append((I_ext, a1, a10, f1, f10))
        log(f"  {I_ext:6.2f} {a1:11.8f} {a10:11.8f} {f1:9.2f} {f10:9.3f} "
            f"{f1 / f10:7.4f}")
    log()
    _amp_dev = max(abs(a1 - a10) / a1 for _, a1, a10, _, _ in tau_rows)
    _ratio_dev = max(abs(f1 / f10 - 10.0) for _, _, _, f1, f10 in tau_rows)
    _in_band_10 = [I for I, _, _, _, f10 in tau_rows if 30 <= f10 <= 100]
    log(f"  READ, not a gate: amplitude agrees to {_amp_dev:.1e} relative across")
    log(f"  tau = 1 and 10 while the frequency ratio is 10 to {_ratio_dev:.1e}.")
    log("  tau divides the whole right-hand side, so this is a rescaling of time")
    log("  and the invariance is a theorem, not a measurement. It is reported as")
    log("  an implementation check that the integrator honours it, and any")
    log("  threshold here would be a free constant several orders above the only")
    log("  values the arithmetic can produce.")
    _grid10 = {I: bands[I][0] / 10.0 for I in bands}
    _any10 = [I for I, f in _grid10.items() if 30 <= f <= 100]
    # The sharper statement, and unlike the tau = 10 row it is not forced by the
    # model's minimum period: for each cycle there is a tau-window that would put
    # IT in band, and those windows do not share a point. So no single membrane
    # constant makes this a gamma oscillator; band membership is per-input.
    _win = {I: (1000.0 / (100.0 * (1000.0 / f)), 1000.0 / (30.0 * (1000.0 / f)))
            for I, f in ((I, bands[I][0]) for I in bands)}
    _lo = max(a for a, b in _win.values())
    _hi = min(b for a, b in _win.values())
    log()
    log("  For each cycle, the tau that would place IT inside 30-100 Hz:")
    for I in _win:
        log(f"    I_ext = {I:4.2f}  tau in [{_win[I][0]:.2f}, {_win[I][1]:.2f}] ms")
    log(f"  intersection over all cycles: "
        f"{'empty' if _lo > _hi else f'[{_lo:.2f}, {_hi:.2f}]'}")
    _pts = sorted([a for a, b in _win.values()] + [b for a, b in _win.values()])
    _best = max(sum(1 for a, b in _win.values() if a <= p <= b) for p in _pts)
    _bestpt = max(_pts, key=lambda p: sum(1 for a, b in _win.values()
                                          if a <= p <= b))
    log(f"  widest overlap: {_best} of {len(_win)} cycles, near "
        f"tau = {_bestpt:.2f} ms")
    gate("G15c no single tau puts EVERY cycle in 30-100 Hz",
         _lo > _hi,
         f"the per-cycle tau windows do not share a point ({_lo:.2f} > {_hi:.2f}), "
         f"so band membership is a property of the operating point and the "
         f"chosen tau together, never of the model. The windows are not empty "
         f"one at a time: every cycle has a tau that puts IT in band, and the "
         f"widest overlap covers {_best} of {len(_win)} of them near "
         f"tau = {_bestpt:.2f} ms, so a tau making this model gamma SOMEWHERE "
         f"exists. What does not exist is one making it gamma throughout. At the textbook tau = 10 ms, "
         f"{len(_any10)} of the {len(_grid10)} cycles are in band; that row alone "
         f"is forced by the minimum period and is reported, not gated (they run "
         f"{min(_grid10.values()):.1f} to {max(_grid10.values()):.1f} Hz there). "
         "The band label on the table above is carried by the stipulated "
         "tau = 1 ms, which neural_gamma_cavity.py:54 sets bare")

    _amp10_bad, _max10_bad = {}, {}
    for I_ext in (1.15, 1.20, 1.50, 2.00, 2.50, 3.00, 3.05):
        def _rhs_bad(t, y, I_, tau):
            E, Iv = y
            return [-E / tau + sg(wEE * E - wEI * Iv + I_),
                    -Iv / tau + sg(wIE * E - wII * Iv)]
        s_ = solve_ivp(_rhs_bad, [0, 20000], [0.3, 0.3], args=(I_ext, 10.0),
                       max_step=0.5, dense_output=True, rtol=1e-9, atol=1e-11)
        t_ = np.linspace(15000, 20000, 200000)
        E_ = s_.sol(t_)[0]
        _amp10_bad[I_ext] = float(E_.max() - E_.min())
        _max10_bad[I_ext] = float(E_.max())
    log()
    log("  The committed two-population Jacobian at neural_gamma_cavity.py:91")
    log("  writes -1/tau + w*f', with the gain term NOT scaled. That is the")
    log("  Jacobian of no Wilson-Cowan model, and it agrees with the correct one")
    log("  only at tau = 1. E and I are FIRING RATES: the sigmoid's whole range")
    log("  is [0, 1], and the correct form keeps them there at every input. The")
    log("  leak-only form does not, which is the unambiguous sign it is not a")
    log("  rescaling of time but different equations:")
    log(f"    {'I_ext':>7s} {'amplitude':>11s} {'max E':>9s}")
    for I_ext in _amp10_bad:
        log(f"    {I_ext:7.2f} {_amp10_bad[I_ext]:11.6f} {_max10_bad[I_ext]:9.3f}")
    _escaped = [I for I, m in _max10_bad.items() if m > 1.0]
    _killed = [I for I, a in _amp10_bad.items() if a < 1e-6]
    gate("G15d the two code paths are different models, not different clocks",
         len(_escaped) > 0,
         f"at tau = 10 ms the leak-only form drives the firing rate to "
         f"{max(_max10_bad.values()):.1f} at I_ext = "
         f"{', '.join(f'{I:g}' for I in _escaped)}, outside the [0, 1] the "
         f"sigmoid can produce, and it kills the limit cycle at "
         f"{len(_killed)} of {len(_amp10_bad)} inputs while keeping a different "
         f"one at I_ext = 2.5. A rescaling of time can do none of those things, "
         "so the two forms are not the same equations. The cycle is not removed "
         "everywhere, which is why the escape from [0, 1] and not the cycle count "
         "is the predicate")
    log()
    gate("G15 the model runs on a limit cycle with no parameter change",
         len(bands) > 0,
         f"limit cycles at I_ext = {', '.join(f'{I:g}' for I in bands)}, at the "
         "same w, alpha and theta the document calls insufficient. The sampled "
         f"cycles run {min(1000.0 / f for f, _ in bands.values()):.1f} to "
         f"{max(1000.0 / f for f, _ in bands.values()):.1f} time constants, "
         "which is the hull of a sample and so a FLOOR on the model's range "
         "rather than the range: the period diverges at both folds. "
         "converting that to Hz needs a tau this model does not contain, so the "
         "BAND is not gated here (see G15c)")

    log()
    from scipy.optimize import fsolve as _fsolve

    def _maxre(v, I_ext):
        E, Iv = v
        dsg_ = lambda x: alpha_ * sg(x) * (1 - sg(x))
        xE, xI = wEE * E - wEI * Iv + I_ext, wIE * E - wII * Iv
        Jw = np.array([[-1 + wEE * dsg_(xE), -wEI * dsg_(xE)],
                       [wIE * dsg_(xI), -1 - wII * dsg_(xI)]])
        return float(np.max(np.linalg.eigvals(Jw).real))

    _CENSUS_GRID = 21

    def equilibria(I_ext, grid=_CENSUS_GRID):
        """EVERY equilibrium, by multi-start, not the one a fixed seed finds.

        A single seed does not track a branch: where three equilibria exist,
        [0.25, 0.25] returns the upper one, which reads as a loss of stability
        that has not happened. Residuals are CHECKED here, which is the defect
        G16 is about.
        """
        found = []
        for a in np.linspace(0, 1, grid):
            for b in np.linspace(0, 1, grid):
                with np.errstate(over='ignore', invalid='ignore'):
                    v = _fsolve(fp_residual, [a, b], args=(I_ext,))
                if (np.linalg.norm(fp_residual(v, I_ext)) < 1e-10
                        and not any(np.allclose(v, u, atol=1e-6) for u in found)):
                    found.append(v)
        return sorted(found, key=lambda v: v[0])

    def interior_stability(I_ext):
        v = _fsolve(fp_residual, [0.25, 0.25], args=(I_ext,))
        return v, _maxre(v, I_ext), float(np.linalg.norm(fp_residual(v, I_ext)))

    def has_cycle(I_ext):
        """Integrate and look, instead of testing a dict key."""
        s_ = solve_ivp(rhs, [0, 2000], [0.3, 0.3], args=(I_ext,), max_step=0.05,
                       dense_output=True, rtol=1e-9, atol=1e-11)
        t_ = np.linspace(1500, 2000, 50000)
        E_ = s_.sol(t_)[0]
        return float(E_.max() - E_.min()) > 1e-4

    log(f"  {'I_ext':>7s} {'interior fixed point':>22s} {'max Re(lambda)':>15s} "
        f"{'cycle?':>8s}")
    log("  " + "-" * 58)
    stab, seed_res = {}, {}
    for I_ext in (0.5, 0.8, 1.10, 1.15, 2.0, 3.0, 3.5):
        v, re_max, res_ = interior_stability(I_ext)
        stab[I_ext], seed_res[I_ext] = re_max, res_
        has = "yes" if has_cycle(I_ext) else "no"
        log(f"  {I_ext:7.2f} {f'({v[0]:.3f}, {v[1]:.3f})':>22s} {re_max:+15.3f} "
            f"{has:>8s}")
    # The error model, stated because a bare threshold is a number and not a
    # law: the residual of a converged Picard/fsolve point is bounded by the
    # solver's own tolerance times the map's local Lipschitz constant, and both
    # sit many orders above eps. What separates a converged row from a failed
    # one here is not a fine cut but a chasm: G16's failures sit at O(0.1), the
    # sigmoid's own range. So the gate is the RATIO between the two populations,
    # which is what makes it a statement rather than a threshold.
    _res_ok = max(seed_res.values())
    gate("G20 the single-seed stability rows ARE fixed points, by a margin that "
         "is a gap and not a cut",
         _res_ok < 1e-8 and _res_ok * 1e6 < 0.1,
         f"worst residual over the seven rows is {_res_ok:.1e}, which is more "
         f"than a million times below the O(0.1) residuals G16 reports for the "
         f"rows that genuinely fail to converge (the sigmoid's range is 1, so "
         f"0.1 is the scale a non-fixed-point sits at). The two populations do "
         f"not overlap, so no choice of threshold inside that gap changes the "
         f"verdict; the check was previously absent from the rows that feed the "
         "stability claim")
    log()
    log("  The rows above track a SEED, not a branch. Where three equilibria")
    log("  exist, [0.25, 0.25] returns the upper one. The census below returns")
    log(f"  all of them, from a {_CENSUS_GRID}x{_CENSUS_GRID} multi-start grid, "
        f"and the branch that was stable below is a different row:")
    log()
    log(f"  {'I_ext':>7s}  {'#eq':>3s}  equilibria, each with max Re(lambda)")
    log("  " + "-" * 66)
    census = {}
    for I_ext in (1.00, 1.10, 1.125, 1.13, 3.06, 3.07):
        eqs = equilibria(I_ext)
        census[I_ext] = [(v, _maxre(v, I_ext)) for v in eqs]
        row = "  ".join(f"({v[0]:.3f},{v[1]:.3f}){r:+.3f}" for v, r in census[I_ext])
        log(f"  {I_ext:7.3f}  {len(eqs):3d}  {row}")
    log()
    log("  NEITHER EDGE IS A HOPF OF THE EQUILIBRIUM THE CYCLE REPLACES. The")
    log("  lower branch is still STABLE at I = 1.125 and is destroyed with the")
    log("  middle saddle by I = 1.13, a saddle-node, which is where the cycle")
    log("  appears. The upper edge is the mirror image: one equilibrium at 3.06,")
    log("  three at 3.07 with the new one stable, and the cycle gone. Both edges")
    log("  are folds, and the period diverging at BOTH is the signature.")
    log()
    log("  Where the stable branch exists, up to the fold between I = 1.125 and")
    log("  1.126, is where the")
    log("  committed script's ~12 Hz ringing comes from. That reading is correct")
    log("  for the quiescent branch and was reported for the whole model.")
    _low_stable = census[1.125][0][1]   # the LOWEST-E branch, which is what the
    # detail names; min() over max-Re would be "the most stable", another row
    _n_eq = [len(census[i]) for i in (1.10, 1.125)]
    gate("G19 the cycle's lower edge is a fold, not a loss of stability",
         _low_stable < 0 and all(n == 3 for n in _n_eq)
         and len(census[1.13]) == 1 and not has_cycle(1.125)
         and has_cycle(1.13),
         f"at I = 1.125 there are {len(census[1.125])} equilibria, the lowest "
         f"still stable at {_low_stable:+.3f} and no cycle from the [0.3, 0.3] basin "
         f"(one trajectory per input, so a coexisting cycle with a small basin "
         f"would be missed); by I = 1.13 only "
         f"{len(census[1.13])} remains and the cycle IS there, which is now "
         f"tested rather than asserted. A Hopf would keep "
         "the equilibrium and change its stability; this destroys it")
    gate("G19b the upper edge is a fold too",
         len(census[3.06]) == 1 and len(census[3.07]) == 3
         and min(r for v, r in census[3.07]) < 0,
         f"{len(census[3.06])} equilibrium at I = 3.06 and "
         f"{len(census[3.07])} at 3.07, the new one stable at "
         f"{min(r for v, r in census[3.07]):+.3f}, with the cycle gone")
    log()
    log("  What the committed script reports instead is a linearisation at the")
    log("  quiescent branch, which is where its iteration converges:")
    v = picard(2.0)
    res2 = float(np.linalg.norm(fp_residual(v, 2.0)))
    log(f"    at I_ext = 2.0 the iteration returns ({v[0]:.3f}, {v[1]:.3f}), "
        f"residual {res2:.4f}")
    gate("G16 the reported fixed point at I_ext = 2.0 is not a fixed point",
         res2 > 0.01,
         "a residual that size against a sigmoid whose whole range is 1 means "
         "the point was never found, so the stability reading taken THERE is "
         "undetermined")
    log()
    log("  Scope that: the iteration does not fail everywhere. It fails inside")
    log("  the window where the equilibrium is unstable, which is where the one")
    log("  measured failure sits.")
    _grid = np.arange(0.0, 10.01, 0.25)
    _res = [(I, float(np.linalg.norm(fp_residual(picard(I), I)))) for I in _grid]
    _bad = [I for I, r in _res if r > 1e-8]
    _crit_rows = [r for I, r in _res if I in (0.0, 0.25, 0.5)]
    log(f"    converges on {len(_res) - len(_bad)} of {len(_res)} grid points; "
        f"fails on I_ext = {min(_bad):g} to {max(_bad):g}")
    log(f"    the three rows that oscillate and produce I_crit (I = 0, 0.25, "
        f"0.5) have residuals {[f'{r:.0e}' for r in _crit_rows]}")
    gate("G16b the Picard failure is local, and NOT what withdraws I_crit",
         len(_bad) < len(_res) / 2 and max(_crit_rows) == 0.0,
         f"{len(_res) - len(_bad)} of {len(_res)} points converge, and the rows "
         "that carry I_crit converge to residual exactly 0.0. I_crit is "
         "withdrawn because 0.00 is the first sampled grid point, not because "
         "the iteration missed it")
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
_pe_fwd = int(np.count_nonzero(Elec[np.ix_(pi, si)]))
_pe_rev = int(np.count_nonzero(Elec[np.ix_(si, pi)]))
log(f"  pharynx <-> soma, electrical   : {elec_cross} "
    f"({_pe_fwd} stored pharynx->soma, {_pe_rev} soma->pharynx, so each entry "
    f"here is one connection rather than half of a reciprocated pair)")
_zr_elec = [i for i in _zero_rows
            if Elec[i, :].any() or Elec[:, i].any()]
_zr_sia = [i for i in _zero_rows if _names[i][:3] in ("SIA", "SIB")]
log(f"  of the {len(_zero_rows)} empty chemical rows, {len(_zr_elec)} have "
    f"electrical entries; the {len(_zr_sia)} SIA/SIB among them carry "
    f"{min(int(Elec[i, :].astype(bool).sum()) for i in _zr_sia)} to "
    f"{max(int(Elec[i, :].astype(bool).sum()) for i in _zr_sia)} each")
gate("G24b an empty chemical row is not the same thing as a neuron leaving the "
     "model",
     len(_zr_elec) > 0 and len(_zr_sia) > 0,
     f"{len(_zr_elec)} of the {len(_zero_rows)} neurons with no outgoing "
     f"chemical synapse do have gap junctions INSIDE the 300, so 'its target is "
     f"muscle' is right for the motor neurons and wrong as a blanket reading. "
     f"The {len(_zr_sia)} SIA/SIB interneurons are the clear case: no chemical "
     f"output at all in this dataset, "
     f"{min(int(Elec[i, :].astype(bool).sum()) for i in _zr_sia)} to "
     f"{max(int(Elec[i, :].astype(bool).sum()) for i in _zr_sia)} electrical "
     "entries each, all of which the analysis discards. A dataset where every "
     "empty row were electrically isolated too would fail this gate")
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
log()
log("  The comparison that decides mixing is against the REAL matrix, not")
log("  between mixing levels: three levels that agree with each other agree")
log("  just as well when none of them has mixed at all. Rewiring LOWERS the")
log("  nullity, so the gap below is the real matrix minus the ensemble.")
_nul_real = N - gf_rank(W_raw, PRIMES[0])
_rng_one = np.random.default_rng(SEED + 2)   # OWN stream, same reason
_one = [N - gf_rank(degree_preserving_swap(W_raw, _rng_one, passes=1)[0],
                    PRIMES[0]) for _ in range(12)]
log(f"    real connectome : {_nul_real}")
log(f"    passes = 1      : {np.mean(_one):.2f}   (already a full sweep over")
log("                       the edge list, so this is NOT an unmixed control;")
log("                       it is the first point at which mixing saturates)")
log(f"    passes = 10     : {means[0]:.2f}")
gate("G18 the ensemble has moved away from the real matrix",
     _nul_real - means[0] > 4 * max(sds),
     f"passes = 10 gives nullity {means[0]:.2f} against the real connectome's "
     f"{_nul_real}, a gap of {_nul_real - means[0]:.2f} against a within-level "
     f"sd of up to {max(sds):.2f}. This is the direction that matters: the null "
     "is a different matrix from the one under test")
gate("G18b more mixing does not move it further",
     max(means) - min(means) < 2 * max(sds),
     f"the three mixing levels differ by {max(means) - min(means):.2f} against "
     f"a within-level sd of up to {max(sds):.2f}; this is a saturation check "
     "and is NOT evidence of mixing on its own")
log()

log("=" * 75)
n_gates = len(gate_names)   # counted from the calls, not from the log text
EXPECTED_GATES = 44
# Nine gates sit inside `try: import scipy...` blocks. Without scipy the run
# would print a smaller total and pass, and the document's gate count would
# quietly stop being true. The count is therefore itself a gate.
if n_gates != EXPECTED_GATES:
    failures.append(f"GATE CENSUS: {n_gates} gates ran, {EXPECTED_GATES} "
                    f"expected. A skipped optional-import block, or the "
                    f"constant needs updating with the new gate.")
log(f"{len(failures)} of {n_gates} gates FAILED" if failures
    else f"ALL {n_gates} GATES PASS")
for f in failures:
    log(f"  FAILED: {f}")
log("=" * 75)
log()
log("WHAT THIS DOES NOT SHOW")
log("-----------------------")
log("It does not test the SIGN-PATTERN half of condition (b) of")
log("PROOF_PALINDROME_NEURAL on a network that could satisfy it; G0b settles")
log("the question here on a count, which is a different and cheaper thing. It does not identify WHICH modes go unpaired, and")
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
