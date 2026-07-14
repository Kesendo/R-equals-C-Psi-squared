r"""grid_proof_sweep.py  (deterministic grid+CRT runner; see experiments/F89_SEED_EXISTENCE_REDUCTION.md).

DETERMINISTIC PROOF that the committed cross_form (simulations/cross_triple_orthogonality.py)
vanishes identically on the variety  V = {cos a1+cos a2+cos a3 = 0} x {cos b1+cos b2+cos b3 = 0},
over Q(i).  This runner drives the vectorized rigorous grid+CRT engine in core_grid.py and
checkpoints every (item, prime) result to simulations/_grid_proof_progress.txt.

======================================================================================= THE PROOF
Coordinates z_k = e^{i a_k}, w_k = e^{i b_k}; constraints Qz = z3^2 + Sz z3 + 1 = 0 and
Qw = w3^2 + Sw w3 + 1 = 0 (Sz = z1+1/z1+z2+1/z2, Sw = w1+1/w1+w2+1/w2).  Work in the field
K = Q(i)(z1,z2,w2)[z3]/(Qz) with w1 a FREE variable.  The proof has four steps; the first two are
exact-symbolic (in residue_assembly_close.py), the last two are this deterministic grid+CRT sweep.

  (1) ELIMINATION  [PROVED, exact].  cross_form == F0(w1) + F1(w1) w3 over K, the 288-term
      w3-elimination pinned to the committed cross_form (err ~ 1e-9 on V).  cross_form = 0 on V
      <=>  F0 == 0 AND F1 == 0 in K(w1).

  (2) SIMPLE POLES  [PROVED, exact].  Within every term the two w1-denominator factors are coprime
      over K (Sylvester resultant != 0), so every finite w1-pole of F0, F1 is SIMPLE.  Distinct
      factors of different terms DO share roots (Sylvester resultant = 0): union-find on the
      shared-root relation groups the 37 w1-carrying factors into 25 COMPONENTS (13 singletons + 12
      resultant-coupled quad-quad pairs).

  (3) NO FINITE POLE  [PROVED here, deterministic].  For each component C, F0's (and F1's) principal
      part at the roots of Pi_C = prod_{f in C} f comes only from the C-incident terms, and vanishes
      <=>  Pi_C | P_C, where P_C = sum_{t incident} Num_t * prod(C-factors except own)(distinct
      co-factors except own) is a DIVISION-FREE structured product.  Pi_C | P_C is tested by reducing
      P_C mod Pi_C in R_p[w1] = (GF(p)[z3]/(Qz))[w1] on a FULL tensor grid of (z1,z2,w2) and checking
      the remainder is 0 (both z3-components), for 17 primes.  25 components x 17 primes.

  (4) WINDOW + ENDPOINTS  [PROVED here, deterministic], contingent on (3).  With no finite pole, the
      w1 polynomial-part window of every summand is a subset of {-1,0,1} (STEP_ENDS_WINDOW, exact),
      so F0 (and F1) = c_-1/w1 + c_0 + c_1 w1 in K.  Evaluating F at three distinct fixed integers
      w1 in {2,3,5} gives three independent linear equations (generalized Vandermonde, det != 0) in
      (c_-1,c_0,c_1); each slice is N(w1=w_k) = sum_t Num_t * prod(all 43 distinct factors except the
      term's own two)|_{w1=w_k}, a division-free structured product in (z1,z2,w2), tested == 0 on the
      full grid for 17 primes.  Three slices x F0, F1 = 6 endpoint items.
      N(2)=N(3)=N(5)=0  =>  c_-1 = c_0 = c_1 = 0  =>  F0 = F1 = 0 in K(w1).

  CONCLUSION.  (1)+(2) exact, (3)+(4) deterministic over 17 primes with CRT (below):
      F0 == 0 and F1 == 0 in K(w1)  =>  cross_form == 0 on V over Q(i).                       QED

============================================================================ RIGOR OF THE GRID + CRT
TENSOR-GRID LEMMA (with skip-and-enlarge).  A Laurent polynomial of z1,z2,w2-degree <= d_v per axis
(after shifting min-exp to 0) that vanishes on a full product grid of d_v+1 DISTINCT GOOD residues
per axis is identically 0.  A grid point where a needed leading coefficient / denominator is a
non-unit (its Qz-norm is 0 in GF(p)) is UNUSABLE and simply not part of the grid; it is replaced by
a fresh good residue (the ring is resampled, over-provisioned by one per axis), which is legitimate
because the lemma only needs SOME full product grid of d_v+1 distinct good residues.  For 30-bit
primes the non-unit locus has density ~1/p, so bad points essentially never occur.

DEGREE BOUNDS (core_grid.span_bounds).  Each certified remainder / numerator is a Laurent polynomial
whose per-axis exponent span is bounded RIGOROUSLY by running the identical eval+reduce over a
tropical span ring (multiply adds exponent intervals, add/subtract unions them, the Sz = z1+1/z1+
z2+1/z2 norm contributes its (-1,1) intervals exactly).  The grid width is that bound + 1 per axis.
(Monomial-leading components -- linear singletons and two of the inc-32 quads -- use the monic
remainder, whose span is not inflated by the pseudo-remainder's leading-coefficient factor, roughly
halving the grid; multi-leading components use the rigorous pseudo-remainder span but reduce via the
cheaper monic algorithm at that width, valid because monic-remainder = 0 <=> pseudo-remainder = 0 at
every good grid point.)

CRT.  Every remainder / numerator is a polynomial over Z[i] after clearing denominators, with a
coefficient-height bound H ~ 2^447.  Vanishing mod each of the 17 primes p (30-bit, p = 1 mod 4) on
a grid that beats its degree forces it = 0 mod p; the product of the 17 primes is 2^510 (exact) > 2H, so
by CRT the coefficients (Gaussian integers of modulus < H) are 0 over Z[i].  --assert re-checks that
prod(primes) > 2H.

======================================================================================= HONEST STATUS
The ENGINE and BOUNDS are complete and validated (bit-exact against the original scalar pseudo-
remainder path on components C06, C07; controls that corrupt one numerator coefficient are caught
at every prime).  The SWEEP is heavy: after the two levers, per-prime cost ranges from ~9 s (cheap
singletons) to ~520 s (the 12 pairs, grid ~470k) to ~3500 s (the one inc-32 multi-lead quad C12,
grid ~3.17M) to ~1000 s (the six 1.2M endpoint slices).  The full sweep is a multi-day run; it is
CHECKPOINTED and RESUMABLE.  Items already proved on full grids over all 17 primes are asserted by
--assert; the giant items (C12 + the pairs + endpoints) may still be in progress -- see the progress
file and this runner's --status for the exact per-(item,prime) state.

Usage:
    python -u simulations/grid_proof_sweep.py            # run/resume the sweep (redirect to a file)
    python -u simulations/grid_proof_sweep.py --status   # print per-item progress, no compute
    python -u simulations/grid_proof_sweep.py --assert   # assert completed items + CRT sufficiency

Authors: Thomas Wicht and Claude, 2026-07-10.
"""
import sys, os, time, math, multiprocessing as mp
sys.path.insert(0, r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations")
import numpy as np
import core_grid as CG
from grid_proof_close import build_structure, component_plan, factor_key
import residue_assembly_close as RA

ROOT = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations"
PROGRESS = os.path.join(ROOT, "_grid_proof_progress.txt")
NPRIMES = 17
W1VALS = [2, 3, 5]   # fixed integer w1 slices for the endpoint window {-1,0,1} (same across primes)
RESERVE_CORES = 4    # logical cores left free for the user; --workers overrides

# controls: one representative per size class, control run EVERY prime
CONTROL_ITEMS = {"C00", "C04", "C06", "C08", "C10", "C13"}

def parse_workers():
    """--workers N  (default: logical cores - RESERVE_CORES, at least 1)."""
    if "--workers" in sys.argv:
        i = sys.argv.index("--workers")
        return max(1, int(sys.argv[i + 1]))
    return max(1, (os.cpu_count() or 1) - RESERVE_CORES)

def load_done():
    done = {}
    if os.path.exists(PROGRESS):
        with open(PROGRESS) as f:
            for ln in f:
                parts = ln.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    done[(parts[0], parts[1])] = parts[2]
    return done

def load_controls():
    ctrl = {}
    if os.path.exists(PROGRESS):
        with open(PROGRESS) as f:
            for ln in f:
                parts = ln.rstrip("\n").split("\t")
                if len(parts) >= 4 and parts[3].startswith("ctrl="):
                    ctrl[(parts[0], parts[1])] = (parts[3] == "ctrl=True")
    return ctrl

def append(item, prime, verdict, control, elapsed):
    with open(PROGRESS, "a") as f:
        f.write(f"{item}\t{prime}\t{verdict}\t{control}\t{elapsed:.2f}\n")
        f.flush(); os.fsync(f.fileno())

# ------------------------------------------------------------------ endpoint slices (window {-1,0,1})
# Given the residue result (F0, F1 have NO finite w1-pole), each of F0, F1 is a Laurent polynomial in
# w1 with powers in {-1, 0, 1} (STEP_ENDS_WINDOW), i.e.  F(w1) = c_-1/w1 + c_0 + c_1 w1.  Evaluating
# F at THREE distinct scalar w1 = w_k forces c_-1 = c_0 = c_1 = 0.  Over K, F_which = sum_t Num_t /
# (den0_t den1_t) with the SAME 43 distinct denominator factors as the residues; at w1 = w_k this is
# a w1-free K-element N(w_k)/D(w_k), D = prod(distinct factors|_{w1=w_k}) != 0 at generic w_k.  The
# certified numerator is  N(w_k) = sum_t Num_t * prod(distinct factors except the term's own two),
# all at w1 = w_k -- a DIVISION-FREE structured product in (z1,z2,w2) over z3, degree ~108 (NOT the
# ~1000 of the pseudo-division form).  We test N(w_k) == 0 on the (z1,z2,w2) grid, all primes.
def build_endpoint_terms(F):
    """Return per which in {0,1}: list of (Num_t, key0, key1) and the distinct-factor dict {key:poly}."""
    out = {}
    for which in (0, 1):
        terms = []; distinct = {}
        for (A, B, dens) in F:
            Num = A if which == 0 else B
            k0, k1 = factor_key(dens[0]), factor_key(dens[1])
            distinct.setdefault(k0, dens[0]); distinct.setdefault(k1, dens[1])
            if not Num: continue
            terms.append((Num, k0, k1))
        out[which] = (terms, distinct)
    return out

def endpoint_span(terms, distinct):
    """SPAN-ring: N = sum Num_t * prod(distinct != own two).  w1 folds to a constant (no z-span)."""
    R = CG.SpanRing()
    dval = {k: _ev0(R, d) for k, d in distinct.items()}
    Ptot = (None, None)  # not needed for span; sum directly
    N = (None, None)
    for (Num, k0, k1) in terms:
        term = _ev0(R, Num)
        for k, dv in dval.items():
            if k == k0 or k == k1: continue
            term = CG.rp_mul(R, term, dv)
        N = (R.fadd(N[0], term[0]), R.fadd(N[1], term[1]))
    b = {0: (0, 0), 1: (0, 0), 4: (0, 0)}; got = False
    for part in N:
        if part is None: continue
        got = True
        for ax in CG._AX: b[ax] = CG._ivun(b[ax], part[ax])
    return b, got

def endpoint_num_array(R, terms, distinct):
    """ARRAY-ring N via Ptot * own^{-1}: N = Ptot * sum Num_t * inv(den0)*inv(den1).  (NUM, ok)."""
    dval = {k: _ev0(R, d) for k, d in distinct.items()}
    Ptot = (np.ones_like(R.Sz), None)
    for dv in dval.values():
        Ptot = CG.rp_mul(R, Ptot, dv)
    dinv = {}
    for k, dv in dval.items():
        inv, ok = CG.rp_inv_array(R, dv)
        if not np.all(ok):
            return None, False
        dinv[k] = inv
    S = (None, None)
    for (Num, k0, k1) in terms:
        term = CG.rp_mul(R, _ev0(R, Num), dinv[k0])
        term = CG.rp_mul(R, term, dinv[k1])
        S = (R.fadd(S[0], term[0]), R.fadd(S[1], term[1]))
    NUM = CG.rp_mul(R, Ptot, S)
    return NUM, True

def _ev0(R, poly):
    """Evaluate a w1-free (or w1-folded) 6-var poly to a single R_p element (a0,a1)."""
    d = CG._eval_poly(R, poly)
    a0 = a1 = None
    for w, ab in d.items():
        a0 = R.fadd(a0, ab[0]); a1 = R.fadd(a1, ab[1])
    return (a0, a1)

def run_endpoint_prime(terms, distinct, W, p, w1val, seedbase=200, do_control=False):
    t0 = time.time()
    gz1, gz2, gw2 = W[0] + 1, W[1] + 1, W[4] + 1
    for attempt in range(6):
        R = CG.build_ring(p, gz1, gz2, gw2, seedbase + 37 * attempt, CG.tonelli_i(p))
        R.w1val = w1val
        NUM, ok = endpoint_num_array(R, terms, distinct)
        if not ok:
            gz1 += 1; gz2 += 1; gw2 += 1; continue
        nz = (not R.fis_zero(NUM[0])) or (not R.fis_zero(NUM[1]))
        control_ok = None
        if do_control:
            t2 = list(terms); Num0, a, b = t2[0]
            corrupted = dict(Num0); k0 = tuple([0] * RA.NVARS)
            corrupted[k0] = corrupted.get(k0, RA.ZERO) + RA.ONE
            t2[0] = (corrupted, a, b)
            NUMc, okc = endpoint_num_array(R, t2, distinct)
            control_ok = bool(okc and ((not R.fis_zero(NUMc[0])) or (not R.fis_zero(NUMc[1]))))
        return (not nz), control_ok, time.time() - t0
    raise RuntimeError("endpoint skip-and-enlarge exceeded attempts")

# ------------------------------------------------------------------ shared plan builder
def build_plan():
    """Structure + per-component plans/widths + endpoint slices + the grid-sized worklist.
    Returns (F, factors, incidence, order, primes, plans, ends, work).  work is a list of
    (grid, item, kind, ref, ctrl) with kind in {'comp','end'} and ref = idx or (which, wk)."""
    F, factors, incidence, comps, edges = build_structure()
    order = sorted(comps, key=lambda c: (len(c), sum(len(incidence[k]) for k in c)))
    primes = CG.gen_primes(NPRIMES)
    plans = {}
    for idx, C in enumerate(order):
        plans[idx] = CG.component_widths(F, factors, incidence, C)   # (monic, W, inc2, distinct)
    eterms = build_endpoint_terms(F)
    ends = {}
    for which in (0, 1):
        terms, distinct = eterms[which]
        b, got = endpoint_span(terms, distinct)
        Wend = {0: b[0][1] - b[0][0], 1: b[1][1] - b[1][0], 4: b[4][1] - b[4][0]}
        for wk in W1VALS:
            ends[(which, wk)] = (terms, distinct, Wend)
    work = []
    for idx, C in enumerate(order):
        monic, W, inc2, distinct = plans[idx]
        g = (W[0] + 1) * (W[1] + 1) * (W[4] + 1)
        item = f"C{idx:02d}"
        work.append((g, item, "comp", idx, item in CONTROL_ITEMS))
    for (which, wk), (terms, distinct, Wend) in ends.items():
        g = (Wend[0] + 1) * (Wend[1] + 1) * (Wend[4] + 1)
        work.append((g, f"E{which}_w{wk}", "end", (which, wk), (wk == W1VALS[0])))
    return F, factors, incidence, order, primes, plans, ends, work

# ------------------------------------------------------------------ multiprocessing worker
# One (item, prime) task per call; structures rebuilt ONCE per process in the initializer (never
# per task).  The worker returns a result line; only the PARENT appends to the progress file.
_W = {}

def _init_worker():
    t0 = time.time()
    F, factors, incidence, order, primes, plans, ends, work = build_plan()
    _W.update(F=F, factors=factors, incidence=incidence, order=order,
              plans=plans, ends=ends, setup=time.time() - t0)

def _worker(task):
    kind, item, ref, p, ctrl = task
    F, factors, incidence = _W["F"], _W["factors"], _W["incidence"]
    try:
        if kind == "comp":
            C = _W["order"][ref]
            monic, W, inc2, distinct = _W["plans"][ref]
            proved, control_ok, secs = CG.run_component_prime(
                F, factors, incidence, C, p, monic, W, inc2, distinct, do_control=ctrl)
        else:
            terms, distinct, Wend = _W["ends"][ref]
            proved, control_ok, secs = run_endpoint_prime(terms, distinct, Wend, p, ref[1], do_control=ctrl)
        verdict = "PROVED" if proved else "FAIL"
        cstr = f"ctrl={control_ok}" if ctrl else "-"
    except Exception as e:
        verdict = f"ERROR:{type(e).__name__}"; cstr = "-"; secs = 0.0
    return (item, p, verdict, cstr, secs, _W["setup"])

# ------------------------------------------------------------------ main sweep (parallel)
def main():
    print("=" * 80); print("DETERMINISTIC GRID SWEEP (25 residue components + 3 endpoint slices x 17 primes)")
    print("=" * 80)
    t0 = time.time()
    F, factors, incidence, order, primes, plans, ends, work = build_plan()
    print(f"[structure] {len(order)} components; primes ({NPRIMES}, 30-bit =1mod4): {primes[0]}..{primes[-1]}", flush=True)

    if "--status" in sys.argv or "--assert" in sys.argv:
        return report(order, F, factors, incidence, primes, with_heights=("--assert" in sys.argv))
    if "--heights" in sys.argv:
        return height_report(F, factors, incidence, order, ends, primes)

    for idx, C in enumerate(order):
        monic, W, inc2, distinct = plans[idx]
        g = (W[0] + 1) * (W[1] + 1) * (W[4] + 1)
        print(f"  C{idx:02d} monic={monic} W={(W[0], W[1], W[4])} grid={g:,}", flush=True)
    for which in (0, 1):
        terms, distinct, Wend = ends[(which, W1VALS[0])]
        print(f"  E{which} terms={len(terms)} distinct={len(distinct)} W={(Wend[0], Wend[1], Wend[4])} "
              f"grid={(Wend[0]+1)*(Wend[1]+1)*(Wend[4]+1):,} (x3 w1-slices)", flush=True)

    done = load_done()
    print(f"[resume] {len(done)} (item,prime) results already in progress file", flush=True)

    # ---- flat (item, prime) task list, ordered LARGEST GRID FIRST so C12 + the endpoint slices
    #      start immediately and the tail is minimized; skip already-checkpointed pairs.
    work_sorted = sorted(work, key=lambda t: -t[0])
    tasks = []
    for g, item, kind, ref, ctrl in work_sorted:
        for p in primes:
            if (item, str(p)) in done: continue
            tasks.append((kind, item, ref, p, ctrl))
    if not tasks:
        print("[worklist] nothing to do -- every (item, prime) already checkpointed.", flush=True)
        return report(order, F, factors, incidence, primes, with_heights=True)

    workers = parse_workers()
    print(f"[worklist] {len(tasks)} (item,prime) tasks remaining; grid sizes "
          f"{work_sorted[0][0]:,} (first) .. {work_sorted[-1][0]:,} (last)", flush=True)
    print(f"[pool] {workers} workers ({os.cpu_count()} logical cores, {RESERVE_CORES} reserved)", flush=True)

    n = 0; setup_reported = False
    with mp.Pool(workers, initializer=_init_worker) as pool:
        for (item, p, verdict, cstr, secs, setup) in pool.imap_unordered(_worker, tasks, chunksize=1):
            append(item, p, verdict, cstr, secs)
            n += 1
            if not setup_reported:
                print(f"[worker-setup] {setup:.1f}s per process (structures rebuilt once, not per task)", flush=True)
                setup_reported = True
            rate = n / (time.time() - t0) * 3600.0
            print(f"  [{n}/{len(tasks)}] {item} p={p} {verdict} {cstr} {secs:.1f}s "
                  f"(~{rate:.1f} tasks/h)", flush=True)

    print(f"[TOTAL] {time.time()-t0:.1f}s for {n} tasks", flush=True)
    report(order, F, factors, incidence, primes, with_heights=True)

# ------------------------------------------------------------------ height bound (owed rigor)
def item_heights(F, factors, incidence, order, ends):
    """Rigorous per-item integer coefficient-height bound H (see core_grid.height_bounds)."""
    H = {}
    for idx, C in enumerate(order):
        H[f"C{idx:02d}"] = CG.component_height(F, factors, incidence, C)
    for (which, wk), (terms, distinct, Wend) in ends.items():
        H[f"E{which}_w{wk}"] = CG.endpoint_height(terms, distinct, wk)
    return H

def height_report(F, factors, incidence, order, ends, primes):
    print("=" * 80); print("COEFFICIENT-HEIGHT BOUND (independent, tight, per item) -- CRT sufficiency")
    print("=" * 80)
    H = item_heights(F, factors, incidence, order, ends)
    prodv = 1
    for p in primes: prodv *= p
    log2prod = math.log2(prodv)
    allp = CG.gen_primes(40)
    print(f"{'item':<9} {'log2(H)':>9}   {'2H<prod?':>8}   primes_needed(prod>2H)")
    Hmax = 0; worst = None; global_ok = True
    for item in sorted(H):
        h = H[item]; Hmax = max(Hmax, h)
        if h == max(H.values()): worst = item
        ok = prodv > 2 * h
        global_ok = global_ok and ok
        need = 1; k = 0
        while need <= 2 * h: need *= allp[k]; k += 1
        print(f"{item:<9} {math.log2(h):>9.1f}   {str(ok):>8}   {k}")
    print("-" * 80)
    print(f"product of {len(primes)} primes = 2^{log2prod:.1f}  (each ~30-bit, p = 1 mod 4)")
    print(f"worst item {worst}: log2(H) = {math.log2(Hmax):.1f}; 2H = 2^{math.log2(2*Hmax):.1f} < 2^{log2prod:.1f}")
    need = 1; k = 0
    while need <= 2 * Hmax: need *= allp[k]; k += 1
    print(f"VERDICT: {'17 primes SUFFICE for all items' if global_ok else 'INSUFFICIENT -- extend primes'}"
          f" (worst item needs {k} primes; {NPRIMES} used)")
    return global_ok, Hmax, k

def assert_realness(F):
    """Soundness prerequisite for the 2H CRT bound.  The certified numerators are built by the
    +,-,* of the reduction from the elimination outputs A (F0/𝔉0-parts), B (F1/𝔉1-parts), and the
    denominator factors of build_F().  Those base coefficients are RATIONAL (imaginary part 0), and
    real coefficients are closed under +,-,*, so every certified remainder and endpoint numerator is
    real.  Hence vanishing mod the 17 primes lifts by ORDINARY integer CRT (prod > 2H), not the
    prod > H^2 a genuinely Gaussian remainder would force.  Asserted exactly here, over every base
    coefficient, so --assert no longer merely assumes realness."""
    n = 0
    for (A, B, dens) in F:
        for poly in (A, B, *dens):
            for k, c in poly.items():
                if c.im != 0:   # explicit raise, not assert: must survive python -O for a soundness gate
                    raise AssertionError(f"non-rational base coefficient at monomial {k}: im = {c.im} "
                                         f"-- the 2H CRT bound would be INVALID (H^2 needed)")
                n += 1
    return n

def report(order, F, factors, incidence, primes, with_heights=False):
    done = load_done()
    items = [f"C{idx:02d}" for idx in range(len(order))] + [f"E{w}_w{wk}" for w in (0, 1) for wk in W1VALS]
    print("-" * 80)
    allok = True
    for item in items:
        res = [done.get((item, str(p))) for p in primes]
        nproved = sum(1 for r in res if r == "PROVED")
        ndone = sum(1 for r in res if r is not None)
        status = "PROVED@all" if nproved == len(primes) else f"{nproved}/{len(primes)} proved, {ndone} done"
        if nproved != len(primes): allok = False
        print(f"  {item}: {status}")
    ctrl = load_controls()
    nctrl = len(ctrl); nctrl_ok = sum(1 for v in ctrl.values() if v)
    print("-" * 80)
    print(f"CONTROLS (corrupt one numerator coeff -> remainder must be nonzero): "
          f"{nctrl_ok}/{nctrl} discriminate" + (" (ALL OK)" if nctrl and nctrl_ok == nctrl else ""))
    if with_heights:
        nreal = assert_realness(F)
        print(f"REALNESS: all {nreal} base coefficients (F0, F1, denominators) rational (im=0) "
              f"-> ordinary integer CRT, prod > 2H is the valid sufficiency bound")
        crt_ok, Hmax, kneed = height_report(F, factors, incidence, order, _ends_from(F), primes)
    else:
        prodv = 1
        for p in primes: prodv *= p
        crt_ok = prodv > (1 << 449)   # coarse gate for --status; --assert recomputes the tight bound
        print(f"CRT: {len(primes)} primes, product ~ 2^{math.log2(prodv):.0f} "
              f"(tight per-item bound checked by --assert / --heights)")
    proof_complete = allok and crt_ok and nctrl and nctrl_ok == nctrl
    print(f"VERDICT: {'PROVED: F0=F1=0 => cross_form == 0 on V over Q(i) (all items, 17 primes, CRT>2H)' if proof_complete else 'INCOMPLETE (see per-item; sweep still running)'}")
    return proof_complete

def _ends_from(F):
    """Rebuild the endpoint slice dict {(which,wk):(terms,distinct,Wend)} for the height report."""
    eterms = build_endpoint_terms(F)
    ends = {}
    for which in (0, 1):
        terms, distinct = eterms[which]
        b, got = endpoint_span(terms, distinct)
        Wend = {0: b[0][1] - b[0][0], 1: b[1][1] - b[1][0], 4: b[4][1] - b[4][0]}
        for wk in W1VALS:
            ends[(which, wk)] = (terms, distinct, Wend)
    return ends

if __name__ == "__main__":
    mp.freeze_support()
    main()
