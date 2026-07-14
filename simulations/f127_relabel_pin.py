r"""A6 -- the w2-relabel pin (hostile-referee patch).

CLAIM UNDER TEST.  The committed structure results of `residue_assembly_close.py`
    STEP_COPRIME   : within every eliminated term the two w1-denominator factors are coprime
                     in the distinguished Laurent variable (=> every pole SIMPLE),
    STEP_ENDS_WINDOW: every summand's distinguished-variable polynomial-part window is a
                     subset of {-1,0,1},
    STEP_POLES     : 43 distinct denominator factors, classified 6 (var-free) + 6 (linear) +
                     31 (quad),
hold VERBATIM when w2 (b-side index 4) is the distinguished Laurent variable instead of w1
(b-side index 3), by the b1 <-> b2 symmetry.  The referee wants this LITERALLY RUN.

MECHANISM (cheapest sound pin, option (i)).  The residue machinery in `residue_assembly_close`
treats variable index W1 = 3 as the distinguished free Laurent variable.  We build the
eliminated form F once, then apply the exponent-index swap sigma : 3 <-> 4 (i.e. w1 <-> w2)
to every polynomial (numerators A, B and both denominator factors) of every term.  In the
swapped data sigma(F), index 3 now carries the ORIGINAL w2 content, so running the SAME
(unmodified) STEP_* functions -- which distinguish index 3 -- is literally the w2-frame
analysis of the original object.

WHY THE SWAP IS SOUND (not a tautological relabel).
  * The w3-elimination `process_term` only ever touches index 2 (z3, via SZ) and index 5
    (w3, via SW = w1+1/w1+w2+1/w2).  SW is INVARIANT under 3<->4 and sigma touches neither
    index 2 nor 5, so sigma commutes with the elimination: sigma(process_term(t)) =
    process_term(sigma(t)).  The swapped dens stay w3-free and Qz-reduced.
  * The b1<->b2 symmetry is verified live: eval_F(sigma(F), vals) == eval_F(F, vals) on V
    (both equal cross_form; the swap is a genuine symmetry of the object, not merely a
    symbol rename), and |eval_F(sigma(F), vals)| ~ 0 on V (the swapped elimination is valid).

If the swapped run produced a DIFFERENT inventory (not 43 = 6+6+31, or a double pole, or a
window outside {-1,0,1}), step 7 of the chain would break.  It does not.

Authors: Thomas Wicht and Claude, 2026-07-14.
"""
import sys
import time

sys.path.insert(0, r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations")

import residue_assembly_close as R
from residue_assembly_close import (
    build_F, STEP_POLES, STEP_COPRIME, STEP_ENDS_WINDOW,
    eval_F, _angle_variety_point, classify_factor,
)
from halfangle_residue_proof import NVARS

import random
from collections import Counter

I3, I4 = 3, 4   # w1, w2 exponent indices


def swap_key(k):
    kk = list(k)
    kk[I3], kk[I4] = kk[I4], kk[I3]
    return tuple(kk)


def swap_poly(poly):
    """Relabel w1 <-> w2 in a sparse Laurent poly (swap exponent indices 3 and 4)."""
    out = {}
    for k, c in poly.items():
        out[swap_key(k)] = c
    return out


def swap_F(F):
    """Apply w1<->w2 to every polynomial of every eliminated term."""
    return [(swap_poly(A), swap_poly(B), [swap_poly(d) for d in dens]) for (A, B, dens) in F]


def swap_vals(vals):
    v = list(vals)
    v[I3], v[I4] = v[I4], v[I3]
    return v


def pin_swap_is_symmetry(F, Fs, ntrials=25, seed=77):
    """eval_F(Fs, vals) == eval_F(F, vals) on V  (b1<->b2 symmetry of the object), and
    |eval_F(Fs, vals)| ~ 0 on V (the swapped elimination is a valid eliminant)."""
    rng = random.Random(seed)
    worst_sym = 0.0        # |eval(Fs, v) - eval(F, v)| ; symmetry of cross_form
    worst_relabel = 0.0    # |eval(Fs, v) - eval(F, swap(v))| ; must be ~0 by construction
    worst_absFs = 0.0      # |eval(Fs, v)| on V ; a wrong elimination is O(1)
    done = 0
    while done < ntrials:
        got = _angle_variety_point(rng)
        if got is None:
            continue
        vals, _, _ = got
        try:
            eF = eval_F(F, vals)
            eFs = eval_F(Fs, vals)
            eFsw = eval_F(F, swap_vals(vals))
        except ZeroDivisionError:
            continue
        worst_sym = max(worst_sym, abs(eFs - eF))
        worst_relabel = max(worst_relabel, abs(eFs - eFsw))
        worst_absFs = max(worst_absFs, abs(eFs))
        done += 1
    return worst_sym, worst_relabel, worst_absFs, done


def inventory_counts(factors):
    return Counter(classify_factor(p) for p in factors.values())


def main():
    print("=" * 78)
    print("A6 RELABEL PIN: STEP_COPRIME / STEP_ENDS_WINDOW / STEP_POLES with w2 distinguished")
    print("=" * 78)
    t0 = time.time()

    print("[build] eliminating w3 (w1-frame) ...")
    tb = time.time()
    F = build_F()
    print(f"  F built: {len(F)} terms ({time.time()-tb:.1f}s)")

    Fs = swap_F(F)
    print(f"[swap] applied w1<->w2 (indices 3<->4) to all {len(Fs)} terms")

    print("[pin] verifying the swap is a genuine symmetry + valid eliminant on V ...")
    ws, wr, wa, nn = pin_swap_is_symmetry(F, Fs)
    print(f"  |eval_F(Fs) - eval_F(F)| on {nn} V-points (b1<->b2 symmetry) max {ws:.2e}")
    print(f"  |eval_F(Fs) - eval_F(F, swap(vals))| (relabel identity)   max {wr:.2e}")
    print(f"  |eval_F(Fs)| on V (buggy elim = O(1))                     max {wa:.2e}")
    assert ws < 1e-6, "b1<->b2 symmetry FAILS -- object not symmetric under the swap"
    assert wr < 1e-9, "relabel identity FAILS -- swap not implemented as a pure index swap"
    assert wa < 1e-5, "swapped elimination does NOT vanish on V -- swap broke the eliminant"

    # ---- reference (w1-frame) inventory, for the side-by-side ----
    print("\n[REFERENCE] w1-frame inventory (unswapped):")
    factors0, incidence0, ptk0 = STEP_POLES(F)
    cls0 = inventory_counts(factors0)

    print("\n[SWAPPED] w2-frame inventory:")
    factors, incidence, per_term_keys = STEP_POLES(Fs)
    cls = inventory_counts(factors)

    # STEP_POLES structure verdict: 43 = 6 + 6 + 31
    nfac = len(factors)
    nfree = cls.get("w1free", 0)
    nlin = cls.get("linear", 0)
    nquad = cls.get("quad", 0)
    print(f"\n[VERDICT STEP_POLES] swapped: {nfac} factors = {nfree} free + {nlin} linear"
          f" + {nquad} quad; reference: {len(factors0)} = {cls0.get('w1free',0)}"
          f" + {cls0.get('linear',0)} + {cls0.get('quad',0)}")
    assert nfac == 43, f"expected 43 distinct factors, got {nfac}"
    assert (nfree, nlin, nquad) == (6, 6, 31), f"expected (6,6,31), got {(nfree, nlin, nquad)}"
    assert dict(cls) == dict(cls0), "swapped inventory class-counts differ from the w1-frame"

    # STEP_COPRIME: all poles simple (raises internally if any term shares a root)
    print("\n[SWAPPED] STEP_COPRIME (all poles simple):")
    STEP_COPRIME(Fs)

    # STEP_ENDS_WINDOW: window subset of {-1,0,1} (asserts internally)
    print("\n[SWAPPED] STEP_ENDS_WINDOW (poly-part window):")
    lo, hi = STEP_ENDS_WINDOW(Fs)
    assert lo >= -1 and hi <= 1, f"window [{lo},{hi}] escapes {{-1,0,1}}"

    print("-" * 78)
    print("A6 PIN OK.  With w2 the distinguished Laurent variable, the committed structure")
    print("holds VERBATIM: 43 factors = 6 + 6 + 31, every pole SIMPLE (STEP_COPRIME), the")
    print(f"poly-part window = [{lo},{hi}] subset of {{-1,0,1}} (STEP_ENDS_WINDOW).  No asymmetry.")
    print(f"[TOTAL] {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
