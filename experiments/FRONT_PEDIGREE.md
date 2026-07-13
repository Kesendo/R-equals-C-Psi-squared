# The Pedigree of the Front: how often a surviving walker was caught, and how long ago

**Date**: 2026-07-14
**Status**: Computationally verified (three methods to machine precision; born as a play session with the renewal cut)
**Script**: `simulations/cone_front_pedigree.py` (output in `simulations/results/cone_defect_arrival/front_pedigree.txt`)
**Question**: the renewal ladder of F126 can be resolved by catch count j, an observable the density-matrix engine hides (a quantum-jump unraveling would sample it as jump counts; the ladder holds it exactly, no sampling). What is the pedigree of the front survivors: how many times was a surviving walker caught on its way, and when was its last rebirth?

## The two answers up front

1. **Front survivors are rebirths, mostly, and lightly caught, always.** At the hardware-flavoured dose (n = 50, γ = 0.05, K = 1.25) only 29% of the surviving front population was never caught; the mean catch count of a survivor is 1.38, against Γt* = 5 for an unconditioned walker: the front selects walkers caught about a quarter as often as the crowd (⟨j⟩/Γt* ≈ 0.28–0.35 across the table), but it is not a purity selection, the majority of what arrives has been caught and re-released at least once.

2. **The rebirths are of every age.** The last catch of the caught front weight is spread broadly over the whole trip (quartiles at 8%, 30%, 66% of the trip before arrival, the 90th percentile at 91%). The halo that rescues the front is not young; it carries collapses from the entire way. This is the fourth follow-up's bulk-dominance finding made visible: the refill integral I₁ draws ~70% of its value from the cone interior, not the caustic sliver, so the re-seeding happens everywhere along the road, not at the edge.

## Method

The F126 renewal ladder (`docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md`), resolved by refill order: S⁽⁰⁾ = |G_{n0}|² (the never-caught term, the coherent front) and S⁽ʲ⁾ = Γ·∫₀ᵗ ds Σ_m |G_{nm}(t−s)|²·S⁽ʲ⁻¹⁾_m(s), so P(j) = S⁽ʲ⁾/Σ_j S⁽ʲ⁾ at the observation point is the catch-count distribution of the arriving population (the universal e^{−Γt} cancels in the ratio). Reading point: site n at the band-edge time t* = n/(2J), interior seed, infinite chain. Three methods agree to machine precision (per-order momentum ladder via the Graf kernel, the implicit full Volterra, a displacement-space Volterra; section [1] of the script, relative deviations ≤ 10⁻¹⁴). The last-catch age density is the final-segment integrand Γ·Σ_m |G_{nm}(τ)|²·S_m(t*−τ).

## Results

**The pedigree table** (script section [2]; θ-free, no threshold anywhere):

| n | γ | K | P(never caught) | ⟨j⟩ | Γt* (free walker) | P(j), j = 0..5 |
|---|---|---|---|---|---|---|
| 20 | 0.05 | 0.50 | 0.529 | 0.67 | 2.0 | 0.529, 0.318, 0.114, 0.031, 0.007, 0.001 |
| 30 | 0.05 | 0.75 | 0.425 | 0.92 | 3.0 | 0.425, 0.336, 0.160, 0.057, 0.016, 0.004 |
| 50 | 0.05 | 1.25 | 0.287 | 1.38 | 5.0 | 0.287, 0.320, 0.216, 0.109, 0.045, 0.016 |
| 30 | 0.10 | 1.50 | 0.159 | 2.10 | 6.0 | 0.159, 0.251, 0.239, 0.170, 0.099, 0.049 |
| 30 | 0.025 | 0.375 | 0.663 | 0.43 | 1.5 | 0.663, 0.262, 0.062, 0.011, 0.002, 0.000 |

**The identity that ties the pedigree to the ceiling story.** P(never caught) = g_coh/g = e^{−(4−A_same)·K} exactly (verified to 10⁻⁹ per row; an identity, not a test): the survival exponent's deficit below the bare rate, the object the third follow-up's ceiling A_∞(γ) bounds, is literally the per-dose rate at which pure pedigree dies out among survivors. What the ceiling keeps of the front's bill is exactly what the pedigree shows as rebirths.

**One reconciliation, and it is the arc's lesson again.** The pedigree's implied A_same (2.86 at n = 30, 3.00 at n = 50, γ = 0.05) at first looked like a handshake break against the committed A_same = 2.7316 (n = 40, section [5] of `cone_front_survival_asymptote.py`). It is not: the two read the same S-ratio at two different times both loosely called t*₀, the band-edge n/(2J) here (where Γt = 4K exactly, the clean pedigree constant) versus the Airy caustic peak (n + 0.809·n^{1/3})/(2J) there. Both methods agree on the underlying object to 2·10⁻⁵ at every point tested; evaluated at each other's time they reproduce each other's numbers to four decimals (script section [1]). The observable is an appointment, for the fourth time in this arc.

**The last-catch age** (script section [3], n = 50, γ = 0.05): quartiles of the last-catch age among the caught front weight at 7.8%, 29.7%, 66.2% of the trip, the 90th percentile at 90.8%; the never-caught fraction 0.287 consistent with the table.

## What this means, and what it does not

The refund of the walk-time arc gains a face. The front that survives the watching is not the untouched remnant of the original wave; at any appreciable dose it is mostly re-released population, of every age, riding the same schedule. The two summary numbers to keep: at K ≈ 1 a surviving front walker was typically caught once or twice (never five times, the free-walker bill), and its last rebirth lies anywhere on the road, spread across the whole trip rather than clustered near arrival (median 30% of the trip before it). Nothing here is a new law: every number is a reading of the proven F126 object, and the one identity stated is algebra. The pedigree observable itself (P(j), the catch-count resolution) is the genuinely new instrument: the ladder representation holds it exactly and in closed form, where the density-matrix engine hides it entirely and a trajectory unraveling could only sample it.

## Links

- The object read here: `docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md` (F126 in `docs/ANALYTICAL_FORMULAS.md`); live witness `inspect --root renewal`; the engine-side cut `compute/MirrorWorld/Renewal.cs`.
- The arc this plays on: `experiments/COUPLING_DEFECT_WALK_TIME_STEP.md` (the ceiling, the bulk-dominance of I₁, the trichotomy of readings).
- The seeing this quantifies: `reflections/ON_THE_REFUND.md` (every observation is a rebirth; the pedigree is that sentence with numbers).
