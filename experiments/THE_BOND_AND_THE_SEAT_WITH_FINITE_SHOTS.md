# A bond and a seat with finite Z-position counts

**Status:** local, ideal-model N=7 simulation at fixed `gamma0=0.05`,
`J_hop=0.075` (`Q=1.5`), 2026-09-23. No F number, hardware flight, or
unrestricted bond-recovery claim.
**Producer:** [finite-shot code](../simulations/handshake_bond_seat_finite_shots.py),
run with `python simulations/handshake_bond_seat_finite_shots.py`.
**Retained run:** [all schedules, Fisher spectra, 13-by-13 confusion matrices,
per-hypothesis Wilson intervals, and paired counts](../simulations/results/handshake_bond_seat_finite_shots_run.txt).
**Controls:** `python -m pytest -q simulations/test_handshake_bond_seat_finite_shots.py simulations/test_handshake_bond_seat_fixed_gamma.py simulations/test_handshake_bond_seat_readout.py`.

The source sweep for this reading found F124's bond-to-mode matrix in the
[F-registry](../docs/ANALYTICAL_FORMULAS.md) and [its proof](../docs/proofs/PROOF_HANDSHAKE_TRANSITION_INVARIANT.md),
and F157's blind-seat count in the registry and [its source](THE_SEAT_THAT_CUTS.md);
the typed
[F124 claim](../compute/RCPsiSquared.Diagnostics/Ptf/BandEdgeTransitionInvariantClaim.cs)
and [F157 claim](../compute/RCPsiSquared.Core/Symmetry/SeatCutBlindnessClaim.cs)
own those distinct statements, not this finite-shot test. The
[experiments](THE_BOND_AND_THE_SEAT_AT_FIXED_GAMMA.md) supplied the fixed-γ
producer, the [Route-B histogram study](ROUTE_B_N4_HISTOGRAM_FILTER.md) and
[ReadoutFisher](../compute/RCPsiSquared.Diagnostics/Foundation/ReadoutFisher.cs)
supplied statistical precedents for other readouts, while the inspected
[hardware flight](IBM_RUN3_PALINDROME.md) and
[Confirmations registry](../simulations/framework/confirmations.py) supplied
no N=7 bond/seat finite-shot measurement. The [glossary](../docs/GLOSSARY.md)
fixes the hopping convention and the F157 seat term; the
[caught-errors ledger](../docs/CAUGHT_ERRORS.md) warns about the XY factor
of two and unequal total dephasing in profile comparisons. The
[OpenArcs registry](../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs)
holds broader inverse questions, with no matching finite-shot conclusion.

In this note, *seat* means the **population readout site**, following the
[coherent bond/seat calculation](THE_BOND_AND_THE_SEAT.md). F157 and the
glossary use *seat* for the **site carrying dephasing**. The readout site 2
and centre dephasing site 3 have different jobs here.

## The experiment that was simulated

The [fixed-positive-gamma response study](THE_BOND_AND_THE_SEAT_AT_FIXED_GAMMA.md)
found six numerical bond-response directions at each off-centre population
readout in its declared continuous-time samples, but left measurement counts
open. This continuation asks a smaller, operational question: if the chain
has no defect or exactly one bond has an additive hopping change of **known
magnitude** `0.0075` (`10% J_hop`), can ideal Z-position counts identify its
bond and unknown sign from the two declared candidate signs?

The prepared state is the nominal highest-energy standing wave `psi1` of the
unperturbed seven-site open chain, held fixed under every defect. In the
single-excitation block,

```text
H_SE = J_hop A_path + sum_b deltaJ_b V_b,
V_b = |b><b+1| + |b+1><b|,
D(rho)_ab = -2(gamma_a+gamma_b) rho_ab  (a != b),
D(rho)_aa = 0.
```

Here `J_hop` is both the matrix hopping and the coefficient of `(XX+YY)/2`
in the full XY book, as in the [coupling convention](../docs/GLOSSARY.md).
The two channel profiles are separately `gamma_a=0.05` at all seven sites
(`uniform`, `sum gamma_a=0.35`) and at site 3 alone (`single_centre`,
`sum gamma_a=0.05`). The per-site `gamma0` stays fixed, while both the total
dephasing budget and its site assignment change. There is no T1,
single-excitation leakage, preparation
error, unknown rate, unknown defect magnitude, or detector mislabeling in
the inference model.

Each time point uses 100,000 independent ideal Z-basis shots, each yielding
exactly one of the seven excitation positions. Ten distinct times cost one
million shots **per experiment**. The site-2 comparator retains only the
Bernoulli count of site-2 outcomes from the **same** seven-category records;
site 3 is a separate reflection control, not the classifier used in the table.
There are 256 independently drawn experiments for each of thirteen fixed
truths, in the order `null,b0-,b0+, ...,b5-,b5+`. Both signs of each one-bond
defect have magnitude `0.0075`. The producer uses NumPy `SeedSequence` with
the root `20260923` and separate profile, schedule, and truth indices;
the retained run records the exact tuple and NumPy/SciPy versions.

For each profile, a deterministic one-for-one exchange starts from
`[2,4,6,8,10,12,15,18,24,30]` and searches the 37 times in the fourfold
subdivided candidate pool. Its objective is the smallest eigenvalue of the
**nominal six-parameter full-position Fisher matrix**. The resulting
schedules are local exchange optima; they are fixed before propagating any
finite defect or drawing counts. Both the original and selected schedules
use one million shots. Site 2 receives the full-position-selected schedule
as a same-schedule comparator, not its own time optimum.

The thirteen candidate probability tables are propagated with the *finite*
`L + deltaJ_b B_b`, not a first-order tangent approximation. For each drawn
count table, the categorical classifier maximizes
`sum_tj n_tj log p_h,tj`; the site-2 classifier maximizes
`sum_t [k_t log p_h,t2 + (100000-k_t) log(1-p_h,t2)]`. A zero-count term is
zero, a positive count on a zero-probability model outcome has score
negative infinity, and exact ties take the first declared hypothesis.

## Observed counts and local information

The following totals are descriptive sums over thirteen **fixed strata**,
each with 256 repeats. They are not pooled-binomial confidence intervals.
The retained run reports every true-versus-predicted confusion row and a
separate two-sided Wilson 95% interval for each stratum and detector.

| Channel profile | Time schedule | Smallest full Fisher eigenvalue | Smallest site-2 Fisher eigenvalue | Full-position correct | Site-2 correct | Full-only / site-2-only correct on paired records |
|---|---|---:|---:|---:|---:|---:|
| uniform | baseline | 190,443 | 0.00160470 | 3328/3328 | 2286/3328 | 1042 / 0 |
| uniform | selected | 344,084 | 0.0000982058 | 3328/3328 | 2476/3328 | 852 / 0 |
| single centre | baseline | 1,144,230 | 0.00429911 | 3328/3328 | 3208/3328 | 120 / 0 |
| single centre | selected | 3,198,093 | 0.00184081 | 3328/3328 | 3137/3328 | 191 / 0 |

The Fisher columns are in inverse squared **additive hopping-rate units**:
each response is differentiated with respect to `deltaJ_b`, so rescaling
the rate unit rescales these eigenvalues. The comparisons hold the same
`gamma0`, `J_hop`, and physical-time convention throughout.

The selected uniform times are
`[8.5,9,9.5,10,10.5,24,25.5,27,28.5,30]`; the selected centre-lit times
are `[10.5,11,11.5,12,19.5,21,22.5,24,25.5,27]`. Every full-position
truth row happened to score 256/256 in these seeded simulations. A 256/256
row has a Wilson 95% lower endpoint of `0.985216`; the run does **not**
establish error-free performance under repeated draws, shifted parameters,
or a real detector. In the uniform baseline run the site-2 null row scored
only 36/256, while its centre-lit counterpart scored 233/256. These are
results for two different specified channel profiles. Because both the total
dephasing budget `sum gamma_a` and placement change, their difference does
not isolate the effect of placement at a matched budget.

The full-position Fisher matrix dominates the site-2 matrix in the
positive-semidefinite sense at the four scheduled fixtures, with computed
minimum eigenvalues of their differences well above the printed roundoff
budgets. The information retained by the full record has a concrete location:
the site-2 bit merges all six outcomes `j != 2` into one "elsewhere" bin.
For two candidate models, categorical KL separation equals site-2 Bernoulli
KL separation plus the probability of "elsewhere" times the KL separation
of the **conditional positions within that bin**. This extra term is
nonnegative and is zero when those conditional position distributions agree.
It is the same bin-merging mechanism tested for a different readout in the
[Route-B histogram study](ROUTE_B_N4_HISTOGRAM_FILTER.md), not an F124 or
F157 consequence. It explains information dominance, not the exact success
count of a finite draw.

The centre-site Fisher rank is three: reflection fixes the centre
population while pairing bonds `b` and `5-b`. Site 2 has no such exact
reflection cap. Yet the site-2 *smallest* Fisher eigenvalue fell under both
full-position-optimized schedules, even as its uniform-profile thirteen-way
score rose; its centre-lit score fell. This is no contradiction: Fisher
describes infinitesimal precision across **six independent continuous bond
parameters**, whereas the classifier chooses among thirteen discrete,
known-magnitude finite alternatives. Neither is a hardware recovery result.

## What the checks establish

The [test controls](../simulations/test_handshake_bond_seat_finite_shots.py)
compare representative finite-defect populations with a separate
matrix-shaped ODE using explicit site-Z operators. That control has no
Kronecker Liouvillian or matrix exponential derivative. Halving only the
producer's channel rate breaks the physical comparison. Hand-computed
three-outcome cases check categorical versus Bernoulli Fisher and
likelihoods, including zero bins and ties. A detector-only outcome
reflection flips the diagnosed end-bond label in a finite-defect
likelihood test, while a physical bond reflection preserves centre-site
probabilities. The producer checks probability normalization, each bond's
zero-sum population slopes, exact per-time shot conservation, all confusion
row totals, Fisher dominance within a stated roundoff budget, and the
centre reflection/rank controls. **No favorable accuracy threshold is a
pass condition.**

This finite-shot run is a local model probe, not a new F theorem or a
hardware confirmation.
