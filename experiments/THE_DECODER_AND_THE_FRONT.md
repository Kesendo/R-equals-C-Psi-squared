# The Decoder and the Front: Which Carrier Confuses

<!-- Keywords: DefectDecoder, walk-time step, arrival-time dictionary, sign-location confusability,
anti-collinear letters, conditioning, carrier choice, bonding mode, site seed, mode participation,
multi-defect decoding, Fisher information, Cramér–Rao bound, shots, equal cost, uniform rescaling, sum rule,
maximum likelihood, de-lossed path -->

**Status:** computed, 2026-10-30. The front's letters, their decoding, the trajectory column of the carrier
table with its conditioning, and in §5 the bond/seat handshake and the window table at Q = 20 are C# readings in `compute/MirrorWorld/WalkTime.cs` (run mode `walk decode`,
pinned in `WalkTimeTests`), and an exact propagation of the same block in
[`simulations/decoder_carrier_dictionary.py`](../simulations/decoder_carrier_dictionary.py) (output
[`decoder_carrier_dictionary.txt`](../simulations/results/decoder_carrier/decoder_carrier_dictionary.txt)) meets
them to three digits (the residual ratios of §1 excepted, which carry the integration step). The linear-letter
column, the mirror-law check, the stationary carriers, the runs past
N = 10, the random carriers of §4, and the rest of §5 (the other handshakes and schedules, the end mixture, the
uniform-rescaling split, the per-bond bounds and the decoder's efficiency) and the ladder of §6 are on the Python
route alone, beside a C# test of §6's exact null. The sum of the letters in §4, the decode at finite shots and the
de-lossed path in §5 are read in [`decoder_letters_sum_and_shots.py`](../simulations/decoder_letters_sum_and_shots.py)
(output [`decoder_letters_sum_and_shots.txt`](../simulations/results/decoder_carrier/decoder_letters_sum_and_shots.txt)),
the sum's two consequences for the bonding mode also in C# tests. Why a carrier confuses is read in part: the
letters of every bond add up to the letter of the common speed-up, which without dephasing a carrier at rest does
not have, so its letters cancel and point against each other on average (§4). Why a carrier on few modes cancels,
and how far its worst pair falls below the mean, is open.

<!-- STAGE 0 SWEEP RECORD (for us, not for the reader: hidden on GitHub, read by every session, agent and hook)

Stores searched 2026-10-30 for carrier choice in the defect reading, a carrier without mirror parity, a bond
defect located by a front, and decoding several defects. docs/ANALYTICAL_FORMULAS.md: F124 holds the band-edge
carrier only ("other carriers give a sum below 2"); its resolution-limit reading treats recovering every bond
from the location dictionary M as an inverse problem, M with the strength column of full rank N minus 1,
smallest singular value the square root of E, typed as BandEdgeResolutionLimitClaim; it names the decoder's 1.5
an alpha-rescaling artifact; F2b and F160 say nothing on carriers or on locating a defect, and F126 (the renewal
representation, PROOF_DEPHASING_FRONT_RENEWAL) is the exact solution of every population this page reads: stated
for a seed site, its derivation carries over to any carrier and any couplings, defects included; joined in
section 6. docs/proofs/:
PROOF_HANDSHAKE_TRANSITION_INVARIANT lines 45 and 58 (interior carriers psi_2, psi_3, definite parity),
PROOF_PTF_CHIRAL_MIRROR_RATE_LAW line 85 (the chiral identity holds with several defects; not a decoding);
nothing on a carrier without parity. experiments/: COUPLING_DEFECT_WALK_TIME_STEP P3' locates the defect from
the front ("the step edge sits at the bond") and contrasts the front with the PTF alpha_i, and never names the
decoder; THE_NUCLEUS_AS_A_KEY reads three couplings of the sodium nucleus back from relaxation rates exactly, a
different inverse problem; nothing decodes several bond defects of one chain. hypotheses/: HANDSHAKE_GEOMETRY
lines 194 to 202 vary the OBSERVABLE with psi_1 fixed (the phase-blindness gate: linear Z and the full Z(t) read
are as anti-collinear as the purity), line 279 "the carrier-parity is load-bearing" for the equivariance sign;
the M4 parity check asks for a mirror-symmetric reading; PERSPECTIVAL_TIME_FIELD tried psi_1 to psi_4, the
product state and the centre seed, and its "Where it degrades" said the centre seed overlaps all seven modes (it
overlaps the four odd ones; corrected there). OpenArcs registry: handshake_decoder (retired; its probe
handshake_carrier_compare.py compared psi_1 with the PairState, f-profile correlation 0.999, both mirror-even)
and f124_inverse_problem_resolution_seam (retired; the de-loss to the signed deviation, on psi_1),
ptf_bonding_class_guard (retired: the Bell+ refusal was a minimizer trap); no arc on carrier choice.
docs/CAUGHT_ERRORS.md lines 586 to 595: the de-loss over-correction ("escapes the anti-collinearity").
fw.Confirmations, docs/GLOSSARY.md: nothing. Typed layer, Core and Diagnostics: DecoderAntiCollinearitySeesawClaim
(cos = 1 minus 2w on mirror pairs, premise the carrier's mirror symmetry), DefectReadingEquivarianceClaim (the
sign for the R-even carrier psi_1), KPartnerSelectionRuleClaim (the bare location couplings of psi_1, rank N
minus 2), BandEdgeResolutionLimitClaim, DefectDecoder (BondingMode hard-coded; DecodeDeviation reads the N = 5
confuser at a squared residual ratio near 516), SingleExcitationPainterDictionary (psi_1 hard-coded), Symphony
(offers a site-0 carrier and Bell+, never used for a dictionary; its comment on localized states was stale for
Bell+). MirrorWorld: WalkTime reads one defect forward, no inverse; Crack and Warble nothing here.
Adjacency checked: WalkTime / COUPLING_DEFECT_WALK_TIME_STEP P3' and DefectDecoder / HANDSHAKE_GEOMETRY locate
one bond defect on one chain and did not cite each other; PERSPECTIVAL_TIME_FIELD's degradation of the rescaling
picture for a localized carrier and this page's collinear linear letters are one finding seen from two sides;
HANDSHAKE's observable sweep and this page's carrier sweep are the two halves of one question; F124's
resolution-limit inverse problem and this page's every-bond decoding are one problem on two carriers. This page
joins them.
For §5 (equal cost), swept again: the typed Fisher apparatus ReadoutFisher and ReadingPowerWitness (FI of a bond
defect against Q, N = 4, psi_1); HANDSHAKE_GEOMETRY's "The resolution law (measured)" (FI about linear in Q) and its
"Identifiability is not FI(Q)" paragraph (rank against precision, the distinction §5 turns into numbers); the Route-B
readout pages ROUTE_B_N4_READOUT_SEARCH, ROUTE_B_N4_HISTOGRAM_FILTER, ROUTE_B_N4_VIRTUAL_READOUT (readout Fisher, no
carrier comparison) and COHERENCE_HORIZON_EP_SENSOR_DEBATE (metrology, parked; nothing on carriers); the arc handshake_decoder's M2c
read-cost test (simulations/handshake_read_cost.py: the dose to the best read falls with Q, not as 2/Q); and the
bond/seat family, experiments/THE_BOND_AND_THE_SEAT.md, THE_BOND_AND_THE_SEAT_AT_FIXED_GAMMA.md and
THE_BOND_AND_THE_SEAT_WITH_FINITE_SHOTS.md (psi_1 at Q = 1.5, the six-bond full-position Fisher matrix, a
thirteen-way classifier at 1e6 shots). None compares carriers at equal cost; all hold psi_1. Adjacency: the
bond/seat Fisher matrix and this page's are one object (reproduced exactly), and that family and the decoder
pages never cited each other; joined here and on the finite-shots page. The bond/seat family's "a uniform change
of all bonds leaves the psi_1 population stationary" is psi_1's worst Fisher direction here, read from the other
side.
For the open items (the sum of the letters, a decode at finite shots, the de-lossed path with other carriers), swept
again: a search for sum rules over docs/, experiments/, hypotheses/, the typed layer and the MirrorWorld base returns
the palindromic sum rule alpha_fast + alpha_slow = 2 Sigma gamma and F49's bond-sum rule, of other objects; PTF's
closure law sums ln alpha over sites, not over bonds; among the hits are also F73's generalisation, the sodium
two-mirror witness, THE_SEAT_THAT_CUTS and F89's pair sum, all of other objects (the list is what was read, not a
closed census). The generator's joint homogeneity of degree one in (H, gamma) is stated in INCOMPLETENESS_PROOF,
Q_BELONGS_TO_NO_SUBSTANCE and CAUGHT_ERRORS; it gives the sum's exact form with dephasing, used in section 4. The
uniform rescaling as psi_1's null: THE_BOND_AND_THE_SEAT ("the ceiling N minus 2 is every bond combination orthogonal
to the uniform change of all bonds"); THE_BOND_AND_THE_SEAT_AT_FIXED_GAMMA (no coherent readout sees the uniform
change, and with gamma fixed the sum of psi_1's six bond slopes at one site is the response to the changed
J/gamma_0). The K-partner column of the location dictionary: KPartnerSelectionRuleClaim and
HANDSHAKE_GEOMETRY ("The consequence is the decoder's rank deficit", naming the confusability the K-partner
near-degeneracy), neither naming the bond vector of the deficit. PROOF_HANDSHAKE_TRANSITION_INVARIANT's two-spaces
caveat: on the full M the floor's bond vector is the staggered (-1)^b, on k = 2..N the K-partner at 0, its vector
unnamed (it is (1, .., 1)). Adjacency: one rank deficit of one square matrix, its channel side the K-partner column
and its bond side the uniform rescaling; HANDSHAKE_GEOMETRY holds the first, THE_BOND_AND_THE_SEAT both without
joining them, and none reads the sum as the source of the opposite signs; joined here and on both. The join exposed
HANDSHAKE_GEOMETRY's sector sentence, which put the K-partner null in the mirror-odd sector by a labelling rule of
simulations/handshake_gram_metric.py; the painted Gram's small mirror-odd eigenvector is the staggered wave, the
mirror-even one the uniform rescaling (repaired there, CAUGHT_ERRORS 2026-10-30). Nothing in the stores reads the sum of the
bond letters for other carriers. Finite shots: THE_BOND_AND_THE_SEAT_WITH_FINITE_SHOTS classifies on psi_1 at 1e6
shots, BandEdgeResolutionLimitClaim reads a conditioning; nothing draws shots for another carrier. The de-lossed path:
DefectDecoder.DecodeDeviation, the bonding mode hard-coded, N up to 5, and HANDSHAKE_GEOMETRY's de-loss update.
-->

The [defect decoder](../hypotheses/HANDSHAKE_GEOMETRY.md) reads where a weakened or strengthened bond sits on an
XY chain under Z dephasing. It starts one excitation in the chain's lowest mode, the bonding mode ψ₁ (the
**carrier**: the single-excitation state the reading is run on), and reads the per-site purity. Its known weakness
is that some of its letters are nearly opposite: at N = 5, weakening bond 3 reads almost like strengthening
bond 1. The [walk-time step](COUPLING_DEFECT_WALK_TIME_STEP.md) reads the same kind of defect with another carrier,
one excitation released at site 0, whose front's arrival is shifted by about −δ/(2J) beyond the bond and not at
all before it. This page reads the front backwards, as a dictionary, and then asks which of the two things that differ between
the two readings sets the decoder's confusion: the carrier or the observable.

Everything here is the single-excitation block of the open chain, coupling J = 1, every coherence charged 4γ, at
γ = 0.05 (Q = J/γ = 20, the decoder's protocol point). A **letter** of bond b is the derivative of a reading with
respect to that bond's δ (J_b = J(1 + δ)): one entry per site, or per site and sample time. Two numbers describe a
dictionary of letters. The **worst signed cosine** is the most negative cosine between two letters; near −1, a
defect on one bond reads like the opposite defect on the other, the decoder's **sign-location confusion**, and
below −0.5 this page calls two letters **nearly opposite**. The **largest cosine** is the most positive one; near
+1, two bonds read alike whatever the sign, a **location blindness** the worst signed cosine does not see. The
**conditioning** is the smallest singular value of the dictionary with every letter normalised to length 1: near 0,
some combination of defects is invisible to the reading (either kind of near-coincidence makes it small); 1 for
orthogonal letters.

## 1. The front, read backwards

The arrival letter of bond b is d(arrival at site i)/dδ_b, by the central difference at δ = ±0.02, the arrival read
as on the walk-time page: the first crossing of 0.2 of the site's clean peak. At N = 5:

| bond | site 1 | site 2 | site 3 | site 4 |
|---|---|---|---|---|
| 0 | −0.347 | −0.458 | −0.545 | −0.528 |
| 1 | +0.013 | −0.460 | −0.553 | −0.547 |
| 2 | 0.000 | +0.062 | −0.554 | −0.547 |
| 3 | 0.000 | −0.001 | +0.163 | −0.528 |

Each letter is a step: near −1/(2J) beyond the bond, zero before it, and a positive bump on the site just before the
bond, larger the farther the bond sits from site 0 (+0.013 at bond 1, +0.163 at bond 3). Nested steps are far from
opposite:

| N | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|
| worst signed cosine | +0.825 | +0.557 | +0.363 | +0.214 | +0.106 | +0.027 | −0.079 | −0.167 |
| conditioning | 0.419 | 0.318 | 0.259 | 0.225 | 0.201 | 0.183 | 0.169 | 0.158 |

against the decoder's painted letters (its per-site α letters, read off propagated runs on ψ₁) at −0.91 to
−0.99 at N = 3 and 5 to 12
([`DecoderAntiCollinearitySeesawClaim`](../compute/RCPsiSquared.Diagnostics/Ptf/DecoderAntiCollinearitySeesawClaim.cs)).
Past N = 10 the bump keeps growing and both numbers keep falling: the worst signed cosine reads −0.300, −0.458 and
−0.544 at N = 12, 16 and 20, the conditioning 0.141, 0.119 and 0.104. The cosine also depends on the threshold
(at N = 12: −0.032 at 0.1 of the peak, −0.300 at 0.2, +0.103 at 0.5). So the front's letters do not stay clear of
the confusion as the chain grows; through N = 10 at this threshold they do.

The dictionary has full rank N − 1 at every N here, and it reads **every bond at once**: least squares over all
letters recovers a defect on every bond together. At N = 5 and N = 8, with δ drawn between −0.02 and +0.02 on each
bond, the largest error is 1.4 · 10⁻⁴ to 6.6 · 10⁻⁴. That error is the second order of the step law and behaves as
one: halving every δ multiplies it by 0.23 to 0.25. [F124](../docs/ANALYTICAL_FORMULAS.md)'s resolution-limit reading
poses the same inverse problem on the bonding mode's location dictionary and reads its conditioning (full rank N − 1
with the strength column, `BandEdgeResolutionLimitClaim`); here the decode is run.

The decoder's own confuser (bond 3 at −0.02 against bond 1 at +0.02, N = 5) the front reads cleanly with the
decoder's single-letter rule: the right bond and sign, squared residual ratios of order 10⁴ and 10⁵ (C# 8093 and
352842, Python 7706 and 315454; the winner's residual is second order in δ and carries the integration step, so the
two routes part in the second digit). This is not the front's advantage. Without shot noise the decoder's de-lossed path,
which reads the signed purity deviation instead of the α stretch, reads the same confuser on ψ₁ as well
(`DefectDecoder.DecodeDeviation`, squared residual ratio near 516, against 1.45 on the α path in the live
`inspect --root decoder`), with its letters still nearly opposite. What nearly opposite letters cost is precision
under shot noise (the finite count of measurements, not the watching γ); for the every-bond decode one measure of it is the conditioning.

## 2. Carrier or observable

The two readings differ in two things at once. The house already varied one of them: with ψ₁ fixed, the linear
⟨Z_i⟩ read and the full ⟨Z_i⟩(t) signal are as anti-collinear as the purity
([HANDSHAKE_GEOMETRY](../hypotheses/HANDSHAKE_GEOMETRY.md), the phase-blindness gate). Here the other is varied.
Two observables are held fixed: the whole population trajectory (every site at t = 0.25, 0.5, …, 20; central
difference at δ = ±0.01 in C#, ±10⁻⁴ in Python) and the decoder's own **linear letter**, one number per site, the
δ → 0 limit of the α fit by which the decoder's painting stage reads each site's purity curve as a stretched copy of
the clean one. The carrier changes:

| carrier | trajectory: worst cos, conditioning | linear letter: worst cos, conditioning |
|---|---|---|
| bonding mode ψ₁ | −0.88 … −0.97, 0.07 … 0.11 | −0.93 … −0.98, 0.0004 … 0.14 |
| mode ψ₂ | −0.90 … −0.99, 0.06 … 0.10 | −0.96 … −0.99, 0.0002 … 0.11 |
| (ψ₁ + ψ₂)/√2 | −0.75 … −0.80, 0.14 … 0.22 | −0.87 … −0.99, 0.005 … 0.03 |
| end pair (\|0⟩ + \|N−1⟩)/√2 | −0.45 … +0.68, 0.45 … 0.63 | −0.08 … +1.00, 0.006 … 0.04 |
| end mixture ½(\|0⟩⟨0\| + \|N−1⟩⟨N−1\|) | −0.14 … +0.19, 0.36 … 0.63 | +0.47 … +0.66, 0.004 … 0.05 |
| centre site (odd N) | −0.26 … +0.40, 0.34 … 0.47 | +0.39 … +1.00, 0.006 … 0.013 |
| site 0 | −0.16 … +0.06, 0.42 … 0.65 | +0.47 … +0.73, 0.020 … 0.040 |

(ranges over N = 5 to 9; every entry per N is in the output file. Below, **localized** means a carrier on one
or two sites: the site seed, the end pair, the end mixture, the centre site.) At N = 12 and 14 the bonding mode's trajectory
reads −0.82 and −0.76 at conditioning 0.07, the site seed's −0.11 and −0.19 at 0.73 and 0.71.

On the trajectory the carrier decides. The modes and the mode superposition read nearly opposite letters and a
dictionary near degeneracy (conditioning 0.06 to 0.22). Every localized carrier reads no nearly opposite pair and a
far better conditioned dictionary (0.34 to 0.65); the site seed's is about four to ten times better conditioned than
the bonding mode's (N = 5 to 9). The house's sweep, with the carrier fixed at ψ₁, found no observable that helps;
this one finds that, for the trajectory, a localized carrier does.

On the linear letter the localized carriers lose the opposite pairs too, but their letters become nearly the same
instead (largest cosine 0.96 to 1.00 at N = 5 to 9): the sign confusion is exchanged for a location blindness, and the
dictionary is conditioned just as poorly as ψ₁'s (0.004 to 0.05, against ψ₁'s 0.0004 to 0.14), better at some N and
worse at others. For these carriers the linear letter is a projection that exists, not a fit the decoder's painting
stage would accept: [PERSPECTIVAL_TIME_FIELD](../hypotheses/PERSPECTIVAL_TIME_FIELD.md) records ("Where it degrades")
that for a localized carrier the rescaling picture behind the α fit does not describe the data, and the near-identical
letters are the same failure read from the dictionary. One stretch per site discards what a localized carrier's
response carries. So for a carrier that holds few modes no observable tried avoids the confusion, and for a localized
one the observable matters as much as the carrier.

At even N the end pair and the end mixture read the same numbers. That is exact for real nearest-neighbour hopping,
with or without bond defects, and no on-site field: the end sites lie on opposite sublattices, the sublattice sign
together with complex conjugation commutes with the generator (it turns H into −H twice and leaves the dephasing
alone) and maps |0⟩⟨N−1| to its negative, so that coherence's contribution to every population has zero real part.

## 3. What does not separate them: mirror parity, stationarity, localization

The decoder's letters obey an exact mirror law under the chain reflection R (site i ↦ N − 1 − i; the spatial
mirror of `DefectReadingEquivarianceClaim`, not the coherence-space mirror of MirrorWorld's `Mirror` or Π). For a
carrier state ρ that R leaves fixed, the letter of bond N − 2 − b is the reflection of the letter of bond b, so the cosine of a mirror pair
is 1 − 2w(b), with w the share of the letter's squared norm in its mirror-odd part
(`DecoderAntiCollinearitySeesawClaim`). On the linear letter it holds here for every such carrier of the table,
the modes (whose states ρ are R-invariant whatever the parity of the wave) and the end pair, end mixture and centre
site, to 4 · 10⁻¹¹ (the finite difference's rounding).

A natural guess is that the parity separates the two groups: the bonding mode keeps the mirror symmetry, the site
seed breaks it. The table rules this out. The end pair, the end mixture and the centre site are R-invariant and read
well; (ψ₁ + ψ₂)/√2 is not R-invariant and reads nearly opposite letters at every N. The mirror law fixes the cosine
of a mirror pair as 1 − 2w and leaves its sign to w; R-invariance is the premise of the law, not the source of the
confusion.

Two further guesses fail as well. Every carrier tested that is stationary under the clean H (ψ₁, ψ₃, ½(P₁ + P₃),
½(P₁ + P₂) at N = 5, 7, 8; P_k the projector on mode k) reads −0.75 to −0.996 on both observables, but a stationary
single-excitation state is a mixture of modes and so R-invariant automatically, which ties this test to the last.
And some moving carriers read nearly opposite letters too, (ψ₁ + ψ₃)/√2 at N = 7 (−0.69 on the trajectory, −0.99 on
the linear letter), while the same carrier at N = 5 reads +0.54 and +0.37. Localization alone does not do it either:
a smooth packet one site wide reads nearly opposite letters in every draw (§4).

## 4. A tendency, not a law: how many modes the carrier holds

The two groups of §2 differ in how many of the chain's modes the carrier holds, counted by the **mode
participation** PR = 1/Σ_k |⟨ψ_k|φ⟩|⁴ (1 for a mode, N for equal weight on all). A mode holds one. A site holds many:
site 0 all N, unevenly (PR 16/3 ≈ 5.3 at N = 7), the centre site the (N + 1)/2 modes ψ₁, ψ₃, … that do not vanish
there, equally (PR 4 at N = 7). Read against PR, over 132 random carriers at N = 7, twelve per group: carriers with Gaussian random weights on m randomly
chosen modes, and Gaussian packets exp(−(i − x₀)²/(4s²)) of width s sites at a random centre x₀:

| carriers | PR (mean) | trajectory: mean, range, below −0.5 | linear letter: mean, range, below −0.5 |
|---|---|---|---|
| on 1 mode | 1.00 | −0.935, −0.990 … −0.882, 12 | −0.951, −0.978 … −0.928, 12 |
| on 2 modes | 1.31 | −0.710, −0.996 … −0.126, 9 | −0.560, −0.999 … +0.156, 7 |
| on 3 modes | 1.42 | −0.742, −0.986 … −0.083, 10 | −0.495, −0.944 … +0.717, 8 |
| on 4 modes | 2.30 | −0.340, −0.735 … +0.277, 5 | +0.125, −0.912 … +0.913, 2 |
| on 5 modes | 2.53 | −0.304, −0.969 … +0.032, 3 | +0.417, −0.985 … +0.907, 1 |
| on 6 modes | 2.96 | −0.133, −0.577 … +0.344, 2 | +0.502, −0.542 … +0.926, 1 |
| on 7 modes | 3.51 | −0.010, −0.483 … +0.240, 0 | +0.645, −0.318 … +0.859, 0 |
| packet, s = 0.5 | 3.18 | −0.252, −0.428 … −0.082, 0 | +0.450, −0.138 … +0.813, 0 |
| packet, s = 1 | 1.81 | −0.713, −0.787 … −0.594, 12 | −0.898, −0.967 … −0.732, 12 |
| packet, s = 2 | 1.22 | −0.746, −0.835 … −0.605, 12 | −0.657, −0.879 … −0.450, 8 |
| packet, s = 4 | 1.17 | −0.643, −0.649 … −0.635, 12 | −0.416, −0.445 … −0.322, 0 |

(worst signed cosine; the last number of each cell counts the draws below −0.5 out of 12.) The correlation of PR with
the worst signed cosine is +0.81 on the trajectory and +0.68 on the linear letter: few modes, nearly opposite
letters; many, not. The trajectory reads the tendency more cleanly. It is a tendency with counterexamples, a
five-mode carrier at −0.969 on the trajectory, and twelve draws per group make it a sighting, not a law.

**What the tendency reads: the letters add up to the common speed-up.** The hopping matrix is the sum of its
bonds, h = Σ_b V_b, so to first order the letters of all bonds add up to the letter U of the uniform rescaling of
every bond (the finite differences meet this to 10⁻⁸ to 10⁻⁶ of the letters' summed norms, their own second-order
error; the C# test gates its δ² law). Without dephasing that rescaling is a rescaling of time, and U = t·dp/dt: the
carrier's own motion. With dephasing the generator is homogeneous of degree one in the couplings and γ together,
so U = t·dp/dt − γ·∂p/∂γ (met to 10⁻⁷ to 10⁻⁶ here). Without dephasing a carrier at rest under H therefore has
U = 0, and its letters add up to nothing. That puts their norm-weighted mean cosine,
(κ²S² − Σ_b|ℓ_b|²)/(S² − Σ_b|ℓ_b|²) with S = Σ_b|ℓ_b| and the **cancellation** κ = |U|/S, at its lowest, κ = 0.
At N = 3 nothing else is free: without dephasing every mode's two letters are exactly opposite (cos −1 to twelve
digits; −0.963 at Q = 20). At N = 4 three letters that close a triangle have cosines fixed by their norms alone,
cos_ab = (|ℓ_c|² − |ℓ_a|² − |ℓ_b|²)/(2|ℓ_a||ℓ_b|), and without dephasing every mode meets that to seven digits (its
mirror pair at −0.920; at Q = 20 the halo moves the pair to −0.829 against the law's −0.904). From N = 5 on the sum,
given the letters' lengths, fixes only their mean cosine, and the bonding mode's norm-weighted mean cosine lies
within 0.016 of −1/(N − 2), the mean of N − 1 letters of one length that add up to nothing (unequal lengths can
only lower it), at N = 5, 7, 9, with or without dephasing.

For the bonding mode the same sum is the weakest direction. Along the norm vector (|ℓ_0|, …, |ℓ_{N−2}|) the
normalised letters add up to U, so without dephasing ψ₁'s dictionary is singular there (conditioning 0, its null
along the norm vector to four digits at N = 5, 7, 9). At Q = 20 the halo lifts it to the conditioning 0.07 to 0.11
of §2, its weakest direction still within 0.997 to 0.999 of the norm vector, while the site seed's weakest
direction lies elsewhere (overlap 0.02 to 0.05, and 0.008 to 0.016 without dephasing). The near-degeneracy of ψ₁'s
row in §2's table, its worst Fisher direction in §5 and the null of its never-caught term in §6 are this one sum in
three coordinates: exact in §6, to 0.997 or better in §2 (against the norm vector) and §5 (against (1, …, 1)). In the bare location dictionary
M[b, k] = ⟨ψ_k|V_b|ψ₁⟩, k = 2..N, it is exact and Q-free: every column sums over the bonds to ⟨ψ_k|h|ψ₁⟩ = 0, so
the one left null vector of that square matrix of rank N − 2 is (1, …, 1) (to 10⁻¹⁵, N = 3..10). The K-partner rule
(`KPartnerSelectionRuleClaim`) is the same deficit seen from the channels: it empties the column k = N, and the sum
names the bond combination left unseen. A rank count alone would allow a dependence among the letters with
coefficients of both signs; that the coefficients are all equal is what makes the letters point against each other
on average.

Dephasing breaks the cancellation only through the halo (§6), the scale of the bonds entering there through J/γ
([THE_BOND_AND_THE_SEAT_AT_FIXED_GAMMA](THE_BOND_AND_THE_SEAT_AT_FIXED_GAMMA.md) sums the bonding mode's six bond
slopes at one site into the response to that changed J/γ₀): at Q = 20 the stationary carriers' κ is 0.03 to 0.06 at N = 5 and 7, against 0.63 to 0.71 for
the site seed (κ = 1 when all letters are parallel and of one sign). The slow beat (ψ₁ + ψ₂)/√2 reads 0.41 and 0.21
at N = 5 and 7, the centre site 0.90 and 0.68. Over the 132 random carriers of the table above, κ tracks the worst
signed cosine more closely than the mode participation: correlation +0.90 against PR's +0.81 at Q = 20, +0.86
against +0.68 without dephasing; κ and PR themselves correlate +0.74 (+0.60), and one set of 132 draws makes
these margins a sighting. On the trajectory at Q = 20 it is a tendency, not a law: worst cosines below −0.5 occur up
to κ = 0.58, and above −0.5 from κ = 0.29. On the linear letter κ tracks both failures, the worst cosine (+0.79
against PR's +0.68) and the largest (+0.65 against +0.40). There the localized carriers' letters barely cancel (κ
0.85 to 1.00 at N = 5 and 7, against 0.04 to 0.16 for the modes and the beat), and letters whose sum nearly reaches
the sum of their lengths are nearly parallel (largest cosine 0.98 to 1.00); why the localized carriers' linear
letters do not cancel is not read. A second guess, that one dominant direction concentrates the letters whatever
their signs and the sum then picks the sign, finds no consistent support: the rank-one share of the normalised
dictionary, which cannot see sign, correlates −0.43 with the worst cosine at Q = 20, the sign the guess predicts,
and +0.31 without dephasing.

What the sum does not fix is how far the worst pair falls below the mean: ψ₁'s worst pair lies at −0.88 to −0.97 on
the trajectory at Q = 20 while its mean lies near −1/(N − 2). The house splits the letter into carrier
channels u_k, whose couplings are closed-form and whose sizes have no closed form (the `handshake_Rk` facet, typed
IsDeadEnd: computable, no closed form), and w, which sets the sign of a mirror pair's cosine, is the balance of
those sizes; whether the worst pair's depth can be read from them is open. Which pair is the worst is not read here,
nor why a carrier on few modes should cancel: κ correlates with PR (+0.74) without following from it, and the beat
(ψ₁ + ψ₂)/√2 confuses at κ = 0.21 to 0.41.

## 5. At equal cost: the Fisher information of the shots

The comparisons above set dictionaries side by side. What a reading costs is counted in **shots**: one shot
prepares the carrier, lets it evolve for a time t, and measures every site in Z; since hopping and Z dephasing keep
one excitation, the shot returns where it is, site j with probability p_j(t). Two carriers are compared at equal
cost when they receive the same shots on the same schedule of times (equal shots per time). The per-shot Fisher
matrix of every bond's δ is F_ab = mean_t Σ_j (∂p_j/∂δ_a)(∂p_j/∂δ_b)/p_j, and its smallest eigenvalue λ_min is the
information in the worst direction when every bond is unknown: after n shots, no unbiased estimate of that
combination of bonds has a variance below 1/(n·λ_min) (the Cramér–Rao bound). The readout and its Fisher matrix
are the house's (`ReadoutFisher`, and the full-position Fisher of
[THE_BOND_AND_THE_SEAT_WITH_FINITE_SHOTS](THE_BOND_AND_THE_SEAT_WITH_FINITE_SHOTS.md)); both are reproduced
exactly: `ReadingPowerWitness`'s 1.5934 at N = 4 (its forward difference at δ = 0.02, so this confirms the
convention), and that page's 190,443 at N = 7 (its baseline schedule and uniform profile, summed over its 10⁶ shots,
in its additive units) and 344,084 (its selected schedule). The cone's RK4 meets 190,443.2487 to its step law
(deviations 0.70, 0.042 and 0.0026 at dt = 0.5, 0.25, 0.125, ratios 16.8 and 16.0 against RK4's 2⁴).

On the decoder's own window, 60 equal-split times on (0, 1/γ] at Q = 20:

| N | bonding mode ψ₁ | site 0 | centre site |
|---|---|---|---|
| 5 | 2.40 · 10⁻³ | 0.260 (108 times) | 0.223 (93 times) |
| 7 | 1.83 · 10⁻³ | 0.217 (119 times) | 0.139 (76 times) |
| 9 | 1.53 · 10⁻³ | 0.185 (121 times) | 0.149 (97 times) |

(λ_min per shot.) From the same shots the site seed and the centre site hold 76 to 121 times the bonding mode's
information in the worst direction, the end mixture 28 to 43 times (6.8, 6.9 and 6.5 · 10⁻²).

Part of that factor is one direction. The bonding mode's worst direction is the uniform rescaling of every
bond (its eigenvector overlaps (1, …, 1)/√(N − 1) to 0.997, 0.998 and 0.998): rescaling every bond rescales H and
leaves ψ₁ an eigenstate, so without dephasing its populations do not move at all
([THE_BOND_AND_THE_SEAT](THE_BOND_AND_THE_SEAT.md): "a uniform change of all bonds leaves the prepared ψ₁
population stationary"), and with it they move only through the changed J/γ
([THE_BOND_AND_THE_SEAT_AT_FIXED_GAMMA](THE_BOND_AND_THE_SEAT_AT_FIXED_GAMMA.md)). A front sees that rescaling as
its speed, and the localized carriers' worst directions lie elsewhere (overlap 0.33 or less). On the directions
orthogonal to the uniform rescaling, the bonding mode reads 2.2 · 10⁻², 7.8 · 10⁻³ and 4.2 · 10⁻³ (9.2, 4.2 and 2.8
times its λ_min), while the localized carriers barely change; the site seed holds 12, 28 and 44 times the bonding
mode there, the centre site 11, 18 and 35 times. So the uniform rescaling carries a factor 9.2 of the 108 at N = 5 (the
remainder 12) and a shrinking share as N grows (4.2 of 119, 2.8 of 121), while the remainder grows with N. Bond by bond the worst variance bound max_b (F⁻¹)_bb,
every other bond unknown, is 33 to 63 times lower for the site seed and the centre site than for ψ₁. The direction
holds under every summary read here; its size depends on which one is read.

Given each carrier the best of the schedules tried (equal splits over γt ∈ (0, 1],
(0, 3] or (0, 10], or all shots at any one time of those grids), the gap widens at Q = 20: ψ₁ reads 3.5, 1.8 and
1.5 · 10⁻³ at N = 5, 7, 9, ψ₂ 3.0, 1.3 and 0.97 · 10⁻³, site 0 0.84, 0.34 and 0.21, the centre site 1.11, 0.79 and
0.64, the end mixture 0.20, 0.12 and 0.065. For the site seed and the centre site the best is a single time: one
snapshot of positions locally reads every bond, given γ and the time, with no outcome to spare (N − 1
independent probabilities for N − 1 bonds). The schedule matters at both watching strengths (at Q = 20 the site
seed loses a factor 16 to 83 between the windows γt ≤ 1 and γt ≤ 10).

At strong watching the decoder's window starves the site seed. At Q = 1.5 it reads 5.9 · 10⁻⁶ at N = 7 and
6.3 · 10⁻⁹ at N = 9 there (by the exact derivative of the propagator; the condition of that Fisher matrix is
10⁸, and the eigenvector route reads 6.5 · 10⁻⁹): its excitation has not reached the far bonds by t = 1/γ. Given the best schedule tried
it reads 8.3, 2.2 and 0.75 · 10⁻² at N = 5, 7, 9, against ψ₁'s 6.6, 2.0 and 1.3 · 10⁻³, and the centre site 0.27,
0.095 and 0.030, while ψ₂ is uneven (8.8 · 10⁻³, 1.8 · 10⁻², 5.7 · 10⁻⁴). On the bond/seat page's baseline
schedule (uniform profile) the centre site holds 34 times ψ₁'s smallest eigenvalue (6,421,954 against 190,443) and
the site seed 0.15 times (29,389); on its selected schedule, which that page's exchange search chose for ψ₁, the
centre site holds 34 times again (11,584,496 against 344,084) and the site seed 0.19 times (64,723).

The decoder's own estimate spends part of what its shots carry. For one known bond, its strength read by the
projection of the linear letters is, linearized, linear in the measured populations, and its variance follows
from the multinomial covariance of the shots. Its **efficiency** is the single-bond Cramér–Rao bound 1/F_bb of the
same shots over that variance (1 when the estimate uses all the information; it cannot exceed 1). On the decoder's
window, mean over bonds, the efficiency on ψ₁ is 0.20, 0.20 and 0.32 at N = 5, 7, 9 at Q = 20, and 0.70, 0.32 and
0.58 at Q = 1.5. At weak watching the decoder's linearized estimate uses a fifth to a third of the information ψ₁'s
shots carry.

**An actual decode.** The bounds can be sampled. With the shots split equally over the decoder's window, draw the
positions and read the defect by maximum likelihood on the exact populations, γ and J known (300 draws per cell).
One defect of 0.02, the decoder's own task, located with the right bond and sign (δ on a grid of 240 values in
[−0.06, 0.06]):

| N, planted defect | carrier | 6,000 shots | 24,000 | 60,000 | 240,000 |
|---|---|---|---|---|---|
| 5, bond 3 at −0.02 | ψ₁ | 0.27 | 0.40 | 0.51 | 0.74 |
| | site 0 | 0.60 | 0.90 | 0.98 | 1.00 |
| | centre site | 0.49 | 0.81 | 0.94 | 1.00 |
| 7, bond 5 at −0.02 | ψ₁ | 0.16 | 0.34 | 0.44 | 0.79 |
| | site 0 | 0.50 | 0.82 | 0.97 | 1.00 |
| | centre site | 0.48 | 0.82 | 0.98 | 1.00 |

With bond 1 at +0.02 instead the rows read alike (ψ₁ 0.50 and 0.49 at 60,000 shots, the localized carriers 0.95 to
0.99). At 60,000 shots the localized carriers locate the defect in 94 to 99 of 100 draws and ψ₁ in 44 to 51; at
240,000 ψ₁ still misses 14 to 32. For the same success ψ₁ needs more than ten times the shots at N = 5 (at 240,000
it reads 0.74 and 0.68, the site seed at 24,000 0.90 and 0.81) and about ten times at N = 7 (0.79 and 0.86 against
0.82 and 0.80). Where the answers go, counted at 60,000 shots: at N = 5 ψ₁'s commonest wrong answer is the
decoder's confuser, a weakened bond 3 read as a strengthened bond 1 in 80 of 300 draws and the reverse in 86; at
N = 7 a weakened bond 5 goes most often to a strengthened bond 2 (63). The localized carriers' few misses there (at
most 19 of 300) land on another bond, the sign kept in all but two.

With every bond unknown (one true configuration per cell, δ drawn in [−0.01, 0.01] per bond; maximum likelihood by
Fisher scoring with step halving, to steps below 10⁻⁹), the site seed meets the Cramér–Rao bound at 60,000 shots
already (the error covariance's trace 0.97 and 0.99 times the bound's, bootstrap 0.87 to 1.07). ψ₁ does not: 71
and 84 of its 300 draws have not settled after 80 steps (its worst direction holds a hundredth of the site seed's
information, §5), and the estimate as run, unsettled draws included, has 1.77 and 1.35 times the bound's variance
(bootstrap 1.39 to 2.16 and 1.15 to 1.54), its spread, 0.12 and 0.13, ten times the defects drawn. At 600,000 shots
every ψ₁ draw settles and meets the bound (1.04 and 0.94, bootstrap 0.83 to 1.18). There the bounds' own spreads put
ψ₁'s error at 8.0 and 7.2 times the site seed's at N = 5 and 7.

**The de-lossed path with other carriers.** `DefectDecoder.DecodeDeviation` reads one number per site, the mean
over the painters' grid (60 times on [0, 1/γ]) of the purity's deviation, against one letter per bond calibrated at
δ = +0.02, and decodes by the single-letter rule. Rebuilt for any carrier, it reads ψ₁'s confuser at N = 5 at the
squared residual ratio 516.52 that the live `inspect --root decoder` reports as 516, and noiselessly it decodes
every defect planted here (two weakenings at −0.02 and one strengthening at +0.01 per N) with every carrier tried.
Its dictionary for the site seed is better conditioned at N = 5 (0.126 against ψ₁'s 0.008, worst cosine −0.58
against −0.95); at N = 4 the site seed's is still better conditioned (0.039 against 0.007) but its worst pair is more
opposite (−0.90 against −0.85), and its noiseless ratios, 15.42 and 29.42, are the lowest read. One time-averaged
number per site keeps little of what a localized carrier's trajectory carries, as with the linear letter of §2.

## 6. Who carries it: the never-caught front and the halo

[F126](../docs/ANALYTICAL_FORMULAS.md) solves exactly the populations this page reads
([PROOF_DEPHASING_FRONT_RENEWAL](../docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md)). Split at the diagonal refill,
every population is a sum over the number j of refills (the times the walker was caught, in the unravelling by
site projectors; [FRONT_PEDIGREE](FRONT_PEDIGREE.md) reads that distribution at the front): the term j = 0 is the
carrier's clean evolution under the universal decay e^{−Γt}, Γ = 4γ (for a site seed, the coherent front), and
everything refilled once or more is the **halo**. The derivation never uses that the couplings are equal or that the
seed is one site, so it carries over to every bond defect and every carrier of this page, the first term becoming
the clean evolution of the carrier's populations. Read on the finite chain through a counting parameter on the
refill (the order j is the coefficient of λ^j when the refill is weighted by λ), the term j = 0 meets that clean
evolution, e^{−Γt}·diag(e^{−iht}ρ₀e^{iht}), to 6 · 10⁻¹⁶, and each order's trace meets the Poisson weight
e^{−Γt}(Γt)^j/j! to 10⁻¹⁵. Each letter splits the same way into a never-caught part and a halo part.

**The bonding mode's never-caught term is blind along the uniform rescaling, exactly.** ψ₁ is an eigenstate of
every rescaled hopping matrix, so its clean evolution does not move when every bond is rescaled together: the
never-caught term's slope along the uniform rescaling is zero (at most 1.2 · 10⁻¹⁰ here, the rounding of a finite
difference at ±10⁻⁶, against 1.2 to 3.3 · 10⁻² for the full populations at Q = 20). Write F⁽⁰⁾ = mean_t Σ_j
∂P⁽⁰⁾∂P⁽⁰⁾ᵀ/P_j, with the full populations below: what the readout would carry if the halo had no slope. It has a
zero eigenvalue (below 10⁻¹⁶, its largest 0.42 to 0.72). Since ψ₁'s worst direction of §5 lies within 0.997 to
0.998 of the uniform rescaling, nearly all of what ψ₁ reads in its worst direction comes through the halo,
walkers caught and released, which carry the scale of the bonds only through J/γ: the
[bond/seat](THE_BOND_AND_THE_SEAT_AT_FIXED_GAMMA.md) page's "changing every hopping changes J/γ₀", split by catch
count.

**The site seed's information lies in both.** At Q = 20 its F⁽⁰⁾ has λ_min 0.141, 0.116 and 0.090 at N = 5, 7, 9:
about half of its full λ_min (0.54, 0.53 and 0.49 of 0.260, 0.217 and 0.185) and 59 to 63 times the bonding mode's
full λ_min. The coherent front by itself already reads every bond; the halo adds the rest (the two parts' cross
terms included, F is not the sum of the parts' matrices). The centre site's F⁽⁰⁾ reads 0.109, 0.086 and 0.079.

Over the decoder's window the never-caught term is a quarter of the population (its trace is e^{−Γt}, whose mean
over γt ∈ (0, 1] is (1 − e⁻⁴)/4 ≈ 0.25, 0.24 on the 60 times), but at Q = 20 it carries 0.53 to 0.77 of the Fisher-weighted letters
(the projection share Σ ∂P⁽⁰⁾·∂P/P over Σ ∂P·∂P/P, a trace over every bond direction, which with the halo's share
sums to 1): ψ₁ 0.73 to 0.77, the site seed 0.53 to 0.61, the centre site 0.60 to 0.70. Weight for weight, the
never-caught term is the more sensitive. At Q = 1.5 the halo carries 0.83 to 0.88
of the letters for every carrier.

## What this does not claim

Shots are drawn only in §5's decode, at Q = 20 on the decoder's window at N = 5 and 7, 300 draws per cell; elsewhere
§5 bounds the shot noise and does not sample it. The readings are first order in δ (the decodes of §1 and §5 and
the de-lossed path excepted), at Q = 20 unless a section says otherwise (§4 also without dephasing, §5 and §6 also
at Q = 1.5), on open chains up to N = 20
(arrival), N = 14 (carriers) and N = 9 (Fisher). The arrival letter is a threshold reading, with the bump and the far end in it,
not the clean step of the walk-time page's N = 60 runs, and its worst cosine depends on the threshold. How the
arrival time, a threshold reading, behaves under shot noise is not computed; the decoder's efficiency is that of its
linearized estimate, and its de-lossed path is compared without noise, on three planted defects per N. The decode
knows γ and J, assumes one defect in its first part, and holds one true configuration per cell in its second. The
sum of §4 accounts for the mean of the cosines and for the bonding mode's weakest direction, not for which pair is
the worst. There is no detector error
in the model. Nothing here was
measured on hardware.

## Reproduce

```bash
dotnet run --project compute/MirrorWorld -- walk decode         # the C# readings, about 10 s
dotnet test compute/MirrorWorld.Tests --filter "FullyQualifiedName~WalkTimeTests"
python simulations/decoder_carrier_dictionary.py                # the exact second route, about 150 s
python simulations/decoder_letters_sum_and_shots.py              # the sum of the letters, the decode at finite shots
```

## Links

- [HANDSHAKE_GEOMETRY](../hypotheses/HANDSHAKE_GEOMETRY.md): the decoder, its anti-collinearity, the mirror
  equivariance, the phase-blindness gate (the observable half of §2).
- [COUPLING_DEFECT_WALK_TIME_STEP](COUPLING_DEFECT_WALK_TIME_STEP.md): the step law this page reads backwards.
- [PROOF_DEPHASING_FRONT_RENEWAL](../docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md) (F126) and
  [FRONT_PEDIGREE](FRONT_PEDIGREE.md): the renewal ladder of §6, the exact solution of every population read here.
- [PERSPECTIVAL_TIME_FIELD](../hypotheses/PERSPECTIVAL_TIME_FIELD.md): the carriers the decoder's painting stage was
  run on, and where its rescaling picture degrades.
- [F124](../docs/ANALYTICAL_FORMULAS.md) and `BandEdgeResolutionLimitClaim`: the every-bond inverse problem on the
  bonding mode.
- [THE_BOND_AND_THE_SEAT_WITH_FINITE_SHOTS](THE_BOND_AND_THE_SEAT_WITH_FINITE_SHOTS.md): the full-position Fisher
  matrix and a finite-shot classifier on ψ₁, the object of §5.
- `compute/RCPsiSquared.Diagnostics/Foundation/DefectDecoder.cs` and `SingleExcitationPainterDictionary.cs`: the
  decoder and its dictionary on the bonding mode.
