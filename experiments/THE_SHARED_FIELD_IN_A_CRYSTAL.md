# The shared field in a crystal: the coherence-variance law meets multiple-quantum NMR data

**Date:** 2026-10-17
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 5.5)
**Verifier:** [`mq_nmr_shared_field.py`](../simulations/mq_nmr_shared_field.py) (parts A to C by default, 30 checks, six of them recorded negatives across the three runs of part C; `--rebirth` adds part D, two more checks, about five minutes)
**Data, not ours:** H. Cho, P. Cappellaro, D. G. Cory, C. Ramanathan, *Decay of highly-correlated spin states in a dipolar-coupled solid*, [arXiv:cond-mat/0608620](https://arxiv.org/abs/cond-mat/0608620) (CaF₂, Fig. 4 and Table 1); H. G. Krojanski, D. Suter, *Decoherence in large NMR quantum registers*, Phys. Rev. A 74, 062319 (2006) ([preprint](https://e3a.physik.tu-dortmund.de/Publications/eprints/Decoherence.pdf)), and Phys. Rev. Lett. 93, 090501 (2004).
**A word fence:** *order* n (or p) is the NMR coherence order, the change of total magnetization Δpopcount; *cluster size* N is the number of correlated spins as the experimenters estimate it; η is the correlation coefficient of the sites' Z fields, as on [the singlet page](SINGLET_UNDER_A_SHARED_FIELD.md); κ is the ratio this page's fit actually determines (defined in part C); *rate* is the inverse of a 1/e time, *width* T_d the standard deviation of a Gaussian decay exp(−t²/2T_d²).

## What this is about

Most of what this repository holds is checked against itself or against qubits we prepare ourselves. On this page one of our formulas meets somebody else's measurement on somebody else's material, a crystal in an NMR lab. The Markovian prediction was written down before any of the data were read; the quasi-static reading was taken up after it failed, when the paper's abstract had already reported the square-root growth and widths linear in the order, so only its shape and constancy criteria were set before the figure was read.

Under Z noise whose sites share a field with correlation η, a coherence |A⟩⟨B| collects the quadratic form 2γ[(1 − η)·h + η·p²], with h its Hamming distance and p its change of popcount ([Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) §2). In NMR, p is the coherence order of spins 1/2 (for a spin 3/2 written in two qubits it is not, [The Sodium Nucleus](THE_SODIUM_NUCLEUS_TWO_MIRRORS.md#the-irreducible-tensors-in-strings)), and multiple-quantum experiments prepare clusters of tens to thousands of correlated nuclear spins and watch coherences of each order die. That is the formula's own object, measured.

Two readings were tested. Read as a decay *rate*, the memoryless Lindblad reading, the formula fails: it gives a rate that grows linearly with the cluster size, and the crystals show roughly its square root; in the five correlation families tried, correlated Markovian Z noise lowers the rate and leaves the growth about linear. Read as the phase *variance* of a field frozen during the decay, the same quadratic form gives the square root, and one constant describes the CaF₂ decay widths across cluster sizes 6 to 96 and orders 2 to 24 to about 5 % rms, its per-size values showing no trend. But the shape it predicts at a fixed size, 1/T_d² affine in n², does not do better than a straight line: unweighted, a line in n fits four of the five sizes slightly better, and the authors' own straight-line fit describes the whole figure slightly better than the law. The quasi-static reading gets the growth with size; the law in the tested form does not get the dependence on order better than an empirical line.

## Abstract

Part A: for Z jumps with a real symmetric Kossakowski matrix a, the generator is diagonal on every |i⟩⟨j| with rate 2vᵀav, v_k = (z_k(i) − z_k(j))/2, met exactly against the full superoperator with integer covariances (residual 0); averaged over the sign patterns of a K-site disagreement it is 2·tr(a), exactly. Part B: read as a Markovian rate, the 1/e decay of a K-spin cluster signal scales as K^0.95 to K^1.20 in five covariance families (K = 7 to 16); by Jensen's inequality correlations can only lower the 1/e rate below the mean 2·tr(a), for a fully shared field to 0.32 to 0.34 of it over K = 7 to 16; Krojanski and Suter measure K^0.48 under free evolution and K^0.43 under decoupling in adamantane. Part C: read as the variance of a quasi-static Gaussian field, the width is T_d ∝ 1/√(N + κ·n²), κ = η/((1 − η)c), h = c·N; on Fig. 4 of Cho et al. one κ fits every point at 4.96 % and 4.86 % rms on two digitizations (κ = 0.324 and 0.319), the per-size κ shows no trend over five sizes (0.23 to 0.38, N = 12 to 96); the authors' own fit (four constants) reaches 4.16 % and 4.04 %, their functional form refit with two 4.52 % and 4.40 %, the Fedorov-Fedichkin two-term form 8.9 % and 8.7 %; unweighted, at fixed size a straight line in n fits better than the law at N = 12, 28, 36 and 66 by 0.3 to 1.6 points of rms, the two split at N = 96, so the shape criterion set in advance is not met, by margins small against the figure's error bars. η itself is not determined: c = 1 reads η = 0.245, c = ½ reads η = 0.139, at identical rms. Part D: an eight-site chain with hopping plus dephasing gives a log-log exponent of 0.48 at one γ, a coincidence; it moves from 0.161 to 0.925 with γ, the rates being an offset plus a slope that approaches 2γ per added spin.

## The sweep: what the repository held

Swept 2026-10-17 by two survey agents and a review agent, stores by name.

- `docs/ANALYTICAL_FORMULAS.md`: no entry names NMR, multiple-quantum coherence, coherence order, spin counting, TROSY or long-lived states; the Hamming-ladder entries and the correlated Z-dephasing residual (lines 120 to 123) are there, none read against NMR.
- `docs/proofs/`: [PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) §2, the source of the quadratic form and its η reading, carried in code by `AbsorptionTheoremClaim` (compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs, line 94).
- `experiments/`, the nearest prior result: [the price-pair flight](PRICE_PAIR_HARDWARE_PREDICTION.md) measured exactly this object on hardware, a per-bond covariance c_ij of the dephasing inside a Gaussian covariance model, with local quasi-static noise (rates adding in quadrature) and a common mode as named decoding rows; a pairwise anti-correlated reading of a quasi-static fit was refuted there and re-attributed to a coherent nearest-neighbour ZZ coupling. [The absorption rung ladder](ABSORPTION_RUNG_LADDER_HARDWARE_PREDICTION.md) runs quasi-static Gaussian and quasi-static common-mode dephasing as named hypotheses. [Singlet Under a Shared Field](SINGLET_UNDER_A_SHARED_FIELD.md) uses the same η for two spins in the Markovian reading and does not say how large η is anywhere; [The Mirror under a Shared Field](THE_MIRROR_UNDER_A_SHARED_FIELD.md), the same Kossakowski structure for the palindrome; [Noise Robustness](NOISE_ROBUSTNESS.md) Q1 and [Decoherence Relativity](DECOHERENCE_RELATIVITY.md), the collective rate 2γ(Δpopcount)² at Δp = 2; [XOR Space](XOR_SPACE.md), the coherences with h = N; [τ_max and the spectral gap](TAU_MAX_SPECTRAL_GAP.md) notes that a claimed NMR validation appears nowhere in the repository. Searches for "coherence order", "multiple-quantum", "spin counting" and "Pines" returned nothing.
- `docs/`: [The Atmosphere and the Cancelled Formulas](../docs/THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS.md), the host of the premise this page measures against, γ₀ as the rate of a field each site has to itself, and where that cancellation stops; [Schrödinger's Cat Translated](../docs/quantum/SCHRODINGERS_CAT_TRANSLATED.md), superdecoherence, the measured η = 1 end in ion traps; [Predictions](../docs/PREDICTIONS.md), a tier header naming NMR as a test bed without an NMR row; [Painter Alternation NMR Bridge](../docs/carbon/PAINTER_ALTERNATION_NMR_BRIDGE.md), a conditional translation that predicts no NMR observable; [Internal and External Observers](../docs/historical/INTERNAL_AND_EXTERNAL_OBSERVERS.md), whose TROSY reading carries a banner; [the literature review](../docs/LITERATURE_REVIEW.md)'s Cheshkov-Sinitsyn entry, spectral palindromes of closed Hamiltonians, a different object.
- The typed layer: `CorrelatedZDephasingPriceClaim` and its live witness; MirrorWorld's `Atmosphere`, whose statement that the Mirror's price survives as the mean is part A's second identity. `fw.Confirmations` and the C# registry: `price_pair_locality_marrakesh_july2026`, the hardware covariance measurement above, and relaxation-rate entries; no NMR entry. The OpenArcs registry, `docs/GLOSSARY.md`, `docs/CAUGHT_ERRORS.md`: nothing on NMR.

Checked for adjacency: the singlet page and the price-pair flight ask this page's question, how much of a field the sites share; the flight measured covariances on hardware, this page reads a ratio off NMR widths in the quasi-static reading, and the singlet page's Markovian η is neither. The price-pair flight is also the fence on part C: a quasi-static fit that reads as a shared field can be a coherent coupling.

**Prior art, read after the derivation.** Cho et al. summarize Palma, Suominen and Ekert (Proc. R. Soc. Lond. A 452, 567 (1996)) as giving decay exp(−f·Γ) for independent and exp(−f²·Γ) for correlated noise, f the Hamming distance; the correlated end of the quadratic form here is p², which equals f² only when every disagreeing bit flips the same way. Fedorov and Fedichkin (2006, as cited by Cho et al.) write the multiple-quantum decay as a two-term mixture p·exp(−αn²t²) + (1 − p)·exp(−(N/2)αt²), with p defined from the couplings and observed to stay constant in adamantane. That a Gaussian decay signals a frozen dispersion of local fields is Abragam's, and the √N growth of a Gaussian rate is the central-limit growth of a sum of independent random phases.

## Part A: the quadratic form, exactly

For jumps on the Z_k with a real symmetric Kossakowski matrix a, D(ρ) = Σ_kl a_kl(Z_k ρ Z_l − ½{Z_l Z_k, ρ}). On |i⟩⟨j| every term is a scalar, and the decay rate is ½·Σ_kl a_kl (z_k(i) − z_k(j))(z_l(i) − z_l(j)) = 2vᵀav (an imaginary antisymmetric part of a would add a frequency shift, outside this page). The verifier builds the full superoperator for integer a (3·I and random B·Bᵀ, N = 3, 4, 5), finds it exactly diagonal and equal to 2vᵀav entry by entry, residual 0. Averaged over the 2^K sign patterns of a K-site disagreement the cross terms cancel and 2vᵀav averages to 2·tr(a) exactly (K = 3, 5, 7). For uniform correlation, a = γ[(1 − η)I + η·11ᵀ] gives 2γ[(1 − η)h + η·p²], p = Σv.

## Part B: the Markovian reading fails

**Written before any multiple-quantum data were read:** in the Markovian reading the 1/e rate of a K-spin cluster signal grows linearly in K, whatever η; a clearly sublinear measured scaling means this family does not describe the experiment.

Krojanski and Suter measure, in adamantane, rate ∝ K^0.48 under free dipolar evolution (Gaussian decays, up to 4899 correlated spins) and K^0.43 under MREV-16 decoupling (exponential decays, which they say may indicate that the interaction can be modeled as a Markovian process). The prediction fails, including in the regime where the Markovian picture is most plausible.

Their own reading of the sublinear growth is that the perturbations on different spins are correlated. Within Markovian Z noise, the families tried do not produce it. The cluster signal, uniform over the sign patterns of K sites, was computed at γ = 1, K = 7 to 16:

| covariance | exponent of the 1/e rate in K | 1/e rate over 2K at K = 16 |
|---|---|---|
| local, η = 0 | 1.000 | 1.000 |
| uniform η = 0.5 | 0.991 | 0.832 |
| fully shared, η = 1 | 0.946 | 0.322 |
| exponential, ρ = 0.9 | 1.197 | 0.723 |
| power law (1 + d)^−½ | 1.087 | 0.846 |

(The alternating family ρ = −0.9 is the exponential one under the gauge v_k → (−1)^k v_k and reads the same.) Since S(t) = E[e^{−rt}] ≥ e^{−E[r]t}, the 1/e rate never exceeds the mean 2·tr(a): correlations lower it, here to as little as a third, but not by a factor that vanishes with K: for the fully shared field the ratio stays between 0.32 and 0.34 over K = 7 to 16, and the exponent 1.197 of the exponential family is a crossover toward linear (its ratio climbs from 0.615 at K = 7 to 0.723 at K = 16). Five families at K ≤ 16 are a survey, not a proof for every covariance, and the uniform average over sign patterns sets h = K for every term and weights the orders binomially, which is not the multiple-quantum state of the experiments. The decoupled K^0.43 is not explained on this page.

## Part C: the quasi-static reading

Let the field be frozen during the decay: site k feels δ_k, Gaussian with covariance Σ ∝ a. The coherence |i⟩⟨j| acquires the phase 2vᵀδ·t, and its average is exp(−2t²·vᵀΣv): a Gaussian in time, exp(−t²/2T_d²) with T_d = 1/(2√(vᵀΣv)), the standard deviation Cho et al. fit to each decay. The same quadratic form is now a variance, and its square root sets the rate. In CaF₂ under free evolution the frozen field is the dipolar local field of the surrounding fluorine spins, a coherent coupling read as a static dispersion: ZZ couplings inside the cluster give no relative phase when both partners flip, the flip-flop part of the dipolar Hamiltonian moves the operator rather than adding a phase (Fedorov and Fedichkin drop it), and the measured T_d is a loss of refocusing after the reversed double-quantum evolution, not the free dephasing of one fixed coherence. The price-pair flight is the house precedent that such a reading can hide a coherent term. The fit below therefore characterizes how the local fields of a cluster co-vary, not an external bath.

For a cluster of N correlated spins in coherence order n, take h = c·N with an unknown constant c. (In Cho et al.'s Eq. 7 a term has n raising, m lowering and q longitudinal factors, so h = n + m and the order is n − m; the paper gives N = n + m + q beside its Table 1. This page writes n for the order.) Substituting into (1 − η)h + η·n² and dividing by (1 − η)c,

    T_d ∝ 1 / √(N + κ·n²),   κ = η / ((1 − η)·c),

so only κ, the shared variance per unit order² over the local variance per correlated spin, is determined by widths alone.

**Written after reading the paper's abstract and before its equations or Fig. 4, in its η form (with c = 1, which is the κ form reparametrized; the criteria below do not depend on c):** one constant and one scale must fit all cluster sizes and orders; it must not drift with N; T_d·√N must depend on n/√N alone; fail if it drifts or if, at fixed N, 1/T_d² is not affine in n². The abstract already states that the decay's standard deviation scales linearly with coherence number; the criterion was set against that, without a numerical threshold.

Cho et al. summarize Fig. 4 by T_d = A·N^−0.57 − B·n·N^−1.21 (A = 51.4 µs, B = 15.6 µs; in their Eq. 10 the exponents are rounded to −½ and −1). Fig. 4 of the arXiv version was digitized twice from one screenshot of the figure: by eye from pixel positions (reading A, 41 points) and by a separate agent through pixel analysis (reading B, 45 points, ±0.02 µs on isolated markers); the two agree to 0.04 µs on shared points, and reading B flags five points where markers merge, hide or meet the axis (N = 28 and 36 at n = 12, N = 36 at n = 16, N = 66 at n = 24, N = 96 at n = 20). The series map to cluster sizes through Table 1: MQ growth times 43.3, 86.8, 130.2, 173.6, 217, 260.4, 303.8 µs give N = 6, 8, 12, 28, 36, 66, 96. The per-size rows below need at least three (κ) or four (law against line) points per size, so N = 6 and 8 enter only the global fits.

| | reading A | reading B |
|---|---|---|
| T_d·√N at the lowest order, N = 6 … 96 | 38.5, 38.2, 36.5, 40.8, 36.0, 37.2, 35.6 | 38.5, 38.2, 36.4, 40.9, 36.1, 37.3, 35.5 |
| one κ for all points | 0.324, rms 4.96 % | 0.319, rms 4.86 % |
| the same fit read as η (c = 1 / c = ½) | 0.245 / 0.139 | 0.242 / 0.138 |
| the authors' own fit (four constants, fitted by them) | rms 4.16 % | rms 4.04 % |
| the same with exponents rounded to −½ and −1 | rms 23.3 % | rms 25.1 % |
| the authors' functional form A/√N − B·n/N, refit (2 parameters) | rms 4.52 % | rms 4.40 % |
| Fedorov-Fedichkin two-term (p fitted freely, α in the scale; T_d where S = e^−½) | rms 8.91 % | rms 8.71 % |
| κ per size, N = 12, 28, 36, 66, 96 | 0.32, 0.38, 0.23, 0.37, 0.25 | 0.33, 0.38, 0.26, 0.37, 0.23 |
| trend d ln κ / d ln N | −0.07 | −0.11 |
| rms at fixed N, law / straight line in n: N = 12, 28, 36, 66, 96 | 1.3/1.0, 1.2/0.7, 1.4/0.7, 4.7/3.2, 1.1/1.2 % | 1.3/1.0, 1.2/0.6, 2.3/1.8, 4.8/3.2, 1.7/1.6 % |

Removing reading B's five flagged points (the verifier's third run) moves κ to 0.302 at 4.51 % rms and leaves every check's verdict as it is.

What holds: the scaling. T_d·√N is nearly constant at the lowest order, and one κ serves all sizes; its per-size values scatter (0.23 to 0.38) without a monotone trend. The two-term form, the theory the experimenters compared with, fits worse. The scaling does not discriminate, though: the authors' straight line has the same collapse in n/√N, and √N alone is the central-limit growth of any sum of random phases.

What does not hold: a better shape than a line. At fixed size and unweighted, the widths sit closer to a straight line in n than to the law's even dependence on n at N = 12, 28, 36 and 66, by 0.3 to 1.6 points of rms; the two split at N = 96 (the law ahead in reading A, the line in reading B), and the authors' own straight-line fit describes the whole figure slightly better than the law. Against the figure's error bars, about 1 % at low order and up to ±0.5 µs at high order, these margins are small; only N = 66 and the ordered residual signs carry real weight. Across sizes the residual of the one-κ fit is also not a pure offset: N = 6, 8 and 28 lie above the curve at every point (by up to 10.5 % at N = 6), N = 96 mostly below, and N = 36 and 66 run in ordered signs.

What part C does not show: a value of η (the same widths read η = 0.24 or 0.14 depending on c); a fit weighted by the error bars, which were not digitized; more than one crystal and one figure.

## Part D: renewal on a small chain

Whether the renewal of [the dephasing front](../docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md), the light catching the excitation and the hopping releasing it, bends the Markovian scaling was checked on an eight-site chain with local Z dephasing, through the fidelity of a string X on sites 0 … K − 1. At γ = 0.2 the log-log exponent is 0.478 with XX + YY hopping and 0.426 with a dipolar-like ZZ − ½(XX + YY) coupling, the measured 0.48 and 0.43 to two digits. It is not a scaling law. With XX + YY the exponent is 0.161, 0.478, 0.785 and 0.925 at γ = 0.05, 0.2, 0.5 and 1, and the rates are the light-free spreading rate plus a slope that approaches 2γ per added spin (0.83, 0.82, 0.94 and 0.98 of 2γ over K = 4 to 7); the hopping and the light add, and the renewal does not pull the growth below linear.

## What this changes

- The Markovian reading, the one the theorem pages use, does not describe how decoherence grows with cluster size in these crystals, and in the five correlation families tried correlations do not rescue the growth law. A page that turns a repository rate into a crystal rate should say which reading it uses.
- The quasi-static reading of the same quadratic form gets the growth with cluster size and carries one constant whose per-size values show no trend; in the tested form it does not describe the dependence on order better than an empirical line. Among the NMR pages the sweep found, this is the first comparison with NMR measurements, and its verdict is mixed.
- The shared fraction η needs more than widths: an independent count of the transverse factors per cluster term, or a direct covariance measurement of the kind the price-pair flight made on hardware.

## Open

- The order dependence. The transverse count h of a cluster term is itself tied to its order (h ≥ |n|), which puts a term linear in n into the variance that the law as tested leaves out. Whether that accounts for the straight lines needs a model of h(n, N) fixed before it meets the data; it has not been tested.
- The decoupled K^0.43 of Krojanski and Suter: exponential in time, sublinear in size, explained by neither reading here.
- Whether the adamantane data, or the CaF₂ time-suspension data (Fig. 7 of Cho et al., the same features with widths about 70 times longer), carry the same κ.

## Links

- [PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) §2, the quadratic form and its η.
- [The price-pair flight](PRICE_PAIR_HARDWARE_PREDICTION.md), the hardware covariance measurement and its refuted shared-field reading.
- [The Mirror under a Shared Field](THE_MIRROR_UNDER_A_SHARED_FIELD.md), [Singlet Under a Shared Field](SINGLET_UNDER_A_SHARED_FIELD.md), [XOR Space](XOR_SPACE.md), [The Atmosphere and the Cancelled Formulas](../docs/THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS.md).
- [The dephasing front renewal](../docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md), part D's hypothesis.
- [Painter Alternation NMR Bridge](../docs/carbon/PAINTER_ALTERNATION_NMR_BRIDGE.md), the conditional NMR translation of the carbon chains.
