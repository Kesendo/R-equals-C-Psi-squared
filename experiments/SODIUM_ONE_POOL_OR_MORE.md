# One Sodium Pool or More: the Mean-Rate Identity against Published Agarose Data

**Date:** 2026-11-14
**Authors:** Thomas Wicht, Claude (Opus 5.5)
**Status:** A test (Tier 4 reading) of published fitted values against an exact relation of the sodium model. The script recomputes every number on this page; the data are read by hand from two open-access papers, a few values off their figures; the error bars are the authors' fit uncertainties, taken as independent, and carry none of their systematic errors.
**Script:** [`sodium_one_pool_or_more.py`](../simulations/sodium_one_pool_or_more.py), output in [`output.txt`](../simulations/results/sodium_one_pool_or_more/output.txt) (about six minutes; `--global` adds a differential-evolution check of the free-and-bound fit)

[The sodium page](THE_SODIUM_NUCLEUS_TWO_MIRRORS.md#what-a-sodium-measurement-can-read) holds one relation between three rates a measurement can read: in a single sodium pool the slow transverse rate is the mean of the two longitudinal rates, whatever the spectral densities, a residual splitting, an offset or a dynamic frequency shift. Its open list asks for that relation against measured rates and for what a voxel with two compartments does to it. This page takes the first published data we found that report all three rates for the same agarose concentrations: 2, 4 and 6 % agarose with 154 mM NaCl at 9.4 T, measured by Reichert, Schepkin, Kleimaier, Zöllner and Schad (2024).

The published fitted values do not fit one pool, and the conflict sits in one place. The measured slow transverse time equals the authors' weighted longitudinal time T₁m within the errors, as one pool with nearly equal J₁ and J₂ would give. Their T₁ fit, with its amplitudes fixed to the one-pool values, reports instead a split of the two longitudinal rates that one pool could only combine with a slow transverse time 14 to 17 % shorter. The mechanisms inside one pool that this page and the agarose page name move the gap the other way, and a signal stored along z by imperfect pulses has the right sign but not the size. Free sodium beside bound sodium does not reproduce the numbers. The authors' own triple-quantum amplitudes, compared across the three samples, rise in a way that one pool with the fitted split cannot give and two bound pools can. What remains is two readings: the fitted split is not one pool's, or there are two pools on scales beyond 10 to 17 µm. We do not know which; the absolute size of the longitudinal triple-quantum signal against one pool's would tell.

<!-- STAGE 0 SWEEP RECORD (for us, not for the reader: hidden on GitHub, read by every session, agent and hook)

One survey agent swept the stores by name for a sodium voxel holding two pools, without and with exchange, and for a check broken by a mixture; the data search afterwards was a narrow search for measured rates, not for theory (Tom, 2026-11-14), and every number on the page was read from the PDFs by hand. Three review rounds with the repository lens repeated the sweep and added the adjacency items below.

- docs/ANALYTICAL_FORMULAS.md: F167 indexes the sodium pages; its scope line says multiple compartments and exchange are outside the model, its open item asks for the mean-rate identity against measured sodium, and it calls the bright pairs' line merging two-site exchange in form. F166 (qubit names): blocks pair about their own centres while one global centre fails. Nothing on mixtures of generators otherwise.
- docs/proofs/: DIRECT_SUM_DECOMPOSITION (sectors palindromic about one centre); nothing on compartments, pools or exchange.
- experiments/: THE_SODIUM_NUCLEUS_TWO_MIRRORS, section "What a sodium measurement can read" (the identity at every splitting; a deviation says more than one environment or motion beyond Redfield), its open item, the triple quantum as the tool sodium MRI uses to tell compartments apart, the triple-quantum path (rank 1 to rank 3, zero where its two rates meet) of which this page's longitudinal argument is the twin, the satellite block in the Bloch-McConnell form of two-site exchange with its coalescence rule, Jaccard, Wimperis and Bodenhausen cited, cross-correlation of the quadrupolar part with rank-1 mechanisms named and not argued here, and an ordered medium not following one spectral-density curve; THE_NUCLEUS_AS_A_KEY (the identity as one of five relations and as a model check that is not new; the single-motion relation); SODIUM_IN_AGAROSE_LAYER_BY_LAYER (the identity applied to another agarose data set as the floor T2 at least T1m/2.5, broken there by a fitted population from below, and used to exclude a motional explanation; exchange without a frequency difference argued there and computed here; diffusion over the measurement far beyond the pores; the gel slowing the water, its mobile populations faster than saline; its candidates exchange between sites of different frequency, rank-1 fields and the second-order shift of slow motion shorten the central coherence, while its pathways stored along z relax at T1 and have this page's sign); THE_SODIUM_CHAIN (coupled nuclei, tensor products; its shared gradient exact results about correlation, not mixture); THE_SHARED_FIELD_IN_A_CRYSTAL (the first test against foreign data; Jensen on a rate mixture). Null results on mixtures of sodium pools: none. No sodium hardware flight exists. In qubit names, F81_VIOLATION_HARDWARE_BRIDGE (a T1 that fluctuates over a job: the unlabeled mixture, fitted as one exponential, read as a false asymptote at 5 to 14 sigma, and a two-rate mixture fitted all three qubits) and F120_MOMENT_TOWER_KINGSTON (a check bent by a comparison across two T1 epochs) are about the same object as this page, a check broken by a hidden mixture; they do not cite each other and neither cites the sodium pages. docs/outbound/NOISE_ASYMMETRY_SYMMETRY_SCALAR (Jensen: the averaged trajectory is slower than the exponential at the mean rate) and docs/quantum/TELEPORTATION_TRANSLATED (the unlabeled mixture read as false temperature) hold the same lesson. These are named here in the record only.
- hypotheses/ and reflections/: nothing on pools or on this identity.
- The OpenArcs registry: nothing on sodium, NMR or pools; the arc outbound_label_adapters holds the heating-leg two-rate mixture named above.
- fw.Confirmations: nothing on NMR; f84_heating_leg_attribution_kingston_july2026 and f120_moment_tower_kingston_june2026 are the hardware records of the mixture-broken checks above.
- docs/GLOSSARY.md: nothing on sodium pools.
- docs/CAUGHT_ERRORS.md: sodium entries on other objects; the entry on an uncoupled pair being a Kronecker sum marks the distinction this page keeps: two nuclei, coupled or not, live on a tensor product, while a mixture of pools, each nucleus in one of them, is a direct sum of generators.
- The typed layer: Spin32QuadrupolarMirrorsClaim and SodiumChainMirrorLadderClaim in Core; SodiumTwoMirrorsWitness gates the triple quantum at the mean longitudinal rate for one nucleus, SodiumChainWitness the chain; TwoByTwoExceptionalPointClaim and its witness carry the two-site-exchange form of the bright pairs; nothing on pools.
- A search of the repository for Reichert, Madelin, IRTQTPPI and TQTPPI: nothing before this page.

Checked for adjacency: the agarose page applies the same identity to a different agarose data set; this page links it, and the open items it answers on the sodium page, the key page, the agarose page and F167 point here.
-->

## The relation

Hubbard's rates for one spin 3/2 under secular Redfield relaxation by fluctuating electric field gradients give the longitudinal magnetization two components, with rates 2J₁ (weight 0.2) and 2J₂ (weight 0.8), and the transverse magnetization two, with rates J₀ + J₁ (weight 0.6) and J₁ + J₂ (weight 0.4), where J₀, J₁, J₂ are the spectral densities at 0, ω₀ and 2ω₀ in the convention of the agarose page; the sodium page and the script write the same rates three times larger (6J₁, 6J₂, 3(J₀ + J₁), 3(J₁ + J₂)). Hence

  1/T₂s = (1/T₁s + 1/T₁f)/2.

Here s and f mean slow and fast, as in the two papers; the agarose page uses the subscript s for short, so its T₁ₛ and T₂ₛ are the fast components. The relation is linear in the densities and holds for any shape of the spectral density. When J₁ ≥ J₂, as one spectral-density curve of any distribution of correlation times gives, the component of weight 0.8 is the slower one, the T₁s of a fit that fixes the weights; an ordered medium need not follow one curve ([the sodium page](THE_SODIUM_NUCLEUS_TWO_MIRRORS.md#beyond-spin-32)). The sodium page shows that the rates do not depend on a residual splitting; [the agarose page](SODIUM_IN_AGAROSE_LAYER_BY_LAYER.md) gives the reason for every static term diagonal in m: the jumps cannot move weight into or out of the central coherence, and such a term at most shifts its frequency (an offset and a dynamic shift do; a first-order splitting does not, since it moves the levels m = ±½ alike). The finite polarization at 9.4 T near room temperature (ħω₀/kT ≈ 2·10⁻⁵) and the terms the secular approximation drops are far below anything read here. The authors' own Eqs. 7 and 8 contain the relation; they do not test it.

Two consequences are used below. A single pool with J₁ ≥ J₂ has T₁m/2.5 ≤ T₂s ≤ T₁m, where T₁m, the authors' Eq. 17, is the mean of the two longitudinal rates weighted 0.8 and 0.2; the lower end is the agarose page's floor, the upper end follows the same way. And the relation's gap splits exactly as

  (1/T₁s + 1/T₁f)/2 − 1/T₂s = (1/T₁m − 1/T₂s) + 0.3·(1/T₁f − 1/T₁s),

a part that compares T₂s with T₁m and a part that is the split of the two longitudinal rates. T₁m is the combination the data fix best: in the first paper's Fig. 5 it stays flat while T₁s and T₁f drift apart as the late points are added (from about 47.5 and 35.8 ms to 49.2 and 30.9 ms), so the split comes from the late tail.

## The data

The first paper (NMR in Biomedicine) gives T₁s and T₁f from an inversion recovery with a triple-quantum filter, in its Table 2 for four phase cycles, fitted with its Eq. 13 jointly on the single- and the triple-quantum signal with the amplitudes held at their one-pool values 0.8 and 0.2. Its free fit (its Table 1), joint with the triple quantum as well, was unstable: at 2 and 4 % it put about 94 % of the single-quantum signal at one time near 43 ms and the other time at the fit's bound, while at 6 % it resolved a split (39 % at 45 ms, 60 % at 29 ms). Its triple-quantum peak (its Fig. 2B) is what the authors take as proof that T₁ is bi-exponential, and their Eq. 16 gives the fitted ratio of its amplitude to the single-quantum amplitude: 5.3, 8.7 and 15.9 % in phase cycle PC2, 3.1, 6.0 and 11.5 % in PC1. Its tabulated T₁m differs from Eq. 17 of its tabulated times by 0.1 to 0.4 ms. Its sampling is described as dense to 65 ms, which its stated step and count would end at 32 ms, and its Fig. 2A runs to about 640 ms; the script samples to 65 ms densely and to 305 ms in all.

The second paper (Magnetic Resonance in Medicine) gives T₂s, T₂f and their amplitudes, in its Table 1, from a triple-quantum sequence with a refocusing pulse at 26 evolution times between 0.4 and 200 ms, "determined in the same way" as its Eq. 11, which we read as the single-quantum signal alone fitted as a biexponential plus an offset; its 0 % sample gives T₂ = 54.19 ± 0.95 ms. Both papers are at 9.4 T with the same concentrations; neither says that the samples or the temperature were the same. The first paper's Fig. 6A shows its own T₂s values, about 1 ms shorter than the second's (44.1, 38.1 and 32.2 ms, read off the figure), and its Fig. 6C the transverse triple-quantum ratio, 23.1, 27.0 and 33.8 % (read off).

## The check

| agarose | T₁s (ms) | T₁f (ms) | T₂s for one pool (ms) | T₂s measured (ms) | gap in the slow transverse rate |
|---|---|---|---|---|---|
| 2 % | 49.22 ± 0.73 | 30.87 ± 0.23 | 37.94 | 45.03 ± 1.46 | 4.15 ± 0.75 s⁻¹, 5.6 σ |
| 4 % | 43.92 ± 0.98 | 25.65 ± 0.31 | 32.39 | 39.17 ± 0.73 | 5.35 ± 0.59 s⁻¹, 9.1 σ |
| 6 % | 37.92 ± 0.90 | 22.90 ± 0.27 | 28.56 | 33.05 ± 0.47 | 4.76 ± 0.59 s⁻¹, 8.1 σ |

The T₁ values are those of PC2. With the T₁s and T₁f errors fully correlated (the papers give no covariance) the gaps stand at 5.4, 7.8 and 6.7 σ; over the four phase cycles they run from 4.0 to 5.4 s⁻¹; with the first paper's own T₂s values from its figure they are 3.7, 4.6 and 4.0 s⁻¹.

The split above shows where the gap sits. The part comparing T₂s with T₁m is +0.5, +0.5 and −0.4 s⁻¹, and with the authors' tabulated T₁m and both errors T₂s − T₁m = +0.9 ± 1.5, +0.5 ± 1.1 and −0.9 ± 0.8 ms, within the errors. The whole gap is the second part, 3.6, 4.9 and 5.2 s⁻¹, from a split of the two longitudinal rates of 12 to 17 s⁻¹. So the measured T₂s sits where one pool with a nearly mono-exponential T₁ would put it, while the fitted T₁ split asks for a T₂s 14 to 17 % shorter. One of the two does not belong to one pool.

## What can carry the conflict

The mechanisms that can add to the relaxation of one pool, as far as this page and the agarose page name them, move the gap the other way, toward a central coherence that decays faster than the mean of the longitudinal rates:

- **Rank-1 fields**, paramagnetic or dipolar to protons: to first order their parts at ω₀ add equally to the central rate and to the mean of the longitudinal rates, and their part at zero frequency adds to the central rate alone. Their cross-correlation with the quadrupolar part, which the sodium page names, is not argued here.
- **Exchange between sites of different central frequency**: in the fast limit it adds p(1 − p)·δ²/k to the central coherence alone, as the agarose page notes.
- **Slow motion beyond the Redfield picture**: the part of the coupling along T₂₀ commutes with the populations and with the central transition however slowly it moves, and the parts sampled at ω₀ and 2ω₀ stay perturbative while ω_Q ≪ ω₀; what is left is the agarose page's fluctuating second-order shift of the central transition, extra decay.
- **Unrefocused field inhomogeneity and diffusion** shorten T₂ as well.

A distribution of correlation times that every nucleus samples keeps the relation. The authors find that correlation times and coupling strengths computed from the T₁ pair and from the T₂ pair disagree, and read that as "a distribution of correlation times and quadrupole strength parameters". Their comparison mixes this relation with [the key's single-motion check](THE_NUCLEUS_AS_A_KEY.md#back-to-the-motion); a distribution breaks the second and cannot break the first. Over the measurement a sodium ion travels micrometres, far beyond the pores, as the agarose page puts it: in a T₁ of 40 ms about √(6Dt) ≈ 17 µm, or √(2Dt) ≈ 10 µm out of a slab (D ≈ 1.2·10⁻⁹ m²/s, recalled). So each ion samples the gel's microstructure, and the distribution is of the kind that keeps the relation. Three readings remain:

- **(i) The fitted split is not one pool's.** The amplitudes are fixed at their one-pool values, the split grows with the late points, and at 2 and 4 % the free fit found the single-quantum recovery close to mono-exponential, which is what T₂s ≈ T₁m asks for.
- **(ii) More than one pool** of nuclei that do not average during the relaxation, which needs regions beyond 10 to 17 µm with different relaxation, a heterogeneity of the gel on that scale.
- **(iii) The measurement.** With imperfect pulses a signal stored along z relaxes at longitudinal rates and pulls an apparent T₂s toward T₁, the sign seen here ([the agarose page](SODIUM_IN_AGAROSE_LAYER_BY_LAYER.md) names such pathways for its echo). The size does not fit. Added to one pool's transverse signal and fitted with Eq. 11, it would have to be 0.64, 0.56 and 0.36 of that signal if it relaxed at T₁s, 1.6 to 1.7 times it at T₁m, and would lift the slow amplitude to 0.55 to 0.76 against the measured 0.40 to 0.41. Stored along z for half the evolution, even the whole slow signal would decay at 42.9, 37.3 and 32.6 ms, 1.5, 2.6 and 1.0 σ short of the measured T₂s.

## Free and bound sodium

The simplest picture of more than one pool is free sodium relaxing like saline, beside one bound pool. The script builds it from the spin-3/2 operators with a pool label beside each nucleus: the free pool in extreme narrowing at T₁ = T₂ = 54.19 ms, the bound pool with J₀ ≥ J₁ ≥ J₂ ≥ 0, exchange with detailed balance. It reads the signals the sequences record, the rank-1 and rank-3 parts after an inversion and after an excitation, and puts them through the authors' fits: Eq. 13 of the first paper for T₁, with the sign between the two longitudinal terms that the physics requires (the printed equation shows a minus) and the triple-quantum part weighted at the reported ratio, and Eq. 11 of the second paper for T₂, on 26 evenly spaced evolution times (their spacing is not stated; geometric spacing is shown beside). Six checks hold: the jump operators, with a Hamiltonian, against the repository's gated sodium generator ([`sodium_two_mirrors.py`](../simulations/sodium_two_mirrors.py)) to 2·10⁻¹⁵, and the voxel's pool block against it to the same; the voxel's signals without exchange against Hubbard's for each pool to 9·10⁻¹⁶; the central-transition cells under exchange, exactly decoupled from the rest, against the closed two-site exchange matrix, to the last digit; a single pool, which comes back from the two fits unchanged; and the fit of Eq. 13 started from swapped times, which still finds them.

- **One pool alone** gives χ² = 26.5, 69.6 and 62.3 for the five independent numbers.
- **Free and bound, without exchange** (four parameters, one degree of freedom): χ² = 15.8, 39.5 and 11.7; with geometric spacing 16.5, 46.0 and 20.2, which makes 2 % the weakest instead of 6 %; halving or doubling the triple-quantum weight changes them by less than one. The free fraction comes out at 4 to 7 %, the model's T₂s stays 3.8, 4.8 and 1.7 σ short and its slow amplitude 1.1, 3.7 and 2.9 σ high. The verdict rests on the free pool's time: at 58 ms the χ² values fall to 12.9, 30.0 and 6.9, at 65 ms to 7.4, 16.5 and 2.4, while the agarose page's mobile populations relax faster than saline, not slower; with the free time fitted under its measured 54.19 ± 0.95 ms they are 15.7, 38.0 and 11.2. Much of the χ² at 4 and 6 % is the slow amplitude, which the paper gives to 0.01 and which pulses and readout move most, and correlations within one fit, which the papers do not give, can move these values substantially either way.
- **With exchange**, profiled over rates from 1 to 1000 s⁻¹: χ² rises with the rate toward the one-pool values, since fast exchange makes one averaged pool.

Under the fit errors taken as independent, and with free sodium like the 0 % sample, the free-and-bound picture fails in all three samples.

## Two bound pools

Two bound pools without exchange, each with J₀ ≥ J₁ ≥ J₂, have seven parameters for five numbers, so the fit does not determine them. What it does constrain is the longitudinal triple-quantum signal, written TQ₁ here (B1 in the script's output) and measured against what one pool with the same fitted times gives (TQ₁ = 1 for one pool). A spread of longitudinal rates across nuclei that are each close to mono-exponential makes little rank-3 longitudinal order. The script searches with TQ₁ held at or above a floor, by a soft penalty and from three starts:

| agarose | TQ₁ free | TQ₁ ≥ 0.25 | TQ₁ ≥ 0.5 | TQ₁ ≥ 0.75 |
|---|---|---|---|---|
| 2 % | χ² 0.3 at TQ₁ 0.001 | 2.1 | 6.4 | 14.0 |
| 4 % | 0.0 at TQ₁ 0.07 | 0.0 | 4.1 | 24.7 |
| 6 % | 0.0 at TQ₁ 0.39 | 0.0 | 0.5 | 11.0 |

The lowest χ² at 2 % comes with both pools mono-exponential and no longitudinal triple quantum at all, which the measured peak excludes if that peak is the pools' own rank-3 order. Within χ² of about 2, two pools allow TQ₁ up to about 0.23 at 2 %, 0.37 at 4 % and 0.54 at 6 %: at most between a quarter and a half of one pool's. The transverse triple quantum of these solutions stays at 0.99 to 1.00 of one pool's, in line with the second paper's agreement between its measured and its theoretical transverse signal. A typical solution at 2 %: 91 % of the sodium with longitudinal times 48.5 and 44.7 ms, 9 % with 28.0 and 16.1 ms. The fit with its amplitudes fixed to 0.8 and 0.2 reads such a pair as T₁s and T₁f of one pool; each pool obeys the relation, and the gap is the reading. The table gives values found by a local search; nudged starts move them by a few units of χ², and the band edges are soft.

## The longitudinal triple quantum across the samples

The two papers already constrain TQ₁ relatively. For one pool the ratio of the longitudinal to the transverse triple-quantum amplitude, each over its own single-quantum amplitude, is the same in every sample: the relaxation times enter only the fitted shapes, and the pulses after the evolution and the acquisition are shared by the two pathways, assuming the pulses act alike in the three samples. Measured (the first paper's Table 2 and Fig. 6C):

| | 2 % | 4 % | 6 % | 6 % over 2 % |
|---|---|---|---|---|
| ratio, PC1 | 0.134 | 0.222 | 0.339 | 2.5 ± 0.7 |
| ratio, PC2 | 0.230 | 0.321 | 0.470 | 2.0 ± 0.3 |

One pool with the fitted split would keep the ratio flat; it rises by a factor of 2 to 2.5, 2.4 and 3.9 σ from flat. Relative to 6 %, the measured pattern 0.39 to 0.49 : 0.65 to 0.68 : 1 lies close to the pattern of the two-pool band edges, 0.43 : 0.69 : 1. This disfavors one pool with the fitted split, which reading (iii) would need; it does not separate readings (i) and (ii), and the absolute scale needs the authors' sequences simulated with their pulses, phase cycles and acquisition, since the phase cycle itself changes the ratio (PC2 against PC1).

## What this is, and what it is not

Not new: Hubbard's rates, from which the relation is one line, and which the authors write themselves; the relation as a model check (the sodium page); the relation applied to published agarose values as a floor, and used there to exclude a motional explanation; diffusion over the measurement far beyond the pores; the second-order shift of slow motion as extra decay (the agarose page); triple-quantum filtering as the way to tell intrinsic bi-exponential relaxation from a superposition of compartments (Jaccard, Wimperis and Bodenhausen, 1986, cited on the sodium page and by the first paper), the sodium page's triple-quantum path, of which the longitudinal argument here is the twin, and the second paper's own modelling of two compartments with their triple-quantum signals; sodium in tissue and gels as more than one pool; free and bound sodium as a model. What this page adds: the relation tested as an equality on three independently reported rates, and the split that places the whole gap in the fitted longitudinal split; the argument that slow motion beyond Redfield, rank-1 fields and the other mechanisms inside one pool named here all move the gap the other way; the size a pathway along z would need; free and bound sodium fitted through the authors' own procedures, where it fails; two bound pools reproducing the published numbers through the fixed amplitudes, with the ceiling on the longitudinal triple quantum they allow; and the triple-quantum ratio across the samples, which disfavors one pool with the fitted split. The test uses fitted numbers, not the measured curves.

## Open

- The absolute size of the longitudinal triple-quantum signal against one pool's: the sequences of both papers simulated with their pulses, phase cycles and acquisition, for one pool with a nearly mono-exponential T₁ and for the two-pool solutions, against the reported amplitudes.
- The measured curves: with them the pool models could be fitted to the data instead of to fitted numbers; the authors have them.
- The authors' values at 21.1 T (the first paper's Fig. 6B); whether this data set, with its two fields, serves as the agarose page's second agarose data set; other data sets reporting T₁s, T₁f and T₂s for the same samples; a direct measurement of the triple-quantum decay, which the model ties to the same slow transverse rate.
- A free pool unlike saline; exchange between two bound pools; the cross-correlation of the quadrupolar part with rank-1 mechanisms.

## Links

- [The Sodium Nucleus and Its Two Mirrors](THE_SODIUM_NUCLEUS_TWO_MIRRORS.md), the relation and why it survives a splitting.
- [The Nucleus as a Key](THE_NUCLEUS_AS_A_KEY.md), the relation as one of five checks, and the single-motion check.
- [Sodium in Agarose, Layer by Layer](SODIUM_IN_AGAROSE_LAYER_BY_LAYER.md), the same relation on another agarose data set, the floor, the pores and the pathways stored along z.
- Reichert S, Schepkin V, Kleimaier D, Zöllner FG, Schad LR. Comparison of triple quantum (TQ) TPPI and inversion recovery TQ TPPI pulse sequences at 9.4 and 21.1 T. *NMR in Biomedicine* 37 (2024) e5106, doi:10.1002/nbm.5106 (open access).
- Reichert S, Schepkin V, Kleimaier D, Zöllner FG, Schad LR. Sodium triple quantum MR signal extraction using a single-pulse sequence with single quantum time efficiency. *Magnetic Resonance in Medicine* (2024), doi:10.1002/mrm.30107 (open access).
