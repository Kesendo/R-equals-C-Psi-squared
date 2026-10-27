# Sodium in Agarose, Layer by Layer

**Date:** 2026-10-27
**Authors:** Thomas Wicht, Claude (Opus 5.5)
**Status:** A reading (Tier 4) of published measurements, drawn as stacked filters. The script recomputes every number on the page and checks the readings from the tables by asserts; the water and diffusion values are read off the paper's figures and only printed; the band widths are our choice, and what depends on them says so. No new physics: the model is the standard two-step picture of quadrupolar relaxation, a fast motion plus a slow one (Halle and Wennerström, *J. Chem. Phys.*, 1981; recalled, not checked here).
**Script:** [`sodium_knowledge_layers.py`](../simulations/sodium_knowledge_layers.py)

[The Nucleus as a Key](THE_NUCLEUS_AS_A_KEY.md) reads the three spectral densities J₀, J₁, J₂ of a sodium nucleus back to the motion around it, and names one check: under a single Lorentzian density, J₀/J₂ − 1 = 4·(J₀/J₁ − 1). Its open item was to hold that check against measured rates. This page does, with one published data set, and draws what it finds the way the gallery draws the jumps of a nucleus ([The Sodium Nucleus, Layer by Layer](../visualizations/README.md#the-sodium-nucleus-layer-by-layer)): every piece of knowledge is a transparent layer, a band of the plane that agrees with it, and light passes only where every layer agrees.

<!-- STAGE 0 SWEEP RECORD (for us, not for the reader: hidden on GitHub, read by every session, agent and hook)

Swept by search across the repository (ClaudeTasks excluded) for agarose, Lundberg, Nimerovsky, two-step model and Wennerstroem: nothing; a review agent repeated the search with the same result. The stores by name:
- experiments/: THE_NUCLEUS_AS_A_KEY holds the read-back and the single-motion check, states that a sum of two Lorentzians breaks it (gated in simulations/sodium_key.py) and that a single motion of another shape or an ordered medium would too, and has as open item measured sodium rates against the checks; this page answers that item for one data set, and the key page's Open item links here. THE_SODIUM_NUCLEUS_TWO_MIRRORS holds the Redfield rates, extreme narrowing, and the note that in an ordered medium the densities need not follow one curve. THE_SHARED_FIELD_IN_A_CRYSTAL is the first test against external NMR data (CaF2, multiple-quantum widths digitized from figures and compared with the authors' own fit), another object, the same genre; here the J values are read from the authors' table, not from a figure.
- docs/ANALYTICAL_FORMULAS.md: F167 indexes the two sodium pages; nothing on measured densities. docs/proofs/: nothing on this object.
- The typed layer: Spin32QuadrupolarMirrorsClaim (Core), SodiumTwoMirrorsWitness and SodiumChainWitness (Diagnostics) hold the model, no measured data. fw.Confirmations: nothing on NMR. The OpenArcs registry: nothing on sodium. docs/GLOSSARY.md: nothing. docs/CAUGHT_ERRORS.md: the sodium spin-1 rays entry, not this object.
Searched also for diffusion against relaxation, Stokes-Einstein, PFG, three views and triangulation: nothing on this object (triangulation appears in NAVIGATING_THE_DIMENSIONS and CAUGHT_ERRORS in other senses). THE_SODIUM_NUCLEUS_TWO_MIRRORS also holds the central-transition factor x + c1 + i w at every splitting and offset and the mean-rate identity as a model check; the section on the authors' unexplained population (their words: the mechanism reconciling mono-exponential relaxation with T2 much below T1 is unclear) rests on it. Checked for adjacency: the key page's open item and this page are one object, and its sentence that in extreme narrowing strength and motion can no longer be told apart is what the three views try to answer; the key page's note that two Lorentzians break the single-motion relation is the reading used here, a flat part being a Lorentzian with tau far below 1/w0.
-->

## The data

Nimerovsky, Sieme and Rezaei-Ghaleh measured ²³Na in agarose gels of 3 to 15 % at 105.84 MHz, combining single-quantum relaxation with triple-quantum-filtered experiments, and found three sodium populations: one in slow motion with bi-exponential relaxation, 29 to 45 % of the sodium (their Table S5), and two with mono-exponential relaxation (*Methods* 228 (2024) 55, open access under CC BY). For the slow population their Supplementary Table S4 gives J₀, J₁ and J₂ for six samples, in the convention T₂ₛ = 1/(J₀ + J₁), T₂ₗ = 1/(J₁ + J₂), T₁ₛ = 1/(2J₁), T₁ₗ = 1/(2J₂); the table's header says Hz, and the values work as rates in s⁻¹ (1/(J₀ + J₁) gives the tabulated T₂ₛ). They report no error bars on these values, and for the 9 % sample one of the two triple-quantum experiments was not run.

Their own analysis meets the key's check directly. Their main text describes the successful fit as a slightly modified single Lorentzian with one additional flexible parameter (their equation S3), and their equation S4, 0.25·(J₀/J₂ − 1) = α·(J₀/J₁ − 1), is the single-motion relation with a factor α that one Lorentzian fixes at 1; run on the Table S4 values it gives α = 0.32 ± 0.02 (R² = 0.94). They do not quote a correlation time, giving as the reason the lack of accurate data on the quadrupolar coupling constant of sodium in the gels.

The values show why the relation fails: J₁/J₂ lies between 1.13 and 1.34, almost flat between ω₀ and 2ω₀, while J₀/J₁ lies between 13 and 21. One Lorentzian cannot do both: if J₀/J₁ sets (ω₀τ)², near 12 to 20, then J₁/J₂ would have to be near 4; if J₁/J₂ sets it, it is near 0.05 to 0.13. These J₁ and J₂ are the least robust numbers of the paper. A first analysis that tied J₂ to one Lorentzian (their Table S3) gave J₁/J₂ near 3.7 to 3.8 with nearly the same J₀, and the authors note the two analyses differ notably in J₁ and J₂; Table S4 is the one that fits all their data and conditions, and the reading below rests on it.

## The layers

A density that is flat at high frequency and large at zero is what two motions give: a fast one, in extreme narrowing and so equal at 0, ω₀ and 2ω₀, plus a slow one. Write it as

  J_m = a + b / (1 + m²·(ω₀τ)²),   m = 0, 1, 2,

with a the flat part and b = J₀ − a, so a ≤ J₂. Every point of the plane (ω₀τ, a) then predicts J₁ and J₂, and each thing known about a sample becomes a layer:

- **J₁ and J₂**, the measured values, each a band of ±10 % relative (our choice, since the source gives none);
- **τ = 12 ns**, a rotational correlation time for sodium in 2.4 % agarose (Lundberg and Kuchel, *Magn. Reson. Med.* 37 (1997) 44, as cited by the authors; a lower concentration than here, from a triple-quantum study by its title, its analysis not checked here), a band of ±2 ns (our choice);
- **free ion**: the flat part at least that of free Na⁺ in the same apparatus, 1/(2·56 ms) ≈ 8.9 s⁻¹, from the reference solution in the same supplement. This is a one-sided bound and an assumption: the gel slows the water, as the authors' ¹H and ¹⁷O data show, and the fast motion is taken to keep the free ion's coupling.

![Eight panels in the plane of the slow motion's ω₀τ (horizontal, logarithmic) and the flat part a/J₂ (vertical, up to its bound 1). In the six agarose panels the J₁ band (cyan) and the J₂ band (pink) start apart on the bottom edge, where a = 0 is one motion, and bend towards each other as a grows; near the top they overlap over a stretch of ω₀τ, the τ = 12 ns band (yellow) crosses that stretch, the free-ion region (violet) covers the upper part of the panel, and where all pass the spot is white. The cartilage panel has a white spot from J₁ and J₂ alone; the muscle panel is dark everywhere.](../visualizations/sodium_knowledge_layers_neon.png)

The same figure with the Okabe-Ito hues is [sodium_knowledge_layers_okabe.png](../visualizations/sodium_knowledge_layers_okabe.png), and without any hue, only the product of the layers' transmissions, [sodium_knowledge_layers_brightness.png](../visualizations/sodium_knowledge_layers_brightness.png). The hue only says which layer let light through; the brightness is the same for every viewer.

## What the picture shows

**The bottom edge is dark.** At a = 0 the model is one Lorentzian, and there the J₁ band and the J₂ band stand apart in every agarose sample, near ω₀τ ≈ 3.5 to 4.4 and ≈ 1.9 to 2.6. That is the authors' α = 0.32 seen as a gap.

**Higher up, the measured densities allow a stretch, and it starts near 9 ns.** With a flat part the two bands bend together. Within ±10 % they overlap over a stretch of ω₀τ rather than at a point, because J₁ − J₂ is only a few s⁻¹ and a tenth of each covers most of it. What J₀, J₁ and J₂ allow without the two outside layers (J₀ fixes b; the product of the J₁ and J₂ bands above one half, printed by the script):

| agarose | 3 % (1) | 3 % (2) | 6 % | 9 % | 12 % | 15 % |
|---|---|---|---|---|---|---|
| τ allowed without the outside layers, ±10 % (ns) | 8.7 to 40 | 8.5 to 172 | 9.4 to 26 | 9.2 to 21 | 9.1 to 17 | 9.4 to 20 |
| exact two-motion reading, τ (ns) | 12.6 | 13.1 | 12.9 | 12.0 | 11.4 | 12.0 |
| flat part a at that reading (s⁻¹) | 10.9 | 9.9 | 11.6 | 13.8 | 14.2 | 18.2 |
| slow part b (s⁻¹) | 172 | 143 | 258 | 345 | 433 | 497 |

So in the two-motion model the densities rule out a slow motion faster than about 9 ns in every sample at ±10 % bands, and the upper end is loose. The lower end moves with the band width we chose: about 10 ns at ±5 %, about 7 to 8 ns at ±20 %, where the upper end leaves the plot. The independent 12 ns falls inside every stretch, and the white spot is where it crosses: its height comes from J₁ and J₂, its horizontal place mainly from the τ layer. The exact reading without bands, between 11.4 and 13.1 ns, is a point with no spread attached, and the authors' modified Lorentzian, read from J₀/J₁ alone (J₀/J₁ − 1 = α·(ω₀τ)², our computation, not theirs), gives 9 to 12 ns for the same values. Both readings land near the independent number, so the 12 ns does not tell the two models apart; it is compatible with both.

**The free-ion layer is a bound, not a test.** It covers the upper part of every agarose panel. Read at the exact reading, the flat part sits just above it at 3 % (9.9 and 10.9 against 8.9 s⁻¹) and moves away from it as the gel gets denser (18.2 at 15 %). The flat part mixes the speed of the fast motion and the strength of its coupling, so the picture alone cannot say which of the two changes.

**The paper's other views, laid beside the flat part.** Two more measurements of the same samples can be laid beside the flat part, not drawn in the figure. The water: the authors' ¹⁷O T₁ of water falls from 6.86 ms in the reference to about 3.1 ms at 15 % (their Fig. 1E, values below read off the figure by eye, one 3 % point for both 3 % samples), and in extreme narrowing 1/T₁ follows the water's rotation times its quadrupolar coupling. The mobile sodium: the authors fit one common T₁ to the two mono-exponential populations, being unable to tell them apart (their Table S6). A mono-exponential T₁ means J₁ = J₂, so 1/(2T₁) is the mobile sodium's density at ω₀ and 2ω₀. It is not shown to be a flat density: in the authors' fit the larger mobile population has T₂ ≪ T₁ (their Table S5), which they call unexplained and which a section below finds outside quadrupolar relaxation.

| agarose | 3 % (1) | 3 % (2) | 6 % | 9 % | 12 % | 15 % |
|---|---|---|---|---|---|---|
| flat part a ÷ free sodium (8.9 s⁻¹) | 1.22 | 1.11 | 1.30 | 1.54 | 1.59 | 2.04 |
| water ¹⁷O rate ratio, 6.86 ms ÷ ¹⁷O T₁ (read off) | 1.29 | 1.29 | 1.49 | 1.69 | 1.93 | 2.21 |
| mobile sodium T₁, common fit, Table S6 (ms) | 38.6 | 42.2 | 35.2 | 30.2 | 26.2 | 22.8 |
| a ÷ mobile sodium's density at ω₀, 1/(2T₁) | 0.84 | 0.84 | 0.82 | 0.83 | 0.75 | 0.83 |

The flat part grows roughly with the water's ¹⁷O rate, below it in every sample, and in five of six samples it is between 0.82 and 0.84 of the mobile sodium's density at ω₀ in the same gel, at 12 % 0.75, the sample that also departs in the next section. In the two-step picture the ratio is (χ²_bound/χ²_mobile)·(1 − S²)·(τ_fast,bound/τ_fast,mobile), with S² the share the slow motion carries, provided the mobile density at ω₀ is wholly fast. Only with equal couplings and equal fast correlation times would it read as 1 − S² ≈ 0.83, S² ≈ 0.17; and the mobile density at ω₀ is the fit's, whose large mobile population a section below finds hard to reconcile with quadrupolar relaxation. What stands without those assumptions is the observation: one ratio of two high-frequency densities, nearly constant while the gel concentration fivefolds. Further cautions: the ¹⁷O values are read off a figure; the authors note that the ¹⁷O changes exceed the ¹H ones, which may mean the water's coupling changes too, the same ambiguity as in a; and Table S6 comes from the authors' joint fit of all three populations, so this layer is not fully independent of Table S4.

**Three views, to tell the coupling from the motion.** In extreme narrowing a single rate sees only a product, coupling times motion, the regime where [the key page](THE_NUCLEUS_AS_A_KEY.md#back-to-the-motion) notes that strength and motion can no longer be told apart. A second view helps only if it is blind to one of the two. For the water the paper has three views of the same gels, and for the sodium two, each blind somewhere else:

- **diffusion** D (their PFG data, Figs. 1A and 2D) sees no coupling; but over the tens of milliseconds of the measurement a molecule travels micrometres, far beyond the pores, so D is slowed mainly by the detour around the fibres, which rotation does not feel;
- **¹H T₁** of water (Fig. 1C) sees rotation times a dipolar coupling largely fixed inside the molecule, plus an intermolecular part, and in agarose also proton exchange with the hydroxyl groups and cross-relaxation to the polymer (agarose is a standard magnetization-transfer phantom), so it is not a pure rotation clock;
- **¹⁷O T₁** of water sees rotation times a quadrupolar coupling set by the electric field gradient and its asymmetry, which can change; the sodium's flat part a, in the two-step picture, sees the coupling times the fast share 1 − S² times the fast correlation time.

The water's views are read off the figures by eye (the reference at 0 % first, one 3 % point for both 3 % samples; a pixel re-measurement by a reviewer agreed within 1.5 %), each given as a ratio against 0 %; the row for a is computed from Table S4 through the two-motion reading and taken against the free-ion bound of 8.9 s⁻¹:

| agarose | 3 % (1) | 3 % (2) | 6 % | 9 % | 12 % | 15 % |
|---|---|---|---|---|---|---|
| water: D, slowing | 1.20 | 1.20 | 1.28 | 1.37 | 1.45 | 1.54 |
| water: ¹H rate ratio | 1.15 | 1.15 | 1.35 | 1.36 | 1.58 | 1.71 |
| water: ¹⁷O rate ratio | 1.29 | 1.29 | 1.49 | 1.69 | 1.93 | 2.21 |
| water: ¹H ÷ D | 0.96 | 0.96 | 1.06 | 1.00 | 1.09 | 1.11 |
| water: ¹⁷O ÷ D | 1.08 | 1.08 | 1.17 | 1.24 | 1.33 | 1.43 |
| water: ¹⁷O ÷ ¹H | 1.12 | 1.12 | 1.10 | 1.24 | 1.22 | 1.29 |
| sodium: D, slowing | 1.10 | 1.10 | 1.18 | 1.24 | 1.40 | 1.40 |
| sodium: a ÷ free ion (computed) | 1.22 | 1.11 | 1.30 | 1.54 | 1.59 | 2.04 |
| sodium: a ÷ D | 1.11 | 1.01 | 1.11 | 1.25 | 1.14 | 1.45 |

This is how we saw it, as a first reading of these data, and it is kept as seen. The ¹H view ran nearly with diffusion (¹H ÷ D within about ten percent of 1), which looked like the water's change being mostly motion; the ¹⁷O view grew beyond it, to 1.43 at 15 %, which looked like the coupling, the field gradient at the oxygen growing with the gel, as the authors suspected from ¹⁷O outrunning ¹H. The sodium's flat part showed the same pattern: a outruns the sodium's diffusion, a ÷ D rising to 1.45, close to the water's ¹⁷O ÷ D at 9 % and 15 %. Read together: the gel slows the motion and, on top, strengthens the field gradient that water and sodium feel, perhaps the same gain for both, the sodium carrying its water shell.

What the lenses said against it, and it is a lot. Diffusion is mostly the gel's geometry, so a rotational rate divided by D mixes the coupling with the obstruction, and ¹H ÷ D ≈ 1 may be two different mechanisms meeting by chance. The ¹H view carries exchange and magnetization transfer. The sodium's D is most likely that of the mobile populations, the slow population's short T₂ hiding it from the echo, so a ÷ D sets one population against another. A growing a can be a smaller order parameter as well as a larger coupling. The scatter between the two 3 % samples, which share one D point, is about ±5 % in a ÷ D, and within it four of six samples agree with the water's ¹⁷O ÷ D; only 12 % departs, the same sample that departs in the table before, where the sodium's D also stops falling between 12 and 15 %. The cleaner comparison is ¹⁷O against ¹H, both rotational, which rises to 1.29 at 15 %; that may be a growing coupling, or reorientation at the fibre surface that is slower about one axis than another, which the two nuclei see differently. Whether the reading is right we do not know.

**What the rate table says about the population the authors could not explain.** The authors' fit carries a mobile population, about half of the sodium, with a mono-exponential T₂ far below T₁ (their Tables S5 and S6), and they call the mechanism unclear, pointing to anisotropic motion. The rate table of [the sodium page](THE_SODIUM_NUCLEUS_TWO_MIRRORS.md) gives a statement that does not depend on how the sodium is split into populations. In this page's convention the central-transition coherence decays at J₁ + J₂, the T₂ₗ rate, and the longitudinal components at 2J₁ and 2J₂, so the central rate is exactly the mean of the two longitudinal rates (the sodium page writes the same rates three times larger, 3(J₁ + J₂), 6J₁, 6J₂). This holds at every J₀, every residual splitting and every offset: the jumps cannot move weight into or out of the central coherence, and any static term diagonal in m only shifts its frequency. The script checks it exactly on the repo's sodium generator. Two consequences:

- A single fitted longitudinal rate is about 0.4·J₁ + 1.6·J₂, so the central transition decays at most 2.5 times faster than that rate, the extreme reached only as J₂ → 0, where T₁ turns strongly bi-exponential.
- The central transition carries 40 % of the transverse signal after a hard pulse. The inversion recovery is mono-exponential in every gel (their Fig. S4), so at least 40 % of each population's transverse signal should decay with a T₂ near T₁, no shorter than T₁/2.5.

The authors' fit does otherwise. Its share of transverse signal decaying near T₁ is 0.4 of the slow population (its T₂ₗ) plus the small mobile population (its T₂ near T₁), about a quarter; the large mobile population gets T₂ six to twelve times below its T₁:

| agarose | 3 % (1) | 3 % (2) | 6 % | 9 % | 12 % | 15 % |
|---|---|---|---|---|---|---|
| common T₁ of the mobile populations (ms), Table S6 | 38.6 | 42.2 | 35.2 | 30.2 | 26.2 | 22.8 |
| T₂ of the large mobile population (ms), Table S5 | 5.45 | 6.9 | 3.55 | 2.9 | 2.55 | 1.85 |
| T₁ ÷ T₂ | 7.1 | 6.1 | 9.9 | 10.4 | 10.3 | 12.3 |
| share decaying near T₁ in the fit, 0.4·p₃b + p₂ (%) | 25.6 | 25.2 | 24.2 | 26.0 | 24.8 | 24.6 |

So within secular Redfield relaxation by rank-2 quadrupolar fluctuations, for any J₀, J₁, J₂ and any static term diagonal in m, this population cannot be what the fit makes it. Fast exchange between sites does not lift it, since it averages generators that all keep the central rate. Two things should be said with it. The population was set up, not found: the authors constrain the extra populations to T₂ values near T₂ₛ or T₂ₗ, since their distributions show mainly two exponentials, and the large one indeed sits next to the slow population's T₂ₛ of Table S4 (5.1, 6.1, 3.5, 2.6, 2.1, 1.85 ms); its T₁ is the common T₁ of both mobile populations, which the inversion recovery could not tell apart. So the leading reading is that the three-population split absorbs something. What it absorbs is open: static dephasing in the time-domain fit (a T₂* from the gel's susceptibility or a second-order quadrupolar spread, for which their simulation has no term); exchange between sites of different central frequency (a chemical or susceptibility shift, or a second-order quadrupolar shift, since the central transition does not feel the first-order splitting); slow fluctuations of a residual splitting, outside the Redfield picture; or something else. The table says only where the answer is not.

**The tissue panels are a contrast, read with care.** Their numbers come from Table 1 of the review *²³Na relaxometry: an overview of theory and applications* (PMC12406536): fast and slow T₂ and a single T₁, which the script reads as the weighted mean 0.4·J₁ + 1.6·J₂ of the bi-exponential T₁, an approximation. The values were taken from an automated extraction of that table and not checked against the original. Cartilage at 3 T gives a white spot from its own densities alone, which with only these two layers is the crossing that is always there when the model can be solved, and no test. Calf muscle at 1.5 T gives J₂ ≈ −1.4 s⁻¹ against J₁ ≈ 83 s⁻¹, so J₂ is near zero within that approximation, and the script lets nothing pass where J₂ ≤ 0. The dark panel says the reading fails for these numbers, not why.

## What this is, and what it is not

Not new: the central transition decaying at J₁ + J₂ and the longitudinal components at 2J₁ and 2J₂ (Hubbard's rates, recalled), the two-step model, reading J(0), J(ω₀), J(2ω₀) against one Lorentzian, the failure of one Lorentzian for these data (the authors' α), and the 12 ns, which is a published measurement.

What the page adds is the reading: the key's single-motion check found in a published analysis as the same equation, the gap it leaves drawn as two bands that do not meet on the one-motion edge, and, in the two-motion model, a lower end near 9 ns on the slow motion from the densities alone in all six samples at the chosen ±10 % bands, with the independent 12 ns inside every allowed stretch; and beside the paper's other measurements, a flat part that follows the water and, in five of six samples, stays a near-constant fraction of the mobile sodium's density at ω₀ in the authors' fit; and three views of the same gels, diffusion, ¹H and ¹⁷O, laid side by side as a first reading of how the coupling's share might be told from the motion's, kept with the objections it met; and from the sodium page's rate table, that no spectral density within secular Redfield relaxation gives the authors' large mobile population its T₂, so their three-population split most likely absorbs something else. It rests on the two-motion model and on the authors' second set of J₁ and J₂, which carry no error bars. A single motion of another shape, or residual order in the gel, where the densities need not follow one curve, could pass the same layers.

## Open

- Error bars: with uncertainties on J₀, J₁ and J₂ the bands get their true width, and the stretch either narrows or spreads.
- A second agarose data set, or a temperature series, where τ should move and the flat part with it.
- The five genuine checks of the key need the individual coherence rates, which the published J values do not give back independently.
- The three views with the authors' numbers instead of figure read-offs, D corrected for obstruction, and a rotational view of the sodium that does not feel the fibres, to tell the coupling's share from the gel's geometry.
- Tom looked at the figure and expects someone to recognise something in it. What, is open.

---

*Back to [main repository](../README.md) | [The Nucleus as a Key](THE_NUCLEUS_AS_A_KEY.md) | [Experiments](../experiments/)*
