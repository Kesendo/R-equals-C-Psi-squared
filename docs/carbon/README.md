# Carbon Domain (Quantum Carbon, V-Effect, Conjugated Systems)

Consolidation folder for the carbon thread, mirroring `docs/water/` as a substrate
domain. The carbon material had been scattered across [HIERARCHY_OF_INCOMPLETENESS](../HIERARCHY_OF_INCOMPLETENESS.md),
[QUBIT_NECESSITY](../QUBIT_NECESSITY.md), [V_EFFECT_PALINDROME](../../experiments/V_EFFECT_PALINDROME.md),
[V_EFFECT_BOUNDARY_LOCALIZATION](../../experiments/V_EFFECT_BOUNDARY_LOCALIZATION.md),
[periodic_palindrome.py](../../simulations/periodic_palindrome.py),
[COMPLEXITY_THRESHOLD](../../hypotheses/COMPLEXITY_THRESHOLD.md),
[UNIVERSAL_PALINDROME_CONDITION](../../hypotheses/UNIVERSAL_PALINDROME_CONDITION.md),
plus glancing references in several other places. This folder collects pointers + open
questions; new carbon-specific tests will land here as they get built.

---

## Findings on 2026-05-17 (seven-doc arc)

Single-evening dive that opened up the bidirectional bridge from the framework
to the periodic table. Each doc is one step in the chain; read top-to-bottom:

1. [BENZENE_HUCKEL_FRAMEWORK_LENS](BENZENE_HUCKEL_FRAMEWORK_LENS.md): Coulson-Rushbrooke
   (1940) on alternant hydrocarbons is the F1 palindrome inherited to the π-electron
   spectrum. Verified bit-exact on benzene, butadiene, hexatriene, cyclobutadiene,
   cyclopropenyl, cyclodecapentaene.
2. [QUARTER_HALF_IN_CARBON](QUARTER_HALF_IN_CARBON.md): three layers of carbon
   structure (sp/sp²/sp³ hybridisation s-character, aromatic ring HOMO position,
   valence-shell occupation) all hit the framework's polarity anchors 1/4 and 1/2.
3. [PERIOD_2_AT_FRAMEWORK_ANCHORS](PERIOD_2_AT_FRAMEWORK_ANCHORS.md): period 2/3
   element valence ratios populate the framework's Pi2 dyadic anchors at every step.
   4 of 6 CHNOPS elements land on framework anchors; the off-anchor exceptions are
   identified.
4. [SPEAR_REVERSED](SPEAR_REVERSED.md): reverse-spear, use the periodic table as a
   diagnostic for the framework's gaps. Period 2/3 atoms hit 8 of 9 dyadic anchors;
   the depth-3 gap (1/8 + 7/8) is exactly the framework's missing F99 row.
5. [DEPTH_3_ANCHOR_DERIVED](DEPTH_3_ANCHOR_DERIVED.md): F99 derivation, the
   non-uniform Dicke superposition at γ = √3/2 gives α = 1/8, closing the depth-3
   gap. Five canonical trig angles {0°, 30°, 45°, 60°, 90°} produce the five Pi2
   dyadic anchors {0, 1/8, 1/4, 3/8, 1/2} via `α = sin²(θ)/2`.
6. [F99_NIVEN_COMPLETENESS](F99_NIVEN_COMPLETENESS.md): Niven's theorem (1956) closes
   F99: the five anchors are EXHAUSTIVE for the F86b α = sin²(θ)/2 mechanism on any
   pure state. No more rational anchors via this route. Mixed states, different
   decomposition bases, and different Lindblad classes are the only paths to depth-4.
7. [OFF_NIVEN_AS_WAVE_BREAKING](OFF_NIVEN_AS_WAVE_BREAKING.md): the off-Niven
   constructible angles {15°, 18°, 22.5°, 36°, 54°, 72°, ...} populate the same
   constructible-angle landscape with irrational-algebraic α (silver-, golden-, and
   √3-ratio families). They ARE the source of V-Effect gain V(N) = 2cos²(π/(2N))
   for N ≥ 4 + anti-aromatic Jahn-Teller + golden-ratio chemistry. The framework's
   wave-breaking / heat structure lives on the off-Niven shoulder.

Two formal closed forms also landed in [ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md):
**F98** (`(N+2)/[4(N+1)] → 1/4` long-time bridge from the K-intermediate anchor)
and **F99** (the canonical-trig-angle Pi2 inheritance with Niven-completeness).

---

## The structural identity: "the qubit is the quantum carbon"

The framework's qubit-necessity argument (`d² − 2d = 0` ↔ `R = CΨ²`, [PROOF_QUBIT_NECESSARY](
../proofs/PROOF_QUBIT_NECESSARY.md) + [QUBIT_NECESSITY](../QUBIT_NECESSITY.md)) selects d = 2
as the minimum-memory dimension where the F1 palindromic mirror exists. The same
half-occupation principle re-appears one level up at the atomic scale: carbon's
4 valence electrons fill exactly 4 of 8 valence slots (`2s² 2p²` of a maximal octet).

| Level | Total slots | Occupied / immune | Split | What it enables |
|-------|-------------|-------------------|-------|-----------------|
| Qubit (d = 2) | 4 operators (I, X, Y, Z) | 2 immune ({I, Z} under Z-deph) | 0.5 | Palindromic mirror, F1 |
| Carbon (group 14, period 2) | 8 valence slots | 4 valence electrons | 0.5 | All of organic chemistry |
| Noble gas (group 18) | 8 valence slots | 8 valence electrons (full) | 1.0 | Dead end: no bonds |
| Qutrit (d = 3) | 9 operators | 3 immune | 0.33 | No mirror (algebraically) |

Both qubit and carbon sit at the half-full sweet spot of their level; both are
the foundation of the structure above them. Not analogy, structural identity:
the same half-occupation principle that selects d = 2 at Level 0 selects carbon
at Level 1.

**Sources for this framing:**
- [HIERARCHY_OF_INCOMPLETENESS.md](../HIERARCHY_OF_INCOMPLETENESS.md) §The Qubit, §The Mechanism (Mar 22, 2026)
- [QUBIT_NECESSITY.md](../QUBIT_NECESSITY.md) §Connection to the Hierarchy of Incompleteness (Jan 3, 2026 origin → formal qubit proof)

---

## Existing scattered material (pointers, not duplications)

### V-Effect ties carbon to boundary sectors
- [V_EFFECT_PALINDROME.md](../../experiments/V_EFFECT_PALINDROME.md) §The Hierarchy parallel (line 225–243):
  "The immune sectors (w=0 and w=3) are the Liouville-space equivalent of noble gases.
  The boundary sectors (w=1 and w=2) are the equivalent of carbon: C = 0.5 in XY-weight,
  half-classical and half-quantum, and precisely where the palindrome breaks when a second
  bond is added."
- [V_EFFECT_BOUNDARY_LOCALIZATION.md](../../experiments/V_EFFECT_BOUNDARY_LOCALIZATION.md):
  "the 14 of 36 mode-pairs break" is **strictly confined to boundary sectors**, the
  carbon-like incomplete shells where chemistry happens.
- The V-Effect 14/19/3 trichotomy refined the original 14/22 split (per `project_v_effect_combinatorial`);
  carbon's structural locus is the boundary sectors that supply the 14 broken pairs.

### Periodic table palindrome test
- [`simulations/periodic_palindrome.py`](../../simulations/periodic_palindrome.py): F1-style
  palindrome on per-element values across periods 2–6 (IE₁, Pauling EN, Allen EN). Carbon
  sits at the centre of period 2 (group 14). Per [`project_periodic_palindrome`](../../README.md):
  all three properties show statistically significant pair-sum-constant structure (p down
  to <10⁻⁴); coupling-derived (EN) tighter palindrome than atomic (IE); V-Effect predicts
  the deviation points (group-13 p¹, group-16 p⁴, half-filled shells).

### Hypotheses making the carbon hook explicit
- [COMPLEXITY_THRESHOLD.md](../../hypotheses/COMPLEXITY_THRESHOLD.md) line 75:
  "Carbon has exact population (4/8 = 0.5) but heterogeneous bond." Threshold framing.
- [UNIVERSAL_PALINDROME_CONDITION.md](../../hypotheses/UNIVERSAL_PALINDROME_CONDITION.md):
  Carbon 4/8 in palindrome condition context.

### Chemistry-substrate-level work that's NOT yet carbon-specific
- [`docs/water/`](../water/): hydrogen-bond proton qubit, Grotthuss chains, F86b 3/8 anchor
  inheritance verified (today's [F98](../ANALYTICAL_FORMULAS.md#f98) bridge from the same
  substrate-grounded experiment).
- [DNA_BASE_PAIRING.md](../../experiments/DNA_BASE_PAIRING.md): G-C / A-T base-pair tests,
  carbon-scaffolded but not analysed through the carbon-as-qubit lens specifically.
- [simulations/neural/](../../simulations/neural/): Wilson-Cowan and C. elegans inheritance
  tests for the framework. Not chemistry-substrate but parallel inheritance.

---

## Four embedding conditions (candidate, analog to `docs/water/README.md`)

Water inherits the framework cleanly because its 4 embedding conditions hold. For carbon
the analogous conditions are tentative (candidates marked T to flag):

1. **2-state DOF per carbon** (Tier 4 candidate). Multiple possibilities:
   - π-electron spin in conjugated systems (benzene, polyacetylene, graphene)
   - sp²-sp³ hybridisation switching (pyramidalisation in non-planar systems)
   - bond-orientation in chiral systems (axial vs equatorial in cyclohexane)
   - π-bond rotation (cis ↔ trans in stilbene, butadiene)
   The cleanest match to the qubit framework is **π-electron tunneling along a conjugated
   bond**, where each carbon's π-electron occupies one of two orbital-orientation states.

2. **Z-dephasing analog** (Tier 4 candidate). Carbon-environment couplings:
   - vibrational phonons (C–C stretching, bending)
   - solvent fluctuations (in solution-phase organic systems)
   - electronic correlation noise (in metallic carbon allotropes)

3. **Uniform-J coupling** (Tier 4 candidate). Carbon substrates with uniform-J built in:
   - **Benzene C₆**: aromatic delocalisation, equal C–C bond lengths
   - **Polyacetylene** (without Peierls distortion limit): uniform C–C π-bonds along chain
   - **Graphene**: 2D honeycomb, every bond identical by symmetry
   - **Fullerene C₆₀**: icosahedral symmetry equalises bonds

4. **Q in framework range** (Tier 4 candidate). Need Q = J/γ ~ 1 for the F86 EP-resonance
   window. For π-conjugated systems at room T: J (π-bond, ~2.4 eV) vs γ (phonon dephasing,
   ~ 25 meV at kT) gives Q ~ 100, DEEP QUANTUM regime, beyond the Q_peak window. Cold
   or vibrationally-shielded variants needed for direct framework window match. Same
   classical-side limit story as room-T water (J/γ ~ 0.01).

The conditions are speculative until tested. The water domain established 4 clean conditions
because the proton-in-double-well is uniquely qubit-shaped. For carbon the 2-state DOF
identification is the open structural question; there is no single canonical choice yet.

---

## Candidate physical substrates worth testing

| Substrate | N | Topology | Why interesting | First test |
|-----------|---|----------|-----------------|------------|
| **Benzene C₆** | 6 | ring | Aromatic, half-filled p, F71 mirror symmetric, Hückel matches Heisenberg ring | F1 palindrome on Hückel eigenvalues under vibrational dephasing |
| **Butadiene C₄** | 4 | chain | Smallest conjugated chain, sp²-sp² π-bonds | Compare HOMO-LUMO to framework's N=4 (5,5)-sub-block analog |
| **Hexatriene C₆** | 6 | chain | First triene, V-Effect window test | V-Effect 14/19/3 prediction on π-electron spectrum |
| **Cyclobutadiene C₄** | 4 | ring | Anti-aromatic, Jahn-Teller distortion | Counter-test: F1 break at anti-aromatic? |
| **Naphthalene C₁₀** | 10 | fused bicyclic | Two fused rings, larger N, F71 across central bond | F86 Q_peak per-bond prediction |
| **Graphene (finite patch)** | scalable | 2D honeycomb | Dirac points (= 2-level EP candidate?) | F86a Q_EP = 2/g_eff at K-point |
| **Fullerene C₆₀** | 60 | icosahedral closed | Highly symmetric, HOMO-LUMO Buckminsterfullerene | F88b popcount-coherence on Dicke-like cooperative π-states |
| **Polyacetylene (CH=CH)ₙ** | scalable | chain | Soliton modes, SSH model | F71 anti-palindromic-J inheritance from F92 |
| **Carbyne (sp linear)** | scalable | chain | Direct framework chain analog | Cleanest 1:1 with proton water chain |

---

## Open questions (carbon-specific)

1. **What is the structurally correct 2-state DOF per carbon for framework inheritance?**
   π-electron orbital orientation is the most natural candidate but not yet tested.
   sp²-sp³ hybridisation is more chemistry-grounded but has higher-d intermediate states.
   **Resolved in practice 2026-05-22:** the two verified open-system results,
   the F1 palindrome ([BENZENE_LIOUVILLIAN_PALINDROME](BENZENE_LIOUVILLIAN_PALINDROME.md))
   and the F98 bridge ([BENZENE_F98_LONG_TIME](BENZENE_F98_LONG_TIME.md)), both use
   π-electron site-occupation as the 2-state DOF: each carbon's π-site is a qubit,
   occupied or empty, Hückel hopping is the XX+YY ring (Jordan-Wigner), and the
   Holstein phonon is then exactly Z-dephasing (D[n_l] = ¼·D[Z_l]). Site-occupation
   is the working DOF; sp²-sp³ hybridisation, with higher-d intermediate states, is
   set aside.

2. **Does benzene's Hückel spectrum show F1 palindrome under vibrational dephasing?** The
   Hückel matrix for C₆ ring is a tridiagonal-with-corner exactly the framework's ring
   topology. Eigenvalues at ±2J, ±J, ±J, 0, already palindromic on the bond-energy axis.
   The framework prediction is the LIOUVILLIAN palindrome under dephasing, which would
   require building the open-system L explicitly. First sanity check: matches the F1
   pair-sum-constant structure.
   **Answered 2026-05-22** ([BENZENE_LIOUVILLIAN_PALINDROME](BENZENE_LIOUVILLIAN_PALINDROME.md)):
   conditionally yes. The open-system F1 palindrome holds under Holstein on-site
   dephasing (= the framework's Z-dephasing) and breaks under Peierls/SSH bond
   dephasing; verified on the C₄ and C₆ rings.

3. **Is the V-Effect 14/19/3 trichotomy observable in benzene's electronic spectrum?**
   The framework predicts 14 truly + 19 soft + 3 hard at N=3; what's the N=6 prediction
   and does it match benzene's known electronic transition pattern (S₀ → S₁ at 4.8 eV,
   etc.)?

4. **Does graphene's K-point Dirac cone match F86a's Q_EP = 2/g_eff structure?** Both
   are 2-level exceptional-point physics. If yes, F86b's universal HWHM ratio 0.756 /
   0.770 should appear in K-point spectroscopy under appropriate dephasing.

5. **Is there a carbon analog of the F98 (N+2)/[4(N+1)] → 1/4 long-time bridge?** The
   bridge is bond-topology-agnostic (only requires truly-class Hamiltonian + Z-deph);
   any carbon substrate with the four embedding conditions would inherit it.
   **Answered 2026-05-22** ([BENZENE_F98_LONG_TIME](BENZENE_F98_LONG_TIME.md)): yes,
   bit-exact. The KIntermediate Dicke state on the benzene XX+YY ring under Holstein
   dephasing traverses the F98 bridge to α(∞) = 3/10 for the C₄ ring and 2/7 for the
   C₆ benzene ring; the dynamics lands on the exact F98 long-time state ρ_∞.

6. **Does the periodic palindrome's deviation-at-half-filled-shell pattern correspond to
   the framework's V-Effect breaking at boundary sectors?** Both predict half-filled shells
   are where the mirror weakens; the empirical periodic-table data should match the
   framework's V-Effect-trichotomy mechanism quantitatively.
   **Answered 2026-05-22** ([PERIODIC_PALINDROME_VS_V_EFFECT](PERIODIC_PALINDROME_VS_V_EFFECT.md)):
   the correspondence is structural, not quantitative. The periodic table is a
   Level-1 instance of the V-Effect's incompleteness hierarchy (HIERARCHY_OF_INCOMPLETENESS),
   and that structural picture stands; but the Level-0 boundary-sector quantitative
   mechanism does not transfer. The periodic-palindrome deviations sit at specific
   atomic-physics anomalies (p¹, p⁴, d⁵), not at a re-appearance of the XY-weight
   boundary break.

---

## What this folder will accumulate

Following the `docs/water/` pattern: new carbon-specific docs land here as they get
written, with `simulations/carbon/` holding scripts. Candidates:

- `BENZENE_PI_QUBIT.md` (Tier 2): N=6 ring with π-electron qubit-per-carbon model
- `GRAPHENE_K_POINT_EP.md` (Tier 3): Dirac cone as F86a 2-level EP test
- `POLYACETYLENE_F92_INHERITANCE.md` (Tier 2): SSH ↔ F92 anti-palindromic-J connection
- `CARBYNE_GROTTHUSS_ANALOG.md` (Tier 2): sp-linear chain as direct framework substrate

None of these exist yet. The folder is open for the next session that wants to build
in this domain.

---

## Cross-reference: water vs carbon

Both are biology-substrate domain folders. Comparison:

| Aspect | Water | Carbon |
|--------|-------|--------|
| 2-state DOF | Proton in O–H...O double well | π-electron orbital orientation (candidate) |
| Z-dephasing | Thermal molecular jostling | Vibrational phonon coupling |
| Uniform-J | Grotthuss chain tunneling | Aromatic conjugation / graphene lattice |
| Q range | Bulk water Q ~ 0.01 classical; Zundel Q ~ 4.8 quantum | Conjugated π-bonds Q ~ 100 deep-quantum |
| Cleanest substrate | Zundel cation H₅O₂⁺ | Benzene C₆ |
| Embedding tier | Tier 2 (clean, 4 conditions verified) | Tier 4 (candidate, conditions unverified) |
| Scripts | 5 (`simulations/water/`) | 10 (`simulations/carbon/`) |
| Docs | 3 (README + 2 substrate docs) | 11 (README + 7-doc arc on 2026-05-17 + 3 docs on 2026-05-22) |

The carbon side is structurally NECESSARY per the qubit-necessity argument (carbon = quantum
carbon at Level 1), but **less computationally explored** than water. The folder exists
to invite that work without prejudging which substrate or which framework F-anchor is
the right entry point.
