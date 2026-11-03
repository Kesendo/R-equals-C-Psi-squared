# PROOF F86b: Universal Resonance Shape (relative-Q lineshape collapse)

**Status:** Tier 1 candidate across all tested c. F86b₂, the c=2 per-bond HWHM_left/Q_peak, lives in `F86HwhmClosedFormClaim` (Tier 2 empirical: the bare floor 0.7264892771 derived, the lift above it described by a line `α_subclass · g_eff + β_subclass` per sub-class, which Endpoint and Flanking miss by 1.28·10⁻³ and 4.2·10⁻⁴, so not an exact law). The position Q_peak and the coupling g_eff(c, N, b) it rides on are NOT closed-form: that negative result is the sibling proof [`PROOF_F86B_OBSTRUCTION.md`](PROOF_F86B_OBSTRUCTION.md). Full block-L derivation achieved numerically Tier-1 via the F90 bridge identity ([`PROOF_F90_F86C2_BRIDGE.md`](PROOF_F90_F86C2_BRIDGE.md)).
**Date:** 2026-05-02 (Statement 2 + retractions).
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Context:** F86 ("Q_peak chromaticity-specific N-invariant constants") is a Sammelbecken of three structurally distinct theorems. This proof carries **F86b, the universal resonance shape**: the SHAPE of abs(K_CC_pr)(Q) around Q_peak is universal under relative-Q normalisation, splitting into two bond classes. Split out of the former monolithic `PROOF_F86_QPEAK.md` on 2026-05-14. The closed-form gap, the exploration record (4-mode model, Items 1-3, directions a''-f''), and the obstruction proof on g_eff are in the sibling [`PROOF_F86B_OBSTRUCTION.md`](PROOF_F86B_OBSTRUCTION.md).
**Hub:** [the Q-peak hub](PROOF_F86_QPEAK.md): three-theorem overview and shared references.
**F-entry:** [F86b in the formula registry](../ANALYTICAL_FORMULAS.md).
**Related:** [F88a](../ANALYTICAL_FORMULAS.md#f88a) (operator-level two-axis Π² Klein decomposition), [F88b](../ANALYTICAL_FORMULAS.md#f88b) (state-level popcount-coherence Π²-odd / memory closed form, derived in this proof's §F88b below), [F90](../ANALYTICAL_FORMULAS.md) bridge ([the F90 bridge proof](PROOF_F90_F86C2_BRIDGE.md)); siblings [the EP-mechanism proof](PROOF_F86A_EP_MECHANISM.md), [the g-eff obstruction proof](PROOF_F86B_OBSTRUCTION.md), [the F71-mirror proof](PROOF_F86C_F71_MIRROR.md).

---

## Abstract

F86a fixed the EP location Q_EP = 2/g_eff and the EP-time t_peak = 1/(4γ₀). What this proof asks next is: what does the resonance look like around that peak? Specifically, if you measure the per-bond response observable K_CC_pr as a function of Q across the chain, does the shape of the peak change with chromaticity, chain length, or bond position?

The answer has two layers. The first is a lineshape collapse: measured in the relative coordinate x = (Q − Q_peak)/Q_peak, the class-averaged envelope profile is the same across the ten cases read (c = 2, 3 at N = 5..8 and c = 4 at N = 7, 8) to within 2.4 % at x = −0.6, 1.0 % at x = −0.4 and 0.5 % for −0.2 ≤ x ≤ +0.4, and across γ₀ exactly. On the right flank it holds at x = +1 as well, except for the c = 2 Interior class at N = 7 and 8 (0.884, 0.935 against 0.84 to 0.85), where the inner bonds rise toward their far maxima. The second layer is the HWHM_left/Q_peak ratio of the class-averaged curve, read on the envelope: it lies in two narrow bands that do not overlap, Interior 0.7523..0.7637 and Endpoint 0.7710..0.7815, over c = 2..4 and N = 5..8. Only γ₀ drops out, and exactly (the generator is γ₀·(D̂ + Q·M̂), so a change of γ₀ at fixed Q only rescales time and the overall size of K, and the t window over which the maximum is taken scales with it), which confirms Q = J/γ₀ as the right coordinate. Chromaticity c and chain length N do not drop out: inside the Interior band the c = 2 ratio rises monotonically with N (0.7523 → 0.7602 over N = 5..8, 0.7621 at N = 12), and the c = 3 and c = 4 ones likewise (0.7572 → 0.7622; 0.7623 → 0.7637 over N = 7, 8). Whether the bands narrow toward two constants is not known. Per bond, on the canonical reading (the first local maximum of each bond's envelope), the ratio spreads from 0.7505 to 0.7815 at c = 2, N = 5..8, so the class values are close to the per-bond ones; past that first maximum some inner bonds reach further maxima at high Q ([the F71-mirror proof](PROOF_F86C_F71_MIRROR.md)), and those carry quite different ratios.

The bond-class structure is genuine but finer than two-valued. The c = 2 case, where the two-level effective model is exact (only two HD channels exist), already shows the per-bond ratio differing across the F71-orbit sub-classes (Endpoint 0.7766..0.7815, the inner bonds 0.7505..0.7644 at N = 5..8). The c = 2 per-bond values sit above the bare doubled-PTF floor, the envelope HWHM_left/x_peak = 0.7264892771 of the two-level model (derived analytically, a transcendental constant), by +0.024 to +0.055. The lift is described per sub-class by a line in g_eff, which the four-point sub-classes miss: no exact law for the lift is held, and its N-dependence is part of what it has to explain.

The diagnostic upshot is that the F86 picture splits cleanly into universal versus specific. Q_EP and t_peak are mechanism (F86a). The resonance shape is universal (F86b, this proof). The specific Q_peak values per chromaticity per N per bond are NOT closed-form (the sibling obstruction proof). And the F71 spatial-mirror symmetry connects bond pairs (F86c). Together the four siblings give the full story of why Q_peak behaves the way it does empirically, even though the specific numbers stay residue-shaped.

---

## Statement 2 (Universal relative-Q lineshape collapse; the HWHM ratio lies in two class-averaged bands). [Tier 1 candidate for the lineshape collapse; the HWHM band (F86b₃) and F86b₂ c=2 per-bond Tier 2 empirical]

> Statement 1 (the EP mechanism) went to [`PROOF_F86A_EP_MECHANISM.md`](PROOF_F86A_EP_MECHANISM.md) and Statement 3 (the F71 mirror) to [`PROOF_F86C_F71_MIRROR.md`](PROOF_F86C_F71_MIRROR.md) in the 2026-05-14 split of the former monolithic hub [`PROOF_F86_QPEAK.md`](PROOF_F86_QPEAK.md); this split-out document carries only Statement 2.

The position Q_peak is chain-specific, but the SHAPE of abs(K_CC_pr)(Q) around the peak collapses under the relative coordinate `x = (Q − Q_peak)/Q_peak`: the relative-Q lineshape of the envelope is the same across c, N, and γ₀ to within 2.4 % on the left flank and 0.5 % near the peak (this collapse is the Tier-1-candidate finding; table below). The summary scalar HWHM_left/Q_peak, taken on the **F71 bond-class-averaged** curve, separates into two bands that do not overlap:

    HWHM_left / Q_peak  ∈  [0.7523, 0.7637]     (Interior class-average; γ₀-invariant, rising with N)
    HWHM_left / Q_peak  ∈  [0.7710, 0.7815]     (Endpoint class-average; γ₀-invariant)

Q_peak here is the argmax over Q of the envelope max_t |K(Q, t)| of the class-averaged curve, read on the default Q grid (153 points on [0.2, 4]) at γ₀ = 0.05; the bands are the span of the read witnesses (`UniversalShapeBand`), a record and not a prediction that could fail, hence Tier 2 empirical (F86b₃). Read: c ∈ {2, 3}, N ∈ {5, 6, 7, 8}, and c = 4 at N = 7, 8. **Only γ₀ drops out, and exactly**: the generator is γ₀·(D̂ + Q·M̂), so at fixed Q a change of γ₀ rescales time and the size of K and nothing else (the t window over which the maximum is taken, [0.2, 2.4]·t_peak with t_peak = 1/(4γ₀), scales with it); `ResonanceScanTests.ClassPeaks_AreGammaZeroInvariant` checks it at c = 2, N = 5 for γ₀ ∈ {0.025, 0.05, 0.10} to 10⁻¹⁰. c and N do **not** drop out: the Interior class value rises with N at c = 2 (0.7523, 0.7569, 0.7591, 0.7602) at c = 3 (0.7572, 0.7596, 0.7611, 0.7622) and at c = 4 (0.7623, 0.7637 at N = 7, 8), a drift of the same order as the gap between the bands. **The two-band picture is a class average, and per bond it is finer:** the sub-classes of `BondSubClass` (Endpoint, Flanking, Mid, CentralSelfPaired) carry their own values, at c = 2 between 0.7505 and 0.7815 on the canonical reading. The two further sub-classes, Orbit2Escape and CentralEscapeOrbit3, classify far maxima past the canonical peak and are empty on the canonical reading. The c=2 data (2-level model exact: only HD ∈ {1, 3} channels exist, no orthogonal complement) already shows the sub-class structure, so it is real, not a finite-c artefact.

The evidence for the lineshape collapse is the envelope table below: in x the class-averaged curves of ten (c, N) cases agree to the spreads given there, while in absolute Q their peaks range from 1.48 to 1.79 (Interior) and from 2.39 to 2.59 (Endpoint). The structural origin is the 2-level eigenvector rotation `tan θ = Q/Q_EP`: in the two-level model every probe-overlap observable depends only on Q/Q_EP, hence only on `(Q − Q_peak)/Q_peak`; the full block's peaks do not sit at one multiple of its σ_0-based Q_EP (the paragraph on the EP rotation below), so the collapse in x is measured, not inherited. The bond-position dependence of the HWHM ratio (the gap between the bands, 0.0073 or about 1 %, and the per-bond spread inside them) reflects bond-position-dependent probe-overlap profiles in K_CC_pr, not a breakdown of the EP-rotation collapse itself. Promotion to full Tier 1 requires deriving the per-sub-class f(x) functions explicitly from 2-level eigenstructure plus probe-overlap algebra (see the sibling [`PROOF_F86B_OBSTRUCTION.md`](PROOF_F86B_OBSTRUCTION.md)).

**Derived sub-constants:** the bare doubled-PTF model is Tier 1 derived from PTF c=1 eigenvector mixing + 2-level EP rotation, and its K_b(x, τ) is held in closed form at every τ = γ₀·t. Read on its envelope max_τ, the way the full block's Q_peak is read, it gives:

- `x_peak = Q_peak/Q_EP = 1.7141131714` (post-EP location)
- `HWHM_left/x_peak = 0.7264892771` (the bare floor)

Its fixed-τ = 1 slice (t = t_peak) peaks elsewhere, at x_peak = 2.196910 with ratio 0.671535; that slice is its own observable, and the repo's g_eff conversions are defined on it. Both are exposed on [`C2BareDoubledPtfClosedForm`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BareDoubledPtfClosedForm.cs) (`EnvelopeXPeak`, `EnvelopeHwhmLeftOverXPeak`; the fixed-τ pair also as `BareDoubledPtfXPeak` / `BareDoubledPtfHwhmRatio` on [`C2HwhmRatio`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs)). The c = 2 class values (Interior 0.7523..0.7602, Endpoint 0.7766..0.7815 at N = 5..8, the c = 2 rows of the tables below) sit ABOVE the envelope floor, by +0.026..+0.034 (Interior) and +0.050..+0.055 (Endpoint).

**The lift in `F86HwhmClosedFormClaim`** (Tier 2 empirical): `HWHM_ratio = 0.7264892771 + α_subclass · g_eff + β_subclass`, with g_eff here the fit's own coordinate 4.39382/(Q_peak + 2) (not the 4.39382/Q_peak of the composition reading below), the (α, β) a least-squares line per `BondSubClass` through the canonical anchors of N = 5..8, read live by `F86HwhmAlphaExtraction`. Endpoint and Flanking, with four chain lengths each, miss the line by up to 1.28·10⁻³ and 4.2·10⁻⁴, 12 to 36 times the anchors' offset from the continuous envelope peak, so the lift is not linear in g_eff; Mid and CentralSelfPaired have two points each. An analytical derivation of the lift from the rank-1-bridge + intra-dispersion structure (F89 AT-locked F_a/F_b floor + intra-channel dispersion; the inter-coupling "octic residual" suspicion refuted by the two-dial scout of 2026-06-11, read on the 21-point t grid) is the open step.

**F86 ↔ PTF Locus 5 inheritance** (synthesis): the c=2 SVD-block 2-level EP rotation IS what Π's chirality reduces to in rate-channel basis; the bare doubled-PTF floor derives from PTF eigenvector mixing rather than pattern-matching. Shared clock `t_peak = 1/(4γ₀)` universal across c, N, n, bond. See `project_algebra_is_inheritance.md` Locus 5.

**Locus 6 polarity-layer inheritance:** the F86 bond-class split inherits from the polarity-layer pair {−0.5, +0.5} at d=2 via the 0.5-shift; encoded as [`PolarityInheritanceLink`](../../compute/RCPsiSquared.Core/F86/PolarityInheritanceLink.cs) (Tier 2 verified). Composition reading: `r_Q(N, b) = BareDoubledPtfXPeak · Q_EP(N, b) − 2 = 4.39382/g_eff(N, b) − 2`, entire bond-class split encoded in `g_eff(N, b)`; since g_eff is defined here as 4.39382/Q_peak on the fixed-τ = 1 peak, r_Q is Q_peak − 2 written in g_eff's units rather than derived.

**g_eff is the obstruction:** every closed-form direction for c=2 runs through `g_eff(N, b)`. The full obstruction analysis (L1–L6, six routes proven blocked) is in [`PROOF_F86B_OBSTRUCTION.md`](PROOF_F86B_OBSTRUCTION.md). Numerical Tier-1 closure of full block-L derivation is achieved via the F90 bridge identity ([`PROOF_F90_F86C2_BRIDGE.md`](PROOF_F90_F86C2_BRIDGE.md), all 22 bonds of N=5..8 agree with F86's scan to 10⁻¹⁰).

---

## Proof of Statement 2 (Universal resonance shape)

### Structural origin: 2-level EP analytics

For two adjacent rate channels (HD = 2k−1, HD = 2k+1) coupled by a bond, the effective 2×2 Liouvillian has diagonal entries (−2γ₀(2k−1), −2γ₀(2k+1)) and same-sign-imaginary off-diagonals. After shifting by the trace midpoint, the dynamics is governed by

    L_eff − (trace/2)·I  =  [ −Δ/2     iJ·g_eff ]
                            [ +iJ·g_eff   +Δ/2  ]      with Δ = 4γ₀

The same-sign-imaginary off-diagonal pattern is the non-Hermitian form that gives this toy 2×2 reduction its genuine EP at the Q_EP = 2/g_eff normalisation. This is "PT-phenomenology-like" and shares the exact Π P-type anticommutation discussed in [the PT-symmetry analysis](../../experiments/PT_SYMMETRY_ANALYSIS.md), distinct from Bender-Boettcher PT (Π is linear; classical PT requires anti-linear operators). After Π²-parity resolution Π phase-normalizes to an involution, but the complete irreducible-sector SRP class remains OPEN. See also [The Fragile Bridge](../../hypotheses/FRAGILE_BRIDGE.md) for the SEPARATE Σγ=0 gain-loss system, whose threshold at two qubits per chain is generically its own EP2 on the real γ axis, a Jordan block at γ_crit where two imaginary eigenvalues meet (K=403 is one finite-offset reading of its K ∝ 1/δ); the full Σγ=N·γ₀ block-L is genuinely non-normal near Q_peak, which is not itself an EP, while separate narrow real-axis seeds are certified at N=5,7,9 (see [`LocalGlobalEpLink`](../../compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs), OpenQuestion). The eigenvalues are

    λ_±  =  ±√((Δ/2)² − J²·g_eff²)  =  ±√(4γ₀² − J²·g_eff²)    (relative to trace/2)

with EP at J·g_eff = 2γ₀, equivalently Q_EP = 2/g_eff.

In the 2-level basis, the eigenvector rotation parameter τ = tanh(θ/2) (hyperbolic for non-Hermitian) satisfies

    τ²  =  (J·g_eff − 2γ₀) / (J·g_eff + 2γ₀)  =  (Q − Q_EP) / (Q + Q_EP)

below the EP, switching to a phase parameterisation above. The probe overlap with eigenvectors thus depends only on the dimensionless ratio Q/Q_EP. Any observable function of probe weight on dressed modes is a function of Q/Q_EP alone. Its peak then sits at a fixed multiple of Q_EP, Q_peak = x_peak·Q_EP (x_peak = 1.7141 for the bare two-level model read on its envelope), so `(Q − Q_peak)/Q_peak` is itself a function of Q/Q_EP: the argument needs Q_peak proportional to Q_EP, not equal to it.

This gives the analytical structural reason for universality: the 2-level EP physics depends only on the dimensionless ratio of detuning (J·g_eff) to gap (Δ/2 = 2γ₀). Specific values of g_eff (chain-N, bond-position, c) shift Q_peak; they don't reshape the resonance. In the full block the proportionality is not observed with F86a's g_eff = σ_0, which is one number per block: Q_peak/Q_EP reads 2.05 to 2.29 for the Interior class (rising with N) and 3.48 to 3.62 for the Endpoint class at c = 2, N = 5..8, against the bare model's 1.714, because the probe does not see the EP partners directly ([the obstruction proof](PROOF_F86B_OBSTRUCTION.md), L3). So this argument explains the form of the collapse within the two-level model; the collapse itself is measured (the table below), and the per-bond positions it does not explain.

The same 2 × 2 shape is written out in other settings of the repository, gathered in [The Double Root](../THE_DOUBLE_ROOT.md#one-2--2-in-many-settings).

### F88b: popcount-coherence Π²-odd / memory closed form (state-level inheritance from F88a, formerly "Structural inheritance from F88")

Popcount-coherence pair states |ψ⟩ = (|p⟩ + |q⟩)/√2 with popcount(p) = n_p, popcount(q) = n_q, HD(p, q) = h have a Π²-odd-fraction-within-memory determined by a precise closed form derived from Krawtchouk polynomial reflection symmetries. The formula generalises across **all** popcount pairs (adjacent n_q = n_p + 1, non-adjacent, and intra-sector n_p = n_q) and HD values, verified bit-exact via [`PopcountCoherencePi2Odd`](../../compute/RCPsiSquared.Core/Symmetry/PopcountCoherencePi2Odd.cs) against [`MemoryAxisRho`](../../compute/RCPsiSquared.Diagnostics/Foundation/MemoryAxisRho.cs) on 213 configurations at N = 2..7 (max deviation 8.88e−16).

**Static fraction** s is HD/bit-position invariant:
- Inter-sector (n_p ≠ n_q): s = 1/(4·C(N, n_p)) + 1/(4·C(N, n_q))
- Intra-sector (n_p = n_q): s = 1/C(N, n)

The kernel of L (= span{P_n} for Heisenberg + Z-dephasing) absorbs only popcount-sector content; off-diagonals project to zero in the kernel.

**α = Π²-odd-fraction of the kernel projection** has three anchor categories, all in closed form, driven by the following Krawtchouk identity (proven below):

  **Lemma.** Σ_s (−1)^s · C(N, s) · K_n(s; N) · K_m(s; N) = 2^N · C(N, n) · [n + m = N].

The lemma follows from Krawtchouk reflection K_n(s; N) = (−1)^s K_{N−n}(s; N) plus orthogonality Σ_s C(N, s) K_a(s; N) K_b(s; N) = 2^N · C(N, a) · δ_{a, b}: substituting gives Σ_s C(N, s) K_{N−n}(s; N) K_m(s; N) = 2^N C(N, m) δ_{N−n, m}. Applying this to E − O := Σ_s (−1)^s · C(N, s) · (A_s + B_s)² (with A_s = K_{n_p}(s; N)/C(N, n_p), B_s = K_{n_q}(s; N)/C(N, n_q)) yields three indicator-weighted contributions:

  E − O = (2^N / C(N, n_p)) · [n_p = N/2] + (2^N / C(N, n_q)) · [n_q = N/2] + (2 · 2^N / C(N, n_q)) · [n_p + n_q = N]

The three anchor categories follow:

- **α = 0** at popcount-mirror n_p + n_q = N (covers inter-mirror n_p ≠ n_q and intra-mirror n_p = n_q = N/2 at even N). Total Π²-odd is forced to vanish by the reflection K_{N − n}(s; N) = (−1)^s K_n(s; N) cancelling between the two sectors P_{n_p}/C(N, n_p) and P_{N − n_p}/C(N, n_q).
- **α = K-intermediate** at K-vanishing (even N with exactly one of {n_p, n_q} equal to N/2, no mirror): K_{N/2}(s; N) = 0 for odd s due to bit-flip symmetry of the half-popcount sector. The value is **α = C(N, N/2) / (2 · (C(N, n_other) + C(N, N/2)))**, where n_other ∈ {n_p, n_q} is the entry not equal to N/2. (Adjacent special case n_other = N/2 ± 1: simplifies to (N + 2)/(4·(N + 1)) using C(N, N/2 ± 1) = C(N, N/2) · N/(N + 2).) Worked values: 3/7 at N = 4 popcount-(0, 2); 3/10 at N = 4 popcount-(1, 2); 10/21 at N = 6 popcount-(0, 3); 5/13 at N = 6 popcount-(1, 3); 2/7 at N = 6 popcount-(2, 3); 5/18 at N = 8 popcount-(3, 4); 3/11 at N = 10 popcount-(4, 5).
- **α = 1/2** generic (none of n_p + n_q = N, n_p = N/2, n_q = N/2 holds). All three indicators vanish, so E − O = 0, hence Σ_{s odd} = Σ_{s even}, hence α = 1/2 exactly. **Proven for all such (N, n_p, n_q).**

**Total Π²-odd of ρ** is HD-dependent:
- HD < N (at least one matching bit): 1/2.
- HD = N (p, q complementary): 0. The off-diagonal Re(|p⟩⟨q|) has only X-and-even-Y-count Pauli strings (Y² = I cancellation kills odd-Y-count terms), and with no matching bits there is no Z-content; all surviving terms are Π²-EVEN. The diagonal residual also vanishes for popcount-(0, N) where each sector has a single basis state.

**Π²-odd / memory closed form:**

  Π²-odd / memory = ┌  0                          if HD = N (Π²-classical)
                    └  (1 / 2 − α · s) / (1 − s)  otherwise

Verified state-level table across all anchor categories:

| N | (n_p, n_q, HD) | category | static s | α | Π²-odd / memory |
|---|----------------|----------|----------|---|-----------------|
| 2 | (0, 2, 2) | HD = N (Bell+ / GHZ_2) | 1/2 | 0 | **0** |
| 2 | (1, 1, 2) | HD = N (Singlet/Triplet) | 1/2 | 0 | **0** |
| 3 | (0, 3, 3) | HD = N (GHZ_3) | 1/2 | 0 | **0** |
| 4 | (0, 4, 4) | HD = N (GHZ_4) | 1/2 | 0 | **0** |
| 5 | (0, 5, 5) | HD = N (GHZ_5) | 1/2 | 0 | **0** |
| 3 | (1, 2, 1) | inter mirror, HD < N | 1/6 | 0 | **3/5 = 0.6** |
| 5 | (2, 3, 1) | inter mirror, HD < N | 1/20 | 0 | **10/19 ≈ 0.5263** |
| 7 | (3, 4, 1) | inter mirror, HD < N | 1/70 | 0 | **35/69 ≈ 0.5072** |
| 5 | (1, 4, 3) | inter mirror, non-adjacent | 1/10 | 0 | **5/9 ≈ 0.5556** |
| 6 | (2, 4, 2) | inter mirror, non-adjacent | 1/30 | 0 | **15/29 ≈ 0.5172** |
| 4 | (2, 2, 2) | intra-mirror N/2 | 1/6 | 0 | **3/5 = 0.6** |
| 4 | (1, 2, 1) | adjacent K-intermediate | 5/48 | 3/10 | **≈ 0.5233** |
| 4 | (0, 2, 2) | non-adjacent K-intermediate | 7/24 | 3/7 | **9/17 ≈ 0.5294** |
| 6 | (0, 3, 3) | non-adjacent K-intermediate | 21/80 | 10/21 | **30/59 ≈ 0.5085** |
| 5 | (1, 2, 1) | inter generic | 3/40 | 1/2 | **1/2 exact** |
| 7 | (2, 3, 1) | inter generic | 2/105 | 1/2 | **1/2 exact** |
| 3 | (1, 1, 2) | intra generic | 1/3 | 1/2 | **1/2 exact** |

This is **F88b's bilinear-apex 1/2 inheriting to the F86 state level**: K_CC_pr's measurement subspace is centred on the framework's half-anchor in the generic case, with structured deviations at specific (N, n_p, n_q, HD) configurations. The HD = N anchor is the **Π²-classical extreme**: GHZ_N, Bell states, and intra-sector all-bits-differ states have zero Π²-odd content. This connects to [F60](../ANALYTICAL_FORMULAS.md) (pair-CΨ = 0 for GHZ_N): the same "classical" classification read from two orthogonal observables (F60 via partial-trace pair-tomography, F88b via Π²-projection on the full state).

The EP-rotation universality of Statement 2 is the dynamic consequence: a Q-resonance profile centred on a half-anchor with a c-and-N-specific shift in the structured-anchor cases (mirror, K-intermediate); the Π²-classical states at HD = N sit outside the K_CC_pr measurement scope (they have no Π²-odd memory content for the J-derivative observable to resonate with).

**Inheritance pattern (state-class extension).** The F88b closed form lifts to multi-state superpositions via two distinct inheritance mechanisms:

- **Popcount-weight invariance**: the kernel projection ρ_d0 of |ψ⟩⟨ψ| depends only on the popcount weights w_n = Σ_{i: popcount(b_i) = n} |c_i|² (kernel = span{P_n}). Any multi-state superposition with the same {w_n} as a pair state inherits the pair-state static-side formula directly. **W states** (|D_1⟩, intra-popcount-1 single sector with w_1 = 1) and the **Bonding-Bell-Pair** ((|0_R, vac⟩ + |1_R, ψ_k⟩)/√2 at N+1 qubits, w_0 = w_2 = 1/2) inherit the pair formula on both static and memory sides.

- **X⊗N-symmetry root**: the HD = N pair-state anchor (GHZ_N, Bell states, intra-sector complement pairs) and the **Dicke-mirror anchor** (|D_n⟩ + |D_{N − n − 1}⟩)/√2 at popcount-mirror 2n+1 = N share a single algebraic mechanism. Both classes are X⊗N-eigenstates; X⊗N · σ_α · X⊗N = (−1)^{bit_b(α)} · σ_α, so X⊗N-eigenstates have ⟨σ_α⟩ = 0 for all Π²-odd σ_α. The two F-axes converge: F60 reads such states as pair-CΨ = 0 (partial-trace blindness), F88b reads them as Π²-EVEN-only (projection blindness). Two structurally-distinct state-class families (HD = N pair vs. Dicke-mirror multi-state) inherit Π²-classicality from a single root symmetry, an instance of the broader F-chain inheritance pattern (cf. `project_algebra_is_inheritance.md` in memory).

The Dicke superposition (|D_n⟩ + |D_{n+1}⟩)/√2 has its own three-anchor structure for total Π²-odd-of-ρ: 0 (Dicke-mirror, 2n+1 = N), 3/8 (Dicke-K-intermediate, even N with n ∈ {N/2 − 1, N/2}), 1/2 (generic).

**Closed-form proof (Tier 1 derived, 2026-05-17):**

    α_total = (1 − γ²)/2     where γ = ⟨ψ|X⊗N|ψ⟩

Proof (4 lines):

1. Decompose ψ = c_+ψ_+ + c_-ψ_- where ψ_± are X⊗N-eigenstates (eigenvalues ±1), with |c_+|² + |c_-|² = 1.
2. ρ = |ψ⟩⟨ψ| splits under X⊗N conjugation into Π²-EVEN diagonal blocks (|c_±|²·|ψ_±⟩⟨ψ_±|) and Π²-ODD off-diagonal blocks (c_+c_-*·|ψ_+⟩⟨ψ_-| + h.c.), since X⊗N · σ_α · X⊗N = (−1)^{bit_b(α)} σ_α.
3. Π²-ODD Frobenius² = 2·|c_+c_-|² = 2·|c_+|²·|c_-|² (the two off-diagonal terms are Hilbert-Schmidt orthogonal).
4. With γ = ⟨ψ|X⊗N|ψ⟩ and the standard decomposition |c_±|² = (1±γ)/2: α_total = 2·(1+γ)(1−γ)/4 = (1−γ²)/2.

The three anchors follow from X⊗N's action on Dicke pairs: X⊗N maps |D_k⟩ → |D_{N−k}⟩, so X⊗N|ψ⟩ = (|D_{N−n}⟩+|D_{N−n−1}⟩)/√2. The overlap γ = (1/2)·(δ_{n,N−n} + δ_{n,N−n−1} + δ_{n+1,N−n} + δ_{n+1,N−n−1}) hits exactly three values depending on which Kronecker deltas align with the (n, n+1) pair:

- **Mirror (N odd, 2n+1=N)**: both δ_{n,N−n−1} and δ_{n+1,N−n} hit → γ = 1 → α = 0. ψ is X⊗N-eigenstate, no Π²-odd content.
- **K-intermediate (N even, n=N/2)**: only δ_{n,N−n} hits → γ = 1/2 → α = 3/8.
- **K-intermediate (N even, n=N/2−1)**: only δ_{n+1,N−n−1} hits → γ = 1/2 → α = 3/8.
- **Generic**: no δ hits → γ = 0 → α = 1/2.

Bit-exact verified for N = 3..7 across all (n, n+1) pairs ([`simulations/f86b_dicke_pi2odd_closed_form.py`](../../simulations/f86b_dicke_pi2odd_closed_form.py)). The proof uses the same orthogonal-symmetry-decomposition mechanism that produced F50's max-spin Dicke endpoint ladder rungs (commit 5523171, [the weight-1 degeneracy proof § max-spin closed-form](PROOF_WEIGHT1_DEGENERACY.md)): symmetry under X⊗N (or equivalently Π = Z⊗N for the F50 case) gives an eigenbasis decomposition that makes the operator algebra computation explicit, replacing the opaque Krawtchouk sum with a clean ±1-eigenspace overlap.

The mechanism extends the F88b algebraic-inheritance principle (`compute/RCPsiSquared.Core/Symmetry/PopcountCoherencePi2Odd.cs`): the **X⊗N-symmetry root** that gives HD = N pair states (GHZ_N, Bell, intra-complements) their Π²-classical content also gives the Dicke superposition its three-anchor structure, but with finer 4-element δ-resolution because the Dicke superposition is a 2-state superposition (not a single state) under X⊗N.

The 4-mode-basis vectors {|c_1⟩, |c_3⟩, |u_0⟩, |v_0⟩} from [`FourModeBasis`](../../compute/RCPsiSquared.Core/Decomposition/FourModeBasis.cs) all report Π²-odd/mem = 0.5000 when embedded as density matrices; per-qubit [`BlochAxisReading`](../../compute/RCPsiSquared.Diagnostics/Foundation/BlochAxisReading.cs) makes the probe-vs-EP-partner orthogonality visible as a **single-body fingerprint**: |c_1⟩ has 1-body Pauli content (X+ per qubit), while |c_3⟩, |u_0⟩, |v_0⟩ have **zero** 1-body Bloch (their content is purely multi-body). The probe-EP-partner orthogonality central to the 4-mode structure manifests at the per-qubit reading.

For the typed-claim lineage in [`Pi2KnowledgeBase`](../../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBase.cs): the F86 layer surfaces three trio claims simultaneously, `HalfAsStructuralFixedPointClaim` (1/2 number-anchor), `BilinearApexClaim` (apex at p = 1/2 of any bilinear), and `NinetyDegreeMirrorMemoryClaim` (F80's i factor is the 90°-rotation that produces the same-sign +iJ·g_eff structure of `L_eff`). The EP-rotation universality is one structural fact viewed from the F86 layer.

### Why HWHM_left is universal but HWHM_right diverges

Pre-EP region (Q < Q_EP): discriminant `4γ₀² − J²·g_eff² > 0`, eigenvalues real, dressed-mode weight rises monotonically. Universal in Q/Q_EP.

Post-EP region (Q > Q_EP): discriminant negative, eigenvalues complex conjugate pair. Probe sits near 99 % on dressed modes (saturated). Long-tail behaviour depends on:
- Higher-channel dressed-mode contributions (chain-N, c specific)
- Time-averaging behaviour of complex-eigenvalue oscillations in K_CC_pr

These chain-specific details enter the post-EP tail. Hence: a pre-EP rise that is nearly universal within each bond class (HWHM_left/Q_peak in the bands 0.7523..0.7637 Interior and 0.7710..0.7815 Endpoint, with a residual N-drift inside the Interior band), and a post-EP tail that is not: at x = +1.0 the envelope reads, at c = 2, N = 5, 6, 7, Interior 0.837 / 0.843 / 0.884 and Endpoint 0.864 / 0.874 / 0.876, so Endpoint lies above Interior at N = 5 and 6 and below it at N = 7, where the Interior class mean carries inner bonds rising toward their far maxima (`ShapeFunctionWitnesses`).

The closed-form derivation of the HWHM_left/Q_peak constants, the exploration record behind it, and the obstruction proof on g_eff are carried in the sibling [`PROOF_F86B_OBSTRUCTION.md`](PROOF_F86B_OBSTRUCTION.md).

---

## Empirical universal-shape data

y = K/|K|max of the class-averaged envelope max_t |K|, evaluated at relative shift x = (Q−Q_peak)/Q_peak, with Q* and |K|max the class peak on [0.2, 4] (the reading of the HWHM band) and the curve read on [0.2, 6] in steps of 0.025 so that x = +1 stays on the grid; linear between grid points, as in the C# `ShapeFunctionWitnesses`. Read by [`simulations/f86_envelope_second_route.py`](../../simulations/f86_envelope_second_route.py) `--lineshape N c`; its Q* equals the class Q_peak of the tables below, and at c = 2, N = 5..7 its x = +1 values equal the live C# ones.

**Interior:**

| x = (Q−Q*)/Q* | c2N5 | c2N6 | c2N7 | c2N8 | c3N5 | c3N6 | c3N7 | c3N8 | c4N7 | c4N8 | range |
|---|---|---|---|---|---|---|---|---|---|---|---|
| −0.60 | 0.722 | 0.729 | 0.733 | 0.734 | 0.730 | 0.733 | 0.736 | 0.737 | 0.738 | 0.740 | 2.4 % |
| −0.40 | 0.900 | 0.903 | 0.905 | 0.906 | 0.904 | 0.905 | 0.906 | 0.907 | 0.907 | 0.908 | 1.0 % |
| −0.20 | 0.980 | 0.981 | 0.981 | 0.981 | 0.981 | 0.981 | 0.981 | 0.982 | 0.982 | 0.982 | 0.2 % |
| 0.00 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0 % (peak) |
| +0.20 | 0.987 | 0.987 | 0.988 | 0.988 | 0.987 | 0.988 | 0.988 | 0.988 | 0.988 | 0.988 | 0.1 % |
| +0.40 | 0.956 | 0.958 | 0.959 | 0.960 | 0.959 | 0.959 | 0.960 | 0.960 | 0.960 | 0.961 | 0.4 % |
| +1.00 | 0.837 | 0.843 | 0.884 | 0.935 | 0.844 | 0.846 | 0.847 | 0.849 | 0.849 | 0.851 | 10.5 % |

**Endpoint:**

| x = (Q−Q*)/Q* | c2N5 | c2N6 | c2N7 | c2N8 | c3N5 | c3N6 | c3N7 | c3N8 | c4N7 | c4N8 | range |
|---|---|---|---|---|---|---|---|---|---|---|---|
| −0.60 | 0.759 | 0.767 | 0.767 | 0.766 | 0.750 | 0.759 | 0.762 | 0.764 | 0.760 | 0.761 | 2.2 % |
| −0.40 | 0.917 | 0.921 | 0.921 | 0.920 | 0.913 | 0.917 | 0.918 | 0.919 | 0.917 | 0.918 | 0.9 % |
| −0.20 | 0.984 | 0.985 | 0.985 | 0.985 | 0.983 | 0.984 | 0.984 | 0.984 | 0.984 | 0.984 | 0.2 % |
| 0.00 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0 % (peak) |
| +0.20 | 0.989 | 0.990 | 0.990 | 0.990 | 0.989 | 0.989 | 0.990 | 0.990 | 0.989 | 0.990 | 0.1 % |
| +0.40 | 0.965 | 0.967 | 0.967 | 0.967 | 0.963 | 0.965 | 0.966 | 0.966 | 0.965 | 0.965 | 0.5 % |
| +1.00 | 0.864 | 0.874 | 0.876 | 0.874 | 0.856 | 0.864 | 0.868 | 0.870 | 0.865 | 0.867 | 2.2 % |

### HWHM_left/Q_peak across all read cases (two bond classes)

The envelope reading: Q_peak is the argmax over Q of max_t |K(Q, t)| of the class-averaged curve, at γ₀ = 0.05 on the default Q grid (153 points on [0.2, 4], parabolic vertex, linear half-maximum), read by `ResonanceScan` and on a second route by [`simulations/f86_envelope_second_route.py`](../../simulations/f86_envelope_second_route.py), the two agreeing to 10⁻¹⁰ (at c = 4, N = 8 to 3·10⁻¹¹ in the ratio and 2·10⁻¹⁰ in Q_peak). Against a continuous search in Q, made per bond at c = 2, N = 5..8 (not on the class-mean curves), the grid reading is off by up to 2.2·10⁻⁴ in Q_peak and 3.5·10⁻⁵ in the ratio.

**Interior bonds:**

| case | HWHM_left/Q_peak | Q_peak |
|------|------------------|--------|
| c=2 N=5 | 0.7523 | 1.4831 |
| c=2 N=6 | 0.7569 | 1.5652 |
| c=2 N=7 | 0.7591 | 1.5975 |
| c=2 N=8 | 0.7602 | 1.6110 |
| c=3 N=5 | 0.7572 | 1.5997 |
| c=3 N=6 | 0.7596 | 1.6789 |
| c=3 N=7 | 0.7611 | 1.7119 |
| c=3 N=8 | 0.7622 | 1.7232 |
| c=4 N=7 | 0.7623 | 1.7605 |
| c=4 N=8 | 0.7637 | 1.7914 |

Interior band: 0.7523..0.7637 (1.5 %). The values rise monotonically with N at each c; the c = 2 value continues to 0.7621 at N = 12.

**Endpoint bonds:**

| case | HWHM_left/Q_peak | Q_peak |
|------|------------------|--------|
| c=2 N=5 | 0.7766 | 2.5171 |
| c=2 N=6 | 0.7815 | 2.5806 |
| c=2 N=7 | 0.7813 | 2.5601 |
| c=2 N=8 | 0.7808 | 2.5431 |
| c=3 N=5 | 0.7710 | 2.3926 |
| c=3 N=6 | 0.7764 | 2.5482 |
| c=3 N=7 | 0.7785 | 2.5837 |
| c=3 N=8 | 0.7792 | 2.5833 |
| c=4 N=7 | 0.7768 | 2.5752 |
| c=4 N=8 | 0.7777 | 2.5908 |

Endpoint band: 0.7710..0.7815 (1.4 %).

γ₀-invariance: exact by the scaling of the generator, γ₀·(D̂ + Q·M̂), with the t window scaling alongside; `ResonanceScanTests.ClassPeaks_AreGammaZeroInvariant` checks both classes at c = 2, N = 5 for γ₀ ∈ {0.025, 0.05, 0.10} to 10⁻¹⁰.

### What this tells us

1. **Two class-averaged bands, separated.** The Interior values lie in 0.7523..0.7637, the Endpoint values in 0.7710..0.7815; the bands do not overlap, the gap between them (0.0073) is smaller than the width of either (0.0113, 0.0105), and much of the Interior width is the systematic rise with N, not random scatter. Per bond the sub-classes carry their own values (at c = 2 on the canonical reading, 0.7505..0.7815); the two bands are the F71-class average.

2. **γ₀ invariance is exact.** It follows from Q = J/γ₀ and a t window that scales with t_peak = 1/(4γ₀): at fixed Q a change of γ₀ rescales time and the size of K. This confirms Q as the inside-observable scale.

3. **c=2 confirms the structure is real, and N does not drop out.** At c=2 the 2-level model is exact (only HD ∈ {1, 3} channels exist, no orthogonal complement), and the per-bond sub-class structure plus the class-averaged split are already present (Interior 0.7523..0.7602, Endpoint 0.7766..0.7815, N = 5..8). So the bond-class distinction lives in the bond-position-dependent probe-overlap profile, not in higher-c orthogonal-complement physics. At c = 2 and c = 3 the Interior value rises with N, by about 0.008 and 0.005 over N = 5..8; whether it settles is not known.

4. **The post-peak tail does not separate the classes.** At x = +1.0 the c = 2 envelope reads Interior 0.837 / 0.843 / 0.884 and Endpoint 0.864 / 0.874 / 0.876 at N = 5, 6, 7: Endpoint above Interior at N = 5 and 6, below it at N = 7, where the Interior class mean carries the inner bonds rising toward their far maxima.

---

## Pointers

**Hub:** [the Q-peak hub](PROOF_F86_QPEAK.md): three-theorem overview and the shared reference list.
**Sibling theorems:** [the EP-mechanism proof](PROOF_F86A_EP_MECHANISM.md) (F86a), [the g-eff obstruction proof](PROOF_F86B_OBSTRUCTION.md) (the closed-form gap, exploration record, and g_eff obstruction proof), [the F71-mirror proof](PROOF_F86C_F71_MIRROR.md) (F86c).
**F90 bridge:** [the F90 bridge proof](PROOF_F90_F86C2_BRIDGE.md): F86 c=2 K_b = F89 path-(N−1) per-bond Hellmann-Feynman; the numerical-Tier-1 route for Direction (b'').
**State-level inheritance:** F88b ([the formula registry](../ANALYTICAL_FORMULAS.md#f88b)) inheriting from F88a, `PopcountCoherencePi2Odd`, `MemoryAxisRho`.
**HWHM closed form (F86b₂):** `F86HwhmClosedFormClaim`, `BondSubClass` in `compute/RCPsiSquared.Core/F86/Item1Derivation/`.
**Chiral classification anchor:** [the PT-symmetry analysis](../../experiments/PT_SYMMETRY_ANALYSIS.md).
**Scripts:** [`simulations/f86_envelope_second_route.py`](../../simulations/f86_envelope_second_route.py) (the envelope Q_peak and HWHM per bond and class, second route of `ResonanceScan`); on the earlier 21-point t-grid reading, [`eq022_b1_step_e_resonance_shape.py`](../../simulations/eq022_b1_step_e_resonance_shape.py) + [`eq022_b1_step_e_inspect.py`](../../simulations/eq022_b1_step_e_inspect.py) (the lineshape table for c=3, c=4 at γ₀=0.05) and [`eq022_b1_step_f_universality_extension.py`](../../simulations/eq022_b1_step_f_universality_extension.py) (c=2 sweep + γ₀ ∈ {0.025, 0.10} check that established the two-bond-class refinement).
**C# OOP layer:** `compute/RCPsiSquared.Core/F86/` carries `UniversalShapeBand`, `UniversalShapePrediction` + `UniversalShapeWitness`, `ShapeFunctionWitnesses`, `C2UniversalShapeDerivation`. CLI: `rcpsi inspect --root f86 --with-measured`.
