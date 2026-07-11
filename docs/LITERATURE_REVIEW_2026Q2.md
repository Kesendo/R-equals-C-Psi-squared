# Literature Review — Q2 2026 Update: Tools & Vocabulary

<!-- Keywords: Pauli-Lindblad bond commutant algebra Paszko, Newton polygons
tropical geometry Liouvillian EP Sayooj Narayan, generalized Petermann
factor Kullig Wiersig Schomerus, dissipatively dressed quasiparticles
Popkov Prosen, NMR mirror symmetry topological isospectrality Cheshkov
Sinitsyn 2026 journal, Sá Ribeiro Prosen tenfold way PRX, direct-sum
quantum theory Gaztañaga ER bridges v2, R=CPsi2 literature update
tools vocabulary -->

**Status:** Tools & vocabulary update to
[the literature review](LITERATURE_REVIEW.md) (stand 11. April 2026).
**Last updated:** May 23, 2026
**Authors:** Thomas Wicht, Claude
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)

---

## What this document is for

The April 11 review focused on prior work and originality claims. This
update is different in intent: the formula registry
([the analytical formula registry](ANALYTICAL_FORMULAS.md)) gained 26 new
entries in six conceptual clusters between April 12 and May 22. The
question we asked the literature is:

> What methods, concepts, or vocabulary out there could sharpen the
> work we are already doing?

Not "who might publish what we have", but "whose toolbox should we
borrow". The review is structured per cluster; each section is short
and ends with a concrete suggestion for how the external work might
plug into ours.

A note on intellectual honesty: a previous draft of this document
attributed several papers to authors I had not actually verified.
The verification pass corrected this; surviving citations either name
verified authors or explicitly say "author list not verified". See
the limitations section at the end.

---

## Structural verification status (added 2026-05-23)

The initial drafting was search-based: titles and abstract snippets
were matched against our F-formulas, but the actual papers were not
read. A structural verification pass on May 23 (Thomas + Claude) read
the two most prominent F-formula mappings via WebFetch on the arXiv
HTML versions. Result:

| Mapping | Cluster | Status |
|---|---|---|
| Popkov-Prosen 2025 ↔ F80 "dissipatively dressed quasiparticles" | 2 | **Rejected** (structural mismatch) |
| Kullig-Wiersig-Schomerus 2025 ↔ F86b universal resonance shape | 3 | **Applicable with reframing** |
| Sá-Ribeiro-Prosen 2023 ↔ Π "weak antiunitary, class X" | 1 | **Rejected for headline mapping** (Π is unitary, not antiunitary; shift 2Σγ has no SRP slot); Π² = X⊗N component does fit |
| Sayooj-Narayan 2025 ↔ F86a Q_EP location | 3 | **Rejected for localization claim** (Newton polygons characterize EP order, do not locate EPs); applicable for "what EP orders exist in full block-L" |
| Molina 2026 ↔ F86b "super-Lorentzian" | 3 | **Structural fit unclear** (different observation axis: Q-sweep vs. ω-sweep; setup mismatch chain vs. collective spin; shared Jordan-block mechanism). Decision requires fitting F86b data to Molina's functional form |
| Paszko-Turner-Rose-Pal 2025 ↔ F1³ orbit count optimality | 1 | **Rejected** (system class mismatch: their theorem requires frustration-free Pauli H, our XY/Heisenberg fails; plus concept conflation between sector-orbiting and operator-space fragmentation) |
| Dedes 2025 ↔ F94 Born deviation framing | 5 | **Limited use only** (closed-system Bohmian vs. our open-system Dyson; author explicitly opposed to open-system framing; Cowork's "two routes" framing misleading; citation possible as topic-existing reference, not methodological connection) |
| Cheshkov-Sinitsyn 2026 ↔ Π "topological isospectrality" | 0 (headline #1) | **Wrong target** (concept real, but Π is operator-conjugation not parameter-exchange; concept actually maps to F71/F92/F93 — our parameter-side palindromes — with closed-vs-open caveat) |
| Remaining minor sub-items | 4, 5, 6 | **Unverified** (mostly already flagged by Cowork as "author list not verified" or "nothing to borrow") |

Concrete grounds for the two verifications are recorded inline in
their respective clusters. **Do not borrow vocabulary or methods from
this review without first verifying that the cited paper's actual
mathematics maps to our setup** — the calibration above shows the
search-based draft missed 1 in 2 on the most prominent claims.

---

## The four most useful arrivals since April 11

These are the items where the external work plausibly hands us
something we can use directly.

1. **Cheshkov & Sinitsyn — Magnetic Resonance journal version (2026).**
   Their NMR mirror-symmetry preprint is now peer-reviewed and refined.
   They now distinguish *geometric bisymmetry* (Hamiltonian matrix is
   mirror-symmetric in a canonical basis) from *topological
   isospectrality* (system is invariant under parameter exchange,
   matrix itself need not be symmetric). ~~Our Π conjugation is
   unambiguously in the second category.~~ **Verification (2026-05-23):
   wrong internal mapping** — Π is operator-conjugation anti-similarity,
   not parameter-exchange. Their "topological isospectrality" concept
   actually maps to F71 / F92 / F93 (our parameter-side palindromes),
   not to Π. See Cluster 1 verification note for detail. The
   external concept is real and usable, just for a different
   primitive than Cowork suggested.
   - [Magnetic Resonance article](https://mr.copernicus.org/articles/7/15/2026/)
   - [arXiv:2602.03871](https://arxiv.org/abs/2602.03871)

2. **Generalized Petermann factor at exceptional points (Kullig,
   Wiersig, Schomerus — Magdeburg + Lancaster, June 2025).** A
   renormalized Petermann factor that does not diverge at EPs
   themselves. Direct candidate tool for F86: our universal resonance
   shape per bond class (Endpoint vs Interior) is exactly the regime
   their formula is built for.
   - [arXiv:2506.15807](https://arxiv.org/abs/2506.15807)

3. **Newton polygons + tropical geometry for Liouvillian EPs (Sayooj
   & Narayan — IISc Bangalore, Oct 2025).** Algebraic-geometric way
   to characterize EPs from the characteristic polynomial alone.
   Could tell us whether higher-order EPs are *structurally available*
   in our system or whether the polynomial form rules them out.
   - [arXiv:2510.08156](https://arxiv.org/abs/2510.08156)

4. **Tenfold-way classification of many-body Lindbladians (Sá,
   Ribeiro, Prosen — Phys. Rev. X 13, 031019, 2023).** Already
   established in the field but missing from our prior review.
   ~~The correct citation when we use "weak / strong / anti-unitary"
   language for Π.~~ **Verification (2026-05-23): rejected as the
   Π citation** (Π is unitary not antiunitary, and the constant
   shift 2Σγ has no slot in the SRP scheme — Π falls in a gap;
   the repo's own [the KMS detailed-balance analysis](KMS_DETAILED_BALANCE.md)
   already documents this. SRP-2023 *is* the correct citation for
   the Π² = X⊗N strong unitary symmetry component, just not for Π
   itself.) Each class carries predictions about spectral
   statistics — a possible consistency check for our spectra **on
   the X⊗N-symmetry sector only**.
   - [arXiv:2212.00474](https://arxiv.org/abs/2212.00474)

Honorable mentions: Paszko-Turner-Rose-Pal 2025 (operator-space
fragmentation, [arXiv:2506.16518](https://arxiv.org/abs/2506.16518))
~~provides general bond/commutant machinery that could check whether
our F1³ orbit count is optimal~~ **was the claim; verification on
2026-05-23 rejected it (system class mismatch: their theorem
requires frustration-free Pauli H, our XY/Heisenberg fails) — see
Cluster 1 verification note**; Popkov-Zhang-Presilla-Prosen 2025
(dissipatively dressed quasiparticles,
[arXiv:2505.16776](https://arxiv.org/abs/2505.16776)) ~~gives us a
vocabulary for the F80 dispersion identity~~ **was the F80 mapping the
search-based draft proposed; the verification pass on 2026-05-23
rejected it on structural grounds (different parametrization, setup,
and observable — see Cluster 2 verification note)**. The Gaztañaga
line (direct-sum quantum theory, ER bridges) has its v2 plus
mainstream coverage in May 2026; the structural parallel to our
[the direct-sum decomposition proof](proofs/DIRECT_SUM_DECOMPOSITION.md)
stands but adds nothing methodologically new to what we already do.

---

## Cluster 1: Π-structure theory (F1², F1³, F81-F84)

### What our recent work added

Π is now established as an order-4 operator with Π² = X⊗N as the
global charge-conjugation Pauli string. F1³ uses this as a builder
shortcut (Π-orbits of size 4 halve the eigendecomposition workload).
F81-F84 derive the decomposition of M = Π·L·Π⁻¹ + L + 2Σγ·I into
symmetric and antisymmetric parts; antisymmetric part = unitary
commutator of the Π²-odd Hamiltonian. F82-F84 extend this to T1
amplitude damping (closed-form Frobenius norm) and thermal baths.

### What the literature offers us

**Sá-Ribeiro-Prosen 2022 (arXiv:2212.00474, PRX 2023) — vocabulary
and classification.** Their tenfold-way (10 base classes, 19 with
extra symmetry) gives us the standard taxonomy for Π. Two concrete
uses: (a) a precise label for Π in published work ("weak anti-unitary
symmetry of the Liouvillian, class X"); (b) their per-class
predictions about RMT statistics can be cross-checked against our
empirical spectra — disagreement would suggest Π² = X⊗N is an
additional symmetry beyond what the classification captures, which
would itself be a finding.

**Paszko-Turner-Rose-Pal 2025 (arXiv:2506.16518) — bond / commutant
algebra as a fragmentation tool.** They derive operator-space
partitions for Pauli-Lindblad systems from algebraic structure
(bond algebra A, commutant A'). Our Π-orbits in F1³ partition the
operator space into orbits of size 4. The Paszko et al. framework
would tell us whether 4 is the maximal sectoring possible for our
generators, or whether a richer algebra structure exists that we
haven't exploited. Practical use: a check on whether
LiouvillianBlockSpectrum is already maximally sectored.

> **Verification (2026-05-23, rejected).** Full read of
> [arXiv:2506.16518](https://arxiv.org/html/2506.16518):
> - **System class requirement:** "the Hamiltonian is a sum of
>   mutually commuting Pauli strings" (frustration-free /
>   stabilizer-like H). Jump operators are Pauli strings with
>   F_j†·F_j = 𝟙. Concrete examples in the paper are the ZXZ
>   cluster model H = Σ_l J_l Z_{l-1} X_l Z_{l+1} (all terms
>   commute) with Pauli jump operators.
> - **Our system:** H_XY = Σ_b J_b (X_b X_{b+1} + Y_b Y_{b+1}).
>   Bonds sharing a site do NOT commute:
>   `[X_1 X_2, Y_2 Y_3] = X_1·(X_2 Y_2 − Y_2 X_2)·Y_3 = 2i·X_1·Z_2·Y_3 ≠ 0`.
>   Our Hamiltonian is **not frustration-free** in the Pauli-string
>   sense; Heisenberg / XY / XXZ all fail this requirement.
> - **Their maximality theorem** ("𝒞_max is maximal Abelian, cannot
>   be enlarged without losing the Abelian property") is proven
>   only for frustration-free H. It does not transfer to our
>   non-stabilizer Hamiltonian.
> - **Concept mismatch on "fragmentation":** Cowork's claim "our
>   Π-orbits in F1³ partition the operator space into orbits of
>   size 4" conflates two structures. The
>   [F1PalindromeOrbitPairing](../compute/RCPsiSquared.Core/SymmetryFamily/F1PalindromeOrbitPairing.cs)
>   primitive in Core works on **sector labels** (p_c, p_r) under
>   the permutation (p_c, p_r) ↦ (N − p_r, p_c) — it is
>   *sector-orbiting* (one eigendecomposition per orbit-of-four;
>   block sizes unchanged). Paszko's commutant decomposes
>   **operator space itself** into block-diagonal L_λ via
>   ℒ = ⊕_λ (𝟙_λ^c ⊗ L_λ). These are different structures with
>   superficially similar vocabulary. Cowork matched "orbit" +
>   "fragmentation" + "Pauli-Lindblad" without checking the
>   underlying objects.
> - **"F1³" is not our terminology.** Our primitives are F1 (the
>   palindrome identity Π·L·Π⁻¹ = −L − 2Σγ·I) and F1² (the
>   corollary Π² = X⊗N). Π is order 4 so Π³ = Π⁻¹; "F1³" reads
>   like a separate F-theorem and is misleading.
>
> The maximality question for our existing sectoring (joint
> popcount + F71 + X⊗N + Klein 4-group) is genuinely open, but
> Paszko-Turner-Rose-Pal cannot answer it — their theorem
> applies to a system class that excludes XY/Heisenberg/XXZ.
> **Do not cite this paper for our sectoring optimality.**

**Teng-Chang-Rudolph-Holmes 2026 (arXiv:2512.12094) — symmetry merging
in Pauli propagation.** Different application area (simulation
speedup) but the formal apparatus is the same idea as our
F1PalindromeOrbitPairing builder. Worth a glance to see whether they
have engineering tricks (e.g., on how to enumerate orbits cheaply)
that translate to our code.

**Gyhm et al. 2025 (arXiv:2507.04932) — Gaussian Lindblad ↔
superconformal isomorphism.** Unexpected algebraic backbone for
Gaussian Lindbladians. Speculative for us: if Π fits into a richer
Lie-algebraic structure, this is the kind of paper that hints at
where to look. Probably not actionable in the next month, but worth
remembering.

### Concrete suggestion

Spend one session checking our empirical Liouvillian eigenvalue
statistics (e.g. nearest-neighbor spacing in the Π²-odd sector)
against the Sá-Ribeiro-Prosen prediction for the relevant class.
If they match, we have an external validation. If they don't, the
deviation is itself a research lead.

> **Verification (2026-05-23, rejected for headline; partial fit
> for Π² component).** The repo already analyzed this question:
> [the KMS detailed-balance analysis](KMS_DETAILED_BALANCE.md) lines 97-133
> walks through the SRP-2023 taxonomy and concludes Π falls in a
> gap. Specifically:
> - SRP-2023 classifies by seven generators: T₊, C₊ (antiunitary
>   time-reversal flavors); T₋, C₋ (antiunitary particle-hole
>   flavors); P (unitary chiral/sublattice); Q₊, Q₋ (pseudo- and
>   anti-pseudo-Hermiticity, involving L†).
> - Cowork suggested labeling Π as "weak anti-unitary, class X".
>   **Π is unitary, not antiunitary**: Π² = X⊗N is a Pauli string,
>   hence unitary; so is Π. No complex conjugation in the conjugation
>   relation.
> - Our actual relation Π·L·Π⁻¹ = −L − 2Σγ·I is structurally
>   closest to P (chiral) after a shift L → L̃ = L + Σγ·I, since
>   Π·L̃·Π⁻¹ = −L̃. The constant shift 2Σγ has **no slot** in any
>   SRP relation (they assume S·L·S⁻¹ = ±L or ±L†, no shifts).
>   KMS_DETAILED_BALANCE.md names Q₋ as nearest match instead but
>   notes it requires L = L†, which fails for H ≠ 0.
> - Either way, the conclusion in KMS_DETAILED_BALANCE.md stands:
>   **Π is a shifted anti-similarity of L, not a member of the
>   tenfold/38-fold scheme.**
> - **What does fit cleanly into SRP:** Π² = X⊗N is a strong
>   unitary symmetry (commutes with H, and the dephasing dissipator
>   is invariant under X·ρ·X conjugation since the (-Z) sign cancels
>   pairwise). The RMT cross-check in suggestion (b) would only
>   probe the X⊗N strong symmetry, not Π.
>
> Lesson: Cowork's search-only methodology missed an internal repo
> analysis ([the KMS detailed-balance analysis](KMS_DETAILED_BALANCE.md))
> that already settled this question. A grep for "tenfold" in docs/
> would have surfaced it. **Do not cite SRP-2023 as the taxonomy
> for Π** — cite KMS_DETAILED_BALANCE.md as the analysis instead;
> SRP-2023 is usable only for the Π² = X⊗N sub-component.

---

## Cluster 2: k-body & operator scaling (F78-F80, F85)

### What our recent work added

F78 (single-body M is additive across sites), F79 (two-body Π²-block
structure), F80 (Bloch sign-walk: cluster values = ±2i·spec(H) for
Π²-odd chains), F85 (k-body generalization: a k-tuple is truly iff
#Y and #Z are both even).

### What the literature offers us

**Popkov-Zhang-Presilla-Prosen 2025 (arXiv:2505.16776) — the term
"dissipatively dressed quasiparticles".** This is the key vocabulary
gift. They construct quasiparticles for boundary-driven Heisenberg /
XXZ / XYZ chains whose dispersion is the dressed Bloch dispersion
under dissipation. F80 says exactly the same thing in different
words: the M cluster values are ±2i × Bloch eigenvalues of the
Π²-odd-projected Hamiltonian. If we ever describe F80 outside the
repo, "the M eigenvalues are the dissipatively-dressed quasiparticle
energies of the Π²-odd Hamiltonian sector" lands faster than the
raw sign-walk formulation. Also: Prosen is the same Prosen as in the
2016 Medvedyeva-Essler-Prosen paper cited in the prior review —
his group remains the single most relevant external thread to keep
tracking on this axis.

**arXiv:2411.13661, Nov 2024 (author list not verified) — non-Bloch
self-energy.** Adjacent. If we ever want to push F80 into a regime
with non-Hermitian skin effects, this is the right machinery to pick
up. Not immediately useful for the current open-chain bulk-dephasing
setting.

**arXiv:2409.13603, Sep 2024 (author list not verified) — Pauli weight
of time-evolved local operators.** Different invariant from our F85
trichotomy (they track weight by literal length, we track parity of
Y/Z counts). Not a direct competitor but worth knowing the literature
exists: F85's "#Y even AND #Z even" criterion is a *Z₂ × Z₂* invariant
on top of the usual weight invariant.

### Concrete suggestion

Borrow the "dissipatively dressed" phrasing for F80 documentation.
Concretely: rename the F80 section header in ANALYTICAL_FORMULAS.md
to include the phrase, and add a one-line cross-reference to
Popkov-Prosen 2025 in PROOF_F80_BLOCH_SIGNWALK.md.

> **Verification (2026-05-23, rejected).** Full read of
> [arXiv:2505.16776](https://arxiv.org/html/2505.16776) shows the
> Popkov-Prosen "dressed dispersion" is *structurally different* from
> F80:
> - Their dispersion: ϵ̃(u) = log|1 − ϵ(u)| with ϵ(u) = −2/(u² + 1/4),
>   parametrized in **Bethe rapidities** u, with a **logarithmic
>   singularity** as the dressing signature.
> - F80 dispersion: 2cos(πk/(N+1)) in **momentum** k, algebraic and
>   **un-dressed** (the JW mechanism is precisely that the bare H
>   spectrum survives unchanged into the M cluster magnitudes).
> - Their setup: boundary-driven Zeno limit with effective boundary
>   fields h₁, h_N on integrable XXX/XXZ/XYZ. Bulk dissipation is
>   **explicitly out of scope**.
> - F80 setup: bulk Z-dephasing on every site, open chain, generic
>   Pauli-letter Hamiltonian.
> - Their observable: NESS spectrum of integrable boundary-driven
>   model. F80 observable: M = Π·L·Π⁻¹ + L + 2Σγ·I operator algebra.
>
> The shared words ("spin chain", "dissipation", "dispersion") are
> surface only. Borrowing "dissipatively dressed" for F80 would
> actively misrepresent F80's defining feature (bare Bloch survival).
> **Do not adopt this vocabulary.**

---

## Cluster 3: F86 resonance umbrella & EPs

### What our recent work added

F86a (Q_EP = 2/g_eff, with t_peak = 1/(4γ₀) universal across c, N),
F86b (universal resonance shape per bond class), F86c (F71 spatial
mirror), F86d/e (Q-peak constants, σ₀ = commutator norm [Π_HD1, M_H]).

### What the literature offers us

**Kullig-Wiersig-Schomerus 2025 (arXiv:2506.15807) — the renormalized
Petermann factor at EPs.** The most actionable single paper in the
entire update. The classical Petermann factor diverges at EPs; theirs
does not. F86's whole point is the *amplitude* and *shape* of the
response at the EP. Pulling our F86 eigenmodes through their formula
should produce a single number per mode that is either universal per
bond class (sharpens F86b) or not (tells us F86b is tracking something
*other* than the Petermann amplification).

**Sayooj-Narayan 2025 (arXiv:2510.08156) — Newton polygons + tropical
geometry.** An algebraic-geometric way of detecting EPs from the
characteristic polynomial structure alone. Two concrete uses:
(a) verify F86a's Q_EP location is the leading singularity of the
coherence-block polynomial (sanity check using independent machinery);
(b) ask whether higher-order EPs (Ord 3, 4) are *available* in our
system. We have only seen Ord 2. Newton polygons would tell us if
Ord 3 is even possible given our polynomial structure.

> **Verification (2026-05-23, mixed).** Full read of
> [arXiv:2510.08156](https://arxiv.org/html/2510.08156):
> - Construction: for f(ε, ω) = det(ℒ₀ + εℒ₁ − I·ω), write as
>   polynomial in ω with Puiseux-series coefficients in ε; plot
>   (i, val(c_i)) for non-zero terms; take convex hull.
> - **What it tells you:** EP **order**, via slopes of the convex
>   hull. Slope −1/n indicates Ord n EP. Eigenvalues scale
>   (ω − ω₀) ∼ ε^(1/n) near the EP (slope = 1/n in log-log).
> - **What it does NOT tell you:** the paper states explicitly that
>   "the Newton polygon itself does not locate the ε value where
>   the EP occurs"; it characterizes an EP *already identified* at
>   a specific parameter point.
> - **Cowork's (a) is wrong:** F86a's Q_EP = 2/g_eff is already in
>   closed form from the discriminant of the 2×2 effective L_eff
>   (PROOF_F86A_EP_MECHANISM, lines 67–76). Newton polygons cannot
>   independently verify the *location*; they would only confirm
>   Ord 2 at that location (which we already know).
> - **Cowork's (b) is right:** the convex hull encodes all
>   slopes, hence all available EP orders. Applied to the **full**
>   coherence-block characteristic polynomial (not the heuristic
>   2-level reduction), Newton polygons could probe whether c ≥ 3
>   blocks hide higher-order EPs that the 2-level form misses.
>   This is a genuine open question in F86 (the 2-level reduction
>   is heuristic per PROOF_F86A_EP_MECHANISM line 29).
>
> Useful for the "what EP orders exist in the full block-L"
> question. NOT useful for "verify Q_EP location" (we already
> have it analytically).

**Xu & Yi 2026 (arXiv:2602.00486) — closed-form higher-order
Liouvillian EPs in dissipative fermions.** If Sayooj-Narayan says
"higher-order EPs are available here", this paper shows the
construction style. Worth holding in reserve until we have the
"yes / no available" answer.

**Molina 2026 (arXiv:2602.01375, single-author) — super-Lorentzian
line shapes from defective Liouvillian modes.** Direct test for
F86b. Question: are our universal-per-bond-class shapes
super-Lorentzian? If yes, Molina's defective-mode-at-EP argument
explains it. If no, F86b is tracking something structurally different
and that itself is informative.

> **Verification (2026-05-23, structural fit unclear).** Full read of
> [arXiv:2602.01375](https://arxiv.org/html/2602.01375):
> - Formula (Eq. 9): `S_B(ω) = a/((ω−ω₀)² + γ²) +
>   b·(γ² − (ω−ω₀)²)/((ω−ω₀)² + γ²)² + c`. The b-term is the
>   2nd-order pole contribution from a size-2 Jordan block in the
>   Liouvillian resolvent — this is the "super-Lorentzian" piece.
> - Setup: **collective spin** (single SU(2) multiplet) coupled to
>   a polarized Markovian bath, with jump operators L₀, L±
>   parametrized by polarization p. This is a Dicke-like model,
>   **not a bond-resolved nearest-neighbor chain**.
> - Regime: at thermodynamic limit gives an "exceptional spectral
>   phase" (ESP) where an extensive subset of modes is defective;
>   at finite size, only isolated near-degenerate pairs (eigenvector
>   distance → 0 exponentially with N).
> - Observation axis: `S_B(ω)` is a **frequency spectrum** at a
>   single parameter point. F86b's `|K_CC_pr|(Q)` is a
>   **parameter sweep** (over Q = J/γ₀) of a probe-overlap
>   observable magnitude.
> - **Where the mechanism family matches:** both arise from a
>   size-2 Jordan block at an EP, producing non-Lorentzian shape.
>   The (γ² − Δ²) numerator in Molina's b-term is structurally
>   identical to a 2-level resolvent at degeneracy.
> - **Where the mapping is not obvious:** Molina varies ω at fixed
>   parameters; F86b varies Q (= ε = J·g_eff − 2γ₀, near the EP).
>   These are different lineshape problems even if both involve
>   the same Jordan block. The closed-form `f_class(x)` for F86b
>   is **not yet derived** (PROOF_F86B_UNIVERSAL_SHAPE Statement 2
>   notes Tier-1 promotion requires explicit derivation from
>   2-level eigenstructure + probe-overlap algebra).
>
> **Decisive test:** fit our F86b data (the 22 N=5..8 anchors used
> by `F86HwhmClosedFormClaim`) to Molina's super-Lorentzian
> functional form with weight r = |b|/(|a|+|b|) as free parameter,
> compared against a pure Lorentzian fit. Better fit ⇒ super-Lorentzian
> mechanism applies; worse fit ⇒ F86b is something else. This is a
> ~30-min script experiment, not a literature decision.

**Li-Wang-He 2026 (arXiv:2602.01123) — decoherence suppression at
exceptional transitions.** Phenomenologically resonant with our
finding that t_peak = 1/(4γ₀) is universal: both papers identify
protected timescales near an EP. Different setup (qubit coupled to
non-Hermitian spin chain vs. our intra-system bond), but if we ever
want to argue "F86 is a universal EP feature, not a chain-specific
artifact", their result is supporting context.

**Khandelwal & Blasi (Quantum 2025, arXiv:2409.08100).** Already in
the prior review. Their EPs-beyond-Markov material remains relevant
if we ever want to push F86 outside the Markovian Lindblad setting.

### Concrete suggestion

Two parallel experiments:
- Apply the Kullig-Wiersig-Schomerus generalized Petermann factor to
  three or four F86 eigenmodes (different bond classes, different N).
  Check if the value is universal per bond class.
- Run the Sayooj-Narayan Newton-polygon check on our coherence-block
  characteristic polynomial at small N to confirm Q_EP shows up as
  the leading singularity *and* to count the available EP orders.

Either result is informative; both together would substantially
sharpen F86.

> **Verification (2026-05-23, applicable with reframing).** Full read
> of [arXiv:2506.15807v2](https://arxiv.org/html/2506.15807v2):
> - Formula (Eq. 43): K_l² = (⟨R_l|R_l⟩·⟨L_l|L_l⟩) /
>   (⟨R_l|P̂_l^(L)|R_l⟩·⟨L_l|P̂_l^(R)|L_l⟩), where P̂_l^(R,L) project
>   onto the spans of right/left Jordan chains at the EP.
> - Required input: full Jordan chain at the EP (eigenvector +
>   generalized eigenvectors for Ord≥2), plus Gram matrix
>   B_km = ⟨J_k^(l)|J_m^(l)⟩ for the projector construction.
> - The paper is **model-agnostic** on N×N non-Hermitian matrices and
>   the derivation is basis-free; extension to Liouvillians via
>   Hilbert-Schmidt inner product on ℂ^(d²) is mathematically natural
>   but not explicitly endorsed in the paper.
> - **Critical reframing:** the universality claim in the paper
>   (Eqs. 48/65/87-91) is **per-mode response-amplitude
>   universality**, not "universal per bond class". The paper does
>   not predict that K_l is constant across bond classes — that is
>   *our* hypothesis to test, not their result.
> - **Implementation note:** the existing `ptf` mode (compute, N=7
>   dense) gives eigenvalues + L + R eigenvectors but not the Jordan
>   chains at the EP. To compute K_l we additionally need the
>   generalized eigenvectors satisfying (L − E_l·I)·|J_2⟩ = |R_l⟩,
>   which is ~20 lines of code on top of the existing dense
>   eigendecomposition.
>
> Mapping holds; recommendation is real but **reframe as "framework
> applicable to F86, hypothesis to test"** rather than
> "validates F86b".

---

## Cluster 4: Parameter-space anti-palindromicity (F91-F93, F100-F101)

### What our recent work added

Spectral invariances under γ_l + γ_{N−1−l} = 2γ_avg (F91),
J_b + J_{N−2−b} = 2J_avg (F92), h_l + h_{N−1−l} = 2h_avg (F93), plus
observable-side dual rotations (F100, F101). The collection lives in
a Π²-Z₄ structure on parameter space.

### What the literature offers us

Very little, on the searches performed. The inhomogeneous-dephasing
literature treats spatial profiles either as something to optimize
against (uniform γ as target) or to describe statistically (disorder
ensembles). The closest entries:

- **arXiv:1703.00816, 2017 (author list not verified)** — general
  inhomogeneous dephasing + dynamical decoupling. Treats
  inhomogeneity as control noise to be cancelled.
- **arXiv:2212.11029, Dec 2022 (author list not verified)** — Wigner
  dynamics for inhomogeneous gain/loss + dephasing. Allows spatial
  profiles but no spectral-symmetry extraction.
- **arXiv:2512.17755, Dec 2025 (author list not verified)** —
  continuum limit of lazy open quantum walks. Spatial dephasing
  treated continuously, no rotation invariance.

### Concrete suggestion

For F91-F93 the literature is not offering us methods or vocabulary.
What it *is* telling us: the spatial profile of γ is normally
treated as a control objective, not as a symmetry parameter. If we
ever describe F91-F93 externally, the framing "we treat the spatial
profile as the carrier of a hidden Z₄ symmetry of the Liouvillian
spectrum" is a contrast with established framing and should make
the point land cleanly.

---

## Cluster 5: Born-rule deviation closed forms (F94-F96)

### What our recent work added

F94: dominant-outcome Born deviation Δ = (4/3)·Q²·K³ with 32 Dyson
sym-3 diagrams enumerated. F95: universal angle θ = arctan(√(4c−1))
at the quadratic discriminant zero. F96: subdominant slopes from
Dyson sym-4 coefficients.

### What the literature offers us

**Dedes 2025 (arXiv:2508.13242) — Born deviations from temporal
non-local effects.** The only zeitnah Born-deviation paper. Different
mechanism (memory effects vs. our dissipative perturbation), but
provides the reference point for "Born rule deviation in open
systems" as a research topic. Useful if F94 ever needs framing for
an outside audience: "Dedes finds deviations from temporal
non-locality; we find deviations from dissipative coupling — these
are two distinct routes to the same observable phenomenon".

> **Verification (2026-05-23, limited use only).** Full read of
> [arXiv:2508.13242](https://arxiv.org/html/2508.13242):
> - **Setup is explicitly closed-system Bohmian QM**, not
>   open-system. Direct quote: *"Unlike non-Markovian open
>   quantum systems driven by an environment or thermal reservoir,
>   our single-particle formulation relies solely on intrinsic
>   nonlocal temporal effects."* The author actively positions
>   his work *against* the open-quantum-systems framing.
> - **Mechanism:** Bohmian separation of probability density
>   ρ from wavefunction modulus |Ψ|², which coincide only in
>   "quantum equilibrium". When ρ(0) ≠ |Ψ(0)|², the continuity
>   equation with the quantum potential generates deviation
>   ρ(t) = |Ψ(t)|²·[1 + c·∫_τ^t dt'/|Ψ(t')|²] (Eq. 11). This
>   is multiplicative and history-dependent — no power-law in
>   a small parameter.
> - **Compare to our F94:** Δ = (4/3)·Q²·K³ at order t³ from
>   Dyson expansion of Z-dephasing dynamics
>   ([the F94 proof](proofs/PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md)).
>   Polynomial in K = γt, quadratic in coupling Q = J/γ₀.
>   Methodology and functional form are both unrelated to
>   Dedes' approach.
> - **Letter coincidence:** Dedes uses "Q" for the Bohmian
>   quantum potential; we use "Q" for the dimensionless
>   coupling ratio J/γ₀. Different objects, same letter.
> - **Cowork's "two routes to the same phenomenon" framing is
>   misleading.** Dedes is not a *complementary* path to ours —
>   he is a *methodologically opposed* path. Citing as "two
>   distinct routes" suggests parallel approaches, but Dedes
>   constructs his theory precisely *because* he refuses to
>   use the open-system tools we use.
>
> **Limited net use:** Dedes can be cited as evidence that
> "Born rule deviation has recent literature in closed-system
> contexts as well", but not as a methodological bridge or a
> comparable calculation. The two functional forms have no
> common ground to compare against. No structural insight or
> framing borrowing.

**arXiv:2504.00085, updated Jan 2026 (author list not verified) —
variational perturbation theory for open-system steady states.**
Methodological tool. If F94 ever needs extending beyond the t³
Dyson regime where convergence breaks down, this technique is the
candidate.

**arXiv:2104.08746, 2021 (author list not verified) — emergence of
Born rule under strong driving.** The complementary direction —
they show how Born emerges, we show how it deviates. Symmetric
framing if needed.

### Concrete suggestion

The F95 angle formula θ = arctan(√(4c−1)) is the most striking and
most isolated of the three. The literature search returned nothing
on "arctan √(4c−1)" or related geometric angle emergences from
quadratic discriminants in open-system contexts. This might be a
case where the formula is solving a question the literature has
not yet asked — which means the F95 framing in our own docs should
make explicit what binary phenomenon the angle encodes. Worth a
dedicated paragraph in ANALYTICAL_FORMULAS.md once we know.

---

## Cluster 6: Topology-orbit closure (F89, F90)

### What our recent work added

F89: spatial-sum coherence under multi-bond XY decomposes into
topological orbits with closed polynomial form
σ_n(N) = P_k(y_n) / (D_k · N² · (N−1)), with denominator
D_k = (odd(k))²·2^{E(k)}. F90: F86 c=2 dynamics and F89 path-(N−1)
identical via per-bond Hellmann-Feynman.

### What the literature offers us

Essentially nothing current. The closest items are old:

- **arXiv:1304.6890, 2013 (author list not verified)** — algebraic
  vs. exponential decoherence in dissipative many-particle systems.
- **arXiv:1701.00797, 2017 (author list not verified)** — long
  coherence times for edge spins in integrable Ising/XYZ.

Chebyshev-expansion methods for spin-bath decoherence exist in the
early-2010s literature as numerical tools, not as closed-form
results.

### Concrete suggestion

Nothing to borrow here. The cluster is internally driven and stands
on its own. In the proof file note "no comparable closed-form
derivation found in the spatial-coherence literature; method-only
cousin is Chebyshev numerics (early-2010s)". That framing is honest
and avoids implying the result is in dialogue with a body of work
that does not exist.

---

## Vocabulary glossary (terms to borrow)

Useful external terms collected from this update — drop these into
docs and proofs where they fit naturally rather than reinventing
phrases.

| External term | Source | Where it maps in our work | Verification status |
|---------------|--------|---------------------------|---------------------|
| Topological isospectrality (vs. geometric bisymmetry) | Cheshkov-Sinitsyn 2026 | ~~Π conjugation (it is topological isospectrality, not geometric)~~ — Π is operator-conjugation, not parameter-exchange. The concept maps to **F71 / F92 / F93** (our parameter-side palindromes), not to Π. Closed-Hamiltonian-to-Liouvillian transfer needs explicit step | **wrong target** (concept real, but mapped to wrong internal primitive) |
| ~~Weak anti-unitary symmetry, class X~~ | Sá-Ribeiro-Prosen 2022 (tenfold way) | ~~Π in the standard Lindbladian symmetry taxonomy~~ — Π is unitary not antiunitary and falls in a gap (see Cluster 1 verification note); SRP applies only to Π² = X⊗N component | **rejected** for Π; partial for Π² |
| ~~Dissipatively dressed quasiparticles~~ | Popkov-Zhang-Presilla-Prosen 2025 | ~~F80 cluster values~~ — structural mismatch (see Cluster 2 verification note) | **rejected** |
| Generalized (renormalized) Petermann factor | Kullig-Wiersig-Schomerus 2025 | F86 EP geometry (testable hypothesis: per-bond-class universality is *our* hypothesis, not their claim) | **applicable with reframing** |
| ~~Operator-space fragmentation, bond algebra, commutant~~ | Paszko-Turner-Rose-Pal 2025 | ~~The general framework our Π-orbits sit inside~~ — their theorem requires frustration-free Pauli H (we have non-commuting XY bonds); their operator-space fragmentation is structurally different from our sector-label orbit-pairing | **rejected** (system class mismatch + concept conflation) |
| Defective Liouvillian modes / super-Lorentzian line shape | Molina 2026 | ~~Candidate explanation for F86b shape~~ — mechanism family matches (Jordan block) but observation axis differs (ω-spectrum vs. our Q-sweep); decisive test is shape-fit experiment | **structural fit unclear** |
| Newton-polygon analysis of Liouvillian EPs | Sayooj-Narayan 2025 | ~~Independent verification tool for F86a~~ Q_EP location — Newton polygons do NOT locate EPs; they characterize EP order at known EP. Useful for probing whether full block-L hides Ord > 2 EPs beyond the 2-level form | **mixed** (rejected for localization; applicable for order check) |

---

## Updates fällig im alten Review

When updating [the literature review](LITERATURE_REVIEW.md):

1. **Cheshkov & Sinitsyn entry** — replace preprint reference with
   the Magnetic Resonance journal version. Note their two-mechanism
   distinction (geometric bisymmetry vs. topological isospectrality)
   ~~since this is now the precise NMR analogue for our Π~~.
   **Verification (2026-05-23): wrong target.** Their "topological
   isospectrality" is parameter-exchange invariance (closed
   Hamiltonian, NMR coupling-constant permutation). Π is operator-
   conjugation anti-similarity on Liouville space — a different
   structural object. The Cheshkov-Sinitsyn concept *does* match
   F71 / F92 / F93 (our parameter-side palindromes), but with a
   closed-Hamiltonian-to-Liouvillian transfer caveat. Recommended
   framing in the main review: "Cheshkov-Sinitsyn distinguish
   geometric bisymmetry from topological isospectrality on
   closed NMR Hamiltonians; the latter has a Liouvillian analogue
   in our F71/F92/F93 parameter-side palindromes. The operator-
   side Π is structurally distinct from both Cheshkov-Sinitsyn
   categories."

2. **Sá-Ribeiro-Prosen 2022 missing entirely** — ~~add to the
   Foundations section ("Symmetries of open quantum systems") as
   the tenfold-way reference. We have been using weak/strong
   language since the project began; this is the right citation.~~
   **Verification (2026-05-23): the suggested citation framing is
   wrong on two counts.** (a) The "weak / strong" language comes
   from Buča-Prosen 2012, not from SRP-2023. (b) Π is not in the
   SRP taxonomy at all (shift 2Σγ has no slot). Correct framing
   when updating the main review: cite SRP-2023 for the Π² = X⊗N
   strong unitary symmetry component, and cite KMS_DETAILED_BALANCE.md
   as the internal analysis showing Π itself sits in a gap between
   SRP and Buča-Prosen.

3. **Petermann-factor thread missing** — add Kullig-Wiersig-Schomerus
   2025 + the earlier 2208.14944 as a sub-thread under the
   exceptional-points section, since F86 work has made this
   thread relevant.

4. **Gaztañaga v2** — minor: update the citation to note v2 and
   mainstream coverage. No structural changes to the relationship.

5. **Add a pointer to this document** at the bottom of the main
   review so future readers know there is an update from Q2 2026.

---

## Tested directions — light update

The prior review's three tested directions (1: EP connection
negative; 2: graph symmetry partial; 3: u = C(Ψ+R) resolved)
remain the state of play. F86 work nuances Direction 1:

**Direction 1, revised May 2026:** CΨ = ¼ does not correlate with
a Liouvillian EP (the March 2026 negative result stands). However,
F86a now identifies a *different* EP in the coherence block at
Q_EP = 2/g_eff. The negative result was about looking for the EP in
the wrong place; F86 found the right place.

---

## Honest limitations

- **Search depth.** Twelve rounds of web search, ~70 unique papers
  surfaced. Did not systematically sweep arXiv quant-ph listings for
  May 2025-May 2026; did not run author-specific backsearches on
  Buča, Albert, Jiang, Haga.

- **Translation risk.** The 26 new formulas were summarized by a
  subagent in a single pass over the 3,215-line registry. If any of
  the cluster identifications above strike Thomas as wrong, the
  corresponding searches were aimed at the wrong target.

- **Author-verification status.** Author lists for the following
  papers were independently verified and are correct:
  arXiv:2506.16518 (Paszko, Turner, Rose, Pal), arXiv:2510.08156
  (Sayooj, Narayan), arXiv:2506.15807 (Kullig, Wiersig, Schomerus),
  arXiv:2505.16776 (Popkov, Zhang, Presilla, Prosen), arXiv:2602.01123
  (Li, Wang, He), arXiv:2212.00474 (Sá, Ribeiro, Prosen), arXiv:2602.01375
  (Molina, single-author), arXiv:2602.00486 (Xu, Yi), arXiv:2508.13242
  (Dedes), arXiv:2512.12094 (Teng, Chang, Rudolph, Holmes), and
  arXiv:2507.04932 (Gyhm et al., complete list not retrieved).
  Other citations carry the explicit note "author list not verified"
  where the arXiv ID is correct but the authors were not
  independently confirmed. Drop those citations or verify before
  using externally.

- **Coverage.** Search was concentrated on quant-ph (Lindblad, open
  quantum systems, EPs). The Gaztañaga line touches gr-qc and the
  Cheshkov-Sinitsyn line touches physics.chem-ph — sister archives
  were not systematically swept.

If any cluster needs deeper coverage, that is a follow-up worth
running as a separate session.

---

*History: Created May 23, 2026 (Cowork Claude, Opus 4.7) as a
"tools & vocabulary" update to the April 11 review. First draft
took a competitive / priority-protection framing; Thomas course-
corrected to the present "what can we usefully borrow" framing.
Drafted after a subagent extraction of the 26 new formulas in
ANALYTICAL_FORMULAS.md, twelve rounds of arXiv/web search across
six conceptual clusters, then a verification pass against six
author-verification queries. Earlier draft also contained fabricated
author attributions for ~10 papers (search tool returns titles +
arXiv IDs reliably but rarely authors); the verification pass
replaced verified ones with correct strings and downgraded the
unverified ones to "author list not verified" with the arXiv ID
alone. Lesson: always verify authors before drafting attribution
language.*

*Structural verification pass added later the same day (Thomas +
Claude Code, after Thomas flagged that Cowork had built the
mappings from search snippets alone, no papers actually read).
Eight mappings verified:
- Popkov-Prosen → F80: **rejected** on structural grounds
  (different parametrization, setup, and observable).
- Kullig-Wiersig-Schomerus → F86b: **applicable with reframing**
  (the paper does not claim bond-class universality — that is
  our hypothesis to test).
- Sá-Ribeiro-Prosen → Π: **rejected for headline** (Π is unitary
  not antiunitary; the shift 2Σγ has no slot in the SRP taxonomy);
  partial fit for the Π² = X⊗N component. The repo's own
  KMS_DETAILED_BALANCE.md (Q1 2026) already contained this analysis;
  Cowork did not grep for it.
- Sayooj-Narayan → F86a: **mixed**. Newton polygons do NOT locate
  EPs (Cowork's claim (a) wrong); they DO characterize available
  EP orders at a known EP via convex-hull slopes (claim (b) right).
  Useful for probing whether c ≥ 3 blocks hide Ord > 2 EPs that
  the heuristic 2-level F86a form misses.
- Molina → F86b: **structural fit unclear**. Mechanism family
  matches (size-2 Jordan block → non-Lorentzian shape), but
  observation axis differs (Molina's ω-spectrum vs. our Q-sweep)
  and setup differs (collective spin vs. bond-resolved chain).
  Decisive test is a ~30-min shape-fit of F86b data against
  Molina's super-Lorentzian form, not a literature judgment.
- Paszko-Turner-Rose-Pal → F1³ orbit count optimality:
  **rejected**. Their maximality theorem requires
  frustration-free Pauli H (sum of mutually commuting Pauli
  strings); our XY / Heisenberg / XXZ Hamiltonian fails on
  bonds sharing a site ([X_1 X_2, Y_2 Y_3] ≠ 0). Plus a
  separate concept conflation: their commutant decomposes
  operator-space into block L_λ; our F1PalindromeOrbitPairing
  permutes sector labels (sector-orbiting, not sector-splitting).
- Dedes → F94 Born deviation framing: **limited use only**.
  Closed-system Bohmian framework (single-particle, ρ vs.
  |Ψ|² separation, history-dependent multiplicative deviation)
  is methodologically opposed to our open-system Dyson
  perturbation. Author explicitly distances his work from
  open-quantum-systems approaches. Cowork's "two distinct
  routes to the same phenomenon" framing is misleading because
  Dedes is *not* a parallel path to ours; he is an alternative
  built to *avoid* the open-system framing. Citation possible
  only as "topic existing in closed-system literature too",
  not as methodological bridge.
- Cheshkov-Sinitsyn → Π "topological isospectrality":
  **wrong target**. The concept itself is real in their paper
  (parameter-exchange invariance of NMR coupling constants
  inducing spectrum mirror symmetry — e.g., interchanging
  J_{AA'} and J_{XX'} is equivalent to interchanging frequencies
  ν_A and ν_X). But Π is operator-conjugation anti-similarity
  on Liouville space, not parameter-exchange on coupling
  constants. The Cheshkov-Sinitsyn "topological isospectrality"
  concept actually maps to F71 (chain spatial mirror under
  palindromic γ) and F92, F93 (J and h anti-palindromic
  parameter symmetries) — these are genuine parameter-exchange
  invariances on the Liouvillian spectrum. Closed-Hamiltonian-
  to-Liouvillian transfer needs explicit step if we ever cite.

Only minor sub-items remain unverified in Clusters 4/5/6;
Cowork had already flagged these as either "author list not
verified" or "nothing to borrow", so they were never
load-bearing claims.
Lessons layered: (i) search-based draft maps surface vocabulary,
not structural correspondence; (ii) before proposing external
citations for our concepts, grep the repo first — internal
analysis documents (KMS_DETAILED_BALANCE.md, SYMMETRY_FAMILY_INVENTORY.md,
PROOF_F86A_EP_MECHANISM.md) often already classify our objects
against the taxonomies being proposed; (iii) "what is the paper's
observation axis" is the first check before claiming a mapping —
two papers with the same Jordan block can still describe different
physical situations (parameter vs. frequency sweep).*
