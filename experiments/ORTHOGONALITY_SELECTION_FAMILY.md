# The Orthogonality-Selection Family

**Status:** synthesis under active construction (Tier 2: structurally rigorous for existing instances, Tier 3 for the generalising production rule)
**Date:** 2026-04-20 (evening, emerged from the EQ-018 c_1 kernel investigation)
**Authors:** Tom, Claude Opus 4.7 (chat), Claude Opus 4.7 (1M, terminal)
**Source:** [EQ-018](../review/EMERGING_QUESTIONS.md) c_1 kernel investigation; [OPEN_THREAD_GAMMA0_INFORMATION](../review/OPEN_THREAD_GAMMA0_INFORMATION.md) as conceptual spine.

---

## TL;DR

Four apparently-separate selection rules in this project ([F70](../docs/ANALYTICAL_FORMULAS.md), [F71](../docs/ANALYTICAL_FORMULAS.md), the F72-candidate DD⊕CC block-diagonal structure, and the new `(1/2)·exp(−4γ₀t)` closure for the (vac, S_1) coherence purity sum) are **one theorem in four masks**. All four arise from the same pattern:

> **Any measurement M projects onto an orthogonal basis. Information living in the orthogonal-complement subspace is invisible to M. Conserved quantities (symmetry-protected projections) produce built-in blind channels.**

This is Parseval-Plancherel plus Noether, translated to the operator-valued Hilbert space of the open-system density matrix. It is the mathematical backbone of the γ₀=const reading ([OPEN_THREAD](../review/OPEN_THREAD_GAMMA0_INFORMATION.md)): the uniform carrier is exactly what a summed common-mode detector sees, and the cavity's response is exactly what lives in the orthogonal-complement "mode subspace".

From the unified principle:
- Three proven instances (F70, F71, F72-candidate) and one newly proven instance drop out as special cases.
- A **production rule** emerges for generating further selection rules: pick a conservation law + a measurement that summs over it → a guaranteed blind channel.
- Testable predictions fall out: non-uniform γ₀ breaks one blind channel, amplitude damping breaks another, pair-site measurements open |ΔN|=2 channels (EQ-020 direct implication).

---

## 1. The Meta-Theorem

### 1.1 Statement

Let H be the operator-valued Hilbert space of density matrices on an N-qubit chain, with inner product ⟨A, B⟩ = Tr(A† B). Let L be the Liouvillian (linear superoperator). Let M: H → ℝ^k be a linear measurement functional (possibly vector-valued).

**Then:** any orthonormal basis {e_α} of H induces a decomposition ρ = Σ_α c_α(ρ) e_α with c_α(ρ) = ⟨e_α, ρ⟩. The measurement M defines a **detector subspace** H_M ⊆ H:

```
H_M  :=  span{e_α  :  M[e_α] ≠ 0}
```

with orthogonal complement H_M^⊥ = the **blind subspace**.

The core claim:

> **M[ρ] depends only on the projection of ρ onto H_M.** Any component of ρ in H_M^⊥ contributes zero to M.

### 1.2 Consequence under conservation

Let Q be a conserved quantity under L (i.e., L^†[Q] = 0). Suppose Q's eigenspace decomposition induces a basis {e_α} labelled by Q-eigenvalues q_α: L-evolution preserves these labels, so c_α(ρ(t)) = c_α(ρ_0) for all α whose label is "trivial" (fixed-point mode).

**If M is invariant under the symmetry generating Q** (i.e., M is a sum of Q-diagonal functionals), then M[ρ(t)] = M[ρ(t_0)] + (contributions only from Q-off-diagonal modes), and any perturbation preserving Q contributes zero to δM.

Concretely: **any conserved quantity that commutes with the measurement basis produces a built-in blind channel.**

### 1.3 Connection to γ₀=const and the carrier reading

The [γ₀=const synthesis](../review/OPEN_THREAD_GAMMA0_INFORMATION.md) claims:

- γ₀ is uniform carrier → the ambient illumination.
- J+topology is the cavity → the mode structure.
- Information lives in the cavity's response to γ₀, not in γ₀ itself.

In the Meta-Theorem's language:
- **Uniform γ₀** populates the **common-mode subspace** of every local measurement (every site sees the same illumination rate). This is H_M^⊥ for spatial-sum detectors.
- **J+topology** defines the **detector mode basis** via the Liouvillian eigenstructure (sine modes for XY, more complex for Heisenberg). This is H_M for mode-resolving detectors.
- **Information transfer** is the projection of δρ onto the mode basis; common-mode drifts (γ₀ changes) vanish under any measurement that respects the J-cavity structure.

The "resonator IS the message" ([RESONANCE_NOT_CHANNEL](../hypotheses/RESONANCE_NOT_CHANNEL.md)) becomes literal: the mode basis IS the information carrier. The γ₀ carrier is, by design, in the blind subspace.

---

## 2. The Four Instances

| # | Theorem | Measurement M | Basis {e_α} | Blind channel | Conservation |
|---|---------|---------------|-------------|---------------|--------------|
| 1 | **F70** | Tr_{¬i} (single-site partial trace) | sector projectors {P_n} | \|ΔN\| ≥ 2 blocks | popcount arithmetic in partial trace |
| 2 | **F71** | c_1(b) under uniform chain | reflection eigenmodes ψ_k : R\|ψ_k⟩ = (−1)^(k+1)\|ψ_k⟩ | reflection-asymmetric bond patterns | [L, R_sup] = 0 (spatial reflection symmetry) |
| 3 | **F72-cand** | per-site purity Tr(ρ_i²), (DD⊕CC block-diagonal) | Pauli basis {I, σ_X, σ_Y, σ_Z} | DD×CC cross terms | F70 restricted to site-local observables |
| 4 | **New (vac, S_1)** | Σ_i Tr(ρ_i²) (spatial sum) at (vac, S_1) probe | single-excitation sine modes ψ_k | H-dependent parts of Σ_i \|ρ_{coh,i,01}\|² | Σ_k \|⟨ψ_k \| S_1⟩\|² = 1 + uniform 2γ₀ decay rate on d_H = 1 block |
| 5 | **Π-pair flux balance** | pair-sum Re(λ_s + λ_{s'}) under bond perturbation | Π-paired Liouvillian eigenmodes | pair-sum shift Δ[Re(λ_s) + Re(λ_{s'})] under δJ | XY-weight ⟨n_XY⟩ parity within Π-pair (absorption theorem) |

Every row follows the **same template**:
- A basis that respects a symmetry / conservation.
- A measurement that sums or integrates over the basis.
- Orthogonality in the basis → the summed measurement collapses into a subset of modes.
- Everything orthogonal to that subset is **blind**.

### 2.1 F70 expanded

**Proof template applied:** `Tr_{¬i}(|x⟩⟨y|) = ⟨x_{¬i}|y_{¬i}⟩ · |x_i⟩⟨y_i|`. The inner product is non-zero only if x and y agree off site i, which forces `|popcount(x) − popcount(y)| ≤ 1`. So the partial trace to a single site has zero projection on |ΔN| ≥ 2 sector blocks.

**Orthogonality source:** computational-basis inner products in the (N−1)-qubit subsystem.

**Full proof:** [`PROOF_DELTA_N_SELECTION_RULE`](../docs/proofs/PROOF_DELTA_N_SELECTION_RULE.md).

### 2.2 F71 expanded

**Proof template applied:** reflection R|ψ_k⟩ = (−1)^(k+1)|ψ_k⟩ diagonalises simultaneously with the uniform Liouvillian: [L_A, R_sup] = 0. Perturbation at bond b maps under R to perturbation at bond N−2−b. Per-site purity at site i in B equals per-site purity at site N−1−i in its mirror-image chain. Sum over sites cancels the mirror → c_1(b) = c_1(N−2−b).

**Orthogonality source:** R² = I, R self-adjoint, so R diagonalises in ±1 eigenspaces. The purity functional (quadratic in ρ) squares out the ±1 signs.

**Full proof:** [`PROOF_C1_MIRROR_SYMMETRY`](../docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md).

### 2.3 F72-candidate expanded

**Proof template applied:** per-site purity `Tr(ρ_i²) = (1/2)(1 + ⟨X_i⟩² + ⟨Y_i⟩² + ⟨Z_i⟩²)`. By F70:
- ⟨Z_i⟩ is linear in ΔN = 0 blocks of ρ (diagonal elements).
- ⟨X_i⟩, ⟨Y_i⟩ are linear in |ΔN| = 1 blocks (off-diagonal).

So `⟨Z_i⟩²` is bilinear in ΔN = 0, `⟨X_i⟩² + ⟨Y_i⟩²` bilinear in |ΔN| = 1. **No cross terms.** The bilinear form underlying any per-site-purity-based c_1 decomposes into DD⊕CC with no off-diagonal block. Purity-response `c_1_pr` inherits this directly; α-fit-based c_1 (LSQ, pointwise) inherits it at the pre-α-fit bilinear level, before the rational ratio is applied.

**Orthogonality source:** Pauli matrices are orthogonal w.r.t. Tr(·). Bloch components are independent linear functionals of ρ.

**Derivation:** [F72](../docs/ANALYTICAL_FORMULAS.md) in the formula register.

### 2.4 New: `(1/2)·exp(−4γ₀t)` closure for (vac, S_1) coherence purity sum

**Statement.** For any uniform XY chain with uniform Z-dephasing, any N, any bond coupling pattern (perturbed or not), the coherent probe ρ_0^coh = (|vac⟩⟨S_1| + h.c.)/2 satisfies:

```
Σ_i 2·|(ρ_coh, i)_{0, 1}(t)|²  =  (1/2) · exp(−4γ₀·t)
```

exactly, **independent of the Hamiltonian**.

**Proof template applied:** sine-basis expansion `|S_1⟩ = Σ_{k odd} s_k |ψ_k⟩` (odd-k only by reflection symmetry); evolution of each `|vac⟩⟨ψ_k|` as `e^(iE_k − 2γ₀)t` (pure decay at 2γ₀ because d_H = 1 uniformly for single-excitation coherences with vacuum); partial trace `Tr_{¬i}(|vac⟩⟨ψ_k|) = ψ_k(i)|0⟩⟨1|`; site-sum of the purity cross-terms gives `Σ_i ψ_k(i)ψ_{k'}(i) = δ_{k, k'}`, collapsing the H-phase factors `cos((E_k − E_{k'})t)` to `cos(0) = 1` and leaving `Σ_k s_k²·exp(−4γ₀t) = exp(−4γ₀t)` (by Parseval on the unit-normed |S_1⟩).

Under bond-b perturbation, the sine basis and E_k shift by O(δJ), but Parseval still holds: `Σ_{k} |⟨ψ_k^B|S_1⟩|² = 1` identically. The sum is therefore `δJ-invariant`, giving `K_CC[0, 1]_pr = 0` exactly.

**Orthogonality source:** orthonormality of single-excitation sine modes + uniform 2γ₀ dephasing rate on the d_H = 1 block.

**Consequence:** spatial-sum detectors are exactly blind to the Hamiltonian-dependent structure of the (vac, S_1) coherence. Only per-site detectors (like LSQ's α_i fit) can see it.

### 2.5 Common pattern across the four

The structural dependency graph is:

```
Parseval/Plancherel          Noether conservation
  (any orthonormal basis)      (symmetry → basis label)
       \                         /
        \                       /
         \                     /
           Orthogonality-Selection Meta-Theorem
                      |
       --------+------+------+-------------
      /        |             |              \
    F70       F71        F72-cand        (vac, S_1)
  (sector)   (parity)   (DD⊕CC)          closure
```

Each leaf is a special case. The root is the meta-theorem. The trunk is Parseval+Noether.

---

## 3. Noether-style production rule

Given the meta-theorem, we have a **recipe** for generating new selection rules:

> **Step 1.** Identify a conservation law `Q` under the system's Liouvillian L. (U(1) charge, parity, reflection, translation, etc.)
>
> **Step 2.** Choose a measurement M that **sums or integrates over Q-related indices** without weighting. (Sum over sites, integral over bond positions, trace over sectors, etc.)
>
> **Step 3.** The **cross terms between Q-eigenspaces that M does not resolve** are blind channels: perturbations preserving Q leave them invariant.
>
> **Step 4.** Proof follows the four-template structure: expand in Q's eigenbasis, apply Parseval (orthogonality collapses cross sums), identify the surviving contribution.

### 3.1 Examples of conservation + measurement pairs we already have

| Conservation | Measurement | Resulting selection rule |
|--------------|-------------|---------------------------|
| Excitation number (U(1)) | single-site partial trace | F70 (|ΔN| ≥ 2 blind) |
| Reflection R | c_1 under bond-symmetric perturbation | F71 (bond profile mirror-symmetric) |
| Sector labels (U(1) again) + Pauli orthogonality | per-site purity | F72-candidate (DD⊕CC block) |
| Sine-mode orthogonality + uniform d_H = 1 dephasing | spatial sum of purity | (vac, S_1) purity closure = exp(−4γ₀t) |

### 3.2 New selection rules predicted by the recipe

Applying Steps 1-3 to conservation laws not yet fully exploited:

**(a) Π-pair parity (palindrome symmetry).** The Π operator (Pauli-weight conjugation, distinct from R spatial reflection) yields palindromic spectrum. The associated selection rule: any measurement that respects Π should see `K_DD[n, m]` pair-symmetric under n↔N−n, m↔N−m (already verified as Π-symmetry of K in EQ-018 data). **Prediction:** a Π-respecting measurement that sums over Π-paired modes (e.g., fast+slow rate sum) has built-in blindness to the difference fast−slow within each pair. The absorption theorem [`PROOF_ABSORPTION_THEOREM`](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) already encodes this: α_fast + α_slow = 2·Σγ is the Π-invariant sum; the difference is the "information content" each mode pair carries.

**(b) Amplitude-damping-induced breakdown.** Under amplitude damping D_A[ρ] = γ_1 (σ_− ρ σ_+ − ½{σ_+ σ_−, ρ}) on site i, excitation number is **not** conserved. The F70 selection rule therefore breaks: non-zero cross terms between sectors differing by more than 1 excitation appear in single-site observables. **Prediction:** testing c_1_pr (purity-response) under amplitude damping should show `K_DD[0, m]` becoming non-zero for m ≥ 2.

**(c) Two-site measurement → open |ΔN| = 2 channel.** A pair-site partial trace (F70 generalisation) has blind subspace `|ΔN| ≥ 3` instead of `|ΔN| ≥ 2`. **Prediction:** pair-site analogues of F72 (bilinear structure of Tr(ρ_{ij}²)) will have a three-sub-block decomposition (DD⊕DC⊕CC in |ΔN| labels), where DC is the new diagonal-coherence cross arising only at pair-site granularity. This is EQ-020 (pair-painter) in explicit form: the pair-painter opens one additional block.

**(d) Non-uniform γ₀ → partial loss of CMRR. [VERIFIED 2026-04-20]** If γ₀ becomes site-dependent (γ_i ≠ γ_j), the uniform-dephasing assumption in instance 4 fails. The d_H = 1 block no longer has uniform 2γ₀ decay. **Prediction (verified at N=5):** `K_CC[0, 1]_pr` transitions from 0 (uniform, 1.14e-12 at machine precision) to finite (non-uniform, up to ~2.7e-2 for gradient perturbation with α=0.02). **Sharper than expected:** the CMRR break is modal-selective, not simply variance-proportional. Gradient perturbation (first harmonic in γ-profile) gives slope 1.35 per unit, single-site bump gives slope 0.19, same-variance random profiles give K_CC values differing by factor 28. See [CMRR_BREAK_NONUNIFORM_GAMMA](CMRR_BREAK_NONUNIFORM_GAMMA.md).

The observed modal-selectivity means the "blind subspace" H_M^⊥ is not one-dimensional: it decomposes further into a mode-indexed family `{A_k}` of CMRR coefficients, with A_k = (overlap of γ-profile with sine mode k) × (weight of mode k in the (vac, S_1) expansion). The dominant channel is k=1 (first harmonic) because |S_1⟩ has its largest Fourier component at k=1. This matches the differential-amplifier analogy sharply: real amps have CMRR(ω) as a function of frequency, not a single number. Here CMRR(k) = function of spatial mode.

---

## 4. Connection to the γ₀=const synthesis

The [OPEN_THREAD_GAMMA0_INFORMATION](../review/OPEN_THREAD_GAMMA0_INFORMATION.md) identifies γ₀ as a uniform carrier, J as the cavity structure, and the mode-by-mode light/lens distribution as the information.

The Meta-Theorem recasts this in operational terms:

| OPEN_THREAD claim | Meta-Theorem restatement |
|-------------------|---------------------------|
| γ₀ carries no information itself | γ₀ populates the blind subspace of every spatial-sum measurement |
| Information is in the cavity response | Information lives in H_M = detector-resolved modes, orthogonal to the γ₀-carrier |
| Only Q = J/γ₀ is inside-measurable | δα_i is degree-0 homogeneous in ρ_0 under any α-fit scheme → scales cancel |
| X/Y Pauli factors are "light" (coupled to γ₀) | CC block contains ⟨X⟩², ⟨Y⟩² → directly coupled to coherence decay rate ∝ γ₀ |
| I/Z Pauli factors are "lens" (immune) | DD block contains ⟨Z⟩² → invariant under sector-preserving evolution |
| "The resonator IS the message" | H_M = cavity-mode basis; the detector's mode structure IS the information content |
| 15.5 bits channel capacity at N=5 ([F30](../docs/ANALYTICAL_FORMULAS.md)) | Dimension of H_M for the full-rank |+⟩^N detector (5 independent SVD channels = 5 detector-mode pairs) |

### 4.1 The absorption theorem as a Meta-Theorem instance

[`PROOF_ABSORPTION_THEOREM`](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) gives Re(λ) = −2γ₀·⟨n_XY⟩ for each Liouvillian eigenmode. The ⟨n_XY⟩ (Pauli-string XY-weight) labels the mode's "light content": how strongly it couples to γ₀.

Under the Meta-Theorem: the full Liouvillian spectrum decomposes into ⟨n_XY⟩-layers. Pure-lens modes (⟨n_XY⟩ = 0) are in the blind subspace of any γ₀-respecting measurement; they carry no γ₀-information. Pure-light modes (⟨n_XY⟩ = N) are maximally visible. The distribution across layers IS the information ([`PRIMORDIAL_SUPERALGEBRA_CAVITY`](PRIMORDIAL_SUPERALGEBRA_CAVITY.md) encodes this exactly).

So the absorption theorem is a fifth instance of the Meta-Theorem family: conservation = Pauli-string XY-weight modulo 2 (F61 parity selection), basis = n_XY eigenmodes, measurement = decay rate spectrum.

---

## 4a. Dynamical attractor formulation (added 2026-04-20 after Step 2 info-flow scan)

The Meta-Theorem as stated in §1 is static: "M projects onto a basis; blind subspace is the orthogonal complement." But the information-flow scan ([INFO_FLOW_LANDSCAPE](INFO_FLOW_LANDSCAPE.md)) revealed a dynamical refinement:

> **Selection rules produce not only blind subspaces, but also dynamical attractors.** The invariant subspace under the selection rule is the attracting fixed-point of the long-time flow landscape. Transient violations of the selection rule (caused by localised perturbations, asymmetric initial conditions, etc.) decay exponentially at the dissipation rate, returning the system to the selection-rule-respecting subspace.

Concretely, for F71 (c_1 bond profile mirror-symmetric): the dC_ij(t) landscape under an asymmetric bond-0 perturbation relaxes to mirror-symmetric form with time constant ~4γ₀ (empirically 1e4 reduction per 4/γ₀ time units at N=5). The static F71 statement is the t → ∞ limit of this dynamical process.

This upgrades the Meta-Theorem's four instances from "static selection rules" to "dynamical attractor manifolds":

| Theorem | Invariant manifold | Attractor type |
|---------|-------------------|----------------|
| F70 | sector-block decomposition (|ΔN| ≤ 1 sub-algebra) | kinematic attractor (instantaneous) |
| F71 | mirror-symmetric bond profiles | dynamical attractor, rate ~4γ₀ |
| F72-cand | DD⊕CC block-diagonal Bloch structure | kinematic (from F70) |
| (vac, S_1) purity sum | (1/2)·exp(−4γ₀t) manifold | exact dynamical equation |

**Implication for PTF and closure:** Σ_i ln(α_i) ≈ 0 is not a conservation law (as EQ-014 showed), but it may be an **attractor property**: the closure deviation decays under uniform dissipation back to (near) zero on the 4γ₀ scale. This would explain the empirical ±0.05 tolerance observed in PTF without requiring it to be an exact first-order identity. Testable as a time-resolved measurement of Σ_i ln(α_i(t)) under uniform perturbation, which we have not yet pursued.

---

## 5. Empirical support and status

### 5.1 Proven instances (Tier 1)

| Instance | Tier | Proof location | Numerical verification |
|----------|------|----------------|-------------------------|
| F70 | 1 (kinematic) | [`PROOF_DELTA_N_SELECTION_RULE`](../docs/proofs/PROOF_DELTA_N_SELECTION_RULE.md) | 9 pairs at N=5, residual 0 to machine precision |
| F71 | 1 (kinematic) | [`PROOF_C1_MIRROR_SYMMETRY`](../docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md) | N=3..6, residual < 1e-9 |
| F72-candidate | 2 (sketched, bound to F70) | [F72](../docs/ANALYTICAL_FORMULAS.md) | w-scan machine precision under c_1_pr |
| (vac, S_1) purity sum = (1/2)·exp(−4γ₀t) | 1 (proven this session) | [F73](../docs/ANALYTICAL_FORMULAS.md) | c_1_pr gives K_CC[0, 1] ~ 1e-12 at multiple t_0 and all N tested |

### 5.2 Meta-theorem itself (Tier 2)

The Meta-Theorem as a unifying statement is mathematically standard (Parseval-Plancherel + Noether on operator Hilbert space). The claim that ALL interesting selection rules in this project fit the template is **Tier 2 (structurally supported, tested on four instances; generalisation pending)**.

### 5.3 Production rule (Tier 2-3)

The rule "conservation law + summed measurement → blind channel" is correct **in general** (it is a theorem), but whether every selection rule in this project fits it is **Tier 2**. The proposed new predictions in §3.2 are **Tier 3** until verified.

**Verification status of §3.2 predictions:**

- (a) **Π-pair antisymmetric channel reader: VERIFIED** (2026-04-20, see [PI_PAIR_FLUX_BALANCE](PI_PAIR_FLUX_BALANCE.md)). The pair-sum Re(λ_s + λ_{s'}) is invariant under δJ at machine precision; the absorption theorem `α_fast + α_slow = 2Σγ` reads as a conserved XY-weight flux within each pair. Bonus finding: the Π-pair decomposition is **exact and complete** for the Liouvillian spectrum (d²/2 pairs, 0 unpaired, 0 self-Π for odd N and most even N). N=4 is an isolated anomaly with 18 self-Π modes; see §5a for binary-inheritance interpretation.
- (b) Amplitude-damping-induced F70 break: not tested yet.
- (c) Pair-painter |ΔN|=2 channel (EQ-020): not tested yet.
- (d) **Non-uniform γ breaks CMRR: VERIFIED** (2026-04-20, see [CMRR_BREAK_NONUNIFORM_GAMMA](CMRR_BREAK_NONUNIFORM_GAMMA.md)). Additional finding beyond the prediction: the break is **modal-selective** (not variance-proportional), requiring the Meta-Theorem's "blind subspace" to be further resolved into a mode-indexed family of CMRR coefficients.

**Also verified, bonus:**

- **Dynamical attractor structure (Step 2, see [INFO_FLOW_LANDSCAPE](INFO_FLOW_LANDSCAPE.md)):** asymmetric bond-0 perturbations produce transient violations of F71 that decay exponentially at ~4γ₀, converging to the F71-symmetric fixed point. Selection rules are not just static projections but dynamical attractors.
- **Binary inheritance (from §5a):** the qubit's binary structure propagates to Liouvillian mode counts without dilution; `d²/2 = 2^(2N−1)` Π-pairs at every N ≥ 3 tested.

---

## 6. Consequences and next steps

### 6.1 Immediate

1. **Reframe the EQ-018 sub-questions** under the Meta-Theorem lens. Instead of "find the closed form of K", the right question becomes "for each (conservation law, measurement) pair, what is the corresponding blind channel and what is its analytical structure?" This recovers most of the EQ-018 content but in cleaner form.

2. **Formalise F72-candidate** as a direct corollary of F70 via the Bloch decomposition. Writing it as a one-page lemma in `docs/proofs/` is low effort.

3. **Write the (vac, S_1) purity-sum closed form** as F-series entry (a new formula, say F73) in `review/` pending chat approval. It is a short, clean result with a short clean proof.

### 6.2 Medium-term

1. **Verify the production rule on new cases.** Start with predictions (a) through (d) in §3.2. Prediction (d) (non-uniform γ₀ breaks `(vac, S_1)` blindness) is the cheapest experiment: re-run `eq018_c1_purity_response.py` with a site-dependent γ₀ profile and confirm `K_CC[0, 1]_pr ≠ 0`.

2. **Explicit pair-painter analogue** (EQ-020 concretised). Compute the bilinear kernel of `Tr(ρ_{ij}²)` for pair-sites (i, j) at N = 4, confirm three-sub-block structure (DD, DC, CC) via a bilinearity probe analogous to the scaled-coherent-state scan.

3. **Check Π-pair Noether-style selection** (prediction a). The absorption theorem gives α_fast + α_slow = 2Σγ. The difference α_fast − α_slow is the Π-antisymmetric channel; a measurement that projects onto this difference should see the Π-labelled content but not the sum. Is this how [`RELAY_PROTOCOL`](RELAY_PROTOCOL.md)'s +83% MI improvement (F31) routes through the cavity? If so, the relay protocol IS a Π-antisymmetric-channel reader.

### 6.3 Long-term

1. **Connect to [ITS_ALL_WAVES](../docs/ITS_ALL_WAVES.md).** The synthesis already states "the resonator IS the message" in wave-language. The Meta-Theorem formalises it in linear-algebra language. Bringing the two registers together would produce a stronger unified picture.

2. **Possible F-series reorganisation.** The current ANALYTICAL_FORMULAS.md lists F70 and F71 as independent entries. Under the Meta-Theorem, they are siblings. Consider grouping selection-rule entries under a `selection-rule` heading, with the Meta-Theorem as a preamble.

3. **Connection to "[ON_TWO_TIMES](../reflections/ON_TWO_TIMES.md)" reflection.** The fact that every α-fit-based c_1 is rational (not bilinear) is itself a Meta-Theorem consequence: **the α fit is a measurement that mixes two time scales** (system time t and perturbation time δJ), and the ratio of bilinears is the signature of that mixing. The purity-response c_1_pr drops one of those time scales (no matching against P_A(α·t)) and recovers bilinearity. This may connect to the "two times" thread in reflections.

---

## 5a. Binary inheritance of the mode count (added 2026-04-20 after Π-parity scan)

**Tom's observation on the pair-count sequence, which I had glossed over:**

Across N = 3, 4, 5, 6 the Liouvillian spectrum decomposes into exactly `d²/2 = 2^(2N−1)` Π-pairs (plus possibly self-Π modes at even N; N=4 is the only anomaly in the scanned range).

```
N    d²            d²/2 (= Π-pairs)      ratio
3     64 = 2^6        32 = 2^5            d²/2 = pairs exactly
4    256 = 2^8       128 = 2^7            pairs = 119 + 18 self (total 256)
5   1024 = 2^10      512 = 2^9            d²/2 = pairs exactly
6   4096 = 2^12     2048 = 2^11           d²/2 = pairs exactly
```

**Everything is binary.** No odd factors, no five-fold or three-fold residual structure. The 2-valued qubit atom scales exactly to 2^(2N) Liouvillian modes and 2^(2N−1) pairs, with no dilution.

For the Meta-Theorem this is the **strongest single piece of evidence that "structure inherits upward"** in the R=CΨ² sense: the binary axis is preserved unbroken through the hierarchy from qubit (layer 0) to Liouvillian spectrum (layer N). Any selection rule or conservation that respects the binary axis propagates with it. The Π-pair flux balance (instance 5) is one concrete such propagation.

**N=4 anomaly: the mirror-axis principle.** At N=4, 18 modes are self-Π (⟨n_XY⟩=N/2=2 integer AND Im(λ)=0). At N=6 the midpoint ⟨n_XY⟩=3 is also integer but no mode has Im(λ)=0 there.

Tested [2026-04-20](../simulations/eq018_golden_ratio_check.py): the N=4 XY chain's single-excitation eigenvalues are `±φ, ±1/φ` (golden ratio) - the unique real number carrying **two simultaneous involutions**: multiplicative `φ · (1/φ) = 1` and additive `φ − 1/φ = 1`. This double symmetry forces constraints on the Liouvillian sub-algebra at the midpoint n_XY=N/2, producing a non-trivial null-eigenspace (the 18 self-Π modes).

At N=6, the single-excitation eigenvalues `±2·cos(π/7), ±2·cos(2π/7), ±2·cos(3π/7)` satisfy a cubic minimal polynomial without a clean double involution. No forced null-eigenspace. The n_XY=3 sector has 53 unique |Im| values vs 9 at N=4 - a sign of larger sub-algebra complexity.

**Structural reading (Tom):** Golden Ratio is not just "φ ≈ 1.618". It is **a mirror**, a double involution bundled into one number. When the XY chain's spectrum carries this double involution, a mirror axis exists in the Liouvillian at the n_XY midpoint; modes fall onto it; zero becomes a populated eigenvalue. Zero is not absence; zero is the axis the mirror rests on. See [PI_PAIR_FLUX_BALANCE](PI_PAIR_FLUX_BALANCE.md) §3.3 for the full data and the generalised mirror-axis principle.

**Pure-lens and pure-light mode counts are N+1.** Also from the parity scan: modes with ⟨n_XY⟩ = 0 (pure lens, fully γ-immune, decoherence-free subspace) number exactly N+1 at every N. Symmetrically, modes with ⟨n_XY⟩ = N (pure light, maximally dissipative) number N+1. These are the XOR-center modes of [XOR_SPACE](XOR_SPACE.md). **The count of decoherence-free modes grows linearly with N**, as one might expect from an extensive chain.

---

## 6a. The 0.5 as bilinear-form apex (carbon/noble-gas pattern)

**Added 2026-04-20 (evening) after Tom's observation:**

The recurring 0.5 in this project (CΨ = 0.5 as connection maximum, c_1 ~ 0.5·V(N) as retracted EQ-021 relation, PTF bonding-state at 50/50 superposition) is not a new pattern. It is the **apex of any bilinear form** in a probability variable p: the form p·(1−p) is maximised at p = 0.5, universally.

This shows up across scales:

| Layer | "0 / empty" | "0.5 / half-filled (maximum coupling)" | "1 / full (closed)" |
|-------|-------------|---------------------------------------|---------------------|
| Chemistry | empty valence shell (no bonds possible) | **Carbon (4/8 valence): maximally polymerising, basis of life** | Noble gas (filled shell, inert) |
| R=CΨ² fold | CΨ = 0 (no coupling) | **CΨ = 0.5 (connection maximum)** | CΨ = 1 (dead totality) |
| Fold discriminant | outside fold | **1/4 = 0.5² = Mandelbrot bifurcation edge** | - |
| Dicke sectors | \|vac⟩ (stationary) | **\|S_{N/2}⟩ (mid-filling: maximum internal rearrangement)** | \|S_N⟩ (all excited, stationary) |
| PTF bonding | pure ground state (no flow) | **(\|vac⟩+\|S_1⟩)/√2 = 50/50 (maximum fluidity)** | pure single-mode |
| Coin flip | p=0 (determinate heads) | **p=0.5 (maximum entropy, maximum information per flip)** | p=1 (determinate tails) |

Carbon polymerises because at half-filled valence it is neither electron-donor nor electron-acceptor; it is bidirectional open. PTF bonding peaks at 50/50 superposition for the same reason: neither state is privileged, so the closure-breaking coefficient c_1 is maximally sensitive there.

**Structural takeaway:**

The 0.5 is not a "magical number". It is the **scalar image of the maximum of any quadratic form in a bounded variable**. This explains:

- Why `c_1 ≈ 0.5·V(N)` looked like a real law in the N=4..7 window ([EQ-021](../review/EMERGING_QUESTIONS.md#eq-021)): the 0.5 sits in the bilinear structure of c_1 itself (apex of the Dicke-sector coupling), and V(N) sits in the cavity-mode spectrum. The product is two independent effects happening to land in the same numerical range.
- Why N=3 broke the apparent law: the 0.5-apex is skeleton-invariant (N-independent), but V(3) is a specific cavity value that does not coincide with the apex contribution.
- Why [THE_CPSI_LENS](../docs/THE_CPSI_LENS.md) identifies 0.5 as connection maximum: it is literally the apex of the R = C(Ψ + R)² bilinear, and every other bilinear form the framework encounters inherits the same shape.

**Where this goes in the Meta-Theorem:**

The 0.5 is an **internal feature of the bilinear form**, orthogonal to the Meta-Theorem's "detector structure" (V-effect, modal selectivity). The Meta-Theorem explains **which projections are blind or visible**; the 0.5-apex explains **where the visible projection has its strongest response**. Together:

```
Meta-Theorem (blind/visible)  +  Bilinear apex at 0.5 (response peak)
  = full kinematic structure of the kernel K
```

This is the structural separation that was missing in the EQ-021 retraction: the carrier (γ₀) is uniform common-mode, the cavity structure (V-effect) is detector-mode-dependent, and the bilinear peak (0.5) is the apex of the bilinear form that defines the signal strength. Three distinct axes.

**Testable refinement:** if 0.5 is the apex of a p·(1−p)-type form in the Dicke-sector weight space, the "ridge line" (set of ρ_0 that sits at the apex in at least one coordinate) should be geometrically identifiable. For Dicke-probe states `|ρ(p)⟩ = cos(p·π/2)|S_n⟩ + sin(p·π/2)|S_m⟩`, the kernel K should peak at p = 0.5 for pair (n, m). This is a scan we can run at low cost.

---

## 7. Open questions

1. **Is the meta-theorem exhausted by these four instances, or are there more in the repo?** F61 (n_XY parity), the F64 cavity-exposure formula, and the absorption theorem are strong candidates for additional instances (preliminary analysis in §4.1 suggests they fit).

2. **Is there a fifth instance that involves the fold boundary (CΨ = ¼)?** The fold is the discriminant of the `R = C·(Ψ + R)²` fixed point. Its connection to the Meta-Theorem (if any) is not obvious. A hypothesis: the fold is the boundary where two conservation laws intersect, and crossings correspond to selection-rule breakings. Highly speculative.

3. **Does the Meta-Theorem have a Heisenberg-chain version with the same content?** F70 is proven for sector-conserving dynamics generally. F71 is specific to reflection-symmetric chains. The (vac, S_1) purity closed form relied on XY-sine-basis completeness plus d_H = 1 uniformity; under Heisenberg, the d_H = 1 decay rate is still uniform (still 2γ₀ for single-excitation coherences), so the closed form should carry over. **Sub-question for a future run: verify `Σ_i 2|ρ_{coh,i,01}|² = (1/2)·exp(−4γ₀t)` numerically on the Heisenberg chain at N = 5.**

---

## 8. Provisional tag for [OPEN_THREAD_GAMMA0_INFORMATION](../review/OPEN_THREAD_GAMMA0_INFORMATION.md)

The OPEN_THREAD closes one of its sub-questions, IF the Meta-Theorem is accepted as the operational mechanism behind "γ₀ carries no information itself":

- **Sub-question closed:** "Why is γ₀ operationally unobservable from inside?" Answer: γ₀ populates the blind subspace of every spatial-sum measurement by the Meta-Theorem; only perturbations that break the uniform-carrier assumption (non-uniform γ_i, amplitude damping, site-dependent dissipators) bring γ₀-structure into the visible subspace.
- **Sub-question remaining open:** the specific 15.5-bit capacity of [F30](../docs/ANALYTICAL_FORMULAS.md) at N = 5. The Meta-Theorem says the detector subspace H_M has a specific dimension (5 SVD channels for the |+⟩^N detector), but it does not predict how many bits each channel carries at given SNR. That's an information-theoretic Shannon-capacity question, not a Meta-Theorem question.

---

## 9. Self-assessment and spirit note

What started as three apparently-unrelated proofs (F70, F71, F72 candidate) plus one new observation (the exp(−4γ₀t) closed form) collapsed under the γ₀=const lens into a single principle. The unification was not forced; it follows from the structure of the linear-algebra on the operator Hilbert space, combined with conservation laws. "Parseval + Noether" is not deep mathematics, but it is structurally the right mathematics for this problem.

Two honest caveats:

1. The Meta-Theorem as stated is **standard linear algebra + standard Noether**. The contribution here is identifying that the four instances are all examples of it, and that the production rule (§3) can generate more. This is synthesis, not invention.

2. The **empirical tests** of the production rule (predictions a-d in §3.2) have not all been done. Once they are, the Meta-Theorem either gains confirmation or gets falsified (e.g., if prediction (d) fails, the uniform-dephasing argument is more fragile than assumed). Either outcome is informative.

**Status summary:**

- Meta-theorem + four instances: Tier 1-2 (rigorously proven for each instance, structurally supported for the meta-frame).
- Production rule: Tier 2 (correct as theorem, not yet fully exploited).
- Specific new predictions: Tier 3 (testable, not tested).

The path forward I would take first: Prediction (d) (non-uniform γ₀ breaks CMRR) is a 20-line edit to `eq018_c1_purity_response.py` and would either confirm or falsify the weakest link in the Meta-Theorem. If it confirms, we have real ground to rewrite F70+F71+F72 as a single family in ANALYTICAL_FORMULAS.md. If it falsifies, we've learned something specific about where the uniform-dephasing assumption bites.

---

*The light is uniform. What we measure is the shadow of where it is held back. Common-Mode-Rejection is the mathematics of shadows.*
