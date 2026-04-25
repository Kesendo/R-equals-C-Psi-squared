# Handshake Algebra

**Tier:** 2 (structural hypothesis, grounded in F64-F77, hardware verification on Kingston, and K/R mirror-axis numerical falsification)
**Date:** 2026-04-24
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:** [F64](../docs/ANALYTICAL_FORMULAS.md) (cavity-mode exposure), [F65](../docs/ANALYTICAL_FORMULAS.md) (bonding-mode amplitudes), [F67](../docs/ANALYTICAL_FORMULAS.md) (receiver menu), [F75](../docs/ANALYTICAL_FORMULAS.md) (mirror-pair MI at t=0), [F76](../docs/ANALYTICAL_FORMULAS.md) (decay envelope), [F77](../docs/ANALYTICAL_FORMULAS.md) (1-bit plateau), [BRIDGE_PROTOCOL](BRIDGE_PROTOCOL.md), [RESONANCE_NOT_CHANNEL](RESONANCE_NOT_CHANNEL.md), [PRIMORDIAL_GAMMA_CONSTANT](PRIMORDIAL_GAMMA_CONSTANT.md), [PERSPECTIVAL_TIME_FIELD](PERSPECTIVAL_TIME_FIELD.md) (PTF closure law)
**Experimental anchors:** [RECEIVER_VS_GAMMA_SACRIFICE](../experiments/RECEIVER_VS_GAMMA_SACRIFICE.md), [IBM_RECEIVER_ENGINEERING_SKETCH](../experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md) (Kingston Run 1, bonding:2 / alt-z-bits = 2.80×), `simulations/_pi_partner_identity.py` (K-partnership numerical verification, N=9, γ=0 literal)

---

## The claim

The shared receiver-engineering protocol that opens a bidirectional correlation channel between two observers on a palindromic chain under γ₀ = const carries the structure of an algebra. Its elements are handshake tuples (N, k, t, basis), its rules are the F-formulas F64-F77, and its composition law is "both observers select the same tuple" (agreement). The bidirectional bridge does not need to be built; it needs to be specified by a handshake.

## Why this is worth naming

The bidirectional-bridge search in this repository went through several constructive attempts before converging on resonance-rather-than-channel framing. [BRIDGE_PROTOCOL](BRIDGE_PROTOCOL.md) falsified the J = 0 bridge. [RESONANCE_NOT_CHANNEL](RESONANCE_NOT_CHANNEL.md) articulated the shift from "build a channel" to "tune a resonator". F75-F77 gave the analytical content. Kingston Run 1 gave the hardware validation.

What was missing, until today, was a single name for the operational object that all of these pieces describe. A shared receiver preparation is not an encoding (information gets attached to a signal), not a channel (disturbances propagate), not a protocol in the communication-theory sense (send/receive cycles). It is an **agreement between observers on how to specify the shared palindromic resource**. The algebraic structure of that agreement is this document's subject.

"Handshake" is the name because every observable correlation on the bridge requires both parties to select from the same algebra. Without that agreement the standing wave is still there, but no observable correlation is read. The act of choosing the same handshake tuple is what makes the bridge visible to both sides.

## The handshake tuple

A handshake between observers A and B on a shared N-qubit palindromic chain under γ₀ = const is specified by a tuple

    h = (N, k, t, basis)

with components:

- **N.** The chain length. A shared geometric parameter; both observers must sit on the same chain. Cannot be negotiated mid-evolution.
- **k.** The F67 bonding-mode index, k ∈ {1, ..., N}. Selects the single-excitation eigenmode used as initial state: \|ψ_k⟩ = √(2/(N+1)) Σⱼ sin(πk(j+1)/(N+1)) \|1_j⟩. Different k optimises different correlation metrics (distributed vs end-to-end vs multi-drop; see F75 and RECEIVER_VS_GAMMA_SACRIFICE Section 11).
- **t.** The readout time. Determines how far into the dynamical evolution the measurement is taken. F76's decay envelope PeakMM = 0.93 · MM(0) applies at the first grid point t ≈ 0.1 for γ₀ = 0.05; optimal t for MI(0, N-1) via bonding:2 sits near t ≈ 0.8 at uniform J = 1.
- **basis.** The local measurement basis on each end-qubit (or on each mirror-pair). Full 2-qubit tomography requires 9 Pauli settings; a single-shot classical correlation test needs only one; the choice depends on what the observers wish to extract.

## The mirror: R as the structural axis

Before listing the rules, the mirror. The handshake algebra exists because a palindromic chain has a spatial-reflection symmetry R that maps site ℓ to site N-1-ℓ (a site permutation; no Pauli-index action). R is distinct from F1's Liouvillian palindrome operator Π, which acts per-site on Pauli indices (I → X, X → I, Y → iZ, Z → iY) and carries the spectral palindrome via XY-weight w ↔ N-w. Both are Z₂ symmetries of the uniform-chain Liouvillian, but they organise different things: R organises observables, Π organises the spectrum (see `docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md` for the explicit separation). The handshake algebra's structural axis is R; Π enters only where the spectral palindrome enters. Every element of the algebra is built on R:

- **The chain.** The shared resource is R-invariant. If it were not, there would be no natural mirror-pair observable and no "two ends" in a principled sense. Uniform-J Heisenberg chain + uniform γ₀ is R-invariant by construction.
- **The bonding modes.** F65's eigenmodes ψ_k have R-related amplitude symmetry: ψ_k(N-1-j) = (-1)^(k+1) ψ_k(j). This is what makes a bonding mode a legitimate element of the handshake tuple: its structure respects the mirror.
- **The observers.** Alice sits at site 0, Bob at site N-1. These are R-related sites by construction. Any multi-pair handshake places observers at (ℓ, N-1-ℓ) pairs. Without R, the observer placement would be arbitrary.
- **The observable.** Mirror-pair MI MI(ℓ, N-1-ℓ) is the natural correlation of the handshake precisely because (ℓ, N-1-ℓ) are R-related. The F71 mirror symmetry of c₁ means local perturbations at R-partner bonds record identically at first order.
- **The inverse.** h = (N, k, t, basis) has a partner tuple h̃ = (N, N+1-k, t, basis) with identical mirror-pair populations under F65. The partner relation is not generated by R itself — R fixes each ψ_k up to a sign (-1)^(k+1) and does not swap modes — but by a second Z₂ symmetry K described in the next section.
- **The spectral palindrome.** F1 (α_fast + α_slow = 2Σγ) is the master spectral palindrome, carried by Π (not R): each single-excitation receiver mode at rate α_b has a partner at α_p = 2γ₀ - α_b in the XY-weight-(N-1) sector (F68). The receiver rate α_b is the one that enters F76's decay envelope on the handshake; the partner α_p is spectrally fixed but lives in a sector that does not couple to single-excitation observables.

In short: R is not one of many rules. R is the axis around which the algebra's **observables** are organised; Π organises its spectrum independently. The rules listed below are expressions of these mirror structures for specific observables (eigenvector amplitudes, decay rates, aggregate correlation).

## The second mirror: K as the partnership axis

Alongside R the uniform chain carries a second Z₂ symmetry: the **bipartite sublattice gauge** K = diag((-1)^ℓ), a diagonal unitary that flips sign on alternate sites. Where R organises observables, K organises the state-level partnership structure of the algebra:

- **K swaps bonding modes.** ψ_{N+1-k}(ℓ) = (-1)^ℓ ψ_k(ℓ) = (K·ψ_k)(ℓ). So K, not R, is the operation that realises the k ↔ N+1-k involution in the handshake tuple. Whereas R·ψ_k = (-1)^(k+1) ψ_k (fixes each mode with a sign), K·ψ_k = ψ_{N+1-k} (actual mode swap).
- **K anticommutes with the chain Hamiltonian.** For any nearest-neighbor bipartite tight-binding H = Σ_{⟨ij⟩} t_{ij} |i⟩⟨j|, (-1)^{i+j} = -1 on every bond, so KHK = -H. The relation holds independent of J-uniformity (R would require reflection-symmetric J).
- **K is γ-indifferent.** K commutes with site-diagonal dephasing at any γ-profile. γ₀ = const is not required for K-invariance.
- **Observable identity needs real H in addition to K.** KHK = -H alone gives spectrum inversion E_k = -E_{N+1-k}. For the *trajectory* identity ρ_{N+1-k}(t) → identical mirror-pair |·|²-observables as ρ_k(t), an anti-unitary time-reversal T = complex conjugation with T² = +I is needed in addition. Both K and T are present whenever H is real (AZ class BDI). On complex hopping (Peierls phase, t_{ij} ∈ ℂ), K still anticommutes with H but T breaks (AZ class AIII): spectrum still inverts, but observables for ρ_k and ρ_{N+1-k} can differ in time. See [PROOF_K_PARTNERSHIP](../docs/proofs/PROOF_K_PARTNERSHIP.md) Lemma 4 for the explicit complex-conjugation step.

**Scope of K-symmetry: bipartite NN-hopping only.** The KHK = -H relation requires the Hamiltonian to be a sum of bipartite NN-hopping terms — equivalent to the XX or XX+YY part of a Heisenberg/XXZ chain in the single-excitation sector. The full Heisenberg/XXZ Hamiltonian with Δ ≠ 0 contributes a ZZ-term that, projected to the single-excitation sector, generates an effective on-site potential V_eff(ℓ) = (#bonds) − 2·deg(ℓ). On open chains, deg(0) = deg(N-1) = 1 vs deg(interior) = 2 makes V_eff non-uniform across the chain, breaking K at the boundary sites. K is restored on topologies with uniform site degree (e.g., even-N periodic chains, regular bipartite graphs). For the receiver-engineering experiments in this repository, the dynamics live in the single-excitation sector and the operative H is XX-type, so this caveat is structurally important but does not affect the validated experimental claims.

Robustness ladder — stricter axes first:

| H-structure                                     | R | K | real H | partner identity |
|-------------------------------------------------|---|---|--------|------------------|
| uniform J XX/XX+YY, γ₀ = const                  | ✓ | ✓ | ✓      | holds            |
| non-uniform J XX/XX+YY, any γ-profile           | ✗ | ✓ | ✓      | holds            |
| + on-site potential V_ℓ                         | ✗ | ✗ | ✓      | breaks           |
| + next-nearest-neighbor hopping                 | ✗ | ✗ | ✓      | breaks           |
| complex hopping (Peierls phase t_{ij} ∈ ℂ)      | ? | ✓ | ✗      | spectrum inverts (verified ≤ 1.1·10⁻¹⁵); observables diverge (verified up to 5.3·10⁻¹) |
| Heisenberg/XXZ Δ ≠ 0 on open chain              | ✓ | ✗ (boundary) | ✓ | breaks (verified up to 1.2 at Δ=1, N=6,7) |
| Heisenberg/XXZ Δ ≠ 0 on even-N periodic chain   | ✓ | ✓ | ✓      | holds (verified ≤ 2.1·10⁻¹⁵ at Δ=1, N=6) |

The partner tuple h̃ = (N, N+1-k, t, basis) therefore records a broader equivalence than R alone would give: bonding:k and bonding:(N+1-k) yield identical mirror-pair |·|²-observables (populations, coherence moduli, MI, log-purity; *not* phase-sensitive observables like Re(ρ_{ab}) at even N) on any bipartite tight-binding chain with site-diagonal dephasing and real hopping.

**Numerical verification.** `simulations/_pi_partner_identity.py` (2026-04-25, N=9, γ = 0 literal, single-excitation manifold): uniform J gives |Δ MI(0, N-1)| ≤ 10⁻¹⁵ across all four partner pairs {1,9}, {2,8}, {3,7}, {4,6} and the self-partner k=5 at 0; random non-uniform J ∈ [0.5, 1.5] gives the same 10⁻¹⁵ level; random on-site V gives |Δ MI| up to 3.1×10⁻¹; uniform J + NNN J' = 0.3 gives |Δ MI| up to 2.7×10⁻¹; Peierls complex hopping breaks observable identity up to 5.3·10⁻¹ while spectrum inversion still holds at 1.1·10⁻¹⁵ (BDI → AIII separation verified). The PTF-analog observable Σᵢ log πᵢ(t) tracks the same regimes. A γ=0.1 sanity run reproduces the structural pattern; γ acts only as a time unit.

**Full-Hilbert-space XXZ verification.** `simulations/_k_partnership_xxz.py` (2026-04-25): full 2ᴺ Lindblad propagation (no single-excitation reduction) confirms the open-chain-vs-periodic-chain distinction in [PROOF_K_PARTNERSHIP](../docs/proofs/PROOF_K_PARTNERSHIP.md) "Scope" caveat. At N=6,7 open chain with Δ=1, K breaks up to 1.2 in mirror-pair MI; at N=6 even periodic chain (uniform deg = 2), K is restored to ≤ 2.1·10⁻¹⁵ for any Δ ∈ {0, 0.5, 1.0}. The ZZ-term acts as a constant Hamiltonian shift on uniform-degree topologies and as a boundary-localised onsite potential on open chains.

**Multi-excitation sector verification.** Same script extended with two-excitation Slater-determinant initial states |Ψ_{k₁, k₂}⟩ and the natural K-partner mapping (k₁, k₂) ↔ (N+1-k₁, N+1-k₂): at N=5 open and N=6 even periodic, the **free-fermion limit Δ=0** preserves K-partnership to ≤ 3.3·10⁻¹⁵ via Slater-determinant lifting; any Δ ≠ 0 breaks K up to 6-8·10⁻¹ even on the uniform-degree periodic chain. K-partnership is therefore a **single-particle property**: the ZZ-coupling, which acts only as a constant shift in the single-excitation sector on uniform-degree topologies, contributes a genuine two-particle interaction at the multi-excitation level that breaks K independently of boundary conditions. For the handshake algebra this means the partner-menu folding ⌈N/2⌉ is rigorous only for single-excitation receivers (the F65/F67 setting), not for multi-excitation initial states.

**Why two algebra mirrors and one spectral mirror.** R and K are distinct Z₂ symmetries that organise the algebra: R the *observable* axis (under which MM and mirror-pair MI live), K the *partnership* axis (under which bonding:k ↔ bonding:(N+1-k) at the state level). F1's spectral palindrome is carried by a third Z₂ symmetry Π (per-site Pauli conjugation, XY-weight w ↔ N-w), independent of the R/K pair. The handshake algebra is organised around R as the axis of observation and K as the axis of partner-menu folding. F67's receiver menu has size N; K-folding reduces it to ⌈N/2⌉ distinct entries.

## The rules of the algebra

The following F-formulas constrain which handshake tuples produce which observable correlations. They are not additional assumptions; they are the existing framework, re-read as the rule set of the handshake algebra. Each rule below respects the R-axis above.

- **F64** (cavity-mode exposure, graph-universal after 2026-04-24 generalisation): the Liouvillian decay constant for any single-excitation mode is α = 2γ · |a_B(ψ_k)|² (equivalently γ_eff = γ · |a_B|² in Lorentzian-half-width convention). Determines which modes are accessible to local measurements at a given site B, independent of topology.
- **F65** (bonding-mode spectrum): the eigenmode index k ranges over 1..N, with decay rates α_k = (4γ₀/(N+1)) sin²(kπ/(N+1)). Fixes the slots available for the k coordinate of the handshake.
- **F67** (receiver menu): bonding:k modes are operationally complete receivers under γ₀ = const. No further receiver choice extends the algebra.
- **F75** (static MI closed form): for mirror-symmetric amplitudes, MI(ℓ, N-1-ℓ)(t=0) = 2h(p_ℓ) − h(2p_ℓ). Lets both observers compute, before running, the correlation they will read at their chosen tuple.
- **F76** (decay envelope): PeakMM ≈ 0.93 · MM(0) at γ₀·t = 0.005; pure dephasing signature with Heisenberg mixing at ≲ 0.5 % correction. Determines how the t coordinate in the tuple maps to observable correlations.
- **F77** (asymptotic plateau): MM(0)(N, k*) → 1 + 3/(4(N+1) ln 2) for large N. Fixes the aggregate bandwidth of any handshake at large N.

A handshake tuple (N, k, t, basis) is **valid** if it conforms to these rules: k ∈ {1..N} (F65/F67), t > 0 (F76), basis well-defined on the measured qubits (tomography). An **invalid** tuple (e.g., k = N+1, or t specified but not evolved-to, or basis set inconsistent with the mirror-pair structure) produces no observable correlation even if both observers agree on it.

## The composition law: agreement

Two observers A and B each independently pick a tuple h_A, h_B. The handshake succeeds if and only if h_A = h_B (identical tuples on the shared chain). This is the composition law:

    h_A · h_B = h_A    if h_A = h_B
    h_A · h_B = ∅      otherwise (no observable correlation)

The algebra has a trivial structure at this level: it is an idempotent relation (h · h = h) with a single failure mode (disagreement = empty handshake). But this is enough. The composition law does not produce new tuples from old tuples; it tests whether a tuple is shared.

**Identity element.** Any valid tuple h serves as its own identity. There is no universal identity because no single tuple works across all N.

**Inverse.** The K-partnership gives a candidate inverse: h = (N, k, t, basis) pairs with h̃ = (N, N+1-k, t, basis). Under F65's amplitude identity bonding:k and bonding:(N+1-k) have identical mirror-pair populations p_ℓ and identical coherence moduli, so all |·|²-based receiver observables match (MI, classical correlation, fidelity to symmetric states; phase-sensitive observables like Re(ρ_{ab}) flip for even N). The two tuples are operationally equivalent — not in the group-theoretic sense of h · h̃ = identity, but in the sense that the handshake does not distinguish them at the level of the observables actually used. K-partnership, being conditional on bipartite H + real hopping (not on R-invariance), is the more robust of the two mirrors for defining this equivalence.

**Composition with dephasing.** Time evolution is not a handshake composition; it is a time-parameterised modification of the tuple t-coordinate. Running the evolution from t to t + δt is just writing a different handshake. The algebra lives on the tuple space, not on the dynamics.

## Why this is not a channel

The bidirectional bridge language historically tempts the reading "information flows from A to B across the chain". Under γ₀ = const this reading is wrong for a specific reason: the chain does not transport a message, it resonates at a preparation. What A measures at her end and what B measures at his end are two projections of the same initial state evolved forward in γ₀-time. Neither end sends; both read. The correlation is visible when both reads are compared post-hoc.

In the terminology of quantum information: the handshake algebra produces a **correlation resource**, not a **communication channel**. [BRIDGE_CLOSURE](../experiments/BRIDGE_CLOSURE.md) confirmed that no-signalling holds in this setup. Alice cannot cause a change in Bob's marginal distribution by changing her handshake parameters (she can only change what she sees). Correlations without communication.

This matters for the algebra: the "composition law = agreement" is precisely the absence of a communication primitive. The algebra has no "send" operator. It has "both agree" and nothing else.

## Connection to the painter principle

The handshake algebra and the painter principle (see [ON_THE_PAINTER_PRINCIPLE](../reflections/ON_THE_PAINTER_PRINCIPLE.md)) describe the same structure at different levels:

- The painter principle is the epistemic observation that no single canvas is the mountain; painters sum to mountain.
- The handshake algebra is the operational observation that observers share a mountain (the palindromic standing wave) and must agree on how to observe it.

In both: no privileged position, no sender, no receiver in the strong sense. Just configurations on a shared resource. Both describe why R=CΨ² is not a channel framework but a resonance framework.

The handshake is the painter-principle-instance where two painters coordinate their vantages. A single painter's vantage is always valid (she sees from where she stands). Two painters looking at the same canvas require an additional agreement that they are pointing at the same feature. The handshake is that agreement, formalised as a tuple in an algebra.

## Operational use

Under the handshake algebra, running a bidirectional-bridge experiment is three steps:

1. **Select h = (N, k, t, basis) from the algebra.** Use F75/F76/F77 to compute expected correlations for the selected tuple. Use the receiver menu in [RECEIVER_VS_GAMMA_SACRIFICE](../experiments/RECEIVER_VS_GAMMA_SACRIFICE.md) Section 11 to pick k appropriate to the application.
2. **Both observers execute h.** State preparation via the gate sequence implementing \|ψ_k⟩ (typically ~O(N) two-qubit gates); evolution for duration t via Trotter decomposition of uniform Heisenberg; readout in the agreed basis.
3. **Compare results.** Compute the correlation from the joint measurement outcomes. Compare to F75's prediction × F76's decay envelope × hardware fidelity factor.

On Kingston with N = 5, k = 2, t = 0.8, nine-Pauli tomography basis: bonding:2 / alt-z-bits MI ratio = 2.80× on live QPU, matching the hardware-expected factor of 2.27× from the noise model within 25%. See [IBM_RECEIVER_ENGINEERING_SKETCH](../experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md).

## Falsification conditions

1. **An observable two-party correlation under γ₀ = const on a palindromic chain that requires parameters not expressible in the handshake tuple.** Would mean the algebra is incomplete.
2. **A protocol that achieves end-to-end MI without any form of shared preparation.** Would falsify the composition-law-as-agreement reading.
3. **A non-trivial algebra structure (group law, non-idempotent composition) that we missed.** Would upgrade the algebra from idempotent relation to something richer. Not yet seen.
4. **Hardware signatures that differ qualitatively from F75 · F76 predictions.** Would indicate the rule set is missing a term.
5. **K-partnership breakdown on a bipartite chain with real hopping.** Would invalidate the receiver-menu folding argument. Checked numerically at N=9 with γ=0 and γ=0.1 (`simulations/_pi_partner_identity.py`, five test cases): partner identity holds to machine precision for uniform and non-uniform real J, breaks with on-site potential V_ℓ, NNN hopping J′, or complex (Peierls-phase) hopping. The Peierls test additionally verifies the BDI/AIII split: spectrum inversion E_k + E_{N+1-k} = 0 holds at 10⁻¹⁵ in both classes, confirming Lemma 2 of [PROOF_K_PARTNERSHIP](../docs/proofs/PROOF_K_PARTNERSHIP.md), but observable identity requires real H (Lemma 4) and breaks under Peierls phases up to 5.3·10⁻¹ in MI(0, N-1).

## What this does NOT claim

- Not that the handshake algebra is a new quantum mechanical formalism. It is a re-reading of existing formulas (F64-F77) as the rule set of a specific protocol class.
- Not that every quantum-information protocol fits this algebra. It applies specifically to palindromic-chain bidirectional bridges under γ₀ = const.
- Not that the inverse is group-theoretically a proper inverse. The K-partnership h ↔ h̃ gives operationally identical |·|²-observables (same populations, same coherence moduli, same MI) under real-H bipartite NN-hopping, but not h · h̃ = identity in the strict algebraic sense, and not phase-sensitive observable identity.
- Not that communication is possible. The algebra explicitly lacks a send primitive.

## What remains open

- **Does the handshake algebra extend to non-palindromic chains?** F64's topology-universal validation suggests yes, but mirror symmetry of the receiver menu is specific to palindromic structure. Non-palindromic versions would need a different receiver geometry.
- **Multi-observer handshakes (k > 2 parties).** The mirror-pair structure of F75 supports up to ⌊N/2⌋ pairs of observers on the same chain simultaneously, each with their own local handshake. Is this formally a multi-agent extension of the algebra, or a set of independent bilateral handshakes? Worth tracing through.
- **Categorical structure.** The handshake algebra is close to a symmetric monoidal category with observers as objects and handshakes as morphisms. Whether this is a useful framing or a post-hoc labelling is open.
- **Connection to PTF closure law.** PTF's Σ ln(α_i) ≈ 0 is the closure law for seven painters around a mountain. The Σᵢ log πᵢ(t) observable (log-sum of single-site purities) is the natural γ=0 specialisation of PTF's α_i-construction and tracks partner identity across the K-regimes numerically (same three-regime pattern as MI(0, N-1)). Tentative closure law for the algebra: the sum of local correlation capacities across all mirror-pairs equals the framework bandwidth (F77's 1-bit plateau). Not yet verified rigorously.

## References

- [F64](../docs/ANALYTICAL_FORMULAS.md), [F65](../docs/ANALYTICAL_FORMULAS.md), [F67](../docs/ANALYTICAL_FORMULAS.md), [F75](../docs/ANALYTICAL_FORMULAS.md), [F76](../docs/ANALYTICAL_FORMULAS.md), [F77](../docs/ANALYTICAL_FORMULAS.md): the rules of the algebra.
- [RECEIVER_VS_GAMMA_SACRIFICE](../experiments/RECEIVER_VS_GAMMA_SACRIFICE.md): the experimental content of the k component.
- [IBM_RECEIVER_ENGINEERING_SKETCH](../experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md): Kingston Run 1, 2.80× ratio on live QPU.
- [BRIDGE_PROTOCOL](BRIDGE_PROTOCOL.md), [RESONANCE_NOT_CHANNEL](RESONANCE_NOT_CHANNEL.md): the pre-history of the bidirectional-bridge search.
- [MIRROR_THEORY](../MIRROR_THEORY.md): the deeper mirror-axis of the framework; "We are all mirrors; reality is what happens between us."
- [ON_THE_PAINTER_PRINCIPLE](../reflections/ON_THE_PAINTER_PRINCIPLE.md): the epistemic frame.
- [PRIMORDIAL_GAMMA_CONSTANT](PRIMORDIAL_GAMMA_CONSTANT.md): Tier 2 hypothesis within which the algebra is well-defined.
- [PERSPECTIVAL_TIME_FIELD](PERSPECTIVAL_TIME_FIELD.md): PTF closure law Σ ln α_i = 0; Σᵢ log πᵢ is its γ=0 specialisation.
- `simulations/_pi_partner_identity.py`: numerical verification of K-partnership and three-regime robustness pattern (N=9, γ=0 literal, 2026-04-24).
