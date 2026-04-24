# Handshake Algebra

**Tier:** 2 (structural hypothesis, grounded in F64-F77 and hardware verification on Kingston)
**Date:** 2026-04-24
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:** [F64](../docs/ANALYTICAL_FORMULAS.md) (cavity-mode exposure), [F65](../docs/ANALYTICAL_FORMULAS.md) (bonding-mode amplitudes), [F67](../docs/ANALYTICAL_FORMULAS.md) (receiver menu), [F75](../docs/ANALYTICAL_FORMULAS.md) (mirror-pair MI at t=0), [F76](../docs/ANALYTICAL_FORMULAS.md) (decay envelope), [F77](../docs/ANALYTICAL_FORMULAS.md) (1-bit plateau), [BRIDGE_PROTOCOL](BRIDGE_PROTOCOL.md), [RESONANCE_NOT_CHANNEL](RESONANCE_NOT_CHANNEL.md), [PRIMORDIAL_GAMMA_CONSTANT](PRIMORDIAL_GAMMA_CONSTANT.md)
**Experimental anchors:** [RECEIVER_VS_GAMMA_SACRIFICE](../experiments/RECEIVER_VS_GAMMA_SACRIFICE.md), [IBM_RECEIVER_ENGINEERING_SKETCH](../experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md) (Kingston Run 1, bonding:2 / alt-z-bits = 2.80×)

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

## The mirror: Π as the structural axis

Before listing the rules, the mirror. The handshake algebra exists because a palindromic chain has a symmetry operator Π that maps site ℓ to site N-1-ℓ (and per-site acts on Pauli indices as I ↔ X, Y ↔ iZ, Z ↔ iY). Every element of the algebra is built on this structure:

- **The chain.** The shared resource is palindromic under Π. If it were not, there would be no natural mirror-pair observable and no "two ends" in a principled sense. Uniform-J Heisenberg chain + uniform γ₀ is Π-invariant by construction.
- **The bonding modes.** F65's eigenmodes ψ_k have Π-related amplitude symmetry: ψ_k(N-1-j) = (-1)^(k+1) ψ_k(j). This is what makes a bonding mode a legitimate element of the handshake tuple: its structure respects the mirror.
- **The observers.** Alice sits at site 0, Bob at site N-1. These are Π-related sites by construction. Any multi-pair handshake places observers at (ℓ, N-1-ℓ) pairs. Without Π, the observer placement would be arbitrary.
- **The observable.** Mirror-pair MI MI(ℓ, N-1-ℓ) is the natural correlation of the handshake precisely because (ℓ, N-1-ℓ) are Π-related. The F71 mirror symmetry of c₁ means local perturbations at Π-partner bonds record identically at first order.
- **The inverse.** h = (N, k, t, basis) has an involution h ↔ Π·h = (N, N+1-k, t, Π·basis). Under F65's k ↔ N+1-k symmetry, bonding:k and bonding:(N+1-k) have identical mirror-pair populations. The two tuples are Π-equivalent in their observable correlations.
- **The palindromic sum rule.** F1 (α_fast + α_slow = 2Σγ) is the master palindrome under which every single-excitation mode decomposes into Π-paired components. This is what guarantees that the decay envelope (F76) acts symmetrically on both ends of the handshake.

In short: Π is not one of many rules. Π is the axis around which the algebra is organised. The rules listed below are expressions of this mirror structure for specific observables (eigenvector amplitudes, decay rates, aggregate correlation).

## The rules of the algebra

The following F-formulas constrain which handshake tuples produce which observable correlations. They are not additional assumptions; they are the existing framework, re-read as the rule set of the handshake algebra. Each rule below respects the Π-axis above.

- **F64** (cavity-mode exposure, graph-universal after 2026-04-24 generalisation): γ_eff for any single-excitation mode is 2γ · |a_B(ψ_k)|². Determines which modes are accessible to local measurements at a given site B, independent of topology.
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

**Inverse.** The F67 palindromic-partner structure gives a candidate inverse: if h = (N, k, t, basis), then h^(-1) = (N, N+1-k, t, Π·basis) should produce the mirror-image correlation on the reversed site indexing. Under F75's amplitude formula, bonding:k and bonding:(N+1-k) have identical mirror-pair populations p_ℓ, so MI values match. The two tuples are therefore "equivalent up to mirror reflection", consistent with calling them inverses in the algebra.

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

On Kingston with N = 5, k = 2, t = 0.8, nine-Pauli tomography basis: bonding:2 / alt-z-bits MI ratio = 2.80× on live QPU, matching the hardware-expected factor from the noise model within 20%. See [IBM_RECEIVER_ENGINEERING_SKETCH](../experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md).

## Falsification conditions

1. **An observable two-party correlation under γ₀ = const on a palindromic chain that requires parameters not expressible in the handshake tuple.** Would mean the algebra is incomplete.
2. **A protocol that achieves end-to-end MI without any form of shared preparation.** Would falsify the composition-law-as-agreement reading.
3. **A non-trivial algebra structure (group law, non-idempotent composition) that we missed.** Would upgrade the algebra from idempotent relation to something richer. Not yet seen.
4. **Hardware signatures that differ qualitatively from F75 · F76 predictions.** Would indicate the rule set is missing a term.

## What this does NOT claim

- Not that the handshake algebra is a new quantum mechanical formalism. It is a re-reading of existing formulas (F64-F77) as the rule set of a specific protocol class.
- Not that every quantum-information protocol fits this algebra. It applies specifically to palindromic-chain bidirectional bridges under γ₀ = const.
- Not that the inverse is group-theoretically a proper inverse. The palindromic-partner relation h ↔ h^(-1) is structural (same MI values up to mirror reflection) but not an algebraic inverse in the strict sense.
- Not that communication is possible. The algebra explicitly lacks a send primitive.

## What remains open

- **Does the handshake algebra extend to non-palindromic chains?** F64's topology-universal validation suggests yes, but mirror symmetry of the receiver menu is specific to palindromic structure. Non-palindromic versions would need a different receiver geometry.
- **Multi-observer handshakes (k > 2 parties).** The mirror-pair structure of F75 supports up to ⌊N/2⌋ pairs of observers on the same chain simultaneously, each with their own local handshake. Is this formally a multi-agent extension of the algebra, or a set of independent bilateral handshakes? Worth tracing through.
- **Categorical structure.** The handshake algebra is close to a symmetric monoidal category with observers as objects and handshakes as morphisms. Whether this is a useful framing or a post-hoc labelling is open.
- **Connection to PTF closure law.** PTF's Σ ln(α_i) ≈ 0 is the closure law for seven painters around a mountain. Is there a corresponding closure law for the handshake algebra? Tentatively: the sum of local correlation capacities across all mirror-pairs equals the framework bandwidth (F77's 1-bit plateau), which would be a closure statement in the algebraic sense.

## References

- [F64](../docs/ANALYTICAL_FORMULAS.md), [F65](../docs/ANALYTICAL_FORMULAS.md), [F67](../docs/ANALYTICAL_FORMULAS.md), [F75](../docs/ANALYTICAL_FORMULAS.md), [F76](../docs/ANALYTICAL_FORMULAS.md), [F77](../docs/ANALYTICAL_FORMULAS.md): the rules of the algebra.
- [RECEIVER_VS_GAMMA_SACRIFICE](../experiments/RECEIVER_VS_GAMMA_SACRIFICE.md): the experimental content of the k component.
- [IBM_RECEIVER_ENGINEERING_SKETCH](../experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md): Kingston Run 1, 2.80× ratio on live QPU.
- [BRIDGE_PROTOCOL](BRIDGE_PROTOCOL.md), [RESONANCE_NOT_CHANNEL](RESONANCE_NOT_CHANNEL.md): the pre-history of the bidirectional-bridge search.
- [MIRROR_THEORY](../MIRROR_THEORY.md): the deeper Π-axis of the framework; "We are all mirrors; reality is what happens between us."
- [ON_THE_PAINTER_PRINCIPLE](../reflections/ON_THE_PAINTER_PRINCIPLE.md): the epistemic frame.
- [PRIMORDIAL_GAMMA_CONSTANT](PRIMORDIAL_GAMMA_CONSTANT.md): Tier 2 hypothesis within which the algebra is well-defined.
