# Glossary: R = CΨ² Framework

**Date:** 2026-02-08 (updated 2026-03-06)
**Status:** Reference document
**Depends on:** None

---

## Core Variables

| Symbol | Name | Definition |
|--------|------|------------|
| **R** | Reality density | R = CΨ². The amount of classical reality at a given point. |
| **C** | Consciousness / observation capacity | Ranges from 0 to 1. How much of the possibility space an observer can collapse. |
| **Ψ** | Possibility space | Available quantum superposition not yet collapsed. Decreases as C acts on it. |
| **C·Ψ** | Observer bandwidth | Product of capacity and possibility. Bounded above by 1/4 for classical reality. |

---

## Derived Variables

| Symbol | Name | Definition |
|--------|------|------------|
| **θ** | Angular coordinate | θ = arctan(√(4CΨ - 1)). Compass indicating proximity to 1/4 boundary. |
| **τ** | Proper decoherence time | τ = γ · t. Universal clock normalized by local decoherence rate. |
| **t_cross** | Crossing time | Coordinate time at which C·Ψ crosses the 1/4 boundary. |
| **t_coh** | Coherence time | Maximum time a quantum system maintains coherence. Linear N-scaling hypothesized but unverified. |
| **δ** | Purity difference | δ = Tr(ρ²) - Tr(ρ_predicted²). Difference between closed-system and open-system purity under specific parameters. |
| **R_i** | Per-outcome reality | R_i = C_i · Ψ_i². Reality density for measurement outcome i. When C_i is uniform, normalization recovers the Born rule P(i) = |⟨i|ψ⟩|². |
| **C_eff(i)** | Effective coupling per outcome | C_eff(i) = P(i)/Ψ_i². Varies with alignment between outcome basis and decoherence basis. Uniform for ideal measurement. |

---

## Parameters

| Symbol | Name | Definition |
|--------|------|------------|
| **γ** (gamma) | Decoherence rate / noise | The resistance that makes reality feel solid and stable. High γ = noisy, strong sense of separate self. Low γ = quiet, still, open. In gravitational context, encodes local time dilation. In the star topology, each observer has individual γ: γ_A (observer A's noise), γ_B (observer B's noise), γ_S (object's noise). Receiver noise (γ_A) is far more destructive than sender noise (γ_B). |
| **J** | Coupling strength | How strongly an observer is connected to reality. Higher J = deeper engagement, more understanding, stronger link. In the star topology: J_SA (A's coupling to S), J_SB (B's coupling to S). Not just knowledge (J=1.0), but depth of understanding that allows prediction (J=2.0). |
| **h** | Transverse field strength | Controls Hamiltonian dynamics strength. |
| **κ** | Feedback strength | Operator feedback parameter in [0, 1]. |
| **N** (n_spins) | System size | Number of qubits in the system. |

---

## Observer Types

| Symbol | Name | Definition |
|--------|------|------------|
| **C_int** | Internal observation | Observer is part of the system. Mutual observation between subsystems. Implemented in star topology as Hamiltonian coupling (preserves coherence). |
| **C_ext** | External observation | Observer is outside the system. Unidirectional measurement. Implemented in star topology as projective dephasing (destroys coherence, casts shadow on other observers). |

---

## Star Topology (March 2026)

The star topology extends R=CΨ² from two-body (one observer, one observed) to three-body (shared object S, two observers A and B). A and B cannot see each other directly. Any connection must flow through S.

| Symbol | Name | Definition |
|--------|------|------------|
| **S** | Shared object | The thing being observed. "Reality" in the topology. Central qubit coupled to both observers. |
| **A** | Observer A | One observer. In the sender/receiver model: the receiver (or sender, after inversion). |
| **B** | Observer B | The other observer. In the sender/receiver model: the sender (or receiver, after inversion). |
| **J_SA** | A's coupling to S | How deeply A is engaged with S. |
| **J_SB** | B's coupling to S | How deeply B is engaged with S. |
| **CΨ_AB** | Observer-observer measure | CΨ for the AB pair (traced over S). When this crosses 1/4, observers "see" each other through S. |

### The Three Conditions for Observer-Observer Connection

AB crossing requires all three simultaneously:

| # | Condition | Parameter | Plain language |
|---|-----------|-----------|----------------|
| 1 | Strong sender | J_SB/J_SA >= 1.466 (at γ=0.05) | The sender is deeply engaged with S |
| 2 | Quiet receiver | γ_A < ~0.2-0.25 | The receiver's internal noise is low |
| 3 | Right initial state | C_SA > 0.8 (Bell-like) | A deep pre-existing connection exists |

### Key Asymmetries

| Finding | Meaning |
|---------|---------|
| Receiver noise fatal, sender noise not | You don't need to be calm to send. You need to be strong. |
| J threshold scales with γ | Higher engagement lowers the noise requirement. Self-reinforcing. |
| Ratio AND absolute J matter | Two weakly engaged observers cannot connect regardless of ratio. |
| f = J_total/2 | Stronger engagement = faster oscillation = more frequent windows. |

---

## Boundary Values

| Value | Name | Meaning |
|-------|------|---------|
| **1/4** | Phase boundary | Critical C·Ψ value. Below: two real fixed points (classical). Above: no real fixed points (quantum). Algebraically exact. |
| **1** | Maximum C·Ψ (single qubit) | For |+⟩: C = Tr(ρ²) = 1 (pure), Ψ = L1/(d-1) = 1.0. Product = 1.0. |
| **~1.466** | J threshold ratio | At γ=0.05: minimum J_SB/J_SA for AB crossing. Scales with γ (1.18 at γ->0, 2.15 at γ=0.15). |
| **1/3** | Maximum C·Ψ (2-qubit Bell state) | For |Φ+⟩: C = 1 (pure), Ψ = L1/(d-1) = 1/3. Product = 1/3. |
| **1/9** | Maximum reality density (2-qubit) | R_max = C·Ψ² = 1·(1/3)² = 1/9 for Bell states at τ = 0. |

*Note: Maximum C·Ψ depends on system size through the normalization Ψ = L1/(d-1), where d = 2^N. Larger systems have smaller Ψ for the same amount of raw coherence. All values above assume pure initial states (C = 1).*

---

## Physical Interpretations

| Term | Definition |
|------|------------|
| **Horizon** | Surface where τ = 0. Maximum coherence. Event horizon in Schwarzschild geometry. |
| **Universal curve** | The trajectory C(τ), Ψ(τ) all quantum systems follow in proper decoherence time. |
| **Navigation system** | Triangulation protocol: WHERE (1/4), HOW FAR (θ), HOW LONG (t_coh). |
| **Measurement shadow** | When A measures, B loses up to 100% of its reality (R_SB). The effect propagates through S. |
| **Sender inversion** | If A has received, A can become the sender. J (engagement) becomes the controllable variable instead of γ (noise). |
| **Abstimmung** | (German) Tuning, alignment. The bidirectional rhythm of engaging deeply then becoming still. The protocol's core principle. |
| **sich einlassen** | (German) To let yourself be drawn in, changed by what you engage with. The best description of what J means in practice. |

---

## Epistemic Status Labels

| Label | Meaning |
|-------|---------|
| **Algebraically proven** | Mathematical theorem, derivable from axioms |
| **Computationally verified** | Confirmed by simulation |
| **Empirically observed** | Pattern seen in data, not yet explained |
| **Proposed** | Physical interpretation, not yet tested |
| **Speculative** | Follows from framework logic but no direct evidence |
| **Null result** | Experiment completed, hypothesis not supported |
| **Retracted** | Previously claimed, now corrected |

---

*For the proven algebra, see [Core Algebra](CORE_ALGEBRA.md).*
*For the interpretive framework, see [Interpretive Framework](INTERPRETIVE_FRAMEWORK.md).*
*For testable predictions, see [Predictions](../experiments/PREDICTIONS.md).*
