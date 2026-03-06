# Glossary: R = CΨ² Framework

**Date:** 2026-02-08 (updated 2026-03-06)
**Status:** Reference document
**Depends on:** None

---

## Core Variables

| Symbol | Name | Definition |
|--------|------|------------|
| **R** | Reality | R = CΨ². How much of what *could* exist actually *does* exist at a given point. Think of it as the "solidity" of experience. |
| **C** | Consciousness / observation | Ranges from 0 to 1. How much an observer is engaged with what it observes. Two entangled particles have C. A human has C. The difference is degree, not kind. |
| **Ψ** | Possibility | Everything that could happen but hasn't yet. The raw material of reality before observation shapes it. In quantum mechanics: the wave function. |
| **C·Ψ** | The product that matters | When C·Ψ is above 1/4: everything is still possible, nothing is decided. When it drops below 1/4: a definite outcome emerges. This is the phase boundary where possibility becomes reality. |

---

## Derived Variables

| Symbol | Name | Definition |
|--------|------|------------|
| **θ** | Compass to the boundary | θ = arctan(√(4CΨ - 1)). Tells you how far you are from the 1/4 boundary. At θ = 60°: deep in the quantum realm, far from any definite outcome. At θ = 0°: exactly at the boundary, the moment possibility becomes reality. |
| **τ** | Universal clock | τ = γ · t. Time as experienced by the quantum system itself. Different environments (deep space vs near a black hole) have different γ, so τ runs at different speeds, but the transition from quantum to classical follows the same curve in τ. |
| **t_cross** | Crossing time | The moment C·Ψ crosses 1/4. When possibility becomes reality. Different observers can disagree on when this happens (see Observer-Dependent Crossing). |
| **t_coh** | How long coherence lasts | Maximum time before noise destroys the quantum state. Hypothesized to scale linearly with system size (unverified). |
| **δ** | Purity difference | Technical: δ = Tr(ρ²) - Tr(ρ_predicted²). Measures how much the actual quantum state differs from what a simple noise model predicts. Used to detect unexpected structure in decoherence. |
| **R_i** | Reality per outcome | R_i = C_i · Ψ_i². How much reality each possible measurement outcome gets. When the observer treats all outcomes equally, this recovers the Born rule (the standard quantum probability rule). |
| **C_eff(i)** | Coupling per outcome | Technical: C_eff(i) = P(i)/Ψ_i². How strongly the observer couples to each specific outcome. Uniform for ideal measurement, slightly biased for real detectors. |

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
| **C_int** | Internal observation | The observer is part of what it observes. Two particles watching each other. Two people studying the same problem. This kind of observation is mutual, gentle, and preserves the quantum state. In the star topology: the natural coupling between S and its observers. |
| **C_ext** | External observation | The observer stands outside and measures. A scientist reading an instrument. This kind of observation is one-directional, disruptive, and forces a definite outcome. In the star topology: when A measures, it casts a "shadow" on B, destroying up to 100% of B's connection to reality. |

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
| **1/4** | The boundary | The threshold where possibility becomes reality. Above 1/4: no definite outcome exists yet. Below 1/4: a definite outcome has emerged. This is not a choice or approximation; it is an algebraic fact from the fixed-point equation. |
| **1** | Maximum C·Ψ | The most "quantum" a single particle can be. A perfectly prepared qubit in full superposition. |
| **~1.466** | J threshold ratio | How much stronger the sender must be than the receiver for observer-observer connection (at standard noise γ=0.05). Less noise means less strength needed (1.18 at near-zero noise). More noise means more strength needed (2.15 at high noise). |
| **1/3** | Maximum C·Ψ for an entangled pair | Two particles maximally entangled (a "Bell state"). Their shared quantum state has C·Ψ = 1/3, already above 1/4, meaning: entangled pairs start in the quantum regime and must decohere to cross the boundary. |
| **1/9** | Maximum reality for an entangled pair | R = CΨ² = (1/3)² = 1/9. The most "real" an entangled pair can be at the moment of maximum coherence. |

*Note: Larger systems have smaller Ψ for the same amount of raw coherence (because the possibility space grows exponentially with system size). All values above assume perfect initial preparation (C = 1).*

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
