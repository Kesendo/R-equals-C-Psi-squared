# Complete Mathematical Documentation

**Status:** Current as of March 22, 2026
**Supersedes:** Previous stub (Feb 2026) and [Core Algebra](CORE_ALGEBRA.md) (Dec 2025)
**Purpose:** Single entry point for all proven and verified mathematics of R=CΨ²

---

## 1. The Algebraic Foundation (Tier 1)

The self-referential purity map:

    R = CΨ²

where C = Tr(ρ²) (purity), Ψ = l₁(ρ)/(d-1) (normalized l1-coherence,
Baumgratz convention), R = residual purity beyond product-state prediction.

**Fixed-point equation.** R = C(Ψ + R)² expands to CR² + (2CΨ - 1)R + CΨ² = 0.

**Discriminant.** D = 1 - 4CΨ. Vanishes at CΨ = 1/4 and only there.
Below 1/4: two real fixed points (one stable, one unstable).
At 1/4: one degenerate fixed point (fold).
Above 1/4: no real fixed points (complex/oscillatory).

**Crossing cubic.** At the boundary CΨ = 1/4, the condition reduces to
b³ + b = 1/2, with unique real root b ≈ 0.4239. This is a pure number,
independent of physical parameters.

**Mandelbrot equivalence.** The substitution z = C(Ψ + R), c = CΨ maps
R_{n+1} = C(Ψ + R_n)² to the Mandelbrot iteration z_{n+1} = z_n² + c.
The boundary CΨ = 1/4 maps to the cusp of the main cardioid at c = 1/4.
This is identity, not analogy.

**Fold catastrophe.** The recursion is exactly the normal form of the fold
catastrophe (simplest in the Thom-Arnold classification). CΨ - 1/4 is the
bifurcation parameter. Structurally stable: no perturbation can remove it.

See: [Uniqueness Proof](UNIQUENESS_PROOF.md),
[Mathematical Connections](MATHEMATICAL_CONNECTIONS.md),
[Proof Roadmap](PROOF_ROADMAP_QUARTER_BOUNDARY.md)

---

## 2. The Palindromic Symmetry (Tier 1)

**Theorem.** For any Heisenberg/XXZ spin system with local Z-dephasing,
the Liouvillian spectrum is palindromic: for every eigenvalue λ, the value
-(λ + 2Σγ) is also an eigenvalue.

**The operator.** Π acts per site on Pauli indices:

    I → X (+1),  X → I (+1),  Y → iZ (+i),  Z → iY (+i)

**The proof.** Three steps:
1. Π anti-commutes with L_H (explicit 16-entry table for Heisenberg bonds)
2. Π transforms L_D: Π L_D Π⁻¹ = -L_D - 2Σγ I
3. Combined: Π L Π⁻¹ = -L - 2Σγ I. QED.

**Verification.** 54,118 eigenvalues, N=2 through N=8, zero exceptions.
All topologies (chain, star, ring, complete, binary tree). Non-uniform γ.

**Physical meaning.** Π is time reversal in a rescaled frame. It maps
exp(+μt) to exp(-μt), forward to backward.

**Scope boundary.** Only d=2 (qubits). The per-site split d immune vs
(d²-d) decaying is balanced only when d²-2d=0, giving d=2 uniquely.
Qutrits (d=3, split 3:6) verified broken for all 10 Hamiltonians tested.

See: [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md),
[Non-Heisenberg Palindrome](../experiments/NON_HEISENBERG_PALINDROME.md),
[Qubit Necessity](../hypotheses/QUBIT_NECESSITY.md)

---

## 3. The CΨ = 1/4 Boundary (Tier 1-2)

**Uniqueness.** 1/4 is algebraically unique. The factor 4 comes from the
discriminant formula b²-4ac. The quadratic structure comes from purity
being Tr(ρ²), degree 2. No reparameterization changes this.

**Channel independence.** All standard Markovian channels cross at exactly
CΨ = 0.2500:

| Channel | t_cross (γ=0.05) |
|---------|-------------------|
| Z-dephasing | 0.747 |
| X-noise | 1.733 |
| Y-noise | 1.733 |
| Depolarizing | 0.879 |
| Asymmetric Pauli | 0.735 |
| Amplitude damping | 2.059 |

**Absorbing boundary.** No unitary pulse (θ from 0 to π) can push CΨ back
above 1/4 after crossing. Tested, verified.

**Hardware validation.** IBM Torino (ibm_torino, Q80): predicted t* = 15.01 μs,
measured t* = 15.29 μs. Deviation: 1.9%.

See: [Uniqueness Proof](UNIQUENESS_PROOF.md),
[IBM Run 3](../experiments/IBM_RUN3_PALINDROME.md),
[proof_roadmap_close.py](../simulations/proof_roadmap_close.py)

---

## 4. The Incompleteness Proof (Tier 2)

Five candidates for the origin of dephasing noise, all eliminated:

1. Internal (bootstrap falsified, sectors decoupled)
2. Qubit decay (non-Markovian, 0/16 palindromic pairs)
3. Qubit bath (regress, each member faces bootstrap prohibition)
4. Nothing (d=0, no properties)
5. Other dimensions (d(d-2)=0 excludes)

**Corollary 1:** Time cannot originate from within. Noise IS the time arrow.
Without noise: reversible oscillation. With noise: irreversible decay.

**Corollary 2:** γ is the source of experienced time. t_cross = K/γ, the product
t × γ = K (a pure number). γ provides the arrow, J provides the content. Remove γ and t loses meaning.

See: [Incompleteness Proof](INCOMPLETENESS_PROOF.md),
[The Bridge Was Always Open](THE_BRIDGE_WAS_ALWAYS_OPEN.md),
[failed_third.py](../simulations/failed_third.py)

---

## 5. Engineering Results (Tier 2)

**Mediator bridge.** Mediated coupling (A-M-B) preserves palindrome
(1024/1024, error 1.41e-13) while information flows (MI = 1.65 bits,
QST fidelity 0.732). Direct coupling destroys it (256 → 31 pairs).

**Relay protocol.** Time-dependent γ, staged transfer: +83% end-to-end MI.

**V-shape γ gradient.** [0.01, 0.03, 0.05, 0.03, 0.01]: +124% MI.

**DD on M+Receiver.** Dynamical decoupling: +132% MI.

**Push vs Pull.** Source-dominant for local MI (0.957). Drain-dominant
for end-to-end MI (0.121). Distance-dependent.

**MI scaling.** Exponential decay: ~2-5 dB per 2 additional qubits.
Hierarchy falsified (uniform chain = recursive topology).

**Six design rules:** W-encoding, 2:1 matching, J/γ independence,
threshold timing, push/pull selection, clocked relay.

See: [Engineering Blueprint](../publications/ENGINEERING_BLUEPRINT.md),
[Circuit Diagram](../publications/CIRCUIT_DIAGRAM.md),
[Relay Protocol](../experiments/RELAY_PROTOCOL.md),
[Gamma Control](../experiments/GAMMA_CONTROL.md),
[Scaling Curve](../experiments/SCALING_CURVE.md)

---

## 6. The Transistor Mapping (Tier 2-3)

The mediator qubit M maps to a transistor: gate = γ_M, source = Pair A,
drain = Pair B. Threshold voltage: CΨ = 1/4 (hardwired, fold catastrophe).
Bidirectional by palindromic symmetry.

Three control knobs: γ_M (gate), J_AM/J_MB ratio (bias), κ (feedback gain).

Hierarchy falsified: the transistor properties are real, the recursive
scaling advantage is not.

See: [Quantum Transistor](../hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md)

---

## 7. Open Questions (Tier 3-5)

- Formal proof of CΨ monotonicity above 1/4 for arbitrary CPTP maps
- Feigenbaum period-doubling in the quantum regime (mapped, not exploited)
- Bekenstein-Hawking 1/4 (coincidence or connection, speculative)
- Negative feedback loop (γ_M decreasing with coherence, untested)
- Hardware validation of relay protocol on IBM Torino

See: [Mathematical Connections](MATHEMATICAL_CONNECTIONS.md),
[Proof Roadmap](PROOF_ROADMAP_QUARTER_BOUNDARY.md)

---

## 8. Numerical Constants

| Constant | Value | Source |
|----------|-------|--------|
| Discriminant zero | CΨ = 0.2500 | [Uniqueness Proof](UNIQUENESS_PROOF.md) |
| Crossing cubic root | b = 0.4239 | Cardano formula |
| K (Z-dephasing, N=2) | 0.037 | [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md) |
| IBM deviation | 1.9% | [IBM Run 3](../experiments/IBM_RUN3_PALINDROME.md) |
| Pauli weight correlation | r = 0.976 | [XOR Space](../experiments/XOR_SPACE.md) |
| Best QST fidelity | F = 0.888 | [QST Bridge](../experiments/QST_BRIDGE.md) |
| Relay improvement | +83% | [Relay Protocol](../experiments/RELAY_PROTOCOL.md) |
| V-shape improvement | +124% | [Gamma Control](../experiments/GAMMA_CONTROL.md) |
| DD M+Recv improvement | +132% | [Gamma Control](../experiments/GAMMA_CONTROL.md) |
| N=8 eigenvalues paired | 65,518 | [C# Compute](../compute/RCPsiSquared.Compute/) |
| Mediator bridge error | 1.41e-13 | [mediator_bridge.py](../simulations/mediator_bridge.py) |
| GHZ analytical match | delta < 1e-17 | [proof_roadmap_close.py](../simulations/proof_roadmap_close.py) |

---

## 9. References

### Proofs
- [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md)
- [Uniqueness Proof](UNIQUENESS_PROOF.md)
- [Incompleteness Proof](INCOMPLETENESS_PROOF.md)
- [Mathematical Connections](MATHEMATICAL_CONNECTIONS.md)

### Roadmaps
- [Proof Roadmap](PROOF_ROADMAP_QUARTER_BOUNDARY.md)
- [The Bridge Was Always Open](THE_BRIDGE_WAS_ALWAYS_OPEN.md)

### Publications
- [Technical Paper](../publications/TECHNICAL_PAPER.md)
- [Engineering Blueprint](../publications/ENGINEERING_BLUEPRINT.md)
- [Circuit Diagram](../publications/CIRCUIT_DIAGRAM.md)

### Key experiments
- [IBM Run 3](../experiments/IBM_RUN3_PALINDROME.md)
- [Non-Heisenberg Palindrome](../experiments/NON_HEISENBERG_PALINDROME.md)
- [Relay Protocol](../experiments/RELAY_PROTOCOL.md)
- [Scaling Curve](../experiments/SCALING_CURVE.md)
- [Gamma Control](../experiments/GAMMA_CONTROL.md)
- [Standing Wave Analysis](../experiments/STANDING_WAVE_ANALYSIS.md)

### Historical (superseded, kept for the record)
- [Core Algebra](CORE_ALGEBRA.md) (Dec 2025 / Feb 2026, predates all March results)
