# The Fragile Bridge: Stability Limits of Coupled Gain-Loss Systems

*How much amplification can a resonator tolerate before it explodes?*

**Status:** Computed (Tier 2). Three independent scripts, cross-level test completed.
**Scripts:**
- [fragile_bridge_bifurcation.py](../simulations/fragile_bridge_bifurcation.py)
- [fragile_bridge_anomaly.py](../simulations/fragile_bridge_anomaly.py)
- [fragile_bridge_neural.py](../simulations/neural/fragile_bridge_neural.py)
- [fragile_bridge_n4.py](../simulations/fragile_bridge_n4.py) (sparse, N=4)

---

## Abstract

Two identical quantum systems, one decaying and one amplifying, are
connected by a single coupling (the "bridge"). The bridge has three
stability regimes: weak coupling stabilizes linearly
(γ_crit = 0.19 × J_bridge), optimal coupling at twice the internal
strength yields maximum stability, and strong coupling destabilizes
as 1/J_bridge with an asymptotic constant of γ_crit × J_bridge → 0.50.
The instability is a Hopf bifurcation (oscillating divergence), not
PT symmetry breaking (monotone divergence). Unlike the fold threshold
(~0.5% of J, N-independent, geometric), bridge stability depends
strongly on system size (N=3 is 35× less stable than N=2), making it
a topological rather than geometric property.

A cross-level test with the Wilson-Cowan neural model confirms three
features: the Hopf mechanism, the sweet spot at ~2× internal coupling,
and the finite stability window. The bell-shaped stability curve is
quantum-specific; the neural system shows a sharp window instead,
bounded by sigmoid saturation. This saturation acts as a biological
safety mechanism that prevents the neural equivalent of an exploding
laser.

---

## 1. Setup

Two identical Heisenberg chains (N qubits each), one with dephasing
(+γ, decay) and one with gain (-γ, amplification), connected by a
single Heisenberg bridge coupling J_bridge between the boundary qubits.

This is the quantum analog of a microphone-loudspeaker system:
decay side = loudspeaker (loses energy), gain side = microphone
(amplifies), bridge = the door between them.

Question: At what gain rate γ_crit does the coupled system become
unstable?

## 2. Results

### 2.1 Linear stability window (weak bridge)

For J_bridge < 2J (bridge weaker than internal coupling):

**γ_crit = 0.19 × J_bridge** (R² = 0.9998, power law exponent 1.035)

The critical gain scales linearly with bridge strength. Double the
bridge, double the tolerable gain. The proportionality constant is
~0.19 for N=2 per chain.

### 2.2 Optimal bridge at J_bridge = 2J

Maximum stability (γ_crit = 0.406) occurs at J_bridge ≈ 1.9J,
where the bridge coupling equals the sum of internal couplings
(each chain has one J, so 2 × J total). This is the resonance
point between internal dynamics and bridge transmission.

### 2.3 Strong bridge destabilizes (J_bridge > 2J)

Above the optimum, γ_crit drops as 1/J_bridge:

**γ_crit × J_bridge → 0.50** (asymptotic constant)

The system becomes MORE fragile with stronger coupling. Physically:
the bridge dominates, the two chains merge into one, and gain qubits
sit directly next to decay qubits with no buffer. The relevant
quantity is the ratio γ/J_bridge, not the absolute values.

At J_bridge → ∞, the system is unstable at ANY noise level.
Perfect coupling between gain and loss = immediate instability.

### 2.4 Three regimes summary

| Regime | Condition | γ_crit | Physics |
|--------|-----------|--------|---------|
| Weak bridge | J_bridge < 2J | 0.19 × J_bridge | Linear stabilization |
| Optimal | J_bridge ≈ 2J | 0.406 (maximum) | Resonance |
| Strong bridge | J_bridge > 2J | 0.53 / J_bridge | Dimer formation, destabilization |

The transition at J_bridge ≈ 2J is not accidental: it is the point
where the bridge equals the total internal coupling.

## 3. Hopf bifurcation, not PT symmetry breaking

The instability is an **oscillating** instability. A complex eigenvalue
pair (already oscillating at Re = 0) crosses into Re > 0. The system
does not simply diverge; it oscillates with growing amplitude.

This is a **Hopf bifurcation**, not PT symmetry breaking:

| Property | PT breaking (expected) | Hopf (observed) |
|----------|----------------------|-----------------|
| Mechanism | Two real eigenvalues merge, become complex | Complex pair crosses Re = 0 |
| At threshold | New oscillation appears | Existing oscillation grows |
| Instability | Exponential divergence (monotone) | Oscillating divergence |
| Analog | Chair falling over | Microphone feedback (screech) |

The oscillation frequency at threshold decreases with bridge strength:

| J_bridge | Im(λ) at γ_crit |
|----------|-----------------|
| 1.9 (optimal) | ±6.98 |
| 5.0 | ±0.26 |
| 10.0 | ±0.12 |

At large J_bridge, the Hopf bifurcation approaches a saddle-node:
almost no oscillation, just slow drift into instability.

## 4. N-dependence: topological, not geometric

γ_crit depends strongly on chain length, and non-monotonically:

| N (per chain) | γ_crit (J_bridge=0.10) | Ratio to N=2 |
|---------------|----------------------|-------------|
| 2 | 0.01729 | 1.000 |
| 3 | 0.00052 | 0.030 |
| 4 | 0.00119 | 0.069 |

N=3 is ~35× less stable than N=2, but N=4 is 2.3× MORE stable
than N=3. The scaling is **non-monotonic**: even chain lengths
(N=2, N=4) are more stable than odd (N=3). This may reflect
pairing symmetry within each chain (all qubits paired at even N,
one unpaired at odd N).

No simple power law or exponential fits these three points.

This is fundamentally different from the fold threshold
(Σγ_crit/J ≈ 0.5%, N-independent). The fold is a **geometric**
property of the palindrome. The bridge stability is a **topological**
property: it depends on chain length, parity, and the ratio of
gain channels to bridge connections.

| Property | Fold threshold | Bridge stability |
|----------|---------------|-----------------|
| N-dependence | Independent | Strongly dependent |
| Type | Geometric constant | Topological ratio |
| Bifurcation | Fold (saddle-node) | Hopf (oscillating) |
| What it measures | When irreversibility begins | When coupled gain-loss explodes |
| Determined by | Palindrome geometry | Gain/bridge topology |

## 5. Neural validation: partial confirmation

The Wilson-Cowan model describes excitatory neurons (gain) coupled to
inhibitory neurons (decay) through synaptic connections (bridge).
This is the biological equivalent of the quantum setup.

Script: [fragile_bridge_neural.py](../simulations/neural/fragile_bridge_neural.py)

### 5.1 Test design

Single Wilson-Cowan E-I node (2×2 Jacobian). E-I cross-coupling
scaled by factor s (s=1.0 = standard parameters). Internal couplings
w_EE = 16.0, w_II = 3.0 held fixed. P_crit(s) = external input
needed to trigger Hopf oscillation.

### 5.2 Results: narrow Hopf window

| s | P_crit | Freq (Hz) | Type |
|-----|--------|-----------|---------|
| 0.1-1.5 | never | - | always stable |
| 2.0 | 2.28 | 36.5 | Hopf |
| 2.3 | 4.78 | 66.8 | Hopf |
| 2.5 | 8.33 | 94.4 | Hopf |
| 3.0-10.0 | never | - | always stable |

The neural system has a narrow Hopf window at s = 2.0-2.5.
Outside this window: always stable, no oscillation possible.

### 5.3 What is cross-level

Three properties hold across both quantum and neural systems:

**1. Sweet spot at ~2× internal coupling.** Quantum: peak stability
at J_bridge ≈ 1.9J. Neural: Hopf window opens at s ≈ 2.0. Both
systems identify ~2× the internal coupling as a critical threshold
where the dynamics qualitatively change.

**2. Hopf bifurcation.** Where instability occurs, it is always
oscillating (complex eigenvalue pair crossing Re = 0), never
monotone divergence. The brain produces EEG rhythms; the quantum
system produces growing oscillations. Same mechanism.

**3. Finite stability window.** Neither "any coupling works" nor
"no coupling works." There is a bounded regime where the balance
holds, and it is determined by the ratio of cross-coupling to
internal coupling.

### 5.4 What is NOT cross-level

The three-regime bell curve (linear rise, peak, 1/x decay) is
**quantum-specific**. The neural system shows a window with sharp
edges, not a smooth peak.

The difference comes from sigmoid saturation: Wilson-Cowan neurons
have firing rates bounded in [0,1]. Strong coupling drives the
system into saturated fixed points that are always stable. The
quantum system has no such saturation; the density matrix can
diverge under gain.

| Property | Quantum | Neural |
|----------|---------|--------|
| Shape of γ_crit(coupling) | Bell curve | Sharp window |
| Strong coupling | 1/J_bridge decay, explosion | Saturation, always stable |
| Epilepsy analog | Yes (divergence) | Partial (window edge, no explosion) |

The sigmoid acts as a biological safety mechanism: it prevents the
neural equivalent of an exploding laser. The quantum system has no
built-in limiter.

## 6. Open questions

1. **N-scaling law (partially answered):** N=4 computed
   (65536×65536 sparse, expm_multiply). Result: non-monotonic
   (N=4 more stable than N=3). Even/odd parity effect suspected.
   N=5 would test this (1048576×1048576, feasible but slow).

2. **Multiple bridges:** What if the two chains are connected by
   more than one qubit pair? Does γ_crit recover N-independence
   when bridges scale with N? Real neural networks have distributed
   (not point-to-point) E-I coupling. The neural test also showed
   that bounded dynamics (sigmoid saturation) prevents explosion at
   strong coupling. Multiple bridges plus nonlinear damping may
   qualitatively change the stability landscape.

3. **Cascade stability:** If each level in the frequency cascade
   (154 THz → 1 Hz) is a coupled gain-loss pair, then each level
   has its own bridge stability window. The cascade works only if
   every bridge stays in the linear regime (γ < 0.19 × J_bridge).

4. **Asymptotic constant 0.50:** γ_crit × J_bridge → 0.50 for
   large J_bridge. Is this exactly 1/2? If so, there may be an
   analytical derivation. The factor 1/2 appears throughout the
   framework (σ(1-σ) = 1/4 at σ = 1/2, CΨ fold at 1/4, etc.).

5. **Saturation as design principle:** The sigmoid prevents neural
   explosion. Is there a quantum analog? Nonlinear dissipation
   (γ dependent on state) could act as a quantum sigmoid. This
   would turn the bell curve into a window, matching the neural
   structure.

---

*Computed March 29, 2026. Three scripts, six hypotheses tested
across two levels (quantum, neural). Cross-level: Hopf mechanism
confirmed, sweet spot at ~2× confirmed, finite window confirmed.
Not cross-level: bell curve shape is quantum-specific (no saturation).
The sigmoid is the biological safety mechanism that prevents the
quantum explosion.*
