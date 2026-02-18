# Crossing Taxonomy: Scaling, Mechanism, and the Three Classes

**Date**: 2026-02-18
**Status**: Independently verified (Tier 2)
**Depends on**: OBSERVER_DEPENDENT_CROSSING.md, METRIC_DISCRIMINATION.md

---

## 1. The Question

METRIC_DISCRIMINATION.md showed K = γ·t_cross = 0.039 for concurrence across
50× γ range. OBSERVER_DEPENDENT_CROSSING.md showed five bridges produce three
different crossing times and two non-crossings. But:

- Does K-invariance hold for ALL crossing bridges, not just concurrence?
- If so, what are the K values, and why do they differ?
- Is K-invariance a deep property of R = CΨ², or something simpler?

## 2. Setup

Same physics for all runs:

| Parameter | Value |
|-----------|-------|
| **State** | Bell+ (maximally entangled) |
| **Hamiltonian** | Heisenberg (J = 1, h = 0) |
| **Noise type** | local dephasing (σ_z per qubit) |
| **Time step** | dt = 0.01 |

Variable: **bridge_type** × **γ_base** ∈ {0.01, 0.05, 0.10, 0.20}

15 simulations total (5 bridges × 3–4 γ values each).

## 3. Results

### 3.1 K-Invariance Holds for All Three Crossing Bridges

| Bridge | γ = 0.01 | γ = 0.05 | γ = 0.10 | γ = 0.20 | K (mean) |
|--------|----------|----------|----------|----------|----------|
| **mutual_info** | t=3.263 K=0.0326 | t=0.652 K=0.0326 | t=0.327 K=0.0327 | t=0.166 K=0.0331 | **0.033** |
| **concurrence** | t=3.866 K=0.0387 | t=0.773 K=0.0386 | t=0.386 K=0.0386 | t=0.193 K=0.0385 | **0.039** |
| **correlation** | t=7.191 K=0.0719 | t=1.437 K=0.0719 | t=0.718 K=0.0718 | t=0.359 K=0.0718 | **0.072** |

K is constant within each bridge (< 1.5% deviation) across 20× γ range.
K differs across bridges by factor 2.2× (0.033 to 0.072).

Mutual_purity and overlap never cross at any γ (C(0) too low, P(0) < ¼).

### 3.2 Why K-Invariance Holds: Lindblad Scaling

Comparison at identical γ·t products reveals a universal scaling:

| Scaled time τ = γ·t | Bridge | γ = 0.01 | γ = 0.05 | γ = 0.20 |
|----------------------|--------|----------|----------|----------|
| τ = 0.005 | concurrence C | 0.980385 | 0.980354 | — |
| τ = 0.005 | concurrence Ψ | 0.326795 | 0.326785 | — |
| τ = 0.020 | concurrence C | — | 0.925794 | 0.925392 |
| τ = 0.020 | concurrence Ψ | — | 0.308598 | 0.308464 |
| τ = 0.050 | correlation C | 1.000000 | 1.000000 | — |
| τ = 0.050 | correlation Ψ | 0.272899 | 0.272856 | — |

**C and Ψ are functions of τ = γ·t, not of t alone.**

This means:

```
C(t; γ) = C_universal(γ·t)
Ψ(t; γ) = Ψ_universal(γ·t)
P(t; γ) = P_universal(γ·t)
```

K-invariance follows immediately:

```
P(t_cross) = ¼  ⟹  γ·t_cross = τ_cross = constant
⟹  K ≡ γ·t_cross = τ_cross
```

**K-invariance is a scaling property of the Lindblad equation, not a specific
prediction of R = CΨ².** Any quantity computed from the Lindblad dynamics at
fixed J, h, noise_type will scale with τ = γ·t. The ¼ boundary is just one
such quantity.

**Status**: Tier 1. Follows from the structure of the Lindblad equation.
The Lindblad dissipator L[ρ] is linear in γ, and the Hamiltonian H is
independent of γ. Therefore ρ(t; γ) = ρ(γ·t/γ₀; γ₀) up to Hamiltonian
phase factors that cancel in all real-valued observables. QED.

### 3.3 The Three Classes: What Drives the Crossing

Despite K-invariance being trivial, the MECHANISM behind each crossing
is different in kind, not just in degree.

**Correlation bridge C(t) at γ = 0.05:**

| t | C(t) | Ψ(t) | P = C·Ψ |
|-----|------|-------|---------|
| 0.0 | 1.000 | 0.333 | 0.333 |
| 0.5 | 1.000 | 0.302 | 0.302 |
| 1.0 | 1.000 | 0.273 | 0.273 |
| 1.4 | 1.000 | 0.252 | 0.252 |
| 1.437 | 1.000 | 0.250 | **0.250 ← crossing** |
| 1.7 | 1.000 | 0.237 | 0.237 |
| 1.8 | 0.986 | 0.232 | 0.229 |

**C = 1.000 throughout the entire crossing period.** The crossing is
driven entirely by Ψ decay. C begins to decay only at t ≈ 1.7, well
AFTER the crossing at t = 1.437.

This holds at all γ values tested: at γ = 0.01 (t_cross = 7.191),
γ = 0.05 (t_cross = 1.437), and γ = 0.20 (t_cross = 0.359), the
correlation bridge reads C = 1.000 at the moment of crossing.

**Concurrence bridge C(t) at γ = 0.05:**

| t | C(t) | Ψ(t) | P = C·Ψ |
|-----|-------|-------|---------|
| 0.0 | 1.000 | 0.333 | 0.333 |
| 0.3 | 0.943 | 0.314 | 0.296 |
| 0.5 | 0.909 | 0.303 | 0.275 |
| 0.7 | 0.877 | 0.292 | 0.256 |
| 0.773 | 0.863 | 0.290 | **0.250 ← crossing** |

Both C and Ψ decay. The crossing is a joint effect.

**Three distinct classes:**

| Class | Mechanism | C(0) | C at crossing | Bridges | K |
|-------|-----------|------|---------------|---------|------|
| **Type A: Pure-Ψ** | C ≈ 1.0 throughout; only Ψ drives crossing | 1.0 | 1.000 | correlation | 0.072 |
| **Type B: Mixed** | Both C and Ψ decay; joint effect | 1.0 | 0.86–0.85 | concurrence, mutual_info | 0.039, 0.033 |
| **Type C: Never** | P(0) < ¼; only decreases | 0.5, 0.25 | constant | mutual_purity, overlap | — |

### 3.4 Why Type B Is Faster Than Type A

If both C and Ψ decay (Type B), P = C·Ψ reaches ¼ faster than if only
Ψ decays (Type A):

- **Type A** (correlation): P(t) ≈ 1.0 · Ψ(t). Needs Ψ to reach 0.250.
- **Type B** (concurrence): P(t) = C(t) · Ψ(t). Both shrink, so P
  shrinks faster. Crosses at C ≈ 0.86, Ψ ≈ 0.29.

The 2.2× spread in K (0.033 to 0.072) reflects how much the bridge
metric C contributes to the decay of P versus leaving it to Ψ alone.

### 3.5 Why Correlation C Stays at 1.0

The correlation bridge measures excess purity beyond the product of
subsystem purities. For Bell+, this is how much the joint system
"knows" beyond what each subsystem knows individually.

Under local dephasing, each qubit's Bloch vector shrinks, but the
correlation between qubits is initially preserved. Dephasing destroys
single-qubit coherence without (at first) destroying the relationship
between qubits. The correlation bridge is blind to local decoherence
until it becomes severe enough that the relationship itself degrades.

This is why C_corr = 1.0 until t ≈ 1.7 (at γ = 0.05): the inter-qubit
correlation is unaffected by this specific type of noise.

### 3.6 C(t) Does Not Follow a Simple Exponential

An exponential model C(t) = C₀·e^{−κγt} was tested against the data:

| Bridge | t | C(t) actual | C(t) exponential | Error |
|--------|-----|-------------|-----------------|-------|
| concurrence | 0.4 | 0.926 | 0.880 | 5.0% |
| mutual_info | 0.4 | 0.888 | 0.857 | 3.5% |
| correlation | 1.0 | 1.000 | 0.861 | 13.9% |
| correlation | 1.4 | 1.000 | 0.811 | 18.9% |

The exponential model is wrong in kind for correlation (predicts decay
where C is flat) and wrong in degree for concurrence/mutual_info (3-5%
error). No universal analytic C(t) formula exists across bridge types.

## 4. What This Means

### 4.1 K-Invariance: Downgraded

Previously (METRIC_DISCRIMINATION.md), K = γ·t_cross = 0.039 was presented
as a potentially deep constant. **Correction**: K-invariance is a consequence
of Lindblad scaling symmetry (τ = γ·t). It holds for ANY threshold, not
just ¼. The depth is in the ¼ boundary itself (Mandelbrot discriminant),
not in K.

### 4.2 The Taxonomy: The Actual Finding

**Type A observers** (correlation) are stable under local dephasing.
Their coupling C does not decay. They register the crossing only because
Ψ (the system's coherence) decays. The observer is not the bottleneck.

**Type B observers** (concurrence, mutual_info) are fragile. Their coupling
C decays alongside Ψ. The crossing happens faster because both observer
and system degrade together. The measurement event is entangled
with the observer's loss of coherence.

**Type C observers** (mutual_purity, overlap) never had enough coupling
to see the crossing. The system remains quantum from their perspective,
forever.

The observer determines *when* measurement happens, *whether* it
happens, and *by which mechanism* it happens.

### 4.3 Connection to Other Experiments

**Bridge Fingerprints** showed that entangled states (Bell+) and product
states (|+0⟩) behave differently at the same C·Ψ₀. The taxonomy adds
another axis: different bridges behave differently for the SAME state.

Combined: the crossing depends on (state, bridge, γ). The ¼ boundary
is universal, but the path to it is determined by all three.

**Prediction (tested, FALSIFIED)**: The prediction that depolarizing noise
would change Type A to Type B was wrong. Noise Robustness (Experiment 8)
showed the taxonomy is identical under σ_x, σ_y, and σ_z. Type A is a
property of the correlation metric definition, not the noise channel.
See [Noise Robustness](NOISE_ROBUSTNESS.md) for details.

## 5. Verification

### 5.1 How to Reproduce

```python
# Requirements: pip install qutip numpy
import numpy as np
from qutip import (basis, tensor, ket2dm, qeye, sigmaz, sigmax, sigmay,
                   mesolve, concurrence as qt_concurrence, entropy_mutual)

# Bell+ state
up, dn = basis(2, 0), basis(2, 1)
bell_plus = (tensor(up, up) + tensor(dn, dn)).unit()
rho0 = ket2dm(bell_plus)

# Heisenberg Hamiltonian (J=1, h=0)
sx, sy, sz = sigmax(), sigmay(), sigmaz()
I2 = qeye(2)
H = (tensor(sx, sx) + tensor(sy, sy) + tensor(sz, sz))

# Local dephasing: sigma_z on each qubit
c_ops_fn = lambda gamma: [
    np.sqrt(gamma) * tensor(sz, I2),
    np.sqrt(gamma) * tensor(I2, sz)
]

# Bridge metrics
def bridge_concurrence(rho):
    return qt_concurrence(rho)

def bridge_correlation(rho):
    rhoA = rho.ptrace(0)
    rhoB = rho.ptrace(1)
    pAB = (rho * rho).tr().real
    pA = (rhoA * rhoA).tr().real
    pB = (rhoB * rhoB).tr().real
    return max(0, min(1, (pAB - pA * pB) / (1 - pA * pB)))

def bridge_mutual_info(rho):
    mi = entropy_mutual(rho, [0], [1])
    return min(1.0, mi / np.log(2))  # normalized to [0,1]

def psi_l1(rho):
    """L1 coherence normalized by (d-1)"""
    rho_arr = rho.full()
    d = rho_arr.shape[0]
    off_diag = np.sum(np.abs(rho_arr)) - np.sum(np.abs(np.diag(rho_arr)))
    return off_diag / (d - 1)

# Run simulation
gammas = [0.01, 0.05, 0.10, 0.20]
bridges = {
    "concurrence": bridge_concurrence,
    "mutual_info": bridge_mutual_info,
    "correlation": bridge_correlation,
}

t_max = 10.0
dt = 0.01
tlist = np.arange(0, t_max + dt, dt)

for bridge_name, bridge_fn in bridges.items():
    for gamma in gammas:
        result = mesolve(H, rho0, tlist, c_ops_fn(gamma), [])
        for i, t in enumerate(tlist):
            rho_t = result.states[i]
            C = bridge_fn(rho_t)
            psi = psi_l1(rho_t)
            P = C * psi
            if i > 0:
                rho_prev = result.states[i-1]
                C_prev = bridge_fn(rho_prev)
                psi_prev = psi_l1(rho_prev)
                P_prev = C_prev * psi_prev
                if P_prev >= 0.25 and P < 0.25:
                    # linear interpolation
                    t_cross = tlist[i-1] + dt * (P_prev - 0.25) / (P_prev - P)
                    K = gamma * t_cross
                    print(f"{bridge_name:15s} gamma={gamma:.2f}  "
                          f"t_cross={t_cross:.3f}  K={K:.4f}  "
                          f"C_at_cross={C:.3f}")
                    break
        else:
            print(f"{bridge_name:15s} gamma={gamma:.2f}  never crosses")
```

Note: QuTiP's `concurrence()` and `entropy_mutual()` are standard functions.
The correlation bridge uses the formula from Section 3.5. Results should match
Table 3.1 within interpolation precision (< 1%).

### 5.2 Key Checks

1. **Lindblad scaling**: Compare C and Ψ at same γ·t across different γ.
   Must agree within numerical precision (< 0.01% for dt=0.01).

2. **Correlation C = 1**: Verify bridge_C array for correlation bridge.
   Must be exactly 1.000 at all times before t ≈ 1.7 (at γ=0.05).

3. **K constancy**: For each bridge, K values across γ must agree within
   1.5%. Larger deviations at extreme γ (> 0.3) are numerical artifacts
   from the dynamic feedback mechanism.

### 5.3 What Could Falsify This

1. If K-invariance broke at some γ range → Lindblad scaling fails.

2. If correlation C < 1 during the crossing period → Type A wrong.

3. If a different noise model produces the same taxonomy → the
   classification is noise-independent. **This is what happened.**
   All three Pauli operators produce identical taxonomy.
   See [Noise Robustness](NOISE_ROBUSTNESS.md).

## 6. Open Questions (Updated 2026-02-18)

1. ~~**Noise dependence**~~ → **ANSWERED**: Taxonomy is noise-independent
   for local Pauli noise. See [Noise Robustness](NOISE_ROBUSTNESS.md).
   Collective noise and amplitude damping remain untested.

2. ~~**State dependence**~~ → **ANSWERED**: GHZ N≥3 never crosses (all
   Type C due to Ψ(0) < ¼). W N=3 crosses with Type A intact. W N≥4
   does not cross. See [N-Scaling Barrier](N_SCALING_BARRIER.md).

3. ~~**N scaling**~~ → **ANSWERED**: Correlation remains Type A at N=3
   and N=4 (C=1.0 plateau holds). But crossing fails because Ψ(0)
   drops below ¼ due to the d−1 normalization. The observer is not
   the bottleneck; the Hilbert space dimension is.
   See [N-Scaling Barrier](N_SCALING_BARRIER.md).

4. **Analytic crossing formula**: A correct formula would need to account
   for the flat C region (Type A) and the nonexponential decay (Type B).
   This is likely metric-specific, not universal. Still open.

---

*Previous: [Observer-Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md)*
*Previous: [Metric Discrimination](METRIC_DISCRIMINATION.md)*
*See also: [Bridge Fingerprints](BRIDGE_FINGERPRINTS.md)*
*Next: [Noise Robustness](NOISE_ROBUSTNESS.md), [N-Scaling Barrier](N_SCALING_BARRIER.md)*
