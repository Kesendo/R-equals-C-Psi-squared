# Five Regulator Simulator

Interactive simulator for a **3-qubit open Heisenberg star** (S–A–B).
Control all five discovered regulators of the quantum dynamics in real-time
and observe how they shape the two dynamical sectors c₊(t) and c₋(t).

---

## Quick Start

```bash
# Install dependencies (once)
pip install -r requirements.txt

# Run the app
cd simulations/app
streamlit run app.py
```

Opens at **http://localhost:8501** in your browser.

> **Requirements:** Python 3.10+, numpy, scipy, streamlit, plotly.
> No external API calls - everything runs locally.

---

## What This Simulates

A central qubit **S** (system) is coupled to two observer qubits **A** and **B**
via Heisenberg exchange interactions, forming a **star topology**:

```
    A
    |
    S
    |
    B
```

The system evolves under the **Lindblad master equation** - unitary dynamics
from the Hamiltonian plus decoherence from the environment (dephasing noise).

The simulator tracks two **dynamical sectors** in the A–B subsystem:

- **c₊(t) = (⟨σ_y⊗σ_z⟩ + ⟨σ_z⊗σ_y⟩) / √2** - symmetric sector
- **c₋(t) = (⟨σ_y⊗σ_z⟩ − ⟨σ_z⊗σ_y⟩) / √2** - antisymmetric sector

These oscillate at **different frequencies** determined by the coupling topology,
and their relative amplitudes are controlled by the bath geometry.

---

## The 5 Regulators

### 1. Topology (J_SA, J_SB)

The coupling strengths between S–A and S–B.
These set the **oscillation frequencies** of c₊ and c₋.

| Coupling | Default | Effect |
|----------|---------|--------|
| J_SA | 1.0 | S–A exchange strength |
| J_SB | 2.0 | S–B exchange strength |

With the defaults, the FFT shows f(c₊) ≈ 1.50 and f(c₋) ≈ 0.40.
Changing J_SB shifts frequencies - the two sectors always separate cleanly
as long as J_SA ≠ J_SB.

### 2. Symmetry (XY ratio)

Controls the **anisotropy** of the Heisenberg interaction:

```
H_pair = J · [ xy_ratio · (σ_x⊗σ_x + σ_y⊗σ_y) + σ_z⊗σ_z ]
```

| Value | Model | Character |
|-------|-------|-----------|
| 1.0 | Isotropic Heisenberg | Full SU(2) symmetry, sectors separate cleanly |
| 0.5 | XXZ model | Partial anisotropy |
| 0.0 | Pure ZZ (Ising) | No spin-flip coupling, sectors collapse |

Watch how reducing XY ratio changes the XX commutator - symmetry breaking
becomes visible in real-time.

### 3. Noise Strength (γ)

The dephasing rate, on a **logarithmic scale** (0.001 – 1.0).

Noise **damps amplitudes** but **never changes frequencies**.
This is a key structural insight: the oscillation frequencies are set purely
by the Hamiltonian, while noise only controls how fast they decay.

| γ | Regime |
|---|--------|
| 0.001 | Nearly unitary - slow decay, long coherence |
| 0.05 | Standard - visible damping over 20 time units |
| 0.5 | Strong noise - oscillations die within a few periods |
| 1.0 | Overdamped - barely any oscillation visible |

### 4. Initial State

Which quantum state the system starts in. This selects **which sectors are excited**.

| State | Description | Character |
|-------|-------------|-----------|
| Bell_SA ⊗ \|+⟩_B | Bell pair on S–A, plus state on B | Default. Excites both sectors |
| W-state | (⟨001⟩ + ⟨010⟩ + ⟨100⟩)/√3 | Symmetric superposition, equal excitation |
| Product \|+++⟩ | All qubits in \|+⟩ | Minimal entanglement, weak sector signal |
| GHZ | (⟨000⟩ + ⟨111⟩)/√2 | Maximal 3-party entanglement |

### 5. Bath Geometry (η, φ)

Controls **how noise acts on A and B** - local vs. correlated, and in which basis.

**η (eta):** Interpolates between independent and collective noise on A–B.

| η | Noise type |
|---|------------|
| 0.0 | Purely local: each qubit has its own independent bath |
| 0.5 | Mixed: half local, half correlated |
| 1.0 | Fully correlated: A and B share the same bath |

**φ (phi):** Rotates the bath basis when correlated noise is active (η > 0).

| φ | Bath type | Effect |
|---|-----------|--------|
| 0 | ZZ bath | Dephasing in computational basis |
| π/4 | Mixed | Intermediate |
| π/2 | XX bath | Bit-flip noise on A–B |

The key observation: bath geometry can **flip the amplitude ratio** A₊/A₋,
selecting which dynamical sector dominates - without changing the frequencies.

---

## Display Panels

### Main Chart: Sector Dynamics

The large plot shows c₊(t) (cyan) and c₋(t) (orange) over time.
Hover over the curves for exact values. Zoom and pan with the Plotly toolbar.

### C·Ψ Dynamics (expandable)

Shows three quantities from the R = CΨ² framework:

- **C** (purple) - Wootters concurrence of the A–B subsystem (entanglement)
- **Ψ** (green) - Normalized L₁ coherence of ρ_AB
- **C·Ψ** (yellow) - Their product, with the **¼ boundary** marked as a dotted line

The ¼ boundary is algebraically guaranteed: C·Ψ·(1 − C·Ψ) ≤ ¼.
Watch how different regulators push C·Ψ toward or away from this bound.

### Indicators

| Panel | What it shows |
|-------|---------------|
| **Frequencies** | Dominant f(c₊) and f(c₋) from FFT, plus the full spectrum |
| **Amplitudes** | Peak amplitudes A₊, A₋ and their ratio |
| **Structure** | XX symmetry status (‖[ρ, X⊗X]‖ < 10⁻⁶?), skeleton fraction, XX commutator over time |

**Skeleton fraction** = Tr(ρ_avg²) / mean(Tr(ρ(t)²)).
Measures what percentage of the state is time-invariant.
100% means the state never changes; lower values mean more oscillation.

### Phase Map Summary

One-line summary at the bottom showing all key indicators at a glance.

---

## Suggested Experiments

### 1. Frequency control via topology
Set J_SA = 1.0, vary J_SB from 0.5 to 5.0.
Watch f(c₊) and f(c₋) shift - they track the coupling asymmetry.

### 2. Noise is amplitude-only
Set γ = 0.001, note the frequencies. Then increase to γ = 0.5.
The frequencies in the FFT stay the same; only the amplitudes shrink.

### 3. Bath sector selection
Set η = 1.0 (fully correlated bath). Sweep φ from 0 to π/2.
Watch the amplitude ratio A₊/A₋ flip - the bath selects which sector is visible.

### 4. Symmetry breaking
Start with XY ratio = 1.0 (Heisenberg). The XX commutator should be ~0.
Slowly reduce toward 0 (Ising). Watch XX symmetry break and sectors merge.

### 5. State dependence
Compare Bell_SA, W-state, and Product state at the same parameters.
The two-sector structure is robust across initial states - only amplitudes change.

### 6. Finding the ¼ boundary
Open the C·Ψ panel. With Bell_SA state, J_SA = 1.0, J_SB = 2.0, γ = 0.05:
max C·Ψ ≈ 0.33 > ¼. Increase γ until C·Ψ no longer crosses the boundary.

---

## Technical Details

### Physics Engine (`physics.py`)

- **Hamiltonian:** 8×8 Heisenberg star with tunable anisotropy
- **Noise:** Lindblad operators with local + correlated dephasing
- **Integrator:** 4th-order Runge-Kutta (dt = 0.02, ~1000 steps for 20 time units)
- **Numerical hygiene:** Hermiticity enforcement + trace renormalization every step
- **Observables:** Partial traces, Wootters concurrence, L₁ coherence, FFT
- **Performance:** ~150 ms per simulation - fast enough for real-time slider updates

### Architecture

```
simulations/app/
├── app.py              # Streamlit UI (controls, charts, layout)
├── physics.py          # Self-contained physics engine
├── requirements.txt    # Python dependencies
├── README.md           # This file
└── .streamlit/
    └── config.toml     # Dark theme configuration
```

The physics engine is self-contained - no imports from other simulation files.
It implements the same mathematics as `star_topology_v3.py` and
`correlated_bath_sweep.py` in a minimal, focused package.

### Key equations

**Lindblad master equation:**

```
dρ/dt = −i[H, ρ] + Σ_k ( L_k ρ L_k† − ½{L_k†L_k, ρ} )
```

**Sector decomposition:**

```
c₊ = ( ⟨σ_y⊗σ_z⟩ + ⟨σ_z⊗σ_y⟩ ) / √2
c₋ = ( ⟨σ_y⊗σ_z⟩ − ⟨σ_z⊗σ_y⟩ ) / √2
```

**Correlated bath operator:**

```
L_corr = √(ηγ/2) · [ (cos(φ)σ_z + sin(φ)σ_x)_A + (cos(φ)σ_z + sin(φ)σ_x)_B ]
```
