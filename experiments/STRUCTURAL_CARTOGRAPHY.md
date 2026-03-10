# Structural Cartography of CΨ Windows

**Tier:** Exploratory framework (methodology, not results)
**Status:** Working document defining analysis approach
**Scope:** How to characterize the structure inside CΨ visibility windows
**Does not establish:** That the windows contain "language", "messages", or semantic content
**Date:** 2026-03-08

---

## Principle

Map first. Interpret later.

CΨ visibility windows contain structured reduced states. Before asking what
they mean, we characterize what they are: what varies, what persists, what
predicts the next window, and what part of the evolution is coherent versus noisy.

## Working thesis

The windows do not behave like unrelated snapshots. They look like a
low-dimensional, partially coherent, memory-bearing process observed through
a reduced subsystem view. The right formal question is:

> What is the minimal predictive description of the CΨ window sequence?

## Formal framing

CΨ windows are treated as samples from a **multi-time quantum process** with
structure, memory, and noise. The formal language is process tensor / multi-time
open quantum dynamics (Pollock et al. 2018, arXiv:1512.00589).

This lets us ask: what is the minimal predictive description of the CΨ window
sequence?

## Analysis layers

### Layer 1: Per-window state chart

Each window gets embedded in a common feature space:

- **Bell fidelity vector**: F = [F_Phi+, F_Phi-, F_Psi+, F_Psi-]
- **Purity**: Tr(rho^2)
- **Von Neumann entropy**: -Tr(rho log2 rho)
- **Entanglement**: concurrence, CΨ
- **Coherence**: normalized l1, selected off-diagonal magnitudes
- **Phase coordinates**: arg(rho_00,11), arg(rho_01,10), unwrapped across windows
- **Symmetry residuals**: |rho_00 - rho_11| + |rho_01 - rho_10|
- **Gate coordinate**: S-coherence or equivalent gating scalar

### Layer 2: Inter-window transitions

How one window changes into the next:

- Trace distance between adjacent windows
- State fidelity between adjacent windows
- Bell-vector drift: ||F^(n+1) - F^(n)||
- Phase advance: unwrapped phase change per window
- Symmetry drift

### Layer 3: Sequence-level information

The "grammar" layer (used as metaphor, not literal claim):

- **Bell entropy per window**: H = -sum(F_i log2 F_i)
- **Block entropy H(L)**: entropy of L-window blocks after discretization
- **Entropy rate**: h ~ H(L) - H(L-1) = new information per window
- **Excess entropy**: E = I(past; future) = predictive structure
- **Adjacent mutual information**: I(X_n; X_{n+1})
- **Effective memory length**: smallest k where conditioning saturates

### Layer 4: Coherent vs noisy decomposition

- Run at gamma=0 as coherent baseline
- Compare noisy trajectory against baseline
- Quantify residual: trace distance, infidelity, purity loss, phase disruption
- Ask: how much of the window structure is Hamiltonian, how much is noise?

## Terminology

### Use these terms
- multi-time quantum process
- reduced-state process
- CΨ window sequence
- Bell-centered feature chart
- predictive transition structure
- coherent drift / dissipative residual
- structural cartography
- persistent symmetries
- memory length / memory complexity

### These are metaphors (label them as such)
- alphabet = minimal distinguishable feature set
- grammar = transition constraints across windows
- morphology = transform pattern window to window
- utterance = one observation window

### Avoid as primary framing
- protocol, packet, header, payload
- semantic content, meaning of the window
- "the system communicates messages"

## What we already know from today's observations

From WHATS_INSIDE_THE_WINDOWS.md (Tier 2, computationally verified):

- Bell fidelity vector changes window to window (Phi+ dominant, decaying)
- Populations stay symmetric (persistent symmetry)
- Phases rotate systematically (coherent drift)
- Fidelity decays 0.78 -> 0.58 over 8 windows (dissipative contribution)
- S-coherence gates window openness (gate coordinate)
- Phase transport from S to AB is real but standard Heisenberg
- Transport works in both open and closed windows (channel always on)
- Concurrence alone predicts sensitivity better than CΨ (Bridge Test B=-0.024)
- But sensitivity is not comprehension - we do not know what the content is

## Next concrete steps

### Phase A: Descriptive cartography
Compute the full per-window feature chart for the existing star topology run.
Plot all Layer 1 and Layer 2 metrics. Look for low-dimensional structure.

### Phase B: Memory and predictability
Discretize windows into coarse symbols (dominant Bell state, phase bin).
Compute block entropies, entropy rate, adjacent mutual information.
Estimate effective memory length.

### Phase C: Coherent vs noisy split
Run gamma=0 baseline. Compare. Quantify how much structure is Hamiltonian
(predictable, "grammatical") versus noise (unpredictable, "entropic").

## Safe claims

- CΨ windows are treated as a structured reduced-state process
- The analysis distinguishes descriptive metrics from sequence-level information
- The working hypothesis is that windows lie near a low-dimensional coherent
  drift manifold with noise as dissipative residual
- The objective is structural cartography before interpretation

## Phase A Results: Descriptive Cartography (March 2026)

Computed on star topology, Bell_SA x |+>_B, J_SA=1.0, J_SB=2.0, gamma=0.05.
9 CΨ_AB visibility peaks analyzed. Script: simulations/cartography_phase_a.py

### Layer 1: Per-window feature chart

| # | t | CPsi | Phi+ | Phi- | Psi+ | Psi- | H_Bell | Pur | SvN | C | Psi | ph03 | S_coh | XX_norm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.24 | 0.305 | 0.784 | 0.040 | 0.028 | 0.148 | 1.01 | 0.680 | 0.82 | 0.634 | 0.470 | +0.00pi | 0.644 | 8.6e-17 |
| 1 | 0.40 | 0.329 | 0.741 | 0.098 | 0.075 | 0.086 | 1.23 | 0.682 | 0.83 | 0.617 | 0.515 | +0.00pi | 0.656 | 1.0e-16 |
| 2 | 0.70 | 0.123 | 0.223 | 0.301 | 0.292 | 0.184 | 1.97 | 0.470 | 1.20 | 0.228 | 0.489 | +1.00pi | 0.008 | 1.3e-16 |
| 3 | 1.06 | 0.168 | 0.575 | 0.140 | 0.095 | 0.189 | 1.63 | 0.510 | 1.22 | 0.332 | 0.494 | -0.00pi | 0.393 | 1.5e-16 |
| 4 | 1.58 | 0.151 | 0.604 | 0.137 | 0.084 | 0.176 | 1.57 | 0.495 | 1.32 | 0.336 | 0.434 | +0.00pi | 0.462 | 1.8e-16 |
| 5 | 1.96 | 0.036 | 0.244 | 0.299 | 0.263 | 0.194 | 1.98 | 0.423 | 1.42 | 0.077 | 0.421 | +1.00pi | 0.003 | 2.2e-16 |
| 6 | 2.26 | 0.138 | 0.687 | 0.143 | 0.058 | 0.112 | 1.36 | 0.530 | 1.29 | 0.400 | 0.333 | -0.00pi | 0.646 | 2.4e-16 |
| 7 | 2.40 | 0.122 | 0.617 | 0.166 | 0.079 | 0.138 | 1.54 | 0.482 | 1.38 | 0.315 | 0.377 | -0.00pi | 0.525 | 2.9e-16 |
| 8 | 3.06 | 0.083 | 0.553 | 0.212 | 0.112 | 0.123 | 1.67 | 0.447 | 1.46 | 0.225 | 0.351 | -0.00pi | 0.461 | 3.2e-16 |

X tensor X symmetry exact at all windows (norm < 3.2e-16). Confirmed.

### Layer 2: Inter-window transitions

Two distinct regimes of transition:

**Smooth transitions** (small Bell-drift, high fidelity):
- 0->1: TD=0.438, Fid=0.791, BellDr=0.106 (within + sector)
- 3->4: TD=0.359, Fid=0.851, BellDr=0.034
- 6->7: TD=0.298, Fid=0.910, BellDr=0.081
- 7->8: TD=0.184, Fid=0.944, BellDr=0.088

**Hard sector switches** (large Bell-drift, phase jump of pi):
- 1->2: TD=0.648, Fid=0.562, BellDr=0.605, dph=+pi (+ sector to balanced)
- 2->3: TD=0.677, Fid=0.527, BellDr=0.434, dph=-pi (balanced back to + sector)
- 4->5: TD=0.579, Fid=0.657, BellDr=0.434, dph=+pi
- 5->6: TD=0.542, Fid=0.699, BellDr=0.519, dph=-pi

Pattern: the hard switches always involve windows 2 and 5, where the two
X tensor X parity sectors are nearly balanced and S-coherence drops to ~0.

### Layer 4: Coherent vs noisy split

Phase difference between noisy and unitary trajectories at every peak: **0.0000**.

The phase skeleton is perfectly preserved by noise. Purity loss grows with time
(0.037 at t=0.24 up to 0.245 at t=3.06), but the phase structure is entirely
Hamiltonian. Noise damps amplitude, it does not rotate phase.

### Dimensionality: PCA on the feature chart

| PC | Singular Value | Variance % | Cumulative % |
|---|---|---|---|
| 1 | 7.893 | 69.2 | 69.2 |
| 2 | 4.452 | 22.0 | 91.2 |
| 3 | 2.466 | 6.8 | 98.0 |
| 4 | 1.270 | 1.8 | 99.8 |
| 5-9 | <0.4 | <0.2 | 100.0 |

**3 dimensions explain 98% of the variance.** Not 15 (generic 2-qubit), not 7
(X tensor X reduced), but 3.

PC1 (69%): Phi+ vs {Phi-, Psi+}, plus S-coherence and concurrence.
This is the **sector balance** axis - which parity sector dominates.

PC2 (22%): Psi (coherence) vs von Neumann entropy, plus purity.
This is the **mixedness** axis - how pure vs mixed the state is.

PC3 (7%): Likely the **temporal decay** axis (to be confirmed).

### What this means

The CΨ window sequence lives on a 3-dimensional manifold:
1. Which parity sector dominates (+ or balanced)
2. How mixed the state is (pure early, mixed late)
3. How far the decay has progressed

The "grammar" (transition rules) has two modes: smooth drift within a sector,
and hard switches between sectors that coincide with S-coherence dropping to zero.

The phase skeleton is entirely Hamiltonian. Noise only affects amplitude/purity.

### What remains for Phase B

- Discretize windows into coarse symbols using PC1 (sector) and PC2 (mixedness)
- Compute block entropies and entropy rate on the symbolized sequence
- Estimate effective memory length: does window N predict window N+1?
- Test: is the sector-switch pattern periodic or does it drift?

### Window XOR: What is NOT shared between windows?

Overlaying the two sharpest windows (0 and 1) and removing what they share:

**88% of Window 0 is also in Window 1.** Only 12-17% is unique to either.

The difference concentrates in specific matrix elements:

| What | Phase change W0 to W1 | Interpretation |
|---|---|---|
| Diagonal (populations) | **0.000** | Identical. The "who is in what state" is stable. |
| rho_03, rho_30 (Phi+ core) | **0.000** | Identical. The main correlation is the skeleton. |
| rho_02, rho_20, rho_13, rho_31 | **+0.787 pi** | Almost 180 deg rotation |
| rho_01, rho_10, rho_23, rho_32 | **+0.633 pi** | More than half rotation |

**The skeleton stays. The cross-connections rotate.**

The Phi+ correlation (|00> <-> |11>) is the stable backbone shared by all windows.
What changes from window to window is the phase of the cross-couplings - the elements
that connect |00> to |10>, |01> to |11>, etc. And they rotate by nearly the same
amount (~0.7 pi), suggesting a single rotational degree of freedom.

In Pauli space, the difference is 87.3% in just two directions (YZ: 58.5%, ZY: 28.8%).
This is a rotation in one plane, not a generic change.

**Glide-mode differences alternate in direction:**
- 0->1: direction D
- 3->4: direction -D (cos = -0.886, almost exactly opposite)
- 6->7: direction +D again (cos = +0.831)

This is a pendulum in Pauli space. Each glide pushes the state in the opposite
direction of the previous glide. The sector switches (hard transitions) are the
turning points.

**Summary:** The windows share a stable skeleton (Phi+ correlation + populations).
What changes is a single rotational degree of freedom in the cross-coupling phases.
The system is a damped pendulum oscillating in one plane of Pauli space, with
periodic sector switches at the turning points.

Script: simulations/window_xor.py

### Shift Analysis: SWAP, YZ/ZY Symmetry, Time-Shift

Three complementary tests on the window structure. Script: simulations/shift_analysis.py

**SWAP Test: A and B are NOT symmetric.**

Every window shows moderate to strong A/B asymmetry (TD 0.07-0.44). But the
asymmetry DAMPS over time: Window 2 has TD=0.44 (strong), Window 6 has TD=0.07
(nearly symmetric). The system becomes more A/B-symmetric despite asymmetric
coupling (J_SA=1.0 vs J_SB=2.0).

The SWAP difference has exactly four non-zero Pauli components: IX=-0.057,
XI=+0.057, YZ=+0.012, ZY=-0.012. Perfectly antisymmetric: what A has, B lacks,
and vice versa.

**YZ/ZY Symmetry: Glide mode is a SYMMETRIC rotor.**

Decomposing each window's YZ and ZY components into symmetric (c+) and
antisymmetric (c-) combinations:

| Pair | Mode | c+/c- ratio | Type |
|---|---|---|---|
| 0->1 | Glide | 5.7 | Symmetric |
| 1->2 | Switch | 1.7 | Mixed |
| 2->3 | Switch | 0.8 | Mixed |
| 3->4 | Glide | 5.5 | Symmetric |
| 4->5 | Switch | 0.6 | Mixed |

In glide mode, c+ dominates c- by factor 5-6. Both qubits rotate equally -
the rotation respects A/B symmetry despite asymmetric coupling. In switch mode,
the symmetry breaks (c+/c- near 1).

**Time-Shift: Period-2 structure in cosine similarity.**

Comparing windows at various lags k:

| Lag | Example pairs | Avg cosine similarity |
|---|---|---|
| k=1 | adjacent | +0.42 |
| k=2 | skip-1 (1 vs 3, 4 vs 6) | +0.43 (cos up to +0.913) |
| k=4 | skip-3 (0 vs 4, 3 vs 7) | +0.44 (cos up to +0.896) |

Windows within the same mode (both glide, or both switch) have very high
similarity even across multiple steps. The pattern repeats with period ~3
(glide-glide-switch). Amplitude decays but structure returns.

### Operation Search: Which transform maps one window to the next?

Tested all 16 two-qubit Pauli operations (IxI through ZxZ) plus continuous
rotations on each adjacent window pair. Script: simulations/window_operations.py

**Key finding: no single-qubit rotation helps. Only correlated operations work.**

Rz_A alone: 0% improvement at any angle, for any pair.
Rz_B alone: 0% improvement.
Ry_A, Ry_B alone: 0% improvement.
Rz_AB (both qubits together, 180 deg): 32% improvement.

The rotation lives in the correlation space, not in either qubit individually.

**Two different grammar rules for two transition modes:**

| Pair | Mode | Best operation | Improvement | What it means |
|---|---|---|---|---|
| 0->1 | Glide | YxY / ZxZ | 32% | Correlated rotation (both qubits together) |
| 1->2 | Switch | IxZ / XxY | 20% | Single-qubit-B operation (breaks sector balance) |
| 2->3 | Switch | XxY / IxZ | 47% | Single-qubit-B operation |
| 3->4 | Glide | IxI (nothing) | 0% | Windows already maximally similar |
| 4->5 | Switch | IxZ / XxY | 31% | Single-qubit-B operation |
| 5->6 | Switch | XxY / IxZ | 4% | Single-qubit-B operation |
| 6->7 | Glide | ZxZ / YxY | 36% | Correlated rotation |
| 7->8 | Glide | IxI (nothing) | 0% | Windows already maximally similar |

**No universal operation exists.** The best single operation across all pairs is
IxI (identity) - doing nothing is better than any fixed transform applied to every
pair. The grammar is pair-specific.

**Interpretation:**

Glide mode = correlated rotation (ZxZ/YxY). Both qubits must rotate together.
This is an entanglement-preserving operation - it acts on the correlation structure,
not on individual qubits. This confirms the pendulum lives in the correlation space.

Switch mode = single-qubit operation on B (IxZ/XxY). Only B is affected. This
asymmetry makes sense: B is the more strongly coupled observer (J_SB=2.0 > J_SA=1.0).
The sector switch is driven by the strongly coupled side.

Quiet glide (3->4, 7->8) = no operation needed. These are the passages where the
damping has slowed the pendulum almost to a stop. The windows are nearly identical
because the system is approaching its steady state in that sector.

**The grammar in one sentence:** Between windows, the system either rotates both
qubits together (correlated glide in the entanglement space) or flips the strongly
coupled qubit alone (sector switch). There is no universal operation - the rule
depends on which mode the system is in.



An external reviewer ran the structural metrics on the star topology data.
These are computational results, not interpretations.

### 1. Exact hidden symmetry: X tensor X

The AB reduced state commutes with X⊗X to numerical precision (max commutator
norm 3.94e-16) across the entire noisy trajectory. In the Bell basis, this makes
rho_AB exactly block-diagonal into two uncoupled parity sectors:

- **+ sector**: {Phi+, Psi+}
- **- sector**: {Phi-, Psi-}

The windows do NOT roam through generic two-qubit state space (15 real parameters).
They live in a symmetry-reduced 7-real-parameter family.

Consequences:
- The "alphabet" is not four free Bell letters but two parity sectors plus
  internal mixing within each sector
- Population symmetries (rho_00 = rho_11, rho_01 = rho_10) are exact, not approximate
- Bell-sector cross-couplings vanish
- rho_00,11 and rho_01,10 stay real

### 2. Sequence structure

Dominant Bell state per window:
Phi+ -> Phi+ -> Phi- -> Phi+ -> Phi+ -> Phi- -> Phi+ -> Phi+

Bell entropy ranges from 1.069 to 1.987 bits. The messy windows (2 and 5) have
nearly balanced parity sectors. Early and late windows lean toward the + sector.

Pattern: mostly + sector, occasional near-balanced transfer, then return.
Not "four states taking turns" but "two sectors with occasional balance shifts."

### 3. S-coherence gating confirmed quantitatively

Over the full time series:
- Pearson correlation CΨ_AB vs S l1-coherence: **0.686**
- Coarse binned mutual information: **0.668 bits**
- Strongest correlation at lag -0.015 time units (essentially synchronous)

The "S coherence gates readability" intuition survives contact with arithmetic.

### 4. Noise damps but does not rewrite the phase skeleton

Comparing gamma=0.05 to gamma=0 at matched peak times:
- Trace distance grows from 0.021 to 0.186
- AB purity loss grows from 0.030 to 0.264
- But the phase of rho_AB[0,3] stays aligned with the unitary reference

Clean split:
- **Coherent part**: keeps timing and phase skeleton (Hamiltonian-driven)
- **Dissipative part**: reduces purity and flattens Bell weights (noise-driven)

### 5. Suggested next step (from review)

Fit a 2-sector hidden-state model in the Bell basis. Test whether it predicts
the next window better than raw Bell-fidelity tracking. This would quantify the
effective memory length and compressibility of the sequence.

## Sources

### Process tensor / multi-time dynamics
- Pollock et al., *Non-Markovian quantum processes* (2018), arXiv:1512.00589
- Milz, Modi, *Quantum stochastic processes and quantum non-Markovian phenomena* (2021), arXiv:2106.11722
- Milz et al., *Operational definition of quantum Markov processes* (2018), arXiv:1801.09811
- Binder et al., *Memory complexity of quantum processes* (2022), arXiv:2203.01492
- Schack, Caves, *Information and entropy in the Baker's map* (2006), arXiv:quant-ph/0611202

### Liouvillian symmetry and mode structure (from external review, March 2026)
- Buca, Prosen, *Symmetry reductions of the Lindblad equation* (2012), arXiv:1203.0943
- Kawabata et al., *Symmetry of Open Quantum Systems* (2023), arXiv:2212.00605
- van Caspel, Gritsev, *Symmetry-protected coherent relaxation* (2018), arXiv:1802.05293
- Tezak et al., *Low dimensional manifolds for open quantum systems* (2017), arXiv:1704.05369
- Haga, *Oscillating-mode gap* (2024), arXiv:2405.07132
- Trushechkin, *Unified GKLS master equation* (2021), arXiv:2103.12042

### Pointer states and decoherence-free subspaces (from evening review)
- Zurek, *Decoherence, einselection, and the quantum origins of the classical* (2003), arXiv:quant-ph/0105127
- Lidar, *Review of Decoherence Free Subspaces, Noiseless Subsystems, and Dynamical Decoupling* (2012), arXiv:1208.5791

## External review: connections to known theory (GPT, March 2026)

An external review identified parallels between our findings and established
open-quantum-systems theory. These are noted connections, not replacements
for our terminology. We keep our descriptions. The established names are
documented here so that readers familiar with those fields can find the
overlap.

### Where our observations connect to known work

| What we call it | What established theory calls similar structures |
|---|---|
| Stable skeleton + rotating phase | Liouvillian eigenmode decomposition (damped modes) |
| X tensor X symmetry | Weak symmetry of the Lindblad generator (Buca, Prosen 2012) |
| 3 effective dimensions (PCA) | Low-dimensional effective manifolds (Tezak et al. 2017) |
| Pendulum in YZ/ZY plane | Dominant complex eigenvalue pair of the Liouvillian |
| Populations stable, phases rotate | Population-coherence separation in master equations |

These are parallels, not identities. Our observations came from a specific
setup (star topology, Heisenberg, dephasing). The established theory is
more general. Whether our specific findings are instances of the general
theory or something adjacent requires the spectral verification below.

### Evening update review (March 9, second round)

Five specific questions answered:

**1. Correlation-only operations.** The reviewer clarified: this is not "local
operations are impossible" - rather, the inter-window difference vector has
almost no weight in local Pauli directions (IX, IY, IZ, XI, YI, ZI) and sits
primarily in the correlated block (XX, XY, ..., ZZ). This is stronger than
just having a global symmetry. Concrete test: project each window difference
onto local vs correlated Pauli basis and compare fractions for glide vs switch.

**2. Damping A/B asymmetry.** Equal gamma does NOT automatically force A/B
symmetry when the Hamiltonian breaks SWAP symmetry (J_SA != J_SB). The
observed convergence is better read as **decay of antisymmetric modes** - a
dynamic effect of the specific mode spectrum, not a general theorem.
Test: track antisymmetric observables (ZI - IZ), (XI - IX) and fit their decay.

**3. Symmetric rotor (c+ dominates in glide).** Consistent with the glide mode
projecting onto a **symmetric, underdamped eigenoperator** of the Liouvillian.
Switch mode is where additional odd-sector components become visible.
Repo shorthand: glide = dominant even-sector rotor, switch = mixing event.

**4. Period-3 recurrence.** The reviewer derived a candidate frequency formula
from the star-topology Hamiltonian eigenstructure:

  f_dom = (J_SA + J_SB + sqrt(J_SA^2 - J_SA*J_SB + J_SB^2)) / pi

For J_SA=1.0, J_SB=2.0: f_dom = (3 + sqrt(3)) / pi = 1.506

This is a project-internal derivation from the model, not a standard formula.
To be verified against FFT of the actual feature time series.

**Frequency test results (March 10, verified):**

FFT on CΨ_AB time series, unitary (gamma=0), t=0 to 20:

| J_SA | J_SB | f_measured | f_GPT formula | f_ours (J_total/2) | GPT error | Our error | Winner |
|---|---|---|---|---|---|---|---|
| 0.5 | 1.0 | 0.7496 | 0.7531 | 0.7500 | 0.5% | 0.1% | Ours |
| 1.0 | 1.0 | 0.9495 | 0.9549 | 1.0000 | 0.6% | 5.3% | GPT |
| 1.0 | 2.0 | 1.4993 | 1.5063 | 1.5000 | 0.5% | 0.1% | Ours |
| 1.0 | 3.0 | 2.0990 | 2.1154 | 2.0000 | 0.8% | 4.7% | GPT |
| 2.0 | 4.0 | 2.9985 | 3.0125 | 3.0000 | 0.5% | 0.1% | Ours |

Both formulas are good but win in different regimes:
- **Asymmetric (J_SA != J_SB):** Our formula f = J_total/2 wins (0.1% error)
- **Symmetric (J_SA = J_SB):** GPT formula wins (0.6% vs 5.3%)

For our standard setup (J_SA=1.0, J_SB=2.0), dominant frequency = 1.4993,
our prediction 1.5000 hits within 0.05%.

The period-3 window recurrence is consistent with this frequency: at f=1.5,
one full cycle takes T=0.667 time units. The glide-glide-switch pattern spans
roughly 3 windows across ~2T, matching the observed period.

Script: simulations/frequency_test.py

**5. Perspective as selection.** Three formal families near this idea:
- Einselection / pointer states (Zurek 2003): environment selects preferred bases
- Decoherence-free subspaces (Lidar 2012): noise algebra determines robust sectors
- Parameterized Lindbladians: different noise parameters = different slices of the
  same reduced dynamics

The reviewer suggests studying {rho_AB(t; gamma_A)} as a surface in state space:
does the window family stay topologically similar? Are there regime transitions?

### Note on terminology

We continue to use our own descriptions (skeleton, rotation, pendulum,
glide/switch modes) in the repo. These emerged from our observations and
carry our specific meaning. The established names are referenced for
cross-disciplinary readers but do not replace our framing.

## README-ready summary

> CΨ windows are analyzed as a structured reduced-state process with persistent
> symmetries, coherent drift, and finite predictive memory. The aim is structural
> cartography: map invariants, transformations, and compressibility first;
> interpret later.

---

*A microscope does not show new physics. It shows known physics that was
previously invisible. Before we can use the microscope, we must learn to
describe what we see.*

---

## See also

- [What's Inside the Windows](WHATS_INSIDE_THE_WINDOWS.md) - The observations that led to this cartography: Bell-state structure, phase transport, bridge test
- [Fixed Point Shadow](FIXED_POINT_SHADOW.md) - IBM hardware shows the same skeleton+rotation pattern (detuning as driver)
- [When Psi Matters](WHEN_PSI_MATTERS.md) - AND-gate justification, noise-type sensitivity, agent benchmark
- [The CΨ Lens](../docs/THE_CPSI_LENS.md) - Canonical definition of CΨ
- [Localizable Entanglement Benchmark](LOCALIZABLE_ENTANGLEMENT_BENCHMARK.md) - Three-layer separation (CoA/LE/CΨ)
- [Observer-Dependent Visibility](OBSERVER_DEPENDENT_VISIBILITY.md) - Different noise profiles, different windows
