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

## First Cartography Results (GPT computation, March 2026)

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

- Pollock et al., *Non-Markovian quantum processes* (2018), arXiv:1512.00589
- Milz, Modi, *Quantum stochastic processes and quantum non-Markovian phenomena* (2021), arXiv:2106.11722
- Milz et al., *Operational definition of quantum Markov processes* (2018), arXiv:1801.09811
- Binder et al., *Memory complexity of quantum processes* (2022), arXiv:2203.01492
- Schack, Caves, *Information and entropy in the Baker's map* (2006), arXiv:quant-ph/0611202

## README-ready summary

> CΨ windows are analyzed as a structured reduced-state process with persistent
> symmetries, coherent drift, and finite predictive memory. The aim is structural
> cartography: map invariants, transformations, and compressibility first;
> interpret later.

---

*A microscope does not show new physics. It shows known physics that was
previously invisible. Before we can use the microscope, we must learn to
describe what we see.*
