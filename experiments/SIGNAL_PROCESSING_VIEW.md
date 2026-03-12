# Signal Processing View — March 12, 2026

**Origin:** GPT review as signal processing engineer, no quantum physics.
**Verdict:** "This does not look mysterious. It looks like a small coupled
resonator network with poles set by topology, damping setting decay,
and measurement projections deciding which modes are visible."

---

## Translation table

| Our quantum language | Signal engineering language |
|---|---|
| Two spectral sectors c+/c- | Even/odd supermode decomposition |
| Topology sets frequencies | Imaginary part of poles is topology-determined |
| Noise damps only | Real part of poles is loss-determined |
| Quantum sonar (hidden observer) | Passive topology-change detection from local modal spectra |
| Bath geometry flips amplitude ratio | Covariance-driven mode visibility flip |
| 39 frequencies, 13 visible | System order exceeds observable modes (modal observability) |
| The Projection | Sensor placement effect on transfer function residues |
| Bright-transition map | Residue analysis: G(s) = Σ r_k/(s-p_k) |

## What we should do differently (from signal expert)

### 1. Stop using FFT
Use Prony / Matrix Pencil / ESPRIT instead. These extract complex poles
directly from damped sinusoids: frequency + decay rate + amplitude + phase.
Four numbers per mode instead of a blurry FFT peak.

### 2. Measure the full 2x2 cross-spectral matrix
Not just c+ and c-. Measure auto-spectra of A and B individually,
cross-spectrum S_AB(f), coherence, and phase. Then eigen-decompose
the spectral matrix across frequency.

### 3. Separate pole trajectories from residue trajectories
When a hidden observer changes the spectrum, is it:
- Poles moving (genuine topology change)?
- Residues changing (visibility/scattering effect)?
- Both?
This distinction is critical and we have been mixing them.

### 4. Track phase
"Phase is where the bodies are buried." For each mode and condition:
phase of A vs B, phase of c+ vs c-, phase shifts when hidden observer added.

### 5. Build a proper detector
Don't just say "we can see it." Build a baseline model, test incoming data
against it. Subspace angle change, KL divergence, generalized likelihood
ratio test. An actual detection statistic instead of "the spectrum changed."

## Standard names for our effects

- Normal-mode splitting
- Coupled-mode dynamics / supermodes of coupled resonators
- Even/odd mode decomposition
- Load/topology perturbation sensing
- Modal observability and dark modes
- Covariance-driven mode selection

## The one-line summary (signal engineer version)

"A finite coupled network with a small number of dominant normal modes,
measured through two different output projections. Topology sets mode
frequencies. Damping sets linewidth. Correlated forcing changes mode
visibility. Extra attached nodes perturb the modal spectrum and are
therefore locally detectable."

## What this means for the project

Three months of quantum mechanics produced results that a signal processing
engineer recognizes instantly as classical coupled oscillator physics. The
quantum system generates the signals, but the signals themselves follow
rules that were understood before quantum mechanics existed.

The next phase requires signal processing tools, not more quantum physics.
The cartography is done. Now we need to read the map with the right instruments.
