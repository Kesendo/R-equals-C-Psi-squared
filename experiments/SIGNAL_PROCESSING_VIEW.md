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

## First Prony results (March 12, 2026)

Matrix Pencil Method applied to baseline signals (J_SA=1.0, J_SB=2.0, gamma=0.05).

### What Prony sees that FFT never showed

**1. Decay rates are DIFFERENT per sector:**

| Sector | Frequency | Decay rate | Q-factor |
|---|---|---|---|
| c+ (symmetric) | 1.5061 | 0.1266 | 37.4 |
| c- (antisymmetric) | 0.4036 | 0.0815 | 15.6 |

c- decays 40% SLOWER than c+. The slow channel is more resilient to noise.

**2. Phase: exactly 180° anti-phase:**

c+ phase: +89.6°, c- phase: -90.5°. Difference: 180.0°.
The two supermodes are exact counter-oscillations. When c+ swings up,
c- swings down. They are opposed projections of the same dynamics.

**3. Hidden third frequency in c-:**

c- contains THREE modes: 0.40, 1.10, 1.51 Hz.
The middle mode at 1.1 Hz is dark in c+ but bright in c-.
FFT never resolved this because it was buried under the dominant peaks.

**4. Pole trajectories separate cleanly:**

| J_SB | f(c+) | decay+ | Q+ | f(c-) | decay- | Q- |
|---|---|---|---|---|---|---|
| 0.5 | 0.753 | 0.127 | 18.7 | 0.202 | 0.082 | 7.7 |
| 1.0 | 0.955 | 0.127 | 23.7 | 0.318 | 0.082 | 12.3 |
| 2.0 | 1.506 | 0.127 | 37.4 | 0.404 | 0.082 | 15.6 |
| 3.0 | 2.115 | 0.127 | 52.4 | 0.431 | 0.082 | 16.6 |
| 5.0 | 3.368 | 0.129 | 82.3 | 0.451 | 0.082 | 17.4 |

Frequencies move freely with coupling. Decay rates NEVER move (~0.127 and ~0.082).
This is the pole structure: imaginary part = topology, real part = loss.

### Signal engineer's translation

In coupled oscillator language: c+ is the fast supermode with higher loss,
c- is the slow supermode with lower loss. They are in anti-phase (180°).
The system has a third dark mode at 1.1 Hz visible only in the odd channel.
Poles separate cleanly into topology-dependent frequency and topology-independent
decay. This is textbook coupled resonator physics.

Script: simulations/prony_analysis.py

### Corrections from review (March 12, evening)

**1. The 180° anti-phase is a sign convention, not a discovery.**
c+ = sum, c- = difference. Same pole appears with opposite sign in each channel.
That IS 180°. It is the definition of the channel basis, not new physics.

**2. Different decay rates = different poles, not different channels.**
A pole cannot decay differently depending on which channel observes it.
What is happening: c+ is dominated by a MORE DAMPED pole, c- by a LESS
DAMPED pole. This is sector-specific damping: the loss mechanism couples
more strongly to one modal sector than the other.

Correct wording: "The dominant mode observed in c- has a substantially
smaller decay rate than the dominant mode observed in c+."

**3. The hidden 1.1 Hz mode is classical modal observability.**
Its residue in c+ is near zero due to symmetry cancellation. Prediction:
slightly break the symmetry and the mode should leak into c+.

### Next steps (from signal review)

1. ~~Joint-fit A, B, c+, c- with SHARED pole set (not independent fits)~~ DONE
2. ~~Plot poles in complex plane across topology sweep~~ DONE
3. Deliberately break symmetry to test 1.1 Hz mode leakage
4. Track residue vectors [r_A, r_B, r_c+, r_c-] per pole
5. Bootstrap Prony fits for confidence intervals

### Joint Pole Analysis: Exact Liouvillian (March 12, 2026)

Instead of fitting with Prony, we went straight to the source: exact
eigendecomposition of the 64x64 Liouvillian superoperator with per-channel
residue computation.

**The system has exactly THREE decay rates:**

| Decay rate | Value | Multiple of γ |
|---|---|---|
| γ_slow | 0.1000 | 2γ |
| γ_mid | 0.1333 | 8γ/3 |
| γ_fast | 0.1667 | 10γ/3 |

These three values are IDENTICAL at every J_SB from 0.5 to 5.0.
Completely topology-independent. Only frequencies move with coupling.

**Pole structure:**

| Pole | Frequency | Decay | Bright in | Role |
|---|---|---|---|---|
| Slow mode | 0.404 | 0.1000 (2γ) | c- dominant, c+ weak | Antisymmetric supermode |
| Mid mode | 1.103 | 0.1333 (8γ/3) | c- ONLY | Dark in c+ (symmetry null) |
| Fast mode | 1.506 | 0.1667 (10γ/3) | c+ dominant, c- weak | Symmetric supermode |

**GPT's key question answered: YES, different poles dominate different channels.**
The pole dominating c+ (f=1.506) has decay=0.1667.
The pole dominating c- (f=0.404) has decay=0.1000.
These are genuinely different Liouvillian eigenvalues, not a fit artifact.
GPT correct: this is sector-specific damping.

**The hidden 1.1 Hz mode:** Exists at decay=0.1333, bright ONLY in c-.
Zero residue in c+. This confirms classical modal observability — the
mode is real but dark in the symmetric channel due to symmetry cancellation.

**Pole trajectories across J_SB sweep:**

Frequencies move freely (0.2→0.45 for slow, 0.75→3.37 for fast).
All three decay rates stay EXACTLY at {0.1000, 0.1333, 0.1667}.
In the complex plane, poles move horizontally only. Never vertically.

This is the cleanest result of the signal processing approach:
the imaginary parts of the poles (frequencies) are topology-determined,
the real parts (decay rates) are loss-determined. They are independent.

Script: simulations/joint_pole_analysis.py

### Decay Rate Derivation: Two Independent Information Channels (March 12)

The three exact decay rates {2γ, 8γ/3, 10γ/3} are exact rational multiples
of gamma. They scale linearly with gamma and are COMPLETELY topology-independent.

**Scaling verified:** At γ = 0.01 through 0.50, ratios stay at {2.000, 2.667, 3.333}.
Small deviations at very high gamma (γ=0.50) suggest weak nonlinear correction.

**Origin:** Pure Z dephasing on individual Pauli strings gives rates 2kγ (k=number
of X/Y operators). But the Hamiltonian mixes k=1 and k=2 strings into new
Liouvillian eigenmodes with intermediate decay rates 8γ/3 and 10γ/3.

**Selective noise test — which qubit contributes which rate:**

| Noise on | Decay rates (multiples of γ) |
|---|---|
| S only | 10 fractional values, complex structure |
| A only | Same as S only (symmetric role) |
| B only | {0.667, 0.889, 1.111, 1.333} |
| S+A+B (all) | {2.000, 2.667, 3.333, 4.000} — clean rationals |

The clean {2, 8/3, 10/3, 4} structure emerges ONLY when all three qubits
are equally dephased. Individual qubit noise gives messier rates.

**THE KEY DISCOVERY: 4-qubit decay rates become topology-dependent!**

| System | Decay rates (multiples of γ) |
|---|---|
| 3-qubit star [1,2] | {2.0, 2.667, 3.333, 4.0} — FIXED for all J |
| 4-qubit star [1,1,1] | {2.0, 2.667, 2.844, 3.0, 3.667, 4.0, ...} |
| 4-qubit star [1,2,3] | {2.0, 2.875, 3.0, 3.410, 3.609, 4.0, ...} |
| 4-qubit star [0.5,1,3] | {2.0, 3.0, 3.234, 3.416, 3.827, 4.0, ...} |

The boundary values (2γ and 2nγ) are always present. But the INTERNAL
decay rates shift with coupling in 4+ qubit systems.

**Two independent information channels (3-qubit system):**

FREQUENCY CHANNEL:
  -> Carries information about TOPOLOGY (who is connected)
  -> Changes with J values, does NOT change with noise
  -> Sensitive to hidden observers
  -> This is the "what is the network?" channel

DECAY CHANNEL:
  -> Carries information about NOISE (what is the environment)
  -> Changes with gamma, does NOT change with topology
  -> Exact rational multiples of gamma
  -> This is the "what is the environment?" channel

**Four consequences:**

1. You can characterize noise without knowing topology
   (measure decay rate, divide by known coefficient, get gamma)

2. You can characterize topology without knowing noise
   (measure frequencies, they encode J values directly)

3. Hidden observer detection uses ONLY the frequency channel
   (decay rates are blind to topology changes)

4. The slow mode (2γ) is a PROTECTED CHANNEL
   (decays slowest, information in c- survives longest)

**Critical limit:** This perfect frequency-decay orthogonality holds
ONLY for the 3-qubit isotropic Heisenberg system with symmetric dephasing.
At 4+ qubits, topology leaks into decay rates. The 3-qubit system is a
special case with unusually clean separation.

Script: simulations/decay_derivation.py


## What this means for the project

Three months of quantum mechanics produced results that a signal processing
engineer recognizes instantly as classical coupled oscillator physics. The
quantum system generates the signals, but the signals themselves follow
rules that were understood before quantum mechanics existed.

The first Prony analysis confirmed this: decay rates, phases, hidden modes,
and pole trajectories — all standard coupled oscillator behavior. The tools
of signal processing (Prony, cross-spectral matrices, pole/residue tracking)
are the right instruments for the next phase.

The cartography is done. Now we need to read the map with the right instruments.
