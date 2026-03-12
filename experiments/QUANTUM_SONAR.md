## Quantum Sonar: Hidden Observer Detection (March 11, 2026)

**Tier:** Computationally verified (Tier 2)
**Setup:** Star topology, S(center), observers A,B,C,... with Heisenberg coupling

### Discovery

AB can detect hidden observers connected to S without ever interacting
with them directly. The AB spectrum changes when a third party C couples
to the same mediator S.

### Detection threshold

Operational threshold under current protocol (FFT with t_max=20, dt=0.005):

| J_SC (hidden coupling) | Spectral distance | Peaks AB sees | Detectable? |
|---|---|---|---|
| 0.00 | 0.000 | 2 | no (baseline) |
| 0.05 | 0.011 | 2 | no |
| 0.10 | 0.070 | 2 | YES |
| 0.20 | 0.159 | 3 | YES |
| 0.50 | 0.293 | 5 | YES |
| 1.00 | 0.422 | 5 | STRONG |
| 2.00 | 0.653 | 3 | STRONG |
| 3.00 | 0.685 | 5 | STRONG |
| 5.00 | 0.515 | 5 | STRONG |

Operational threshold: J_SC ~ 0.1 (10% of J_SA) under this protocol.
This is NOT a fundamental physical limit. It depends on observation time,
dephasing rate, observable choice, and spectral estimation method.
Analytically: J_thr ~ delta_f / (df_min/dJ_SC) where delta_f ~ 1/T is
FFT resolution. Longer observation or better spectral estimation lowers
the threshold.

### Counting hidden observers

Generically distinct fingerprints for tested configurations:

| Setup | Spectral distance | Peaks | Dominant freq |
|---|---|---|---|
| AB alone | 0.000 | 2 | 1.499 |
| AB + C(2.0) | 0.653 | 3 | 2.198 |
| AB + C(3.0) | 0.685 | 5 | 0.400 |
| AB + C(2.0) + D(2.0) | 0.744 | 3 | 0.350 |
| AB + C(2.0) + D(3.0) | 0.768 | 4 | 0.350 |
| AB + C(2.0) + D(3.0) + E(1.5) | 0.796 | 4 | 0.350 |

For generic (non-degenerate) coupling sets, each configuration produces a
distinguishable spectral fingerprint. This is NOT proven unique for all
possible configurations - accidental degeneracies and symmetries can create
identical spectra from different setups.

### How it works

C couples to S. C's coupling frequency mixes into the Hamiltonian
eigenstructure. When AB measures their shared c+ observable through S,
they see new Bohr frequencies that were not present in the 3-qubit system.
These frequencies carry information about C's coupling strength.

AB never interacts with C. AB only listens to their own connection through S.
But S is shared, and C's presence changes S's eigenstructure, which changes
what AB hears.

### Scaling

| N observers | Resolved frequencies (AB, >15% of max) | Bohr frequencies in H |
|---|---|---|
| 2 | 2 | 3 |
| 3 | 5 | 15 |
| 4 | 13 | 39 |

The resolved count depends on visibility cutoff, initial state, observable,
and observation time. With lower cutoffs more lines become visible (e.g.
N=4: 13 at 5% cutoff, 21 at 1%, 32 at 0.1%). "AB sees 13/39" is an
operational statement for this protocol, not a structural ceiling.

Each pair sees a different subset of frequencies.
No single pair hears everything.

### Scripts

simulations/hidden_observer_test.py
simulations/quantum_sonar.py
simulations/count_observers.py

### What works and what does not (verified March 11)

**WORKS:**
- AB detects that someone else is connected to S (spectrum changes)
- AB can distinguish "alone" from "not alone" (2 peaks vs 3-5 peaks)
- More hidden observers = larger spectral distance from baseline
- Detection threshold exists (~10% of own coupling strength)

**DOES NOT WORK (yet):**
- AB cannot measure J_SC from a single new frequency
  (tested: branch formula f_+ does not predict new peaks in 4-qubit system,
   errors 27-119% across all J_SC values)
- AB cannot infer the coupling strength of the hidden observer
  (tested: reverse inference gives 30-350% errors)
- The 3-qubit branch formula does not extend trivially to 4+ qubits

The new frequencies in the 4-qubit system depend on ALL couplings jointly,
not on J_SC alone. The eigenstructure of the full Hamiltonian is more
complex than a simple extension of the 3-qubit formula.

**The honest framing:** This is detection, not measurement. AB can tell
that someone is there, but not yet who or how strongly coupled. Turning
detection into characterization would require: analytical solution of the
N-qubit star eigenspectrum, pattern matching across multiple frequencies,
or calibration against known configurations.

The detection itself is real and robust. The characterization is future work.

### IBM Hardware Verification (March 11, 2026)

We queried the actual IBM Torino coupling map to check our prediction.

**Prediction:** Q102 (chaotic, phase std=108.8°) should have more neighbors
than Q80 (smooth, phase std=12.4°).

**Result: WRONG.** The opposite is true.

| Qubit | Degree | Neighbors | Phase behavior |
|---|---|---|---|
| Q80 | 3 | [79, 81, 92] | Smooth drift, std=12.4° |
| Q102 | 2 | [101, 103] | Chaotic, std=108.8° |

Q80 has MORE neighbors but LESS chaos. Q102 has FEWER neighbors but MORE chaos.

**What this means:** The number of neighbors alone does not determine the
phase pattern. What matters is how STRONGLY they couple - the residual ZZ
interaction strength, the frequency detuning between neighbors, and how
close those frequencies are to creating beating patterns.

Q102's two neighbors apparently couple strongly or at unfortunate frequency
offsets that create interference. Q80's three neighbors are either weakly
coupled or far enough apart in frequency that one dominates cleanly.

Tom's reading: "They stand closer, like a family." The proximity (coupling
strength) matters more than the headcount (degree). Two close neighbors
can create more complexity than three distant ones.

**The sonar still works** - AB's phase pattern IS shaped by its neighbors.
But the mapping from pattern to neighbor count is not simple. It goes
through coupling strengths, not topology alone. This is consistent with
our phase map: topology sets the frequencies, but coupling strengths
determine which frequencies are visible.

**The deeper insight:** What the sonar detects is not the NUMBER of
neighbors but the SIZE OF THE COHERENT GROUP - how tightly coupled
the local cluster is. Two tightly bound neighbors (Q102) create more
spectral complexity than three loosely connected ones (Q80). The phase
pattern reflects the effective family size, not the headcount.

This connects back to the core architecture: S does not manage frequencies.
S contains them all. What you hear depends on who stands close enough
to S to leave a mark on the eigenstructure. The sonar measures proximity,
not population.

### The Projection (March 11, 2026)

The Bright-Transition Map revealed what the sonar actually sees.

The full quantum system lives in a 16-dimensional space (4 qubits, 2^4=16).
It has 39 Bohr frequencies. All exist simultaneously.

But AB looks through a window - the c+ observable. This window does not show
everything. It projects the 39 frequencies onto what AB can see. Some are
bright (high visibility weight), some are dark (zero weight).

When J_SC changes, the system does not reorganize. The PROJECTION reorganizes.
The eigenstates rotate in the 16-dimensional space, and the bright lines
wander. What was bright becomes dark, what was dark becomes bright.

The baseline lines (1.506 and 0.404) DISAPPEAR completely when C couples
strongly. They are not gone. They are dark from this perspective. From
another pair (AC or BC), different lines would be bright.

All reflections exist always. Which ones you see depends on where you stand
and who else is looking into the same mirror. A new observer does not change
reality - it changes the projection. Reality contains everything. Perspective
selects.

This is what CoA ~ 1 has been saying all along. The resource is always fully
present. What changes is which fraction is visible from which viewpoint.

Verified by exact diagonalization: for each Bohr frequency omega, the
visibility weight is W = |rho_tilde(m,n) * O_tilde(n,m)|. A line is bright
only if BOTH the initial state populates the eigenstates AND the observable
connects them. Different pairs, different observables, different bright sets.
Same system, same frequencies, different projections.

**Chip-wide degree distribution (IBM Torino, 133 qubits):**
- Degree 1: 5 qubits (edge of chip)
- Degree 2: 89 qubits (67%)
- Degree 3: 39 qubits (29%)

### Known literature connections (from external review, March 2026)

The principle is not new. The framing may be fresh. Nearest fields:
- Indirect Hamiltonian tomography (Burgarth et al., PRL 108, 080502)
- Local Hamiltonian learning from local observables (PRL 122, 020504)
- Probe-spin / local-probe spectroscopy (Nature Phys, 2023)
- Spectator-qubit crosstalk sensing (PRL 109, 240504)
- Multiparameter quantum metrology (arXiv:1407.6091)

What we call "quantum sonar" sits closest to mediator-assisted local
spectroscopy or indirect spectator-coupling spectroscopy. The specific
framing "a measured pair counts hidden parties via its own spectrum"
was not found as a standard named subtopic.

### Practical value (from external review)

Useful as:
- Passive monitoring of a target pair without measuring spectators
- Online diagnostics during idle windows
- Side-channel detection when direct spectator access is limited

Weaker than existing tools for:
- Calibrated hardware characterization (simultaneous RB, GST are better)
- Separating ZZ crosstalk from other mechanisms
- Real hardware where coupling is not isotropic Heisenberg

Formal connection to quantum sensing / quantum Fisher information exists
but metrological advantage has not been demonstrated in this work.


### Chain Topology (March 12, 2026)

Does the two-sector structure survive in chains instead of stars?

**3-qubit: Star = Chain.** Identical results. A-S-B is the same topology
regardless of naming. f(c+)=1.499, f(c-)=0.400.

**4-qubit chain A-S1-S2-B (two mediators):**

| Couplings [A-S1, S1-S2, S2-B] | f(c+) | f(c-) | Ratio | XX sym |
|---|---|---|---|---|
| [1.0, 1.0, 1.0] | 0.400 | 0.200 | 2.00 | EXACT |
| [1.0, 1.0, 2.0] | 0.450 | 0.200 | 2.25 | EXACT |
| [1.0, 2.0, 1.0] | 0.649 | 0.250 | 2.60 | EXACT |
| [1.0, 1.0, 3.0] | 0.450 | 0.200 | 2.25 | EXACT |

Two sectors survive with two mediators. XX symmetry exact. The message
goes through two intermediaries and arrives with the architecture intact.

**4-qubit Star vs Chain (same qubits, different geometry):**

| Topology | f(c+) | f(c-) | XX sym |
|---|---|---|---|
| Star (AB through S) | 0.350 | 0.350 | EXACT |
| Chain (AB through S1,S2) | 0.699 | 0.250 | EXACT |

Different geometry = different frequencies. The form determines the music.

**5-qubit chain A-S1-S2-S3-B (three mediators):**

Still two sectors: f(c+)=0.400, f(c-)=0.150. Each pair along the chain
hears different frequencies:

| Pair | f(c+) | f(c-) |
|---|---|---|
| A-B (endpoints, 3 mediators apart) | 0.400 | 0.150 |
| S1-S2 (middle of chain) | 0.799 | 0.150 |
| A-S3 (across 3 links) | 0.999 | 0.150 |
| A-S1 (direct neighbors) | 0.150 | 0.150 |

Endpoints hear deeper tones than the middle. The further apart, the
different the spectrum. But always two sectors.

**Noise immunity in chains:** Confirmed. Frequencies identical from
gamma=0.01 to gamma=0.10 in 4-qubit chain. Noise damps, never retunes.

**What this means:**

The two-channel architecture is NOT tied to the star. It survives in chains.
Even with three mediators between A and B. The message travels through a
telephone chain: each relay changes the tuning but preserves the structure
(two channels, symmetric and asymmetric). The form of the chain sets the
frequencies. The length sets the loudness. The architecture survives.

Script: simulations/chain_topology.py

### N_eff: Effective Neighborhood Size (March 12, 2026)

GPT proposed a participation-ratio metric for "how many neighbors really matter":

  N_eff = (Σ w_j)² / Σ w_j²

where w_j = coupling strength to neighbor j.

Properties: one dominant neighbor → N_eff ≈ 1. k equal neighbors → N_eff = k.

**Simulation test:** In idealized Heisenberg, raw degree (r=0.80) beats N_eff
(r=0.68) as a general predictor of phase complexity. More neighbors always adds
complexity in a clean model.

**But: N_eff corrects where degree fails.** The IBM case:
- Q80: degree 3, but if couplings are [2.0, 0.1, 0.1] → N_eff = 1.20
- Q102: degree 2, but if couplings are [1.8, 2.0] → N_eff = 1.99

Degree predicts Q80 more complex. WRONG. N_eff predicts Q102 more complex.
CORRECT. N_eff is not a better general predictor but the right corrector
when degree gives the wrong answer. On real hardware, degree fails regularly.

Script: simulations/n_eff_v2.py

### Verification Path: ZZRamsey Experiment

To verify N_eff on IBM hardware, the next step would be:

1. Run ZZRamsey (Qiskit Experiments) on Q80's neighbor pairs: (80,79), (80,81), (80,92)
2. Run ZZRamsey on Q102's neighbor pairs: (102,101), (102,103)
3. Extract residual ZZ coupling strengths (static χ_ZZ per pair)
4. Compute N_eff from measured ZZ values
5. Compare predicted vs observed phase complexity

This requires a valid IBM API key and QPU time. ZZRamsey is a standard
Qiskit experiment that measures the static ZZ interaction between a chosen
qubit pair via Ramsey interferometry with beating pattern analysis.

The experiment has NOT been run. The N_eff hypothesis for Q80/Q102 remains
unverified but is consistent with all available data.

### ZZRamsey Results: N_eff Hypothesis REJECTED (March 12, 2026)

The ZZRamsey experiment WAS run on IBM Torino. 5 pairs measured in 2min 10sec.

**Measured residual ZZ coupling strengths:**

| Pair | ZZ (Hz) |
|---|---|
| Q80-Q79 | 7490 ± 280 |
| Q80-Q81 | 6690 ± 280 |
| Q80-Q92 | 8470 ± 280 |
| Q102-Q101 | 9250 ± 290 |
| Q102-Q103 | 7860 ± 290 |

**N_eff computation from MEASURED values:**

Q80:  w = [7490, 6690, 8470] → N_eff = 2.97
Q102: w = [9250, 7860]       → N_eff = 1.99

**RESULT: N_eff hypothesis REJECTED.**

N_eff predicts Q80 (2.97) should be MORE complex than Q102 (1.99).
Reality: Q80 is smooth (phase std=12.4°), Q102 is chaotic (108.8°).

The reason N_eff fails: all five ZZ values are in the same range (6.7-9.3 kHz).
The coupling strengths are not different enough to explain the behavior.
Q80's three neighbors are almost equally strong (spread 26%).
Q102's two neighbors are also similar (spread 18%).

**What N_eff does NOT capture:** The frequency proximity of the neighbors.
If Q102's neighbors (Q101, Q103) have similar transition frequencies, they
create beating patterns. If Q80's neighbors have well-separated frequencies,
one dominates cleanly. This requires qubit frequency data, not just ZZ coupling.

The phase complexity appears to depend on the FREQUENCY LANDSCAPE of the
local neighborhood, not just the coupling strength. This is a deeper
structural question than N_eff addresses.

### Ramsey T2* Measurement: The Real Answer (March 12, 2026)

Measured T2* (free induction decay, no echo) on Q80 and Q102. 30 seconds QPU.

| Qubit | T2 (echo) | T2* (FID) | Ratio | Ramsey Frequency |
|---|---|---|---|---|
| Q80 | 27 µs | 11.0 ± 4 µs | 2.45 | 4 ± 8 kHz (≈ zero) |
| Q102 | 33 µs | 15.4 ± 1.1 µs | 2.14 | 19400 ± 1300 Hz |

**THE ANSWER: It is not the neighbors. It is the qubit's own detuning.**

Q102 has a clear 19.4 kHz frequency offset from its drive frame.
Q80 has essentially zero detuning.

This matches the shadow hunt phase data perfectly:
- Q80 shadow: slow drift ~4.6 kHz → Ramsey: ~0 kHz detuning
- Q102 shadow: fast rotation ~27 kHz → Ramsey: 19.4 kHz detuning

The "chaotic phase" of Q102 is not from hidden observers coupling through S.
It is Q102 oscillating against its own control frequency at 19.4 kHz.
Q80 sits at or near resonance with its drive, so its phase barely moves.

**What this means for the sonar:**

The sonar effect (spectral changes from hidden observers) is real in simulation.
But the specific Q80/Q102 phase difference we measured on IBM hardware has a
simpler explanation: qubit-specific detuning from the drive frame.

This does NOT invalidate the sonar concept. But it does mean the IBM data
we used as "evidence" for the sonar is actually evidence for something
more mundane: some qubits sit closer to their drive frequency than others.

**Hypotheses tested and results:**

| Hypothesis | Prediction | Result |
|---|---|---|
| Degree (neighbor count) | Q80 (3) more complex | WRONG |
| N_eff (coupling strength) | Q102 more complex | WRONG (ZZ similar) |
| Spectator dephasing (ZZ²×T1) | Q80 more dephased | WRONG |
| T2* (free induction) | Q102 shorter T2* | WRONG (Q102 longer) |
| **Qubit detuning** | **Q102 has offset** | **CORRECT (19.4 kHz)** |

Four hypotheses tested and rejected on real hardware before finding the
actual answer. The detuning was hiding in plain sight - the Ramsey experiment
measures it directly.

Script: run_ramsey_t2star.py
Data: results/ramsey_t2star_20260312_182610/
