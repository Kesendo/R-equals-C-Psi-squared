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
