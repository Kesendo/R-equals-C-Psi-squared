## Quantum Sonar: Hidden Observer Detection (March 11, 2026)

**Tier:** Computationally verified (Tier 2)
**Setup:** Star topology, S(center), observers A,B,C,... with Heisenberg coupling

### Discovery

AB can detect hidden observers connected to S without ever interacting
with them directly. The AB spectrum changes when a third party C couples
to the same mediator S.

### Detection threshold

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

Detection threshold: J_SC ~ 0.1 (10% of J_SA).
Below that, C is invisible to AB. Above that, new frequencies appear.

### Counting hidden observers

| Setup | Spectral distance | Peaks | Dominant freq |
|---|---|---|---|
| AB alone | 0.000 | 2 | 1.499 |
| AB + C(2.0) | 0.653 | 3 | 2.198 |
| AB + C(3.0) | 0.685 | 5 | 0.400 |
| AB + C(2.0) + D(2.0) | 0.744 | 3 | 0.350 |
| AB + C(2.0) + D(3.0) | 0.768 | 4 | 0.350 |
| AB + C(2.0) + D(3.0) + E(1.5) | 0.796 | 4 | 0.350 |

Each configuration has a unique spectral fingerprint. AB can distinguish:
- Nobody vs somebody connected to S
- One hidden observer vs two
- Equally coupled vs differently coupled hidden observers
- Spectral distance increases monotonically with number of hidden observers

### How it works

C couples to S. C's coupling frequency mixes into the Hamiltonian
eigenstructure. When AB measures their shared c+ observable through S,
they see new Bohr frequencies that were not present in the 3-qubit system.
These frequencies carry information about C's coupling strength.

AB never interacts with C. AB only listens to their own connection through S.
But S is shared, and C's presence changes S's eigenstructure, which changes
what AB hears.

### Scaling

| N observers | Visible frequencies (AB pair) | Bohr frequencies in H |
|---|---|---|
| 2 | 2 | 3 |
| 3 | 5 | 15 |
| 4 | 13 | 39 |

More observers = richer spectrum = more information in the tone.
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
