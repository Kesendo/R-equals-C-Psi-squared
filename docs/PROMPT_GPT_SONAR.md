Attach the repo as ZIP. This is a new discovery from tonight.

---

## What we found (March 11, 2026)

### Quantum Sonar: Hidden Observer Detection

We asked a simple question: if A and B are connected through a mediator S,
can AB detect that a third party C has also connected to S — without ever
interacting with C directly?

**Answer: yes.**

When C couples to S (even at just 10% of A's coupling strength), new
frequencies appear in AB's c+ spectrum. Each configuration of hidden
observers produces a unique spectral fingerprint. AB can:

- Detect that someone else is connected to S (threshold J_SC ~ 0.1)
- Count how many hidden observers there are (spectral distance increases monotonically)
- Distinguish differently coupled observers (each J_SC produces different frequencies)

The mechanism: C changes the Hamiltonian eigenstructure of the full system.
Even though AB only measures their own pair through S, S's eigenmodes carry
the signature of everyone connected to it. More observers = more Bohr
frequencies = richer spectrum.

### Scaling

| N observers on S | Frequencies AB sees | Total Bohr frequencies |
|---|---|---|
| 2 | 2 | 3 |
| 3 | 5 | 15 |
| 4 | 13 | 39 |

Each pair (AB, AC, BC) sees a DIFFERENT subset of frequencies.
No single pair hears the complete spectrum.

### Also verified: Branch formula

Your formula f_+/- = (1 + J_SB +/- sqrt(1 - J_SB + J_SB^2)) / pi
is confirmed across all sweep points. It beats our simpler f = J_total/2
at large J_SB (0.7% error vs 10.4% at J_SB=5.0). The slow branch f_-
is weakly dispersive as you predicted.

---

## Questions for you

### 1. Does "spectral detection of hidden parties" exist in quantum information?

We are describing a scenario where two parties (AB) detect a third party (C)
purely through changes in the spectral structure of their shared mediator.
No direct interaction, no measurement of C, no classical communication.
Just listening to the overtones of their own quantum connection.

Is this known? Does it have a name? The closest things we can think of:
- Quantum spectroscopy of many-body systems
- Spectator qubit effects in superconducting hardware (ZZ crosstalk)
- Quantum sensing / quantum metrology

But we have not seen it framed as "counting hidden observers from the
spectrum of a pair."

### 2. Is the detection threshold meaningful?

We find J_SC ~ 0.1 (10% of J_SA) as the threshold where C becomes visible
to AB. Is this threshold a property of FFT resolution, or is there a
fundamental limit from the Hamiltonian eigenstructure? Can you derive
the detection threshold analytically from the spectral gap?

### 3. What is the theoretical maximum of distinguishable configurations?

With N=4 observers we get 13 visible frequencies in AB. But there are
39 Bohr frequencies in the Hamiltonian. AB sees 13/39 = 33%. Is there
a formula for how many frequencies a given pair can see as a function
of N? Is this related to the partial trace structure?

### 4. Practical relevance: crosstalk detection

In real quantum hardware, spectator qubits cause ZZ crosstalk that
degrades gate fidelity. Our result suggests you could detect and
characterize spectator coupling purely from the spectral fingerprint
of a target pair — without measuring the spectators directly.

Is this useful? Is it better or worse than existing crosstalk
characterization methods (simultaneous RB, correlated error tomography)?

### 5. Connection to quantum sensing

Quantum sensing uses entangled probes to detect external perturbations.
Our setup is different: the "probe" is the AB pair, the "perturbation"
is C coupling to the shared mediator. But the principle is similar —
read changes in the probe to infer properties of the environment.

Is there a formal connection to quantum sensing / quantum metrology
frameworks? Could the spectral fingerprint approach give metrological
advantage (Heisenberg scaling or similar)?

---

As always: honest assessment. If this is a trivial consequence of
"more qubits = more eigenvalues" and everyone already knows it,
say so. If there is something genuinely new here, say that too.

Key files:
- experiments/QUANTUM_SONAR.md
- simulations/hidden_observer_test.py
- simulations/quantum_sonar.py
- simulations/count_observers.py
