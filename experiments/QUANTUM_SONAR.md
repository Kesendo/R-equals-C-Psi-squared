# Quantum Sonar: Passive Detection of Hidden Observers Through Spectral Shifts

<!-- Keywords: quantum sonar hidden observer detection, spectral shift passive
quantum monitoring, Bohr frequency eigenstructure reorganization, star topology
mediator coupling detection, correlated bath sector selection, chain topology
quantum sensing, IBM hardware qubit detuning resolution, quantum side-channel
detection Hamiltonian learning, R=CPsi2 quantum sonar -->

**Status:** Computationally verified (Tier 2); IBM hardware: not yet demonstrated
**Date:** March 12, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Star Topology Observers](STAR_TOPOLOGY_OBSERVERS.md), [Structural Cartography](STRUCTURAL_CARTOGRAPHY.md)

---

## Abstract

A qubit pair AB can detect the presence of a hidden observer C coupled to
their shared mediator S, without direct interaction. In a star topology
(S at center, A and B as outer qubits), coupling a new qubit C to S
reorganizes the system's eigenstructure: new Bohr frequencies (the
discrete oscillation frequencies determined by energy-level differences)
appear and existing ones shift. The detection threshold is J_SC ~ 0.1 (10% of J_SA)
under the current FFT protocol (t_max=20, dt=0.005). Fingerprints are
generically distinct for different coupling configurations. The two-sector
frequency structure survives in chains (tested up to 5 qubits) and is
completely immune to correlated noise, though bath geometry selects which
amplitude sector dominates. IBM hardware investigation (Q80, Q102 on Torino)
traced apparent phase anomalies to qubit-specific detuning (19.4 kHz for
Q102), not the sonar effect. The sonar mechanism is verified in simulation
but remains to be demonstrated on quantum hardware.

---

## What this document is about

Can a pair of qubits detect a hidden third qubit that they cannot see
directly? Yes, in simulation: when a new qubit couples to their shared
mediator, the frequencies visible to the original pair shift and new ones
appear. This is quantum sonar: passive detection through spectral
fingerprinting. The detection threshold is about 10% of the main
coupling strength. An IBM hardware investigation traced apparent
anomalies to qubit-specific frequency offsets, not the sonar effect.
The mechanism is verified computationally but not yet demonstrated on
real quantum hardware.

---

## The Effect (simulation, verified)

In a star topology (S = central mediator, A and B = outer qubits), AB can
detect hidden observers connected to S without direct interaction. When a
new qubit C couples to S, the spectrum visible to AB changes: new Bohr
frequencies appear, existing ones shift. This is exact eigenstructure
reorganization, not an artifact.

Operational detection threshold: J_SC ~ 0.1 (10% of J_SA) under current
FFT protocol (t_max=20, dt=0.005). Not a fundamental limit; depends on
observation time, dephasing, and spectral estimation method.

Fingerprints are generically distinct for generic coupling sets.
Not proven unique for all configurations (symmetries, accidental degeneracies).

## The Projection (key insight)

All Bohr frequencies exist simultaneously. AB sees only those where BOTH
the initial state populates the eigenstates AND the observable connects them:

  W = |rho_tilde(m,n) * O_tilde(n,m)|

A new observer does not change reality; it rotates the eigenstates in the
full space, which changes which lines are bright from AB's perspective.
Different pairs, different observables, different projections. Same system.

Verified by exact diagonalization (bright_transition_map.py).

## Chain topology (March 12)

Two-sector structure survives in chains, not just stars. Tested 4-qubit
(two mediators) and 5-qubit (three mediators). XX exact, noise immune,
all configurations. Each position hears different frequencies.
Architecture is topology-independent.

## IBM hardware investigation (March 12)

We tried to find the sonar effect on IBM Torino by comparing Q80 (smooth
phase, std=12.4°) vs Q102 (chaotic phase, std=108.8°).

**Five hypotheses tested, four rejected:**

| Hypothesis | Prediction | Result |
|---|---|---|
| Degree (neighbor count) | Q80 (3) more complex | WRONG |
| N_eff (coupling strength) | Q102 more complex | WRONG (ZZ all similar) |
| Spectator dephasing (ZZ²×T1) | Q80 more dephased | WRONG |
| T2* (free induction) | Q102 shorter T2* | WRONG (Q102 longer) |
| **Qubit detuning** | **Q102 has offset** | **CORRECT (19.4 kHz)** |

The Q80/Q102 difference is qubit-specific detuning from the drive frame,
not neighbor effects. Q102 oscillates at 19.4 kHz against its control
frequency. Q80 sits near resonance (~0 kHz).

**Measured ZZ couplings (ZZRamsey, 2min 10sec QPU):**
Q80: [7490, 6690, 8470] Hz. Q102: [9250, 7860] Hz. All similar range.

**The sonar effect is real in simulation. The IBM data does not demonstrate it.**

## Literature connections (from external review)

Nearest fields: indirect Hamiltonian tomography (Burgarth, PRL 108),
local Hamiltonian learning (PRL 122), probe-spin spectroscopy (Nature Phys 2023),
spectator-qubit crosstalk sensing (PRL 109), multiparameter quantum metrology.

Useful as passive monitoring or side-channel detection.
Weaker than simultaneous RB (Randomized Benchmarking) or GST (Gate Set Tomography) for calibrated hardware characterization.

## Correlated Bath Sweep (March 12, 2026)

GPT's experiment: test if correlated noise (shared bath on A+B) breaks
the "noise only damps" rule. Two knobs: eta (local→correlated) and
phi (ZZ bath→XX bath).

**Frequencies: COMPLETELY IMMUNE.** f(c+)=1.499 and f(c-)=0.400 in every
single configuration. Local, correlated ZZ, correlated XX, all mixtures.

**Amplitude ratio flips:**

| Bath geometry | A+/A- | Dominant sector |
|---|---|---|
| Local (baseline) | 1.22 | c+ dominates |
| Correlated ZZ (eta=1, phi=0) | 1.09 | c+ still |
| Correlated mixed (eta=1, phi=0.5) | 0.73 | c- takes over |
| Correlated XX (eta=1, phi=1) | 0.46 | c- dominates |

**XX symmetry:** Breaks at mixed correlated bath (phi=0.1-0.9) but
survives at both extremes (pure ZZ AND pure XX). Mixing breaks it.

**Result: GPT's prediction #1 confirmed.** Bath geometry is a sector
selector. Frequencies never move. But the bath chooses which channel
is louder. The phase map needs a fifth role:

1. Topology sets frequencies
2. Symmetry cleans sectors
3. Noise damps amplitude
4. Initial state selects visibility
5. **Bath geometry selects which sector dominates**

Script: simulations/correlated_bath_sweep.py

## Scripts

- simulations/hidden_observer_test.py: detection proof
- simulations/quantum_sonar.py: threshold sweep
- simulations/count_observers.py: scaling with N
- simulations/bright_transition_map.py: exact diagonalization + visibility
- simulations/chain_topology.py: chain vs star
- simulations/n_eff_v2.py: N_eff metric test

## IBM experiment scripts (in AIEvolution/experiments/ibm_quantum_tomography/)

- run_zzramsey.py: ZZ coupling measurement
- run_ramsey_t2star.py: T2* and Ramsey frequency
- query_coupling_map.py: topology query
- Data: results/zz_ramsey_*/ and results/ramsey_t2star_*/
