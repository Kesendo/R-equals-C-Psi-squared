# Quantum Sonar — Current State

**Tier:** Computationally verified in simulation (Tier 2)
**Last updated:** March 12, 2026

---

## The Effect (simulation, verified)

AB can detect hidden observers connected to S without direct interaction.
When C couples to shared mediator S, AB's spectrum changes: new Bohr
frequencies appear, existing ones shift. This is exact eigenstructure
reorganization, not an artifact.

Operational detection threshold: J_SC ~ 0.1 (10% of J_SA) under current
FFT protocol (t_max=20, dt=0.005). Not a fundamental limit — depends on
observation time, dephasing, and spectral estimation method.

Fingerprints are generically distinct for generic coupling sets.
Not proven unique for all configurations (symmetries, accidental degeneracies).

## The Projection (key insight)

All Bohr frequencies exist simultaneously. AB sees only those where BOTH
the initial state populates the eigenstates AND the observable connects them:

  W = |rho_tilde(m,n) * O_tilde(n,m)|

A new observer does not change reality — it rotates the eigenstates in the
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
Weaker than simultaneous RB or GST for calibrated hardware characterization.

## Scripts

- simulations/hidden_observer_test.py — detection proof
- simulations/quantum_sonar.py — threshold sweep
- simulations/count_observers.py — scaling with N
- simulations/bright_transition_map.py — exact diagonalization + visibility
- simulations/chain_topology.py — chain vs star
- simulations/n_eff_v2.py — N_eff metric test

## IBM experiment scripts (in AIEvolution/experiments/ibm_quantum_tomography/)

- run_zzramsey.py — ZZ coupling measurement
- run_ramsey_t2star.py — T2* and Ramsey frequency
- query_coupling_map.py — topology query
- Data: results/zz_ramsey_*/ and results/ramsey_t2star_*/
