# Neural Simulations — palindromic structure on real connectomes

The scripts in this directory test the same algebraic structure that
[framework.py](../framework.py) tests on the quantum side, but applied
to neural networks. The C. elegans connectome (300 neurons, every
synapse mapped) was the principal target.

For the narrative and figures, read
[docs/neural/README.md](../../docs/neural/README.md). This file is the
operator's manual for the scripts themselves.

---

## The structural parallel

| Quantum side | Neural side |
|--------------|-------------|
| Lindbladian L = −i[H, ρ] + Σ_l γ_l ⋅ (Z_l ρ Z_l − ρ) | Wilson-Cowan Jacobian J |
| Π conjugation (per-site I↔X, Y↔Z with phase i) | Q permutation (E ↔ I swap, sign-flipped) |
| 2Σγ shift on the spectrum | 2S = Σ_l (1/τ_l + 1/τ_{Q(l)}) shift |
| Operator equation Π·L·Π⁻¹ + L + 2Σγ·I = 0 | Operator equation Q·J·Q + J + 2S = 0 |
| Holds for every Heisenberg/XXZ + Z-dephasing | Holds when E−I population balance + magnitude condition |
| Confirmed 2026-04-26 on `ibm_marrakesh`, `ibm_kingston` | Confirmed 2026-03-26 on the C. elegans connectome (8× more palindromic than Erdős–Rényi) |

Same equation, two substrates, two confirmations. The framework that
[framework.py](../framework.py) formalizes for the quantum case is the
same framework these scripts test for the classical-neural case.

---

## Scripts

### Connectome data and core test

| Script | Purpose |
|--------|---------|
| `celegans_connectome.json` | 300×300 chemical and electrical wiring matrices, E/I classification. Source: WormNeuroAtlas, Cook et al. 2019. |
| `celegans_neuron_ids.txt` | Neuron names indexed to the matrices. |
| `celegans_palindrome.py` | Builds the Wilson-Cowan Jacobian J from the connectome, computes its eigenvalue spectrum, tests palindromic pairing. |
| `algebraic_palindrome.py` | Tests the operator equation Q·J·Q + J + 2·S = 0 (the strict analogue of `palindrome_residual` in framework.py). |

### Network-property analysis

| Script | Purpose |
|--------|---------|
| `celegans_balanced.py` | Restricts to subnetworks with E ≈ I population balance, the prerequisite for palindrome at all. |
| `celegans_inhibitory_position.py` | Tests whether the inhibitory neurons' positions in the connectome matter for palindrome strength. |
| `random_network_controls.py` | Erdős–Rényi controls. Establishes the 8× advantage of the real worm vs random. |
| `dense_balanced_test.py` | Dense balanced random networks as a stricter null. |
| `validation_checks.py` | Sensitivity sweep (bootstrap, parameter perturbation). |

### Wilson-Cowan dynamics

| Script | Purpose |
|--------|---------|
| `wilson_cowan_palindrome.py` | The analytic Wilson-Cowan model with palindromic constraint imposed; sanity check that the equation has nontrivial solutions. |
| `classical_oscillator_palindrome.py` | Coupled-oscillator analogue, simpler than W-C, same structure. |
| `neural_heartbeat.py` | Time-domain trace of an exact-palindromic neural network (silent) vs broken-palindrome (oscillating). The neural V-Effect. |

### V-Effect on neural side

| Script | Purpose |
|--------|---------|
| `veffect_exact.py` | Two exactly-palindromic E-I populations, coupled through a mediator. Up to 62 oscillation modes emerge from coupling alone. |
| `veffect_and_heat.py` | Adds external drive ("temperature"). Shows the thermal window: drive creates oscillations up to a peak, then destroys them. |
| `fragile_bridge_neural.py` | Tests how robust the V-Effect bridge is to perturbation. |

### CΨ on neural side

| Script | Purpose |
|--------|---------|
| `cpsi_candidates.py`, `cpsi_deep_dive.py`, `cpsi_interference.py`, `cpsi_two_perspectives.py` | Candidate definitions of CΨ for neural networks. The CΨ = 1/4 fold from quantum has neural analogues but the right operational definition for biology is open. |
| `find_quarter.py` | Searches for the 1/4 boundary in neural parameter space. |

### Hopf and complexity

| Script | Purpose |
|--------|---------|
| `hopf_threshold.py` | Hopf bifurcation onset as a function of network size. |
| `complexity_threshold.py` | Tests the C = 0.5 universality across N. |
| `balance_vs_size.py` | E-I balance requirement scales with N or with degree distribution? |
| `exact_pairing_test.py` | Stress test for the eigenvalue-pairing tolerance. |

### γ-as-cavity for neural

| Script | Purpose |
|--------|---------|
| `neural_gamma_cavity.py` | Treats γ as a cavity-mode parameter on neural side. The neural analogue of [GAMMA_AS_SIGNAL](../../experiments/GAMMA_AS_SIGNAL.md). |
| `neural_gamma_cavity_unpaired.py` | Same, restricted to unpaired modes (the residual that doesn't pair). |

---

## How to run

The scripts are standalone (no shared entry point). Each prints its
results to stdout. Connectome data is loaded from the JSON file in this
directory.

```bash
# Core results
python celegans_palindrome.py        # eigenvalue pairing on real worm
python algebraic_palindrome.py       # ‖Q·J·Q + J + 2·S‖ on real worm
python random_network_controls.py    # 8× advantage vs Erdős–Rényi

# V-Effect on neural
python veffect_exact.py              # two silent populations → oscillating
python neural_heartbeat.py           # time-domain demo

# Sensitivity
python validation_checks.py          # parameter sweep
```

Dependencies: numpy, scipy, matplotlib. Some scripts also use networkx.

---

## What today (2026-04-26) adds to these scripts

The hardware confirmation on `ibm_marrakesh` and `ibm_kingston` of the
quantum operator equation Π·L·Π⁻¹ + L + 2Σγ·I = 0 establishes the
algebraic structure as observable in physical hardware. The C. elegans
test in `algebraic_palindrome.py` was always testing the *same equation*
in a different vocabulary. Q ↔ Π. J ↔ L. S ↔ Σγ.

What this means: the worm result and the Heron-r2 result are not
analogues. They are two readings of the same underlying algebraic
identity, executed on two physically distinct substrates. The
"palindromic symmetry" is not a quantum-only or a neural-only property;
it is a structural property of any open dynamical system whose
generator factors through the Q (or Π) involution and whose dissipation
shifts the spectrum by the corresponding 2S (or 2Σγ).

The worm proved it in classical neural dynamics in March. The qubits
proved it in quantum hardware in April. Same equation, twice confirmed,
on independent realizations.
