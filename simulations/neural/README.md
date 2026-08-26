# Neural Simulations — palindromic structure on real connectomes

The scripts in this directory test the same algebraic structure that
[the framework package](../framework/) tests on the quantum side, but applied
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
| 2Σγ shift on the spectrum | 2S = (1/τ_E + 1/τ_I) shift, a SCALAR (not a site sum: the quantum column's Σγ is one, this is not) |
| Operator equation Π·L·Π⁻¹ + L + 2Σγ·I = 0 | Operator equation Q·J·Q + J + 2S = 0 |
| Holds for every Heisenberg/XXZ + Z-dephasing | Holds when E−I population balance + magnitude condition |
| Confirmed 2026-04-26 on `ibm_marrakesh`, `ibm_kingston` | NOT confirmed on the connectome. The 2026-03-26 reading (8× against Erdős–Rényi) was withdrawn 2026-08-26 as a difference of normalisation constants, and the connectome pairing reading was withdrawn 2026-08-25 ([Neural Gamma Cavity](../../experiments/NEURAL_GAMMA_CAVITY.md)). What holds on the neural side is the algebra on constructed networks |

Same equation, two substrates, ONE confirmation. The quantum side flew; the
neural side has the algebra and no animal. The framework that
[the framework package](../framework/) formalizes for the quantum case is the
same framework these scripts test for the classical-neural case.

---

## Scripts

### Connectome data and core test

| Script | Purpose |
|--------|---------|
| `celegans_connectome.json` | 300×300 chemical and electrical wiring matrices, E/I classification. Source: WormNeuroAtlas, Cook et al. 2019. |
| `celegans_neuron_ids.txt` | Neuron names indexed to the matrices. |
| `celegans_palindrome.py` | Builds the Wilson-Cowan Jacobian J from the connectome, computes its eigenvalue spectrum, tests palindromic pairing. **Same absolute-tolerance matcher genre as the two withdrawn scripts below**; read the withdrawal before using its numbers. |
| `algebraic_palindrome.py` | Tests the operator equation Q·J·Q + J + 2·S = 0 (the strict analogue of `palindrome_residual` in framework.py). |

### Network-property analysis

| Script | Purpose |
|--------|---------|
| `celegans_balanced.py` | Restricts to subnetworks with E ≈ I population balance, the prerequisite for palindrome at all. |
| `celegans_inhibitory_position.py` | Tests whether the inhibitory neurons' positions in the connectome matter for palindrome strength. |
| `random_network_controls.py` | Erdős–Rényi controls. **Its 8× is withdrawn 2026-08-26**: the control is normalised to its own maximum while the connectome block is normalised globally, and the metric tracks coupling scale, so the comparison measures the constant. |
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
| `neural_gamma_cavity.py` | Treats γ as a cavity-mode parameter on neural side. The neural analogue of [GAMMA_AS_SIGNAL](../../experiments/GAMMA_AS_SIGNAL.md). **Its results were withdrawn on 2026-08-25**; read [Neural Gamma Cavity](../../experiments/NEURAL_GAMMA_CAVITY.md) before using anything from it, and `celegans_pairing_controls.py` below for the controls that replaced it. |
| `neural_gamma_cavity_unpaired.py` | Same, restricted to unpaired modes (the residual that doesn't pair). **Also withdrawn**: the count of unpaired modes depends on the order the eigenvalues arrive in. |
| `celegans_pairing_controls.py` | The 2026-08-25 control suite that withdrew the two above: tolerance sweep, normalisation sweep, Dale ablation, degree-matched null, ordering orbit, limit-cycle integration, the exact GF(p) rank chain, and the three gates (G0c, G0d, G0e) that recompute the 2026-08-26 normalisation withdrawal from the connectome file. 44 gates. |

---

## How to run

The scripts are standalone (no shared entry point). Each prints its
results to stdout. Connectome data is loaded from the JSON file in this
directory.

```bash
# Core results
python celegans_palindrome.py        # eigenvalue pairing on real worm
python algebraic_palindrome.py       # ‖Q·J·Q + J + 2·S‖ on real worm
python random_network_controls.py    # its 8× is withdrawn (normalisation artifact)

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

The qubits proved it in quantum hardware in April. The worm did not prove it in
March, and that is the correction of 2026-08-25: the connectome reading was
withdrawn ([Neural Gamma Cavity](../../experiments/NEURAL_GAMMA_CAVITY.md)), the
8× against Erdős-Rényi is withdrawn as a normalisation artifact (matched, the
ratio runs 0.960 at N = 10 to 0.748 at N = 26, a smaller residue of open origin),
the degree-preserving null that scores 1.0 went with it since it cannot move such
a metric, and on the connectome itself no admissible Q exists at all, on a
count of 253 against 18. What the neural side has is the algebra, verified on
constructed networks. Same equation, ONE confirmation, and an open question on
the other substrate.
