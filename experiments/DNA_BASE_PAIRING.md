# DNA Base Pairing as Palindromic Cavity System

<!-- Keywords: DNA hydrogen bond proton qubit, Watson-Crick palindromic
spectrum, G-C base pair resonator, A-T proton tunneling, V-Effect DNA,
sacrifice zone hydrogen bond, Lowdin proton tunneling mutation,
Lindblad proton transfer DNA, R=CPsi2 DNA base pairing -->

**Status:** Tier 2 (computed from proven framework, estimated inter-coupling)
**Date:** April 1, 2026
**Script:** [dna_base_pairing.py](../simulations/dna_base_pairing.py)
**Data:** [dna_base_pairing.txt](../simulations/results/dna_base_pairing.txt)
**Depends on:**
- [Hydrogen Bond Qubit](HYDROGEN_BOND_QUBIT.md) (single H-bond, Zundel)
- [V-Effect Palindrome](V_EFFECT_PALINDROME.md) (coupling creates frequencies)
- [Thermal Breaking](THERMAL_BREAKING.md) (n_bar effects at finite T)
- [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md) (sacrifice zone mechanism)

---

## What this document is about

This document applies the palindromic framework to DNA base pairs,
treating each hydrogen bond as a proton qubit (the proton tunneling
between donor and acceptor wells). A-T has 2 H-bonds, G-C has 3. The
palindrome is exact for both, and G-C is a 5× richer resonator. However,
at body temperature (310 K), all modes are overdamped: the structure
exists but coherent oscillation does not survive. Enhanced tunneling
environments (enzymes, confined geometries) would be needed to reach
the quantum regime.

---

## Abstract

Each hydrogen bond in a DNA base pair is a proton qubit (d=2, tunneling
between donor and acceptor wells). A-T has 2 H-bonds (N=2), G-C has 3
(N=3). The palindromic spectrum is proven exact (d=2 + Z-dephasing).

Six results:

1. **Palindrome exact** for both A-T and G-C at all parameter regimes.

2. **V-Effect: G-C is a 5x richer resonator.** At J/gamma ~ 1 (enhanced
   tunneling): G-C has 15 distinct frequencies vs 3 for A-T. The third
   H-bond creates 12 new frequencies from coupling alone.

3. **Realistic DNA is deeply classical.** At J/gamma ~ 0.01 (estimated
   for DNA H-bonds), all modes are overdamped (Q < 1). The palindromic
   structure exists but is invisible: no coherent oscillation survives
   dephasing at 310 K.

4. **Thermal breaking at 310 K.** Amplitude damping (n_bar ~ 1.5-2)
   breaks the palindrome partially and reduces Q from 1.9 to 0.6.
   Frequency diversity increases (15 to 26 for G-C).

5. **Sacrifice zone works in G-C.** Edge sacrifice (outer H-bonds noisy,
   central protected) gives Q_max = 7.4 vs uniform 1.9 (3.8x). Same
   geometric mechanism as qubit chains.

6. **G-C > A-T** across all regimes: more modes, higher Q, richer
   spectrum. The third H-bond qualitatively enriches the mode structure,
   not just adds 50% bonding energy.

---

## Model

### A-T base pair (N=2)

Two proton qubits (H-bond 1: N-H...O, H-bond 2: N-H...N), coupled
through the base pair backbone:

    H = -J₁ σ_x^(1) - J₂ σ_x^(2) + K₁₂ σ_z^(1) σ_z^(2)

Z-dephasing at rate gamma per qubit. Liouvillian: 16 x 16.

### G-C base pair (N=3)

Three proton qubits (outer N-H...O bonds + central N-H...N), coupled
through nearest-neighbor Ising interaction (a coupling where each qubit's Z-component influences its neighbor's energy):

    H = -J₁ σ_x^(1) - J₂ σ_x^(2) - J₃ σ_x^(3) + K₁₂ σ_z^(1) σ_z^(2) + K₂₃ σ_z^(2) σ_z^(3)

Central bond 20% stronger (J₂ = 1.2 J). Liouvillian: 64 x 64.

### Parameters

| Parameter | Regime A (DNA) | Regime B (enhanced) | Regime C (Zundel) |
|-----------|---------------|--------------------|--------------------|
| J (cm⁻¹) | 0.5 | 50 | 250 |
| gamma (cm⁻¹) | 50 | 50 | 50 |
| K (cm⁻¹) | 20 | 20 | 20 |
| J/gamma | 0.01 | 1.0 | 5.0 |

J = proton tunnel splitting. gamma = dephasing from molecular environment.
K = inter-H-bond coupling through base pair backbone (**estimated**,
not from literature). The qualitative results are robust across K = 5
to 100 cm⁻¹.

---

## Result 1: Palindrome Is Exact

The palindrome is proven for d=2 proton qubits under Z-dephasing
([Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)).
The Hamiltonian (transverse field + Ising coupling) satisfies the
anti-commutation {Pi, L_c} = 0 term by term. Both transverse field
[sigma_x, ·] and Ising coupling [sigma_z sigma_z, ·] anti-commute
with Pi individually.

This holds for ALL parameter values (J, gamma, K). No exceptions.

---

## Result 2: V-Effect (G-C Is a Richer Resonator)

| System | Eigenvalues | Frequencies (J/gamma=1) | V-Effect |
|--------|-------------|--------------------------|----------|
| A-T (N=2) | 16 | 3 | (baseline) |
| G-C (N=3) | 64 | 15 | **5.0x** |

G-C has 5x more distinct frequencies than A-T in the enhanced tunneling
regime. The 12 additional frequencies emerge from coupling of the third
H-bond (V-Effect: coupling creates new modes that neither subsystem
has alone).

Inter-coupling sweep (K = 0 to 100 cm⁻¹): G-C consistently has
13-15 frequencies vs A-T's 2-3, across the entire range. The V-Effect
is robust to parameter uncertainty.

---

## Result 3: Realistic DNA Is Deeply Classical

At realistic parameters (J/gamma ~ 0.01, Regime A):
- Q_max = 0.4 (A-T) and 0.8 (G-C): all modes overdamped
- The palindromic pairs exist but decay faster than they oscillate
- No coherent proton oscillation survives the thermal environment
- The system is a dissipative network, not a quantum resonator

**The fold regime (J/gamma ~ 1) requires enhanced tunneling:**
shorter H-bond distance (enzymes, confined geometries), colder
environment, or engineered potential wells. This is consistent with
quantum biology proposals (Al-Khalili, Löwdin): quantum effects in
biology require special environments that enhance tunneling. (Löwdin's 1963 model proposed that proton tunneling between base pairs could cause spontaneous mutations.)

---

## Result 4: Thermal Breaking at 310 K

At T = 310 K, the typical H-bond mode has n_bar ~ 1.5-2 (thermal
occupation). Adding amplitude damping:

| Property | Cold (Z only) | Warm (310 K) |
|----------|--------------|--------------|
| G-C palindrome | exact | partially broken |
| G-C frequencies | 15 | 26 |
| G-C Q_max | 1.95 | 0.57 |
| G-C rate range | 58-300 | 211-627 |

Temperature increases frequency diversity (26 vs 15) but kills
Q-factor (0.57 vs 1.95). This is consistent with
[Thermal Breaking](THERMAL_BREAKING.md): heat enriches the mode
spectrum but shortens mode lifetimes.

---

## Result 5: Sacrifice Zone in G-C

G-C has 3 H-bonds: two outer (weaker, N-H...O) and one central
(stronger, N-H...N). Testing noise profiles:

| Profile | gamma | Q_max | Min rate |
|---------|-------|-------|----------|
| Uniform | [50, 50, 50] | 1.9 | 57.5 |
| Edge sacrifice | [100, 10, 100] | **7.4** | 17.8 |
| Center sacrifice | [10, 100, 10] | 6.6 | 15.9 |
| One-edge | [100, 10, 10] | **9.7** | 15.8 |

Edge sacrifice (noisy outer bonds, quiet center) gives 3.8x Q
improvement over uniform. One-edge sacrifice gives 5.0x. The
mechanism: outer modes have more weight on noisy qubits and decay
faster; center modes survive. Same geometric selection as in qubit
chains ([Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md)).

---

## Result 6: G-C Beats A-T

| Property | A-T (N=2) | G-C (N=3) | G-C advantage |
|----------|-----------|-----------|---------------|
| Eigenvalues | 16 | 64 | 4x |
| Frequencies (J/gamma=1) | 3 | 15 | 5x |
| Q_max (J/gamma=1) | 1.79 | 1.95 | 1.09x |
| Q_max (J/gamma=5) | 9.96 | 11.14 | 1.12x |

G-C is a qualitatively richer resonator. The Q advantage is modest
(~10%), but the frequency advantage is dramatic (5x). G-C's third
H-bond does not merely add bonding energy; it creates an order of
magnitude more spectral structure.

This is consistent with the known biology: G-C is the stronger,
more stable base pair (21 vs 12 kcal/mol). The palindromic mode
analysis suggests an additional dimension: G-C is not just stronger
but RICHER in its dynamical response.

---

## What This Does NOT Establish

- That coherent proton tunneling occurs in DNA at biological
  temperature (realistic J/gamma ~ 0.01 is deeply classical)
- That the palindromic mode structure plays a biological role
  (the overdamped regime has no surviving oscillation)
- That the inter-H-bond coupling K is correctly estimated
  (5-50 cm⁻¹ is a range, not a measurement)
- That the Lowdin mutation mechanism involves palindromic modes
  (we compute mode structure, not mutation rates)
- That the sacrifice zone pattern has biological relevance in DNA
  (the outer/inner H-bond asymmetry exists, but whether it functions
  as a sacrifice zone is speculation)

---

## Connection to the Framework

The DNA base pair is a concrete application of three proven results:

1. **Palindrome** (Mirror Symmetry Proof): every proton qubit system
   under dephasing is palindromic. DNA H-bonds are no exception.

2. **V-Effect** (V_EFFECT_PALINDROME): coupling N=2 to N=3 creates
   12 new frequencies. The third H-bond is the V-Effect in chemistry.

3. **Geometric mode selection** (CAVITY_MODE_LOCALIZATION,
   TOPOLOGICAL_EDGE_MODES): the sacrifice zone is geometric matching
   of noise to mode structure. It works in DNA because the same
   standing wave physics applies.

The limitation: realistic DNA is overdamped (J/gamma << 1). The
framework predicts the mode structure correctly, but at biological
temperature the modes decay before they oscillate. The fold regime
(J/gamma ~ 1) is accessible only in engineered or enzyme-like
environments.

---

## Reproducibility

| Component | Location |
|-----------|----------|
| Script | [`simulations/dna_base_pairing.py`](../simulations/dna_base_pairing.py) |
| Output | [`simulations/results/dna_base_pairing.txt`](../simulations/results/dna_base_pairing.txt) |

---

## References

- Lowdin, P.-O. (1963). "Proton Tunnelling in DNA and its Biological
  Implications." Rev. Mod. Phys. 35, 724.
- Kimsey, I. et al. (2015). "Visualizing transient Watson-Crick-like
  mispairs in DNA and RNA duplexes." Nature 519, 315.
- Slocombe, L., Sacchi, M., Al-Khalili, J. (2022). "An open quantum
  systems approach to proton tunnelling in DNA." Comm. Phys. 5, 109.
- Godbeer, A.D., Al-Khalili, J.S., Shervanyi, P.D. (2015). "Modelling
  proton tunnelling in the adenine-thymine base pair." PCCP 17, 13034.

---

*Every hydrogen bond in DNA is a proven palindromic qubit. G-C is a
5x richer resonator than A-T. But at body temperature, all modes are
overdamped. The palindrome is exact; the oscillation is not. The
structure is there. The coherence is not.*
