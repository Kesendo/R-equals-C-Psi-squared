# Energy Partition

**Status:** Quantum spectral diagnostics; biological interpretation unestablished
**Authors:** Thomas Wicht, Claude (Anthropic)
last refreshed 2026-09-05 (the change history lives in git)

The experiments compare sums of eigenvalue frequencies and decay rates in
Heisenberg qubit chains. These sums have units of inverse time; they are not
physical energies or mode populations. A thermal-bath parameter changes the
quantum spectrum. It has no calibrated mapping to neural drive or metabolism.

The evidence sweep checked [the F-registry](../docs/ANALYTICAL_FORMULAS.md)
(F1 pairing, F8's range/centre ratio, conditional neural F36/F37),
[docs/proofs](../docs/proofs/PROOF_ABSORPTION_THEOREM.md),
[the neural proofs](../docs/neural/proofs/PROOF_VEFFECT_MECHANISM.md), and
[experiments](../experiments/THERMAL_BREAKING.md): quantum rate identities,
thermal protocols and neural mechanism constraints. Hardware-flight searches
and [fw.Confirmations](../simulations/framework/confirmations.py) supplied no
neural hardware confirmation. [GLOSSARY](../docs/GLOSSARY.md),
[OpenArcs](../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs), and
[CAUGHT_ERRORS](../docs/CAUGHT_ERRORS.md) supplied scope distinctions and
instrument failures. The [current neural account](../docs/neural/README.md)
owns the biological support null and constructed counterexamples.

## 1. Quantum pairing and the filtered spectrum

[energy_partition.py](../simulations/energy_partition.py) and Part A of
[thermal_emergence.py](../simulations/thermal_emergence.py) use open
Heisenberg chains with uniform J=1 and Z-dephasing γ=0.1, N=2…5.
They discard eigenvalues with |λ|≤10⁻⁸ before searching for the reflected
partner −λ−2Nγ with a relative tolerance of 10⁻⁶.

The [full-spectrum quantum theorem](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)
pairs the zero modes with roots at −2Nγ. Removing the zero roots leaves
their partners unmatched in the filtered list. It does not break the
generator's palindromic symmetry. The scripts' nearest-partner search also
does not independently test multiplicity.

| N | Matched fraction after zero removal | Unmatched decay rate | Matched mean decay rate | Ratio |
|---|---:|---:|---:|---:|
| 2 | 76.9% | 0.4000 | 0.2000 | 2.0 |
| 3 | 93.3% | 0.6000 | 0.3000 | 2.0 |
| 4 | 98.0% | 0.8000 | 0.4000 | 2.0 |
| 5 | 99.4% | 1.0000 | 0.5000 | 2.0 |

The unmatched roots in these runs are real. All resolved oscillatory
content remains in the matched list. This is a statement about the selected
quantum family and filtering convention, not a necessary condition for
oscillation in arbitrary open systems.

[F8](../docs/ANALYTICAL_FORMULAS.md#f8-2-universal-decay-law-tier-1-corollary-of-absorption-theorem)
identifies the ratio as the full decay interval's width 2Nγ divided by
its centre Nγ. For γ>0 this is 2. At γ=0 that ratio is undefined.
Neither this ratio nor spectral pairing labels a mode as signal or noise;
an observed response also depends on its preparation and readout overlaps.
There is no neural 2× law implied by this calculation.

## 2. Quantum thermal-bath census

Part C of [thermal_emergence.py](../simulations/thermal_emergence.py) sets
H=0 and uses independent local raising/lowering jumps. This uncoupled
generator has no oscillatory eigenvalues in the tested thermal-occupation
sweep. That result concerns this bath and H=0.

Part D compares Z-dephasing and a thermal bath at N=3, J=1. Its diagnostics
are F=Σ|Im λ| and D=Σ|Re λ| on roots with |λ|>10⁻⁸; an oscillatory root
has |Im λ|>10⁻⁸. F and D are spectral sums, not thermodynamic energies.

The pure Z-dephasing rows show that even when the maximum frequency stays
fixed, the number of oscillatory roots and their frequency sum can change:

| γ | Oscillatory roots | Maximum abs(Im λ) | F |
|---:|---:|---:|---:|
| 0.01 | 40 | 6.0 | 160.0 |
| 1.00 | 40 | 6.0 | 148.3 |
| 10.00 | 32 | 6.0 | 119.5 |

The combined bath holds Z-dephasing γ=0.1 and thermal rate Γ=0.1 fixed,
with local jump amplitudes √(Γ(1+n_bar)) for σ₋ and √(Γ n_bar) for σ₊:

| n_bar | Oscillatory roots | F | D | F/D |
|---:|---:|---:|---:|---:|
| 0.00 | 40 | 159.9 | 28.8 | 5.55 |
| 1.00 | 42 | 159.8 | 48.0 | 3.33 |
| 5.00 | 42 | 158.4 | 124.8 | 1.27 |
| 10.00 | 42 | 154.2 | 220.8 | 0.70 |

The bath is present even at n_bar=0 through the σ₋ jump. Raising its
occupation produces two additional resolved oscillatory roots in these
rows, while F decreases and D increases. The table brackets F/D=1
between n_bar=5 and 10; it defines no universal optimal thermal window.
Eigenvalue counts alone establish neither population of those modes nor
conversion of heat into a measured oscillatory signal. The larger
[Thermal Breaking experiment](../experiments/THERMAL_BREAKING.md) is a
separate protocol with its own channels, frequency census and Q diagnostics.

## 3. What a neural comparison would require

The [neural coupling/drive report](../docs/neural/V_EFFECT_NEURAL.md)
measures frequency bins on specified synthetic Jacobians. Its external
input P changes the sigmoid operating point and effective row gains.
P is not n_bar, temperature, ATP production or metabolic power. E/I balance
does not identify a thermal window.

[F36/F37](../docs/neural/proofs/PROOF_PALINDROME_NEURAL.md) require one
involutive Q and one scalar s satisfying both the diagonal and coupling
conditions. They constrain complex partner sums, not class-mean decay
ratios or the protection of oscillation. Exact constructed palindromes
can oscillate and can be unstable; the full committed C. elegans chemical
model fails the support condition. No biological neural network in the
repository is known to pass F36.

A comparison needs a specified neural generator and operating point,
declared observables and units, and controls that can reject the proposed
relation. A spectral census must survive frequency-resolution checks;
an equilibrium-stability claim requires converged equilibria. A sustained
response needs preparation, readout and time-domain evidence. The
[mechanism constraints](../docs/neural/proofs/PROOF_VEFFECT_MECHANISM.md)
provide these gates; no thermal, metabolic or life mechanism is established.

## 4. Open questions

- **The two extra quantum roots:** Track their eigenvalues and invariant
  subspaces through the combined-bath sweep; determine their contributions
  to a specified response. Mode counts alone do not identify their structure.
- **Other quantum channels and Hamiltonians:** Test the full spectrum under
  the actual mirror theorem's assumptions before comparing filtered lists.
  The range/centre identity supplies no theorem for arbitrary channels.
- **Neural comparison:** Define a spectral or response diagnostic on a
  specified converged neural model and test it against matched controls.
  A physical calibration is required before relating its input to a bath.
- **Connection to the quarter:** A spectral F/D crossover alone cannot
  specify state-dependent CΨ. Choose the same generator, state, time and
  normalization and compute both independently. The algebraic fold of
  R=C(Ψ+R)² in [the quarter proof roadmap](../docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md)
  does not supply this missing identification.

## Reproduce

From the repository root in PowerShell:

```powershell
$env:PYTHONIOENCODING = 'utf-8'
python simulations/energy_partition.py
python simulations/thermal_emergence.py
python simulations/neural/neural_translation_gate.py
python -m pytest simulations/neural/tests/ -q
```

The two quantum scripts print their historical diagnostic labels to stdout.
Interpret those labels using the definitions above. For current neural
producers and controls, use [the neural operator's manual](../simulations/neural/README.md).
