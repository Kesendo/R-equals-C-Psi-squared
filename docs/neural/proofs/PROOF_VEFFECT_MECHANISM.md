# Neural V-effect: falsifiers and constraints on a mechanism

Last refreshed: 2026-09-05 (change history lives in git).

No V-effect mechanism is established. The [experiment report](../V_EFFECT_NEURAL.md)
owns the current coupling and external-drive frequency censuses. This page
owns their interpretation constraints: exact palindrome does not force real,
stable or silent modes; the odd-seat coupled matrix fails F36 before coupling;
frequency bins and unconverged endpoints cannot supply a mechanism.

The named stores checked are [ANALYTICAL_FORMULAS](../../ANALYTICAL_FORMULAS.md)
(F36/F37 give conditional complex pairing; F137 gives the trace constraint),
[docs/proofs](../../proofs/MIRROR_SYMMETRY_PROOF.md) (quantum conjugation,
with [F134](../../proofs/PROOF_F134_TWO_ROW_REFLECTION_LAW.md) owning a separate
two-reflection shape), and [docs/neural/proofs](PROOF_PALINDROME_NEURAL.md)
(entry conditions and multiplicity-preserving transport).
[Experiments](../../../experiments/NEURAL_GAMMA_CAVITY.md) and its
[committed controls](../../../simulations/results/celegans_pairing_controls.txt)
supply the connectome support null and instrument checks; the
[neural clock record](../../../experiments/NEURAL_CLOCK_TWO_HANDS.md) supplies
the trace reading, not fixed individual rates. Hardware-flight searches,
including [the IBM synthesis](../../../experiments/IBM_HARDWARE_SYNTHESIS.md),
and [fw.Confirmations](../../../simulations/framework/confirmations.py) return
no neural hardware confirmation. [GLOSSARY](../../GLOSSARY.md) separates
evidence grades; [OpenArcs](../../../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs)
records conditional F36/F37 carriage with substrate questions open;
[CAUGHT_ERRORS](../../CAUGHT_ERRORS.md) records diagonal, normalization and
matcher failure shapes. The current [README](../README.md),
[algebra page](../ALGEBRAIC_PALINDROME_NEURAL.md),
[shared primitives](../../../simulations/neural/neural_palindrome.py),
[translation gate](../../../simulations/neural/neural_translation_gate.py),
[tests](../../../simulations/neural/tests/test_neural_palindrome.py), and four
producers linked below distinguish algebra, census, normalized crossing
arithmetic and non-equilibrium endpoints. The
[cavity result file](../../../simulations/results/neural_gamma_cavity.txt)
rejects its own headline pairing interpretation; it supplies no mechanism.

## The exact implication and its counterexamples

For an involutive permutation Q and one scalar s,

```
QJQ + J + 2sI = 0
    ⇒  λ ↦ −λ−2s preserves the full complex eigenvalue multiset.
```

[F36/F37 and the proof](PROOF_PALINDROME_NEURAL.md) preserve algebraic
multiplicity and transport generalized eigenspaces. The map fixes the single
spectral point −s; it exchanges −s+iω and −s−iω, not a vertical line
pointwise. For real s>0, strict asymptotic stability requires Re λ<0 for
every root; pairing then puts them in the open strip −2s<Re λ<0. The
identity supplies neither inequality.

The canonical gate's exactly representable two-seat example is

```
J = [[−0.5, −0.25], [0.25, −0.25]], Q = (0 1), s = 0.375,
QJQ + J + 2sI = 0,
λ = −0.375 ± (√3/8)i.
```

Its nonreal spectrum disproves “exact palindrome implies silence.” The
current `exact_ensemble_census` uses N=10, 5E+5I, τ_E=5, τ_I=10, density
parameter 0.3, seeds 0…199 and the paired-weight constructor in the report.
It has no nonlinear fixed-point solver or frequency-bin grid. A scalar
pass means relative Frobenius residual <10⁻¹³; complex means some
|Im λ|>10⁻⁸; unstable means max Re λ>0, without an additional tolerance:

| Linear α | Scalar passes | Complex spectra | Unstable spectra |
|---:|---:|---:|---:|
| 0.5 | 200/200 | 24 | 0 |
| 1.5 | 200/200 | 110 | 1 |
| 3.0 | 200/200 | 149 | 15 |
| 5.0 | 200/200 | 159 | 24 |
| 10.0 | 200/200 | 167 | 45 |

The [canonical gate](../../../simulations/neural/neural_translation_gate.py)
reproduces these rows and includes multiplicity, imaginary-part and wrong-Q
negative controls. These are finite constructed-network counts, not
population estimates. A silent seed supplies no silence theorem.

## A fixed coordinate is not a centre eigenmode

Writing J=D+W_eff with W_eff zero-diagonal, F36 requires
d_i+d_Q(i)+2s=0 and W_eff[Q(i),Q(j)]=−W_eff[i,j]. At Q(i)=i, the diagonal
equation only gives d_i=−s. The incident weights must still meet the
off-diagonal condition; the coordinate vector e_i need not be an eigenmode.
For example, with Q=(0 1)(2), s=3/8, this matrix satisfies the full identity:

```
J = [[−3/8,  1/4,  1/8],
     [−1/4, −3/8, −1/8],
     [ 1/8, −1/8, −3/8]].
```

Yet J e₂=(1/8,−1/8,−3/8)ᵀ: the fixed coordinate is not a centre eigenmode.
This is an algebraic counterexample, not a Dale-network model. An
odd-dimensional matrix satisfying the *full* identity does have a root at
−s: noncentral roots pair, leaving odd multiplicity at the centre. This
uses the full spectral identity, not a fixed coordinate or leak alone.

The actual odd-seat matrix in
[`veffect_exact.py`](../../../simulations/neural/veffect_exact.py) has
N_c=21 or 41, built from N=10 or 20 constituents, seeds 42 and 99,
τ_E/τ_I=5/10, α=0.5 and density parameter 0.3. Its fixed mediator has
d_M=−0.2 and s=0.15, giving diagonal defect 2d_M+2s=−0.1 at every coupling.
At c=0 its scalar residual is 0.132127 or 0.095731 respectively. These
belong to the odd-seat matrices, not their exact balanced constituents.
The per-seat fitted off-diagonal instrument gives zero there by removing
the diagonal condition; the current producer uses the full scalar test.
No alternative scalar repairs it: paired E/I seats require s=0.15 while
the fixed mediator requires s=0.2.

The two-master explanation also fails as an operator statement. If Q_A
is the A-swap extended by identity on B, the B-block of
Q_A J Q_A+J+2sI is 2(J_B+sI), already nonzero without a bridge. The script
uses the combined A/B swap with a fixed mediator; it does not construct
two separately valid equations whose competition explains the census.

## What the instrument can establish

Equal frequency tolerances do not make activity binning and pair-sum
binning commute. Sub-threshold frequencies can shift a sum across a bin
edge; distinct exact frequencies can share one activity bin. For K *exact*
positive activity frequencies, pair sums have at most K(K+1)/2 positive
sums, K(K−1)/2 positive differences and K singles from adding a real mode:
at most K(K+1). The rounded K_act cannot replace that exact K.

At N=20, c=0.05 the producer returns (K_act,K_corr)=(7,62) at ε=10⁻⁶,
(8,72) at ε/4, and (7,64) for Jᵀ at the original ε. Transposition
preserves the exact spectrum: these are numerical census limits, not a
violation of the exact-frequency ceiling. The N=10 table survives the
report's two refinements and transpose control, still as a finite census.

The [drive producer](../../../simulations/neural/veffect_and_heat.py) uses
one random Dale network, N=50, seed 42, density 0.3, τ_E/τ_I=5/10,
α=0.3 and the report's sigmoid/converged fixed-point protocol. It is
separate from both the exact ensemble and the odd-seat coupling sweep.
K_corr=124 at P=4 is a bin count at ε=10⁻⁴, not a crossing value; at ε/4
it becomes 286. P is external input, not temperature, heat or metabolism.
No biological, life or Yerkes-Dodson mechanism follows from this model.

Drive changes g_i in W_eff=α diag(g_i/τ_i)W. For the separate constructed
N=20, seed-42 sigmoid sweep at α=0.5, r_s=0.0443 at P=3.5 but
K_act=K_corr=0 at ε=10⁻⁶. Together with the nonreal exact example, this
disproves the candidate that exact pairing enforces silence and breaking
it supplies resolved oscillation. The drive producer's decay groups use
an imaginary-part threshold, not palindrome partners. Their class-mean
ratios establish no 2× law; F37 constrains a partner sum, not that ratio.

## A bounded iteration cannot locate a Hopf bifurcation

[`find_quarter.py`](../../../simulations/neural/find_quarter.py), Idea 6,
uses N=200, 100E+100I, seed 42, density 0.3, τ_E/τ_I=5/10, P=3 and the
report's E/I sigmoid parameters. Starting at x_i=0.3 it permits 3000
synchronous updates with step tolerance 10⁻¹³. Fresh equation residuals
at four returned endpoints are:

| α | Equation residual max abs(x−S(αWx+3)) |
|---:|---:|
| 5.87 | 0.6206 |
| 6.00 | 0.5162 |
| 7.00 | 0.9323 |
| 10.00 | 0.9981 |

These are not equilibria to the stated tolerance. Their Jacobian eigenvalues
support no equilibrium-stability or Hopf verdict. A Hopf candidate first
needs a converged equilibrium branch, then a nonzero imaginary pair crossing
the imaginary axis transversely and the relevant nondegeneracy checks.
Failure of fixed-point iteration does not itself establish instability
of the differential equation.

The sigmoid identity S(θ)(1−S(θ))=1/4 only gives maximal normalized slope;
the actual slope is a/4. Likewise
[`cpsi_two_perspectives.py`](../../../simulations/neural/cpsi_two_perspectives.py)
reads normalized squared E/I amplitudes p_E+p_I=1. Equality forces
p_E=p_I=1/2 and product 1/4 by arithmetic. Its alternate normalization
divides by the initial total squared norm, so that sum generally changes
with time. Neither is a quantum CΨ observable or neural stability boundary.
Its real-part greedy match does not test complex pairing or Q transport.

## Quantum shapes for the next gates

The eleven-axis rulebook generates questions only; no axis is evidence.
Observation fixes units, resolution and operating-point protocol here;
contract/residue separates F36 from each seed's census. The other concrete
translations are below. Each needs its additional neural hypothesis and
a falsifier from below. The open rows are executable gate recipes, not
claims that the proposed experiments are implemented or passed.

| Quantum owner | Neural candidate | Extra hypothesis | Executable falsifier | Current grade |
|---|---|---|---|---|
| [F89 reflection seat](../../../experiments/F89_PATH_K_DIABOLIC.md) | A fixed mediator carries a centre mode | Full scalar identity and a specified candidate vector | Evaluate `scalar_center_residual` including c=0; compute `(J+sI)@e_M`. The three-seat matrix above separates the claims | Exact diagonal constraint; present bridge rejected; coordinate interpretation false |
| [Rate/phase question](../../../reflections/ON_THE_NINETY_DEGREE_GAMMA.md), [F137](../../ANALYTICAL_FORMULAS.md) | Coupling moves phases while preserving rates | Specify mean versus individual rate; keep D and units fixed | Compare `trace(J)/n` and the full complex eigenvalue multiset across c. A moving Re λ rejects individual-rate invariance | Mean trace invariant only; individual-rate claim unsupported |
| [One ledger](../../../reflections/ON_THE_ONE_DIAGONAL.md), [shadow/source](../../proofs/PROOF_F86A_EP_MECHANISM.md) | A specific bridge or readout produces the census | Specify x₀, readout C and modes contributing to `C exp(Jt)x₀` | Compute with `scipy.linalg.expm`; remove a named bridge weight at fixed normalization and test its predicted contribution; repeat binning under transpose/refinement | Declared-resolution census; source mechanism open |
| [Quantum topology](../../THE_STAR_FROZEN_SEAM.md) | Directed support determines a nonoscillatory neural survivor | Match leak, gain, source signs and normalization; test any proposed Q separately | Build matched directed supports, compare leading eigenvalues/subspaces and test the proposed mode's null action directly | Open; no quantum graph or Laplacian law imported |
| [Circle and run](../../../reflections/ON_LEAVING_THE_CIRCLE.md) | Spectral transport predicts a transient | Fixed autonomous J, converged operating point, specified perturbation/readout and time interval | For exact F36 compare `expm(J*t)@Q` with `exp(-2*s*t)*Q@expm(-J*t)`; separately shrink nonlinear perturbations and compare responses | Static Q transport gated; observed-dynamics landing open |
| [One-way absence](../../proofs/PROOF_CODIM1_BY_ADDITIVITY.md) | Reject inadmissible support before fitting spectra | Nonzero gain, positive scales, Dale source signs and bijective Q | Count nonempty E/I source columns; unequal counts reject. Equal counts must proceed to paired support, scalar diagonal and magnitudes | Necessary-condition rejection; a passing count proves no identity |
| [F134 reflections](../../proofs/PROOF_F134_TWO_ROW_REFLECTION_LAW.md), [GammaFold](../../../compute/MirrorWorld/GammaFold.cs) | Two maps on a neural parameter family yield a shift and trajectory factor | Explicit maps/domains; generator difference must be scalar times I | Construct both transformed matrices, test every entry against the proposed shift, then propagators against its exponential factor; mutate one leak | Open family-level candidate; two unequal centres on the same finite J excluded by trace |

These gates carry repo-native questions into the neural model without
promoting a quantum analogy to evidence. No replacement mechanism is
established.

## Rerun the constraints

From the repository root; commands print to stdout. The
[report](../V_EFFECT_NEURAL.md#reproduce-the-census) supplies coupling,
drive, transpose and frequency-refinement reproductions.

POSIX:

```bash
export PYTHONIOENCODING=utf-8
python simulations/neural/neural_translation_gate.py
python -m pytest simulations/neural/tests/ -q
python simulations/neural/find_quarter.py
python simulations/neural/cpsi_two_perspectives.py
```

PowerShell:

```powershell
$env:PYTHONIOENCODING = 'utf-8'
python simulations/neural/neural_translation_gate.py
python -m pytest simulations/neural/tests/ -q
python simulations/neural/find_quarter.py
python simulations/neural/cpsi_two_perspectives.py
```
