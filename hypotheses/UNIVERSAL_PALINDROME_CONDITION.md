# The Universal Palindrome Condition

**Status:** Tier 4 cross-domain hypothesis; Tier 1 conditional algebra underneath
**Authors:** Thomas Wicht, Claude (Anthropic)
last refreshed 2026-09-05 (the change history lives in git)

One operator identity can organize spectra in different mathematical models.
The quantum theorem supplies a physical instance under its stated assumptions;
the neural theorem supplies conditions and synthetic constructions. No biological
neural network in this repository is known to satisfy those conditions. Whether
the identity has a biological realization, and whether such a realization is
inherited from quantum dynamics, are separate open questions.

The stores checked are [the F-registry](../docs/ANALYTICAL_FORMULAS.md)
(F1, F36/F37 and the trace distinction F137),
[docs/proofs](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) (the quantum owner),
[the neural proofs](../docs/neural/proofs/PROOF_PALINDROME_NEURAL.md)
(scalar conditions and spectral transport), and
[experiments](../experiments/NEURAL_GAMMA_CAVITY.md) (the connectome support null).
The hardware-flight search, including [the IBM synthesis](../experiments/IBM_HARDWARE_SYNTHESIS.md),
and [fw.Confirmations](../simulations/framework/confirmations.py) supplied no
neural hardware confirmation. [GLOSSARY](../docs/GLOSSARY.md) separates objects
and interpretive readings; [OpenArcs](../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs)
records conditional F36/F37 carriage and open substrate questions;
[CAUGHT_ERRORS](../docs/CAUGHT_ERRORS.md) records diagonal, normalization and
matcher failures. The [current neural account](../docs/neural/README.md),
[Python gate](../simulations/neural/neural_translation_gate.py), and
[frequency producers](../docs/neural/V_EFFECT_NEURAL.md#reproduce-the-census)
are the evidence anchors below.

## The exact identity and its reach

For a finite matrix X, an invertible linear map Q and one scalar s,

```
Q X Q⁻¹ = −X − 2s I
```

makes X similar to −X−2sI. Its full complex eigenvalue multiset therefore
pairs under λ ↦ −λ−2s, including algebraic multiplicity. This implication
is algebra. Finding a physical generator and a Q satisfying it is the
substantive task.

The [quantum owner F1](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) gives
ΠLΠ⁻¹ = −L−2Σγ I for its Heisenberg/XXZ family under local Z-dephasing.
Its Π need not be involutive. This theorem is not a statement about arbitrary
Hamiltonians, channels, atoms or living systems.

The [neural F36 theorem](../docs/neural/proofs/PROOF_PALINDROME_NEURAL.md)
uses an involutive permutation Q²=I and J=D+W_eff, with D=diag(d_i)
and W_eff having zero diagonal. The identity holds if and only if

```
d_i + d_Q(i) + 2s = 0                     for every i,
W_eff[Q(i),Q(j)] + W_eff[i,j] = 0          for every i ≠ j.
```

With d_i=−1/τ_i and Q exchanging the two E/I time-constant populations,
s=(1/τ_E+1/τ_I)/2. For that s, unequal time constants require opposite-type
exchange. Equal time constants make the diagonal condition automatic for any
involutive Q; the coupling condition still has to hold. A fixed seat requires
d_i=−s as well as the coupling condition.

For W_eff=α diag(1/τ_i)W and nonzero α, with column j the source,
the coupling condition is W[Q(i),Q(j)]=−(τ_Q(i)/τ_i)W[i,j]. At α=0
the effective coupling vanishes, so this W-only condition is not necessary.
Nonuniform sigmoid gains belong in W_eff at the specified operating point.
Dale's Law supplies source signs, not paired support or scaled magnitudes.

| Component | Quantum owner | Neural candidate |
|---|---|---|
| Generator | Liouvillian L in F1's family | Specified Jacobian J |
| Mirror | Pauli-space conjugation Π | Involutive permutation Q |
| Coupling relation | Conjugation of the commutator in F1 | Additional signed support and magnitude relation on W_eff |
| Diagonal relation | Dephasing shift 2Σγ | One paired diagonal sum 2s |
| Spectral consequence | λ+λ′=−2Σγ | μ+μ′=−2s, conditional on F36 |
| Mode statement | Transport under the owner's conjugation | Q transports generalized eigenspaces; no dominance percentage follows |
| Current grade | Proven within the quantum model's scope | Constructed instances; biological landing missing |

For F36, Jv=μv implies J(Qv)=(−μ−2s)Qv. At degeneracy compare
invariant subspaces rather than arbitrary eigensolver basis vectors.
An E/I coordinate swap exchanges their squared weights, but does not
require either population to dominate a mode. Pairing alone implies
neither realness, silence nor stability.

## The biological gate is currently negative

The [committed C. elegans chemical matrix](../simulations/neural/celegans_connectome.json),
under its stored Dale labels, has **271 nonempty output rows: 253 E and 18 I**.
A sign-reversing support permutation needs a bijection between the nonempty
E and I outputs. These unequal counts rule it out for the full matrix with
nonzero gain and positive rate scales, before fitting magnitudes.
[G0b and its output](../simulations/results/celegans_pairing_controls.txt)
own this obstruction. The stored matrix uses source rows; the common-gain
formula above uses source columns.

This does not decide every selected subnetwork or every neural model.
A balanced subnetwork must still pass both F36 conditions on its effective
Jacobian. A mean pair sum fixed by trace, a permissive spectral matcher or
a coupling-norm residual cannot certify that identity. The
[matched audit](../docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md#3-empirical-c-elegans-null)
explains why normalization and a null capable of changing the measured
quantity matter. E/I count balance alone is not a mechanism.

## Coupling and drive: finite censuses, open mechanism

The [neural V-effect](../docs/neural/V_EFFECT_NEURAL.md) is an operational
frequency census on specified synthetic protocols. It counts rounded bins
of |Im λ| and |Im(λ_i+λ_j)|, with an explicit cutoff ε. It does not measure
the birth of persistent nonlinear oscillations or biological rhythms.

| Protocol | Current reading | Limit of the reading |
|---|---|---|
| Exact N=10 linear ensemble, 200 seeds, τ_E/τ_I=5/10, density 0.3, α=0.5 | 24 draws have some abs(Im λ)>10⁻⁸ while all pass scalar residual threshold 10⁻¹³ | Exact palindrome does not enforce silence |
| Coupled N=20 seed-42/99 constituents, α=0.5, same time constants, ε=10⁻⁶ | 48 correlation bins at c=0.01; 62 at c=0.05 | Odd mediator violates F36 even at c=0; counts depend on resolution and backend |
| Random Dale N=50, 25E+25I, seed 42, density 0.3, α=0.3, same time constants, P=4 | 124 correlation bins at ε=10⁻⁴; 286 at ε/4; 403 at ε/16 | P is external input; these are not constants or a temperature calibration |

The linked report specifies constructors, normalization, grids, sigmoid,
solver and backend. The coupled bridge also need not preserve Dale source
signs. Neither its counts nor the drive sweep establish a universal coupling
window, a universal silent–active–silent shape, or a neural V-effect mechanism.
The quantum [V-effect experiment](../experiments/V_EFFECT_PALINDROME.md)
remains a separate model and observable; its interpretation does not prove
the neural mechanism. There is no established 2× neural decay law.
External P is not heat, metabolism or life.

The [mechanism constraints](../docs/neural/proofs/PROOF_VEFFECT_MECHANISM.md)
also show why unconverged fixed-point iterations cannot locate a Hopf
bifurcation or decide equilibrium stability. Those require a converged
equilibrium branch and additional dynamical checks.

## What the quarter does and does not transfer

For nonzero total squared activity, define
p_E=||x_E||²/(||x_E||²+||x_I||²) and p_I=1−p_E. Then

```
p_E p_I = 1/4 − (p_E−1/2)² ≤ 1/4.
```

Equality means equal squared amplitudes. This is bookkeeping, not the
quantum CΨ=1/4 boundary or evidence for a shared physical mechanism.
Likewise a logistic sigmoid has S(1−S)≤1/4; with steepness a its maximum
slope is a/4. Neither identity locates a neural stability transition.
See the [normalization analysis](../docs/neural/proofs/PROOF_VEFFECT_MECHANISM.md)
and the [quantum quarter's scope](../docs/Q_BELONGS_TO_NO_SUBSTANCE.md).

## The hypothesis and the next gates

Can one spectral organization recur across levels through a physical
reduction preserving its conjugation? The identity alone forces no hierarchy,
creation of a next level or persistence. Each proposed step needs a landing.

| Quantum owner | Neural candidate | Extra hypothesis | Gate | Current grade |
|---|---|---|---|---|
| F1 conjugation | F36 on a biological J | A circuit, operating point, Q and s satisfy both entry conditions | Test support, diagonal and magnitudes with scalar_center_residual; then full-complex multiset pairing | Constructed algebra passes; full-connectome support rejects; biological instance open |
| F1 mode transport | A mirrored neural response | Fixed autonomous J, specified perturbation/readout and validity of linearization | Compare expm(J*t)@Q with exp(-2*s*t)*Q@expm(-J*t), then observed responses | Static transport gated; biological dynamics open |
| Quantum coupled-system experiment | Coupling-dependent neural census | A readout and causal bridge contribution beyond binning | Refine ε, transpose J, vary seeds and remove named bridge entries at fixed normalization | Finite census; mechanism open |
| Quantum-to-effective description | Inheritance into a neural model | An explicit reduction taking the quantum generator/conjugation into J and Q | Evaluate proposed generator and conjugation intertwining equations on the same model | No inheritance derivation supplied |

The first two rows can start from the [shared Python primitives](../simulations/neural/neural_palindrome.py)
and [translation gate](../simulations/neural/neural_translation_gate.py).
The latter two are gate recipes, not implemented results. Another candidate
substrate needs its own generator and linear conjugation; two named
populations or an intuitive swap do not establish the equation.

Failure of a biological candidate rejects that landing, not the conditional
theorem. An inheritance claim needs both a biological instance and a
derivation connecting the levels. At present we have the algebra and the question.
