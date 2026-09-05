# The conditional algebraic palindrome in neural networks

Last refreshed: 2026-09-05 (change history lives in git).

There are four distinct layers here: an exact conditional theorem, networks
constructed to satisfy it, a null on the committed C. elegans chemical matrix,
and a gate for the transport of modes. The theorem transfers an operator
relation from the quantum side. It does not establish that a brain has that
relation.

The named stores checked are [ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md)
(F36/F37 give the neural conditions and pairing; F137 gives the trace warning;
F157 gives a quantum blindness candidate),
[docs/proofs](../proofs/MIRROR_SYMMETRY_PROOF.md) (F1's conjugation and
[F134's two-row reflection](../proofs/PROOF_F134_TWO_ROW_REFLECTION_LAW.md)),
[docs/neural/proofs](proofs/PROOF_PALINDROME_NEURAL.md) (the conditional
derivation and the [V-effect counterexamples](proofs/PROOF_VEFFECT_MECHANISM.md)),
and [experiments](../../experiments/NEURAL_GAMMA_CAVITY.md) (the connectome null,
matcher controls and [neural trace reading](../../experiments/NEURAL_CLOCK_TWO_HANDS.md)).
The hardware-flight search and [fw.Confirmations](../../simulations/framework/confirmations.py)
returned no neural hardware confirmation. [GLOSSARY](../GLOSSARY.md) separates
evidence grades; [OpenArcs](../../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs)
records conditional F36/F37 carriage without closing the substrate questions;
[CAUGHT_ERRORS](../CAUGHT_ERRORS.md) records the support, scalar-diagonal and
normalisation hazards. The executable sources inspected are
[neural_palindrome.py](../../simulations/neural/neural_palindrome.py), its
[canonical gate](../../simulations/neural/neural_translation_gate.py) and
[tests](../../simulations/neural/tests/test_neural_palindrome.py), the
[MirrorWorld owner](../../compute/MirrorWorld/NeuralPalindrome.cs) and
[README](../../compute/MirrorWorld/README.md), and the
[C. elegans controls](../../simulations/neural/celegans_pairing_controls.py)
with their [stored output](../../simulations/results/celegans_pairing_controls.txt).
The [V-effect page](V_EFFECT_NEURAL.md) and current
[coupling](../../simulations/neural/veffect_exact.py) and
[drive](../../simulations/neural/veffect_and_heat.py) producers supply synthetic
censuses, not a biological mechanism.

## 1. Exact conditional theorem

For n ≥ 1, let J ∈ ℂⁿˣⁿ, J = D + W_eff, where D = diag(d_i) and W_eff has
zero diagonal. Let Q be an involutive permutation matrix (Q² = I), with
Q[i,Q(i)] = 1, and let s ∈ ℂ be one scalar. Then

```
Q J Q⁻¹ = −J − 2s I
```

holds **if and only if** both conditions hold:

```
d_i + d_Q(i) + 2s = 0                       for every i,
W_eff[Q(i),Q(j)] + W_eff[i,j] = 0           for every i ≠ j.
```

A permutation preserves the diagonal/off-diagonal split, so the two parts
must vanish separately. The centre is **−s**, the partner shift **−2s**.
The identity implies trace(J)/n = −s, but a correct trace does not imply the
identity or a paired spectrum. See
[F36](../ANALYTICAL_FORMULAS.md#f36-neural-palindrome-condition-tier-1-derived-algebra),
[F37](../ANALYTICAL_FORMULAS.md#f37-neural-eigenvalue-pairing-tier-1-from-f36),
and the [full proof](proofs/PROOF_PALINDROME_NEURAL.md).

For the linear neural model dx/dt = Jx, x is a deviation from a specified
operating point. With no self-coupling, d_i = −1/τ_i and τ_i > 0. If each
seat has τ_E or τ_I and Q swaps the types, then
s = ½(1/τ_E + 1/τ_I). At unequal time constants this choice of s requires
opposite-type exchange. At equal time constants the diagonal condition holds
for every involutive Q. A fixed seat Q(i) = i requires d_i = −s; it must
also satisfy every off-diagonal condition. A fixed seat is a coordinate of Q,
not automatically an eigenmode at the spectral centre.

In the common-gain convention W[i,j] is the connection **from j to i**, and
W_eff = α diag(1/τ_i)W with W[i,i] = 0. If α ≠ 0 the coupling condition is

```
W[Q(i),Q(j)] = −(τ_Q(i)/τ_i) W[i,j].
```

The nonzero-gain hypothesis matters: at α = 0, W_eff = 0 and this W-only
condition is unnecessary. In a sigmoid linearisation, row-dependent slopes
belong in W_eff. The W-only ratio above assumes a common slope α. Any
self-coupling belongs in D and changes the diagonal test.

For real Dale weights, the source column fixes the sign of every nonzero
outgoing connection. A type-swapping Q gives opposite signs **where both
partner entries exist**. It does not produce Q-symmetric support or the
required scaled magnitudes. The chemical-connectome file instead stores
presynaptic sources in rows: transpose to the column-source convention before
applying a row-time-constant formula. Transposition preserves the support
question, but row scaling becomes column scaling under transpose.

## 2. Constructed exact networks

The Python constructor builds balanced E/I populations, an involution and
partner weights satisfying the two conditions. This is a construction of
admissible matrices, not an estimate of how frequently biology supplies them.
The canonical gate includes

```
J = [[−0.5, −0.25], [0.25, −0.25]],  Q = (0 1),  s = 0.375,
QJQ + J + 2sI = 0,
λ = −0.375 ± (√3/8)i.
```

All entries are dyadic; the entrywise residual is exactly zero in the run.
The [MirrorWorld helper](../../compute/MirrorWorld/NeuralPalindrome.cs) checks
the real matrix entrywise without an eigensolver, and its
[tests](../../compute/MirrorWorld.Tests/NeuralPalindromeTests.cs) and `neural`
run expose this example and the fixed-seat counterexample
J = diag(−0.5, −0.25, −0.5), Q = (0 1)(2), s = 0.375. The latter has zero
off-diagonal residual but maximum absolute scalar-identity residual 0.25.

The [current 200-seed census](README.md#what-has-been-tested) is N = 10,
τ_E = 5, τ_I = 10, density 0.3, seeds 0…199. At α = 0.5, 1.5, 3, 5 and 10,
all 200 matrices per row pass the scalar residual < 10⁻¹³; the oscillatory
counts are 24, 110, 149, 159 and 167, and the unstable counts 0, 1, 15, 24 and
45. Oscillatory means some |Im λ| > 10⁻⁸; unstable means max Re λ > 0.
These numerical census predicates do not replace the exact entry conditions.

Exact palindrome alone implies neither real eigenvalues, stability, silence,
nor a biological mechanism. In particular, the nonreal two-seat example is
already a counterexample to “palindrome implies silence.” For real s > 0,
strict asymptotic stability requires every Re λ < 0; paired roots then lie in
the open strip −2s < Re λ < 0. Pairing itself does not put them there.

## 3. Empirical C. elegans null

No biological neural network in the repository is known to satisfy F36.
On the [committed chemical data](../../simulations/neural/celegans_connectome.json),
G0b of the [connectome audit](../../simulations/neural/celegans_pairing_controls.py)
counts **271 nonempty outgoing rows: 253 excitatory, 18 inhibitory** under
the stored labels. A sign-reversing Q must biject those two nonempty sets.
The unequal counts rule out a qualifying permutation before magnitudes are
tested, for the full Dale-signed chemical matrix with nonzero gain and
positive rate scales. The [output](../../simulations/results/celegans_pairing_controls.txt)
and [event record](../../experiments/NEURAL_GAMMA_CAVITY.md) hold the result.
An equal count would only remove this obstruction; it would not prove F36.
Selected subnetworks and other models require their own support and rate tests.

The fitted off-diagonal instrument has a different scope from F36. Write
R_off = QW_effQ + W_eff, with W_eff zero-diagonal. If the nonzero supports of
W_eff and QW_effQ are disjoint, orthogonality in Frobenius norm gives

```
||R_off||_F / ||J||_F = √2 ||W_eff||_F / ||J||_F       (J ≠ 0).
```

Here the instrument reads coupling magnitude. Off this disjoint-support set
it can register the placement of partner edges. An empty coupling can score
zero without testing a wiring hypothesis; a per-seat fitted centre can also
hide a failing scalar diagonal condition.

The stored matched-normalisation audit uses 200 balanced blocks per size,
τ_E/τ_I = 10/20, α = 0.3, and `RandomState(trial + 100)`. Both the connectome
block and its random control are divided by their own maximum absolute
weight (empty blocks remain zero):

| Block size | Ratio of mean residuals, connectome/control | Connectome blocks with disjoint partner support |
|---|---:|---:|
| 10 | 0.9602 | 198/200 |
| 20 | 0.8414 | 184/200 |
| 26 | 0.7484 | 177/200 |

These are [G0d/G0e's stored measurements](../../simulations/results/celegans_pairing_controls.txt).
The larger blocks retain a smaller residual than the control; the comparison
does not establish a wiring mechanism. Most blocks meet the norm-collapse
condition, while a minority do not. A control that keeps each weight in its
own row also preserves the row-scaled coupling norm, so equality under that
rewiring is not independent evidence for pairing on the collapse set.

## 4. Gated mode transport

F36 gives JQ = −QJ − 2sQ, hence the exact relation

```
Jv = λv  ⇒  J(Qv) = (−λ − 2s)(Qv).
(J + (λ + 2s)I)^k Q = (−1)^k Q (J − λI)^k.
```

The second line transports generalized eigenspaces as well as ordinary
eigenvectors, including defective eigenvalues. Similarity preserves algebraic
multiplicity. The map fixes only the spectral point λ = −s; it sends
−s + iω to −s − iω, rather than fixing an entire vertical line pointwise.

If Q exchanges E and I, it is norm-preserving and exchanges the squared E/I
amplitudes of v and its **Q-transported** partner. This statement is about
that vector and its image, not an eigensolver's arbitrary choice of basis
inside a repeated eigenvalue. At a degenerate eigenvalue use the whole
invariant subspace.

The [canonical gate](../../simulations/neural/neural_translation_gate.py) calls
`partner_subspace_error` in [neural_palindrome.py](../../simulations/neural/neural_palindrome.py).
It forms eigenvalue clusters by connected components at complex distance
10⁻⁷, uses ordered complex Schur vectors for their invariant subspaces, and
reports the largest sine of principal angles between Q U_C and U_(−C−2s).
It requires matching algebraic dimensions and raises on a missing partner.
This is basis-invariant within the selected subspaces; the cluster tolerance
is a numerical resolution choice, not proof of exact degeneracy.

The gate accepts the constructed network below transport sine 10⁻⁸ and the
fully degenerate example below 10⁻¹². A wrong permutation on a matrix whose
spectrum still pairs gives sine 1 and is rejected. The focused tests also
cover defective clusters. Separately, `spectral_pairing_error` assigns the
entire complex eigenvalue multiset to −λ−2s by minimum total cost and reports
the largest assigned distance. It preserves multiplicity and imaginary parts;
it is not a bottleneck-optimal matcher and does not establish Q transport by
itself.

## Translation roads

The repo's eleven-axis rulebook is used here only to generate candidate
directions; an axis is not evidence. The full walk asks about observation,
mirror, darkness, rates versus phases, the one ledger, topology, contract
versus residue, circle versus run, one-way absence, fixed seats, and two
mirrors. The concrete landings below are separately graded. None transfers a
quantum graph or Laplacian result to a neural network without the stated
support, sign, permutation and dynamics hypotheses.

| Quantum owner → neural candidate | Extra neural hypothesis | From-below gate | Current landing / grade |
|---|---|---|---|
| [F1 mirror](../proofs/MIRROR_SYMMETRY_PROOF.md) → an E/I pair with one fixed spectral span (axis 2) | One involutive permutation; both scalar-diagonal and effective-weight conditions, not Dale signs alone | Full scalar residual, full complex assignment, and Q-subspace transport; mutate a leak or Q | F36/F37 conditional theorem; constructed gates pass; full connectome support null |
| [F89 reflection-seat question](../../experiments/F89_PATH_K_DIABOLIC.md) → a mediator fixed by Q (axis 10) | The neural fixed seat needs d_i = −s and all incident paired weights; no quantum reality conclusion is inherited | Include the mediator diagonal at zero coupling before scanning coupling; the dyadic fixed-seat control rejects it | Exact neural diagonal constraint; the current coupled construction fails it |
| [Clock and trace](../../experiments/NEURAL_CLOCK_TWO_HANDS.md), with [F137](../ANALYTICAL_FORMULAS.md) → invariant spectral mean versus moving individual rates/frequencies (axes 1, 4) | Fixed D, zero-diagonal coupling, fixed time unit and operating-point convention | Record trace(J)/n alongside each Re λ and Im λ; vary coupling with D held fixed | Trace invariant is algebraic; invariance of individual rates is not implied |
| [F157 BlindSeat](../../compute/MirrorWorld/BlindSeat.cs) → a mode invisible to a chosen neural readout C (axes 3, 5) | Specify dx/dt = Jx and y = Cx; readout has no feedback into this J. The quantum real-symmetric H/gcd formula is not assumed for directed J | Build the observability rows C, CJ, …, CJ^(n−1); test a candidate in their joint kernel and propagate it; perturb C so the proposed blindness can fail | Open neural direction; invisibility to C would not by itself imply slow decay or stability |
| [Topology-dependent quantum survivor](../THE_STAR_FROZEN_SEAM.md) → a topology-dependent neural slow mode (axis 6) | Fix the neural leak, gain, Dale assignment and normalization; specify a family of directed supports and any proposed Q | Compare leading eigenvalues and invariant subspaces across supports with matched scales; directly check any claimed null action on a mode | Open neural direction; no inherited chain, star or Laplacian law |
| [Mirror versus the particular response](../../reflections/ON_WHAT_CLOSES_ONLY_WITHOUT_US.md) → separate F36 from one seed's census (axis 7) | State which matrices satisfy the identity and which seed/drive protocol supplies the observation | Keep the operator test beside the 200-seed census; a nonreal exact example defeats a universal silence claim | Conditional algebra and finite census are distinct; no universal V-effect mechanism |
| [Circle and run](../../reflections/ON_LEAVING_THE_CIRCLE.md) → compare a paired spectrum with an actual neural trajectory (axis 8) | Fixed autonomous linearization, specified initial perturbation and readout, justified time interval; nonlinear departures require a separate model | Propagate x and Qx, measure the chosen readout, and test whether observed frequencies and transient amplitudes agree with that J | Static transport is gated above; a dynamical or nonlinear biological landing remains open |
| [One-way absence certificate](../proofs/PROOF_CODIM1_BY_ADDITIVITY.md) → cheap support rejection before spectral fitting (axis 9) | Dale sign by source, nonzero gain and positive scales, a bijective Q | Count nonempty E/I outputs; if counts agree, continue to paired support and magnitudes rather than declare success | G0b certifies absence for the committed full chemical matrix; passing the count alone certifies no presence |
| [F134's two reflections](../proofs/PROOF_F134_TWO_ROW_REFLECTION_LAW.md) and [GammaFold](../../compute/MirrorWorld/GammaFold.cs) → a parameter-family transformation and a time-dependent scalar factor (axis 11) | Name two actual maps on specified neural matrices/parameters and their domains; an affine rate shift must be a scalar generator shift | Construct both transformed matrices; test their entrywise identities and then propagators with the predicted factor; a same-matrix pair of different centres fails the trace condition | Open direction. One finite J has only s = −trace(J)/n; two unequal centres cannot both be its scalar palindrome. No two-mirror explanation of the V-effect is established |

The local-readout and parameter-family rows are candidate experiments, not
new neural laws. Observation provides no additional neural measurement-cost
claim here: the drive P, leak times and frequency resolution must first have
declared units. The concrete obligation is to keep those readings distinct.
Similarly, the quantum threshold CΨ = 1/4 has no neural threshold supplied by
F36/F37. A neural bifurcation requires its own operating-point and eigenvalue
crossing analysis.

## Reproduce and continue

Run the [README commands](README.md#how-to-rerun) from the repository root for
the canonical Python gate, focused Python tests, MirrorWorld neural tests and
`neural` run. Both POSIX and PowerShell UTF-8 forms are provided there.
For connectome work, start with the [data](../../simulations/neural/celegans_connectome.json)
and [support/normalisation audit](../../simulations/neural/celegans_pairing_controls.py),
whose full run rewrites [its result file](../../simulations/results/celegans_pairing_controls.txt).
For dynamics, use the [coupling/drive account](V_EFFECT_NEURAL.md) and its
current producers with explicit frequency resolution and fixed-point residuals.
The V-effect remains a coupling/drive census with mechanism open, not a 2× law.
