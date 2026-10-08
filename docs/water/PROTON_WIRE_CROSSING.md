# The Proton Wire Crossing: what the framework's grading IS in water

**Status:** Tier 2 (computed from proven framework)
**Date:** 2026-07-31
**Authors:** Thomas Wicht, Claude
**Gate:** [`simulations/water/proton_wire_crossing.py`](../../simulations/water/proton_wire_crossing.py), 55 checks
**Companion gates:** [`simulations/f98_scope.py`](../../simulations/f98_scope.py), 30 checks; [`simulations/water/cube_on_water.py`](../../simulations/water/cube_on_water.py), the cube on the wire, 35 checks

Every result in this folder is graded by one number: the popcount `Ŵ = Σ_l (I − Z_l)/2`.
The F4 kernel is `span(P_0, …, P_N)`, one projector per popcount. F98's long-time value
is a ratio of binomials in the popcount. F86b's static Dicke anchor instead has its own
state premise: a specified two-popcount Dicke superposition. F88b likewise reads its
specified states through popcount sectors. These individual state premises neither make
the F86b superposition a popcount eigenstate nor transfer automatically through `[Ŵ, H] = 0`.

The [carbon crossing pass](../carbon/BENZENE_THREE_DEPHASE_LETTERS.md) had to ask what
that operator is in its substrate, and the answer there was easy: the carbon qubit is
the occupation of a π site, so `Ŵ` is the π-electron count, and the filter
`[Ŵ, H] = 0` is the statement that a molecule does not spontaneously ionize.

Water cannot copy that. Its qubit is the **position** of one proton inside its own
hydrogen bond, |L⟩ donor or |R⟩ acceptor ([the hydrogen bond as a qubit](HYDROGEN_BOND_QUBIT.md)).
The proton is always there. Nothing is being counted. So what is `Ŵ` here, and does
anything in this folder still mean what it says?

---

## The answer: popcount is the wire's dipole moment

Set `s_l = 1` when the proton in bond `l` sits on its right-hand oxygen. A wire of N
bonds has N+1 oxygens. Against the neutral reference `s = 0…0`, in which each water
keeps its own two covalent hydrogens, the wire charges are

```
q_0 = −s_0        q_j = s_{j−1} − s_j  (0 < j < N)        q_N = +s_{N−1}
```

These sum to zero in every configuration, so the wire is neutral throughout and the
dipole moment `μ = Σ_j j·q_j` does not depend on where the origin is put. Telescoping
the sum gives

```
μ  =  Σ_l s_l  =  popcount
```

identically, in units of `e` times the O···O spacing. Verified on every basis state at
N = 3, 4, 5, 6: `max|μ − popcount| = 0`, and total charge exactly 0 in all 2^N
configurations.

So the framework's grading is not meaningless in water, and it is not the same thing it
was in carbon. It is **the number of protons that have crossed their own bond**, which
is the total charge displaced along the wire summed over its bonds.

It is *not* a flux through a single cut. The charge lying to the right of the cut at
bond `k` is `Σ_{j>k} q_j = s_k`, the state of that one bond and nothing more (checked
on every basis state at N = 3, 4, 5, zero violations). The popcount is the sum of those
N single-bond readings, so `dμ/dt` is the sum of the N bond currents rather than "the"
proton current.

| | carbon | water wire |
|---|---|---|
| qubit | occupation of a π site | position of a proton in its bond |
| `Ŵ = Σ_l (I − Z_l)/2` | π-electron number | dipole moment = displaced charge |
| `[Ŵ, H] = 0` means | the molecule does not ionize | **the total dipole is fixed** |

### What it is not: the defect count

The obvious guess for a proton wire's currency is the number of ionic defects, and the
obvious operator for it is the domain-wall count `Σ_b (I − Z_l Z_{l+1})/2`, since a
charged oxygen looks like a domain wall in the proton-position string. That operator is
not the charged-oxygen count. It is blind to both termini, where a wire of N bonds keeps
two of its N+1 oxygens: it disagrees with the true count on 12 of 16 basis states at
N = 4 and 24 of 32 at N = 5, and `s = 1111`, a proton pushed the full length of the
wire, has `q = (−1, 0, 0, 0, +1)`, an OH⁻ at one end and an H₃O⁺ at the other, and zero
domain walls. It also counts `|q|`, so it cannot tell H₃O⁺ from OH⁻ at all.

### The bias field is the same operator

`Σ_l Z_l = N·I − 2·Ŵ`, verified to machine zero. The double-well bias `Δ·Σ_l Z_l`,
which [the hydrogen bond as a qubit](HYDROGEN_BOND_QUBIT.md) introduces as the asymmetry
of an unequal hydrogen bond, is therefore exactly the wire's dipole coupled to a uniform
field along it. This matters below, because that term is the one that breaks the
palindrome.

---

## What the grading buys, and what it costs

Which of this folder's operators leave the dipole alone? Measured as
`‖[Ŵ, H]‖_F / ‖H‖_F`, at N = 4 and N = 5:

| Hamiltonian | N = 4 | N = 5 | moves the dipole? |
|---|---|---|---|
| Heisenberg `J(XX+YY+ZZ)` | 0 | 0 | no |
| Ising `K Σ ZZ` | 0 | 0 | no |
| bias `Δ Σ Z` | 0 | 0 | no |
| tunneling `−J Σ X` | 1.000000 | 1.000000 | **yes** |
| TFI `−J Σ X + K Σ ZZ` | 0.917663 | 0.912871 | **yes** |

The ratio 1.000000 on the tunneling term is an identity, not a saturated bound:
`[Ŵ, X_l] = −i Y_l` and `‖Σ Y_l‖_F = ‖Σ X_l‖_F`, so the ratio is 1 by construction.
The XX chain reaches 1.414214, so 1 is not a ceiling.

The split is clean and it is the wrong way round for the folder. Every Hamiltonian that
conserves the framework's grading holds the wire's total dipole fixed. The one operator
that changes it, the tunneling term, is the elementary proton hop. And the model this
folder itself calls "the physical proton model", the transverse-field Ising chain of
[the proton water chain](PROTON_WATER_CHAIN.md), is built on that term.

This shows up in the kernel. Under the same Z-dephasing bath:

| | N = 3 | N = 4 | N = 5 |
|---|---|---|---|
| Heisenberg, `dim ker L` | 4 | 5 | 6 |
| TFI (physical proton model), `dim ker L` | 1 | 1 | 1 |

The Heisenberg column is F4's `N+1`, and its kernel basis is the popcount projectors.
The physical model has a one-dimensional kernel: the maximally mixed state and nothing
else. Reading `Ŵ` as the dipole says why. F4's kernel is one stationary mode per value
of the dipole, which is a conservation law only as long as the dipole is fixed.

**The fixed-dipole/dynamical results in this folder are the F4 kernel and F98's
`α(∞) = (N+2)/[4(N+1)]`; their transfer to the wire requires `[Ŵ, H] = 0`.**
F86b's static Dicke anchor remains under its own two-popcount Dicke-superposition
premise, and F88b's memory split under its individual state premise. Neither is an
automatic fixed-dipole transfer through `[Ŵ, H] = 0`.

Fixed dipole is not zero motion, and the difference is worth being exact about. Under
`XX+YY` protons do cross their bonds; they cross in correlated pairs whose displacements
cancel. Acting on `s = [1,1,0,0]`, the swap across bonds 1 and 2 shifts the oxygen
charges by `Δq = (0, +1, −2, +1, 0)`: two protons move simultaneously in opposite
directions in two different double wells, and the dipole is unchanged. What the
framework's grading forbids is not motion but a net displacement of charge along the
wire. Note also that the sectors are the level sets of `μ`, not `μ = 0`.

That is the inheritance boundary this pass was looking for. It is not a defect in those
results; it is what they are about.

---

## What crosses anyway: the F1 palindrome

The palindrome does not care about the grading. Under Z-dephasing with `Σγ = Nγ`, all
norms Frobenius:

| Hamiltonian | N = 3 | N = 4 | N = 5 |
|---|---|---|---|
| tunneling `−J Σ X` | 6.48e-15 | 1.44e-14 | 3.07e-14 |
| Ising `K Σ ZZ` | 6.48e-15 | 1.44e-14 | 3.07e-14 |
| TFI (physical proton model) | 6.48e-15 | 1.44e-14 | 3.07e-14 |

`‖M‖_F` at machine zero, the physical proton model included. This is F87's criterion
doing its work, and it is insensitive to whether the Hamiltonian moves the dipole.

The bias breaks it, and breaks it exactly linearly:

```
Δ = 0.0  →  ‖M‖_F = 0            Δ = 0.3  →  ‖M‖_F = 19.2
Δ = 0.1  →  ‖M‖_F = 6.4          Δ = 1.0  →  ‖M‖_F = 64.0
```

ratio 3.000000000 against the predicted 3. And no arrangement of unequal bonds repairs
it: the linear map from a per-site bias profile `(δ_0, …, δ_{N−1})` to the residual M is
injective at N = 3, 4, 5, so no nonzero profile lies in a kernel because there is no
kernel. Its singular values are moreover all equal, to `2^(N+1)` exactly (16, 32, 64),
which says the extra thing that the map is a scaled isometry: the size of the break
depends only on `‖δ‖`, never on how the bias is distributed.

The residual M is measured against one fixed mirror. [F158](../proofs/PROOF_PALINDROME_TWO_END_COUNT.md)
answers for all of them, and the answer is the same. With Z dephasing on every site at
positive rates and every tunnelling amplitude `J_l ≠ 0`, the far kernel
`𝒲 = {U : [H, U] = 0, U anticommutes with every Z_l}` of
`H = −Σ J_l X_l + Σ K_ab Z_a Z_b + Σ δ_l Z_l` is spanned by `X^⊗N` when there is no bias
and is empty as soon as one `δ_l ≠ 0`, whatever the couplings `K_ab`. An operator that
anticommutes with every `Z_l` is `D·X^⊗N` with D diagonal; the tunnelling terms make D
the same on every pair of basis states one flip apart, and with every `J_l ≠ 0` the flips
connect all of them, so `D = c·I`; then `[Σ δ_l Z_l, X^⊗N] = 2(Σ δ_l Z_l)·X^⊗N` forces every
`δ_l = 0`. So by F158 no operator carries the palindrome under any bias profile: not a
weighted sum of strings, not a rotated frame (the ones
[the palindrome as a colouring](../../experiments/THE_PALINDROME_AS_A_COLOURING.md) builds), and
not for a profile that sums to zero (W7, exact at N = 3 and 4: dimension 1 without bias, 0
under a uniform bias and under the profile `(δ, −δ, …)`).

Two things this measurement does **not** say, both worth stating because the first
reading of it got them wrong:

- **The bias is not the only breaker, and this is F87, not a new law.** The criterion is
  F87's: a term survives when `#Y` and `#Z` are **separately** even. At N = 4, standalone
  and with no tunneling term present, `Σ X_l` gives `‖M‖_F = 0`, `Σ Y_l` gives `64.0`,
  `Σ Z_l` gives `64.0`, `Σ (XZ + ZX)` gives `78.384`, and `Σ (YZ + ZY)` gives `110.851`.
  That last one is the case that matters: `#Y` and `#Z` are each odd while their sum is
  even, so a criterion phrased on the combined parity would let it through, and it is the
  repo's canonical soft case. Note that `Σ Y_l` breaks on its own, not through
  interference with the tunneling term; X is the only single-site axis that survives
  Z-dephasing.
- **`‖M‖_F` is not a severity meter.** At the same `Δ = 0.3` both models give
  `‖M‖_F = 19.2`, but the damage differs. On Heisenberg the decay-**rate** palindrome
  survives exactly (pairing error 2.22e-14) and only the frequencies are lost. On the
  physical proton model the rates go too (pairing error 4.08e-01). That is F87's own
  soft/hard split showing up on this inventory: the operator residual cannot see it, the
  spectrum can.

---

## What a field along the wire reads

A field along the wire couples to its dipole. With `P = Ŵ − N/2 = −½ Σ_l Z_l`, the popcount
dipole of the first section measured from its middle, the pulsed generator is
`H(E) = H + E·P`, and the dephasing stays `Z_l` on every site. Three things follow, each
exact; the derivations are local and hold at every N, and gate W8 checks them
symbolically in the couplings, rates and field at N = 3 (the analysis was first worked
out by a second model, Codex).

- **Reversing the field reverses the response, without bias.** `X^⊗N` sends `H(E)` to
  `H(−E)` and `P` to `−P`, and leaves the Z dephasing alone, so from any preparation it
  leaves unchanged, `⟨P⟩_E(t) = −⟨P⟩_{−E}(t)` at all times. That is F131's mirror
  order-sorting with the mirror `X^⊗N`, sighted for a longitudinal field in
  [the h thread](../../experiments/LATTICE_H_THREAD.md) §1 (`X^N·H(h)·X^N = H(−h)`), here on
  the wire; its §3 sweeps the four readout cells of the same mirror. The response is not trivially zero: from `|+⟩^⊗N` its curvature at the start
  is `P''(0) = −E·Σ_l J_l`.
- **The total dipole can hide a bias.** On a wire that is mirror-symmetric, in its
  tunnelling, couplings and rates, a bias that is odd under the mirror
  (`δ_l = −δ_{N−1−l}`) keeps the reversal parity of the total dipole: `S = R·X^⊗N`, with R
  the site reversal, again sends `H(E)` to `H(−E)` and `P` to `−P` and leaves the
  dephasing invariant. The far kernel is empty all the same (previous section), so the
  parity of the total dipole is no test for the absence of a bias.
- **Each coordinate reads its own bias.** From `|+⟩^⊗N`, with `p_l = −Z_l/2`,
  `p_l''(0)` at `+E` plus `p_l''(0)` at `−E` equals `4·J_l·δ_l`, whatever the couplings,
  rates and field. The total dipole sees only the sum `4·Σ_l J_l δ_l`, which the
  mirror-odd bias makes zero while each end still shows it.

What a real wire's spectroscopy resolves, and whether single coordinates are
addressable at all, is not claimed here.

## The cube on the wire

[The letter cube](../THE_ONE_SQUARE.md#7-the-cube-three-squares-at-once) places a Pauli string at
(k_Z, k_X, k_Y), the number of its letters that anticommute with Z, X and Y, and reads a
dephasing along a letter as a height along that axis. On the wire the dephasing letter is
Z, the environment reading whether each proton sits on its left or its right oxygen, so a
coherence `|i⟩⟨j|` between two proton configurations, a sum of strings that all share one
k_Z, sits at `k_Z = popcount(i ⊕ j)`, the number of bonds l in which the two configurations
disagree, and decays at `−2γ·k_Z`. That is MirrorWorld's disagreement count, the one thing
its `Pair` lets the watching read, and it is not a dipole difference: the dipole
`μ = popcount` is diagonal and lives on the face `k_Z = 0`, and `|01⟩⟨10|` has the same μ on
both sides while sitting at `k_Z = 2`. The tunnelling `−J X_l` has the coordinates (1, 0, 1),
odd k_Z and even k_X; the Ising pair term `K Z_a Z_b` has (0, 2, 2) and keeps k_Z; the bias
`Δ Z_l` has (0, 1, 1), odd k_X. [The table of moves](../THE_ONE_SQUARE.md#9-the-moves-and-the-two-questions)
asks every move on the cube two questions, whether it commutes with L (a copy, an invariant
subspace) and whether it is the one-sided action `ρ ↦ ρ·F` of an element F that commutes
with H and anticommutes with every jump, a lit element (a reflection,
`S·L·S⁻¹ = −L† − 2σ` with `σ = Σγ`, the far end of
[F158](../proofs/PROOF_PALINDROME_TWO_END_COUNT.md)). Laid on this folder, the two
hydrogen-bond-qubit models and the wire, three of the table's moves are sighted and no new
one; the moves of the Z side, the dark shift by `Z^⊗N`, the centre line, the quarter-turn
and the turn about Z, are broken by the tunnelling field, and the bond-and-field row and
the plane are what the coordinates above and the decay rate read. Gate
[`simulations/water/cube_on_water.py`](../../simulations/water/cube_on_water.py), 35 checks,
every exact row compared to 0 and every one with a control that must fail.

- **Row 1, the corner shift `ρ ↦ ρ·F`, on the two hydrogen-bond-qubit models.** At `Δ = 0`
  every term of [the N = 2 and N = 4 models](HYDROGEN_BOND_QUBIT.md) has even k_X, so
  `[H, X^⊗N] = 0`, `X^⊗N` is lit, and the shift by it is a reflection: `S·L·S⁻¹ = −L† − 2σ`
  entry for entry, in floats at dyadic and at fifty generic couplings alike, since the
  global flip keeps every summation order. That identity is the pairing those models show.
  The bias `Δ = 1/8` breaks the commutator and the reflection exactly, and the spectrum
  stops pairing.
- **The turn row, on the same models.** At `Δ = 0` the conjugation
  `Ad_{X^⊗N}: ρ ↦ X^⊗N ρ X^⊗N`, which is `(−1)^{k_X}` on strings, commutes with L exactly: a
  copy and never a reflection, since a conjugation negates no dissipator. It is the first
  reading of [the field section](#what-a-field-along-the-wire-reads) at zero field, where the
  reversal `H(E) ↦ H(−E)` becomes a commutant of L. The bias breaks it.
- **The half-turn row, on the wire: a copy without a reflection.** The second reading of
  the field section, `S = Rev·X^⊗N` on a mirror-symmetric wire under a mirror-odd bias
  `δ_l = −δ_{N−1−l}`, at zero field. Rev is the site reversal (not the table's R, which is
  the shift by `X^⊗N`), and `X^⊗N` is the per-site half-turn about the axis X ⊥ Z, so S is
  the table's half-turn composed with a site permutation. `Ad_S` commutes with L, exactly in
  rationals at one generic rational point, while the far kernel `𝒲` is empty, as the argument
  in [what crosses anyway](#what-crosses-anyway-the-f1-palindrome) says for every nonzero bias
  on a wire with every `J_l ≠ 0` and every site dephased at a positive rate and the gate reads
  by rank at the mirror-odd profile, and the spectrum does not pair, so no reflection of any
  kind exists. In floats the copy's residual is a function of the order in which H's terms
  are summed and of nothing else in the physics, the case CLAUDE.md's no-rounding rule calls
  reading rather than gating: summed site by site it is a few eps on H's diagonal, the
  bias and ZZ terms that Rev reverses; summed with each term's Rev image first it is 0.0 at
  every one of two hundred generic draws. The shape is the one the table's half-turn row
  records for the Néel split of the block-spectrum pairing pass
  ([CAUGHT_ERRORS](../CAUGHT_ERRORS.md), the entry on what the pairing presumed): there too
  the chain reflection composed with `X^⊗N`, this same `Rev·X^⊗N`, survives as the copy after
  the element that paired the spectrum is gone; the two Hamiltonians differ, the surviving
  operator is the same. Under the odd bias Rev alone is no copy, it sends the bias to its
  mirror image, so the `X^⊗N` is needed; a mirror-even bias `[δ, 0, δ]` breaks this copy
  while Rev alone remains one; without a bias `X^⊗N` is lit again and the spectrum pairs.

The sweep for this section: [the reflection on the one square](../../reflections/ON_THE_ONE_SQUARE.md)
already says in plain words that the wire adds no row and only a reading of the Z height,
and this section is its gated form; the table of moves holds the three moves and not the
wire; MirrorWorld's `Pair` holds the Z height as the disagreement count, in that word; the
formula registry holds the water pages under F3, F6, F86b and F98 and no cube entry for
them; the glossary holds the cube with its two questions and k_Z as the XY-weight; the
OpenArcs registry holds the polarity cube, the Q-of-water audit, the Grotthuss wire among
the incompleteness survivors, the crossing gate as a bit-exact leaf and a parked water
prose note, no move; the Confirmations registry holds no water row; CAUGHT_ERRORS holds the
Néel entry, the water scripts' vectorisation entry, the pairing-assignment entry on the
proton chain's two ladders and the σ⁻ label entry; the typed layer names water in three
inheritance claims and the incompleteness witness, none a move. Checked for adjacency: the
three moves against the table, where each stood without the wire; the first against
[the hydrogen bond as a qubit](HYDROGEN_BOND_QUBIT.md), which reads the same pairing as
exact; the third against the Néel entry, the same surviving operator, and against
[the h thread](../../experiments/LATTICE_H_THREAD.md), where `X^⊗N` is F131's third sighted
mirror after the site reversal, the two factors of S; the Z height against `Pair`. Nothing
on the wire asks a question the chain had not asked; what the wire adds to the cube is a
physical name for MirrorWorld's disagreement count: the number of bonds in which two proton
configurations disagree.

---

## Scope, stated plainly

This is a model of a **water wire**: single-file, 2-coordinate oxygens, one proton per
O···O linkage. That is the geometry of water in a carbon nanotube, in gramicidin A, in
aquaporin. It is not bulk water, where oxygen coordination is about 3.5, nor ice Ih,
where it is 4. On a branching oxygen the charge is a popcount over the incident bonds
and the 1D dipole argument above does not apply as written.

The sharper limit is that **this state space cannot hold an excess proton at all**. It
has N protons in N bonds in every one of its 2^N configurations, so the wire is neutral
by construction, which is exactly what makes the dipole origin-free. Grotthuss transport
moves an excess H₃O⁺ along a wire, and a wire carrying one is charged. That configuration
is outside the model, not merely unlikely inside it. What the dipole measures here is the
polarization of a neutral wire, not the position of a charge carrier.

For the same reason `−J Σ X_l` should be read as the elementary hop and not as the
Grotthuss mechanism. The current literature picture of that mechanism, as concerted
bursts along a directed wire gated by hydrogen-bond exchanges in the second solvation
shell, is **not verified here against primary sources** and is recorded only to mark
where the model stops; a 2-coordinate 1D wire has no second shell in any case. And
nothing here computes a transport rate. As [the proton water chain](PROTON_WATER_CHAIN.md)
already says of itself: we compute structure, not transport.

---

## What this pass changed elsewhere

Asking what carries F98 in water turned up a premise that was simply the wrong one.
The F98 entry and four docs downstream of it stated the result for "any truly-class
(F87) Hamiltonian", citing F4, which is proven for the dephased **Heisenberg**
Liouvillian and uses `[H, Ŵ] = 0` explicitly. Measured at N = 4, the palindrome class
turns out to do no work at all for this asymptote and `Ŵ`-conservation to do all of it:

- `H = Σ_a X_a X_{a+1}` is truly-class by F87's own criterion and does **not** conserve
  `Ŵ`. It sends the K-intermediate Dicke state to `I/2^N` and gives `α(∞) = 0`.
- The soft DM chain `Σ(XY − YX)`, the non-truly `Σ(XX+YY) + Σ_l h_l Z_l` with random
  `h`, and even `H = 0` all conserve `Ŵ` and all land on `α(∞) = 0.3000000000` exactly.

Two further premises in the same sentence turned out to be unnecessary. Connectedness
is not needed, and the mechanism it was cited for is not the mechanism: a Heisenberg
chain on the disconnected bond set {(0,1), (2,3)} has `dim ker L = 9`, not `N+1 = 5`,
and still gives F98 exactly. Uniform γ is not needed either. What is actually doing the
work is that dephasing removes the Z-basis coherences, `[H, Ŵ] = 0` keeps populations
inside their popcount sectors, and this particular initial state already has a uniform
diagonal within each of its two sectors. That last part is a property of the state, not
of H, so the result does not extend to arbitrary initial states: off a sector-uniform
one, H must also mix within a sector, and `Σ ZZ` does not.

The premise is now stated as `[H, Ŵ] = 0` in [the formula registry](../ANALYTICAL_FORMULAS.md)
with a Valid-for / Breaks-for pair, and in the four downstream docs, one committed
script and one typed claim that carried the old one. Every verified instance was already
inside the corrected premise, so no number moved.

---

## Open

1. **The Zundel parameterization.** The older parameterization remains unresolved: its
   energy scale may be a H₅O₂⁺ shared-proton stretch fundamental rather than a tunnelling
   coupling. [The hydrogen bond as a qubit](HYDROGEN_BOND_QUBIT.md) now withholds Zundel
   `Q` and crossing derivations. A primary-source reading is required before assigning a
   tunnelling coupling or making any replacement numerical claim.
2. **The branching oxygen.** The dipole identity is oriented-1D. Bulk water and ice need
   the charge read as a popcount over incident bonds, which is a different operator.
3. **The charged wire.** Adding an excess proton leaves this state space. Whether the
   framework's grading survives that extension, and what it becomes there, is untouched.
4. **A displaced-charge observable.** An infrared or terahertz response along a confined
   wire couples to `Σ_l Z_l`. The model side is in "What a field along the wire reads":
   the reversal parity, a bias it can hide, and the local readout that shows it. Which
   spectroscopy resolves this on a real wire, with what bath and preparation, is open.
