# Benzene Hückel through the Framework Lens

**Date:** 2026-05-17 evening
**Authors:** Tom + Claude
**Status:** Tier 2 (structural observation, not yet a new theorem); first substantive
carbon-domain reading.
**Script:** [`simulations/carbon/benzene_huckel_framework_lens.py`](../../simulations/carbon/benzene_huckel_framework_lens.py)
**Tested:** benzene C₆, cyclodecapentaene C₁₀, butadiene C₄ chain, hexatriene C₆
chain, cyclobutadiene C₄ ring, cyclopropenyl cation C₃ ring.

---

## Structural siblings on distinct objects

Coulson and Rushbrooke (1940) proved that ALTERNANT hydrocarbons, those whose
carbon framework is bipartite, have molecular orbital spectra that are
**palindromic around α** (the on-site Coulomb integral): for every MO at energy
α + x there is a corresponding MO at α − x. Every chemist learns this in their
first MO theory course.

The R=CΨ² framework's F1 palindrome theorem (proven 2026,
`docs/proofs/MIRROR_SYMMETRY_PROOF.md`) says that the Liouvillian spectrum of a
Heisenberg-family Hamiltonian under Z-dephasing is closed under
λ → −λ − 2·Σγ, palindromic around the centre −Σγ, on any graph and at any N.
F98 supplies a distinct long-time result for its specified KIntermediate state.

**Coulson-Rushbrooke and F1 are structural siblings, not one theorem on two
physical levels.** Each has a Z₂-based pairing implication and a centre, but
Coulson-Rushbrooke acts on a single-particle Hückel matrix whereas F1 acts on a
Liouvillian operator space. Their formal framework partners are also distinct:
Coulson-Rushbrooke's is the bipartite sublattice gauge K, while F1 uses Π.

C-R's K trigger fails on an odd ring; F1's Π is topology-blind within its stated
Hamiltonian-and-channel scope. [Why the centre row is graded lower than the
pair-sum row](#why-the-centre-row-is-graded-lower-than-the-pair-sum-row) makes
that distinction explicit.

| Theorem | Date | Level | Centre | Trigger | Breaks at |
|---------|------|-------|--------|---------|-----------|
| Coulson-Rushbrooke | 1940 | Carbon Level 1 (MO) | α (on-site Coulomb) | bipartite C-framework | non-bipartite ring (e.g. C₃) |
| F1 (R=CΨ²) | 2026 | Qubit Level 0 (Liouvillian) | −Σγ | Π conjugation on truly-class H + Z-deph | F1-Brecher (T1, depolarising, a field violating [F138](../ANALYTICAL_FORMULAS.md)'s clause 2) |

The common structural shape is a Z₂ pairing with a centre. C-R's trigger is the
bipartite graph; F1's is the selected Hamiltonian and Z-dephasing channel. The
following calculations preserve that distinction rather than treating a material
mapping as established.

---

## Numerical observations

### Bipartite cases (palindrome agrees at machine precision)

The Coulson-Rushbrooke statement is exact for the stated bipartite Hückel
matrices. The producer's floating-point eigensolver reproduces its pair sums at
machine precision (pair-sum deviation < 1e-15). Energies in units of β:

| System | Spectrum (sorted) | Pair check |
|--------|-------------------|------------|
| Benzene C₆ ring | −2, −1, −1, +1, +1, +2 | 3 pairs at 0 = 2α ✓ |
| Cyclodecapentaene C₁₀ ring | −2, −1.618, −1.618, −0.618, −0.618, +0.618, +0.618, +1.618, +1.618, +2 | 5 pairs at 0 ✓ |
| Butadiene C₄ chain | −1.618, −0.618, +0.618, +1.618 | 2 pairs at 0 ✓ |
| Hexatriene C₆ chain | −1.802, −1.247, −0.445, +0.445, +1.247, +1.802 | 3 pairs at 0 ✓ |
| Cyclobutadiene C₄ ring | −2, 0, 0, +2 | 2 pairs at 0 ✓ |

The 1.618 / 0.618 pattern in C₁₀ is the golden ratio φ and its reciprocal; both
appear in cos(πk/(N+1)) at N = 10 / 4 etc, the same OBC sine-mode dispersion the
framework uses in `XyJordanWignerModes`.

### Non-bipartite counter-case (palindrome breaks)

**Cyclopropenyl cation C₃⁺** (3-ring, odd cycle, non-bipartite, 2π aromatic per
4n+2 with n=0): spectrum {α + 2β, α − β, α − β}. Pair sum (α + 2β) + (α − β) = 2α + β,
**deviation = β from 2α, palindrome explicitly violated.** The framework counterpart
is K, the bipartite sublattice gauge, which breaks on an odd cycle for the same reason
([K partnership](../proofs/PROOF_K_PARTNERSHIP.md)). It is *not* an F1 Brecher: F1's
palindrome is topology-blind and holds on this ring, which is the scope difference the
[centre section](#why-the-centre-row-is-graded-lower-than-the-pair-sum-row) reads.

The C₃ graph input carries the conventional `4n+2` occupation label while its
non-bipartite Hückel matrix fails the K pairing condition. This selected
single-particle calculation distinguishes the occupation label from the
palindrome condition. It does not establish aromaticity, stability, or any other
material outcome for a molecule; such a translation needs its own specified
degree of freedom, Hamiltonian, bath, preparation, and measurement.

---

## Selected 4n+2 vs 4n occupation labels, not a material mechanism

For the stipulated Hückel C₆ and C₄ graph inputs conventionally labelled benzene
and cyclobutadiene, both spectra are palindromic. Their selected filling lines
occupy different positions relative to the palindrome centre:

- **4n+2 label**: the selected HOMO sits BELOW `α` at an isolated paired energy;
  the LUMO sits ABOVE `α` at the palindrome partner. For the C₆ input with `β < 0`:
  HOMO at α + β = α − |β| (below α), LUMO at α − β = α + |β| (above α),
  gap = 2|β|.
- **4n label**: the selected HOMO and LUMO are degenerate at `α`, the palindrome
  centre. For the C₄ input: HOMO = LUMO = 0β, gap = 0β.

Thus the stated Hückel graphs distinguish a filling line that cuts a paired-MO
gap from one that meets a degenerate pair at the centre. They do not by themselves
establish aromatic or antiaromatic stability, Jahn-Teller behavior, bond distortion,
or any material mechanism.

**Conditional framework-lens question:** after an MO basis, Liouvillian
construction, and channel have been selected, could the two selected filling-line
positions be read through a Klein-4 character constraint? The table's `α` to
`−Σγ` Tier-2 analogue is not a material identification, and no current producer
constructs the required MO-basis Liouvillians. This remains an open model proposal.

---

## Framework-vocabulary translation

| Hückel | Framework | Status |
|--------|-----------|--------|
| α (on-site Coulomb integral, ≈ −11.4 eV) | −Σγ analog (palindrome centre) | Tier 2 structural identification |
| β (resonance integral, ≈ −2.4 eV) | J (framework coupling) | Tier 2 structural identification |
| bipartite carbon framework | K H K = −H, the sublattice gauge (single-particle, Δ = 0) | Tier 1 algebraic match |
| Coulson-Rushbrooke pair sum 2α | F1 pair sum −2Σγ | Tier 1 algebraic match |
| C₃ odd-ring palindrome break | K breaks on an odd cycle (not F1, which holds there) | Tier 2: the same break, gate-verified, the identification not derived |
| what the palindrome needs at all | bipartiteness for C-R and K; truly-class letters for F1's operator identity, with F87's soft class holding the spectrum alone | Tier 2: the same role, two different conditions |
| half-filled p-shell on C₆ ring | F86b KIntermediate Dicke (n ∈ {2, 3}) | Tier 3 candidate |
| 4n+2 vs 4n aromaticity | Klein-4-group HOMO character | Tier 4 candidate (open test) |
| benzene HOMO at −β | F86b Dicke superposition γ = 1/2 anchor | Tier 4 candidate |

The `β → J` row is a Tier-2 structural translation inside a selected
Hückel/framework model. Choosing the site occupation
`n = (I − Z)/2`, the XX+YY Hamiltonian, and a local-density jump gives a
conditional model in which `D[n] = D[Z]/4`. It does not choose a material-carbon
degree of freedom, a material β-to-J convention, a bath channel or rate, or a
carbon `γ`, `T₂`, or `Q`; [the README's model contract](README.md#conditional-c4-and-c6-working-model)
and [Q audit](../Q_BELONGS_TO_NO_SUBSTANCE.md) keep those assignments open.

---

## Why the centre row is graded lower than the pair-sum row

The table above calls the pair sums a Tier 1 algebraic match and the centres only a
Tier 2 structural identification. That reads like an understatement, because the two
centres are the same formula: α = tr(H)/N on the Hückel side, since the adjacency
matrix has a zero diagonal, and −Σγ = tr(L)/4ᴺ on the framework side, from
Tr(L_H) = 0 for any Hermitian H together with Tr(L_D) = −γN·d²
([Absorption Theorem §4.4](../proofs/PROOF_ABSORPTION_THEOREM.md); that section works
in decay rates, where the same centre reads +Σγ).

It is not an understatement. For the molecules in this document the centre is the half
of the palindrome that is free.

The [universal palindrome condition](../../hypotheses/UNIVERSAL_PALINDROME_CONDITION.md)
splits a generator against a conjugation Q into two independent requirements: the
coupling part must be mirror-odd, Q A Q⁻¹ = −A, and the bath part must pair to a
constant, B + Q B Q⁻¹ = −2c·I, which is what puts the centre at −c. For Hückel,
A is the off-diagonal part β·A_graph and B the on-site diagonal; the standard Q is
the sublattice sign flip. In an all-carbon molecule B = α·I and −c = α. Measured against that Q
([`huckel_palindrome_conditions.py`](../../simulations/carbon/huckel_palindrome_conditions.py)):

| System | B-condition (fixes the centre) | A-condition (mirror-odd) | palindrome |
|--------|-------------------------------|--------------------------|------------|
| C₆, C₄, C₁₀ rings, butadiene, hexatriene | 0.0 exactly | 0.0 exactly | holds, ~1 eps |
| C₃, C₅, C₇ rings | 0.0 exactly | 6.788 = 2√2·\|β\| | broken by 4.800 / 2.967 / 2.136 eV |
| pyridine, α_N = α + 0.5β | 2.191 | 0.0 exactly | broken by 0.018 eV |
| push-pull 4-chain | 3.578 | 0.0 exactly | **holds**, ~1 eps |

Three of the rows behave; the fourth is the one to read first, because it says what the
B-residual does *not* mean.

**The odd rings break the A-condition and leave the centre untouched.** Their spectra
still have mean α: cyclopropenyl is {α + 2β, α − β, α − β}, mean α, extreme pair sum
2α + β. The row is not quoting a colouring that does not exist, either. B = α·I is a
multiple of the identity, so Q B Q⁻¹ = B for *every* invertible Q and the centre
condition collapses to 2α·I = −2c·I with no reference to the graph: that is an
identity, not a measurement. The A-residual is the one that depends on Q, so 6.788 is
quoted as the *minimum* over all 2ᴺ sign patterns, which an odd cycle attains at one
frustrated edge (it forces an odd number, and 2β at two symmetric positions gives
2√2·\|β\| = 6.788 eV, ring-size-independent). The palindrome break is not: it falls
off with ring size, 4.800 / 2.967 / 2.136 / 1.667 / 1.366 eV at C₃ to C₁₁, and the
script gates each against the closed form α + 2β·cos(2πk/n) rather than against a
threshold. Those figures are the largest deviation over all pairs, which on an odd ring
is a degenerate level paired with itself, 2\|β\| = 4.800 at C₃; the extreme-pair
deviation quoted in the counter-case section above is the smaller \|β\| = 2.400.

**The push-pull chain says the B-condition is Q-relative.** Give a 4-chain a linear
gradient of on-site energies, α + 0.5β at one end down to α − 0.5β at the other.
Against the sublattice flip the B-residual is 3.578 and the naive reading calls the
palindrome broken. It is exact. The rescuing conjugation is Q = P·K, the chain
reflection composed with the sublattice flip, and against *it* both residuals are 0.0.
What the B-condition asks is that the diagonal be **odd under Q about the centre**;
"constant diagonal" is only what that reduces to when Q happens to be diagonal, as the
sublattice flip is. An earlier draft of this section stated the reduction as the
condition, which is the same error as reading one colouring's A-residual as the
A-condition.

**Pyridine really does break it.** The A-condition is exactly satisfied, the graph is
still bipartite, and the palindrome goes anyway, by 0.018 eV, because α_N = α + h_N·β
lifts one entry in a way no Q can absorb. That last part needs no search: X is
diagonalizable, so a Q with Q X Q⁻¹ = −X − 2c·I exists exactly when
Spec(X) = −2c − Spec(X), and summing that identity pins −c at the spectral mean. The
pair deviation *is* the obstruction, and 0.018 ≠ 0 rules out every invertible Q at
once. So the centre condition is a live condition in general, and it is free
in all-carbon Hückel for one reason only: B = α·I, whatever Q. That is why "the centre
is tr/dim" is, *for hydrocarbons*, a restatement of "the diagonal is α·I and the
adjacency has no self-loops", holding for alternant and non-alternant alike. What
Coulson and Rushbrooke proved is the pairing, and on top of the free centre the pairing
needs the A-condition.

The two breaks are the two ways to lose bipartiteness that
[F103 §7.2](../proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) names on the framework
side, an odd cycle and a lifted diagonal, read here on the molecular site graph rather
than on F103's basis-state graph.

The framework side splits the same way and says so:
[§4.4](../proofs/PROOF_ABSORPTION_THEOREM.md) gets the spectral mean from the traces
alone, with no symmetry needed, while the pairwise statement, that a mode and its
partner average the centre, needs F1. The two coincide exactly when the spectrum is
already palindromic.

### The scopes do not match, and that is the useful part

Take the two breaks to the framework and they land differently, which is what kills the
reading that C-R and F1 are "the same theorem at two physical levels".

The **odd cycle** has no F1 counterpart at all. F1 is topology-blind, as its registry
entry already records ("any graph; any N; non-uniform γ per qubit"): its palindrome
pairs every eigenvalue on the C₃ and C₅ rings, the graphs on which Coulson-Rushbrooke
fails. The two live on different objects, an N-site
single-particle matrix against a 4ᴺ operator space, and the graph enters them
differently: it is C-R's whole trigger and, for F1, only where the bonds sit.

The **lifted diagonal** does have one, and it is not the γ profile. A per-site γ is
a lift in the dissipator and F1 survives it; the counterpart of pyridine's lift is an
on-site field term in H. What such a field must satisfy is
[F138](../ANALYTICAL_FORMULAS.md)'s clause 2, a single common axis orthogonal to every
dephasing axis present. F138 is where that law lives, and it records clause 2 as
spot-checked rather than swept, so it is open work and not a settled result to lean
on.

So Coulson-Rushbrooke's partner in the framework is not F1 but
**[K](../proofs/PROOF_K_PARTNERSHIP.md)**, the bipartite sublattice gauge
K = diag((−1)^ℓ) with K H K = −H, which is the A-condition itself, and which breaks on
odd cycles for the same reason C-R does (K_PARTNERSHIP: "odd-N periodic is
non-bipartite: the wrap-around bond connects same-sublattice sites, so K breaks there
independently of Δ; this is a pure topology effect"). Mind K's own scope before
carrying it further than Hückel needs: K H K = −H holds for the single-particle
hopping problem, and the same proof records that a ZZ term survives the conjugation, so
K's chiral structure breaks at any Δ ≠ 0. Hückel is Δ = 0 and single-particle, which is
why the identification is clean here. K is a Hilbert-space involution on 2ᴺ; F1's Π
acts on the 4ᴺ operator space; they are independent maps on distinct objects. What C-R
and F1 *do* share is the implication,
which Open Question 7 of the universal palindrome condition derives once,
substrate-free, from the two sub-conditions alone. Shared implication, different
ingredients, different scopes.

---

## Open questions (this folder, next sessions)

1. **Does the selected C₄/C₆ XX+YY + local-Z model satisfy F1?**
   [The selected Liouvillian comparison](BENZENE_LIOUVILLIAN_PALINDROME.md)
   has the F1 palindrome for its all-site local-Z jump and gives the listed finite
   residuals for its selected bond jump. If the model selects
   `n_l = (I − Z_l)/2` and the density jump, `D[n_l] = ¼·D[Z_l]`; that identity
   does not classify a material bath.

2. **Where is F98 instantiated?**
   F98 requires its KIntermediate initial state, a magnetization-conserving
   Hamiltonian, and Z-dephasing on every site. [The selected C₄/C₆ ring
   instances](BENZENE_F98_LONG_TIME.md) give α(∞) = 3/10 and 2/7, respectively;
   they do not assign a material carbon bath or Q.

3. **Klein-4-group character of HOMO at palindrome centre as aromaticity criterion?**
   Apply `KleinFourGroupSelfPairedRefinement` to the benzene + cyclobutadiene
   Liouvillians; check the Klein character of the HOMO at α; predict + verify the
   4n+2 vs 4n distinction.

4. **Does the V-Effect 14/19/3 trichotomy show up in benzene's electronic spectrum?**
   The framework's V-Effect splits N=3 Pauli-pair bilinears 14 hard + 19 soft + 3
   truly. For N=6 (benzene): what's the prediction, and does it match benzene's
   known electronic-transition pattern?

5. **The Frost-circle construction (1953) and the framework's mode dispersion.**
   Frost circle is the geometric mnemonic that Hückel ring eigenvalues are projections
   of inscribed N-gon vertices on a circle of radius 2β at centre α. The framework's
   `XyJordanWignerModes` uses the same cos(πk/(N+1)) dispersion at chain BC.
   Frost-cyclic and framework-OBC are two boundary-condition choices on the same
   algebraic structure. [The selected-model clock comparison](FROST_CIRCLE_AS_THE_CLOCK_FACE.md)
   supplies its own XX+YY/Z-dephasing frequencies, lifetimes, and Q* values; it
   does not assign them to carbon material parameters.

---

## What the structural comparison supplies

Coulson-Rushbrooke fixes the single-particle Hückel/K statement, and F1 fixes the
Liouvillian/Π statement. Their shared pairing shape is useful for organizing
selected model translations, but evidence for either statement does not validate
the other or select a material-carbon embedding. The next material-facing step is
therefore to choose the degree of freedom, coupling convention, bath channel, and
bath rate before asking whether a carbon experiment realizes one of these models.

---

## Anchor

- Script: [`simulations/carbon/benzene_huckel_framework_lens.py`](../../simulations/carbon/benzene_huckel_framework_lens.py)
- Parent docs: [`docs/carbon/README.md`](README.md), [`docs/water/README.md`](../water/README.md)
  (sister substrate)
- Framework anchors: [F1 palindrome](../ANALYTICAL_FORMULAS.md#f1), [F86b 3/8 Dicke
  anchor](../ANALYTICAL_FORMULAS.md#f86), [F87 trichotomy](../ANALYTICAL_FORMULAS.md#f87),
  [F98 (N+2)/[4(N+1)] bridge](../ANALYTICAL_FORMULAS.md#f98),
  `compute/RCPsiSquared.Core/BlockSpectrum/KleinFourGroupSelfPairedRefinement.cs`
- Literature: Coulson + Rushbrooke (1940) "Note on the method of molecular orbitals",
  Proc. Camb. Phil. Soc. 36, 193; Frost (1953) "Frost circle"; Hückel (1931).
