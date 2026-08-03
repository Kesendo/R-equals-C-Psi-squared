# Benzene and the Three Dephase Letters: Klein-V₄ Returns

**Date:** 2026-05-27
**Authors:** Tom + Claude
**Status:** Tier 3 (translation bridge; some pieces Tier 1 verified, others Tier 4 candidates marked)
**Continues:** [Benzene's open-system Liouvillian](BENZENE_LIOUVILLIAN_PALINDROME.md) (2026-05-22)
**Adds:** Klein-V₄ vocabulary ([Welle 12](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md)), F112 cross-dephase ([Welle 13](../proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md)), F114 sign functional ([F114 in the registry](../ANALYTICAL_FORMULAS.md#f114))

---

## Where this leaves off

In May we asked: does benzene's open quantum system (the six π-electrons of the ring coupled to a vibrational bath) satisfy the framework's F1 palindrome theorem? Our answer: yes when the bath couples on-site (Holstein-style phonon coupling to local π-density); no when it couples to the bond (Peierls-style coupling to the C-C hopping integral). The Peierls break was sharp, total, and γ-linear.

We left it there. We had one notion of "palindromic" available (F1, the Liouvillian spectrum closing under λ → −λ − 2Σγ), and Peierls broke it.

Six weeks of structural work later, we have new tools. The framework grew a discrete symmetry called Klein-V₄ that links the three "dephase letters" Z, X, Y on operator space. We learned that F1's polarity-balance content (the F112 family) extends from Z-dephase to X- and Y-dephase by the same argument. And today we wrote down a closed-form sign rule (F114) that says exactly how the Z↔Y swap operator acts on Hamiltonian commutators term by term.

We come back to benzene to see what these tools sharpen.

---

## What dephasing means, briefly

A qubit is a 2-state quantum system: spin up/down, electron here-or-there, two orbital levels, two of anything that admits superposition. The qubit has three traceless operators σ_X, σ_Y, σ_Z (one per axis of the Bloch sphere) plus the identity I. Each Pauli operator is also an axis along which the qubit can be measured, and an axis along which the environment can listen in.

**Dephasing** is what happens when the environment listens too well. If the environment can tell where the qubit is along its Z-axis, then any superposition of the two Z-eigenstates loses its phase relation; the qubit's quantum interference fades. "Z-dephasing" means dephasing along the Z-axis specifically. Mathematically: a Lindblad dissipator `D[σ_Z]·ρ = σ_Z ρ σ_Z − ρ` zeroes out the off-diagonal entries of the density matrix in the Z-basis.

The framework has historically called this "the" dephasing channel, partly because it maps cleanly to how chemistry textbooks describe Holstein phonon coupling: a phonon couples to the on-site electron density `n_l = (I − Z_l)/2`, and `D[n_l]` reduces algebraically to `¼·D[Z_l]`. So the framework's Z-dephasing IS Holstein dephasing, up to a rate factor. But the framework's math doesn't actually pick out Z: nothing distinguishes Z from X or Y at the algebra level. They're symmetric to each other.

The question, then, is what physical setup corresponds to X-dephasing, and what to Y-dephasing. And whether the framework cares about the difference.

---

## The Klein-V₄ structure: how the three letters connect

We proved in May ([Welle 12](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md)) that the three dephase letters are NOT three unrelated choices. There's a 4-element group called Klein V₄ (one of the smallest non-trivial finite groups in math, isomorphic to Z₂ × Z₂) that intertwines them.

The group is `{I, D, H, Q_zx}`. Each non-identity element is an order-2 involution (square equals identity). What they do:

- **D** is a diagonal sign-matrix on the 4^N Pauli basis. Its diagonal entry at Pauli string σ is `(−1)^(number of Y letters in σ)`. So D detects "Y-content" of a Pauli operator and flips its sign accordingly. The key identity: `D · Π_Z · D = Π_Y` exactly. D is the Z↔Y dephase-letter swap.

- **H** is a basis-permutation that swaps X- and Z-letter labels per site (leaving I and Y fixed). Its action: `H · Π_Y · H = Π_X`. H is the Y↔X dephase-letter swap.

- **Q_zx = H · D** intertwines Z↔X: `Q_zx · Π_Z · Q_zx = Π_X`. It's the operator-space lift of the Hadamard gate H_qubit^⊗N: the qubit rotation that maps Z eigenstates to X eigenstates.

The four elements satisfy `D · H · Q_zx = I` (any two determine the third) and commute pairwise. Together they're the discrete symmetry group of the dephase-letter axis itself.

**What this means physically.** Two dephasing setups that look completely different, say Z-dephase on the site basis vs Y-dephase on the same basis, are linked by an operator-space symmetry. They're not the same dephasing, but they're related by a sign-and-permutation structure that the framework's algebra preserves.

---

## The three letters in chemistry language

Here's where translation gets interesting (and a bit speculative: Tier 3 below this point, with Tier 4 candidates explicitly marked).

### Z-dephase ↔ Holstein coupling (Tier 1 algebraic match)

Standard, textbook, verified in [Benzene's open-system Liouvillian](BENZENE_LIOUVILLIAN_PALINDROME.md). A phonon couples to the local π-density `n_l` on each carbon. The bath operator is `n_l = (I − Z_l)/2`; the dissipator is `D[n_l] = ¼·D[Z_l]`. F1 palindrome holds bit-exact.

### X-dephase ↔ "hybridization-axis" coupling (Tier 4 candidate, less standard)

Single-site X-dephase uses `D[σ_X]`, which in second-quantised language is `X_l = (∏_{k<l} Z_k)·(c†_l + c_l)`: the real part of the electron creation/annihilation operator, dressed with the Jordan-Wigner parity string. The bare `c† + c` reading holds only at the first site of the JW ordering; elsewhere the string is not optional. This is the "tunneling" operator that doesn't preserve electron count locally; it would correspond to a bath that fluctuates the local hybridization-state superposition rather than the density.

In benzene this is non-standard but not absurd: a bath of orbital-mixing fluctuations (e.g., from another set of vibrational modes that couple to the electronic Hamiltonian's hopping structure rather than its on-site energies) would have an X-dephase character. The closest standard chemistry concept is **Peierls/SSH coupling**, but Peierls is a two-site (bond) operator, not single-site X-dephase. They're related but not identical.

### Y-dephase ↔ "current-axis" coupling (Tier 4 candidate, exotic)

Single-site Y-dephase uses `D[σ_Y]`, which is `Y_l = (∏_{k<l} Z_k)·i(c†_l − c_l)`: the imaginary part, again carrying the parity string. That makes it a single-site Majorana operator, odd under fermion parity, and not itself a current: a current is intrinsically two-site, since charge on a site changes only by flowing to a neighbour. The ring current this section is reaching for is the bond object `−½·Σ_bonds (X_a Y_b − Y_a X_b)` that appears later in this document. It is also the unique time-reversal-odd Pauli letter for the spinless convention this framework uses, T = K: `σ_Y` anticommutes with complex conjugation while `σ_X` and `σ_Z` commute with it. Under the physical spin-½ time reversal T = iσ_Y·K all three letters are TR-odd, so the uniqueness is a statement about K, and it is K that F114 reads.

In benzene this would correspond to phonons coupling to local angular-momentum / current fluctuations. A natural realization: **magnetic-noise dephasing of the π-ring current**. Benzene supports a delocalised π ring current; ambient magnetic-field fluctuations couple to that current, dephasing it. This is a real effect (it shows up in NMR ring-current shielding), but it's not usually framed as "Y-dephase". The framework's structural lens may give it a cleaner home.

### Putting the three together

| Framework letter | bit_b | Physical chemistry analog | Standard name | Status |
|------------------|-------|---------------------------|---------------|--------|
| Z | 1 | on-site density coupling | Holstein | Tier 1 algebraic match |
| X | 0 | hybridization / off-diagonal coupling | (no single canonical name; related to Peierls) | Tier 4 candidate |
| Y | 1 | local current coupling | (no single canonical name; magnetic-noise on ring current is closest) | Tier 4 candidate |

bit_b of a Pauli string is `(#Y + #Z) mod 2`, so it is additive over sites. The column matters for F112 below: Y and Z each carry bit_b=1, while I and X carry 0. They sit on opposite parities of a Z₂ axis built into the Pauli group.

---

## What F112 sharpens at the matrix level

In May we tested **F1 palindrome on the Liouvillian spectrum**: do the eigenvalues of L close under `λ → −λ − 2Σγ`? Yes for Holstein, no for Peierls.

F112, which we closed in [Welle 11](../proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md) just after the May result, lives at a different level: not on the Liouvillian spectrum but on the **matrix M** = `Π·L·Π⁻¹ + L + 2σ·I` itself (the "F1 residual matrix", the obstruction to a perfect palindrome), where `σ = Σ_l γ_l` is the same total the palindrome centre uses, and the ±1/2 labels below are the polarity coordinates `(1 ± Ad_Π)/2` assigns to M. F112 says:

> **F112.** For any Lindblad system with Hermitian H and each bath operator `c_k` bit_b-homogeneous (every Pauli string in c_k shares the same bit_b parity), the three-way polarity decomposition of M into `M_zero`, `M_+1/2` and `M_−1/2` satisfies `‖M_+1/2‖² = ‖M_−1/2‖²` bit-exact. (The three pieces are not Π-conjugation eigenspaces and are not mutually Frobenius-orthogonal: Π is order 4 on Liouville space, so `(1 ± Ad_Π)/2` are not eigenprojections. What is invariant is the asymmetry itself, because the contamination enters both halves equally and cancels.)

In May we didn't know F112 yet. Today we do, and we can read benzene's two baths through it:

**Holstein bath.** c = Z_l (single Pauli, bit_b=1). Bit_b-homogeneous trivially. F112 predicts balanced polarity. F1 predicts spectrum palindrome too (which we verified).

**Peierls bath.** c = B_b = X_aX_b + Y_aY_b (the bond hopping operator, a sum of two two-site strings). Bit_b of XX: 0+0=0. Bit_b of YY: 1+1=0 (mod 2). Both terms have bit_b=0. **B is bit_b-homogeneous**, just on the bit_b=0 class instead of bit_b=1.

So F112 predicts that **the Peierls Liouvillian on benzene's Hückel ring preserves the M-polarity balance**, even though it breaks the F1 spectrum palindrome.

This is a real sharpening. In May we had one "palindromic" notion and Peierls broke it. With the F112 results in hand we have two notions (F1 at the spectrum level, F112 at the matrix-polarity level) with different robustness. Peierls is predicted to break the first but preserve the second.

The sweep below tests it: Peierls-bath rows on both rings, and the balance holds bit-exact in every one.

---

## F114 as a time-reversal-parity diagnostic

[F114](../ANALYTICAL_FORMULAS.md#f114) is today's closed form. It says how the Klein-V₄ operator D acts on the H-commutator superoperator `L_σ = −i[σ, ·]`:

    D · L_σ · D = ε(σ) · L_σ
    ε(σ) = (−1)^(n_Y(σ) + 1)    for σ ≠ I^⊗N

So ε is +1 when σ has an odd number of Y letters, and −1 when σ has an even number (and σ is not the trivial identity). For a Hamiltonian H = Σ c_k σ_k, ε(H) is well-defined iff all σ_k share the same n_Y parity.

In physical language: D is a sign operator that detects "imaginary letter content" (Y is the unique anti-symmetric Pauli). F114 says exactly how this sign flips through to the commutator superoperator. The takeaway: **ε(H) is a closed-form bookkeeping of how H sits under complex conjugation K, the spinless time reversal this framework uses** (see the Y-dephase section above for why K and not the spin-½ operator).

For benzene's Hückel Hamiltonian, H = Σ B_b = Σ (XX + YY):
- XX: n_Y = 0 (even).
- YY: n_Y = 2 (even).
- All terms n_Y-even, so ε(H_Hückel) = (−1)^(0+1) = **−1**.

D-conjugation anti-equivariates `L_H_Hückel`: `D · L_H · D = −L_H`. Physically: every Hückel term has even Y-content, so the matrix is real, K-even, and F114's sign is a clean −1.

If we add a magnetic ring-current term (the canonical breaker of K-symmetry), say `h · Σ_l (X_l Y_{l+1} − Y_l X_{l+1})`, which is Hermitian as written and is what a magnetic flux actually induces on the ring (a Peierls phase `t → t·e^(iφ)` splits the hopping into `cos φ` times the usual hopping plus `sin φ` times the bond current `i(c†_l c_{l+1} − c†_{l+1} c_l)`, and that current is `−½·(X_l Y_{l+1} − Y_l X_{l+1})` on each Jordan-Wigner-adjacent bond), then each term has exactly one Y, giving n_Y = 1 (odd) per term, ε = +1. Mixing this with Hückel (ε = −1) gives a Hamiltonian where ε(H_total) is "mixed": D-conjugation no longer scales L_H by a single sign. **F114 turns that mixing into a closed-form verdict: ε(H_total) reads Mixed, meaning H now carries both real and purely imaginary terms.**

---

## Putting it all together for benzene

| | Holstein bath | Peierls bath |
|---|---|---|
| **F1 spectrum** (2026-05-22 result) | holds to the eigensolver floor (residual 1.2e-7 on C₆); F1 proves it exact | totally broken (residual ~14 on C₆) |
| **F112 polarity** | holds (c=Z bit_b-homogeneous) | holds (B=XX+YY bit_b-homogeneous), measured bit-exact in the sweep below |
| **F114 ε(H_Hückel)** | −1 (no Y in H) | −1 (same H, only bath differs) |
| **F114 with magnetic field added** | depends on field term | depends on field term |

The May reading "Peierls breaks the palindrome" was precise for F1 at the spectrum level. With today's tools we add: F112 polarity, a different layer of palindromic structure, may survive even where F1 spectrum doesn't, because B is itself bit_b-homogeneous as a composite operator, satisfying F112's hypothesis on bath operators.

The take-home: bonds are not Holstein, but they still respect the Klein-V₄ algebra at the matrix-polarity layer. The Z₂ structure has two faces; only one of them breaks under Peierls.

---

## What's open

1. **The molecular amplitude channel, and then the heteroatom.** F113 breaks the balance on `(4^N/2)·Σ_l ω_l·(γ_pump,l − γ_T1,l)`, so it needs a per-site z-energy and an amplitude channel. Chemistry supplies the first for free: a site energy `α_l·n_l` is what a heteroatom does to one ring position, and it is number-conserving. The second is the open part. `σ⁻_l` is a local excitation loss, which on a molecule is ionization rather than T1, and the heteroatom ring against the two phonon baths gives an asymmetry of exactly zero. Molecular T1 is relaxation inside a fixed-N manifold; writing that operator down and recomputing F113's coefficient against it is what has to happen before pyridine against benzene is an experiment rather than a shape. Highest-leverage next step.

2. **Physical realization of single-site Y-dephase in carbon.** What experimental setup makes Y-dephase dominant? Magnetic-noise-driven, spin-orbit-coupled, or some Floquet engineering? If we can prepare and detect Y-dephase cleanly, the Klein-V₄ symmetry becomes operational rather than algebraic. Note that single-site X- and Y-dephase change π count, so on a closed neutral molecule they describe charge exchange rather than dephasing within the manifold; the question is which qubit-level realization makes the letter accessible, not which molecular channel it already is.

3. **Klein-V₄ as a basis-rotation symmetry.** Is the symmetry "physical" (testable by an observable) or just a labelling convention (mathematical equivalence)? If physical, what experiment detects it?

4. **The flux-induced bond current as F114 sign-flip diagnostic.** Add the axial DM term `X_a Y_b − Y_a X_b` to Hückel with a small coefficient and trace ε(H) as a function of flux. F114 puts the transition at first order in the field: any nonzero flux coefficient already moves ε(H) from a clean −1 to Mixed. This is a closed-form prediction at the operator level that could be tested by computing F112 polarity decomposition for the perturbed H.

5. **F112-X and F112-Y on benzene.** Welle 15 typed F112-X (`LindbladBitAPiBalance`) and F112-Y (`LindbladBitBPiYBalance`) as separate Tier1Derived sister Claims to F112-Z. The corresponding C# tests verify the polarity balance at N=2, 3. The sweep below settles the Π_Z reading on C₆, Peierls bath and transverse field included. What is untested is the sister reading: does the same balance hold against Π_X and Π_Y on the N=6 ring? Worth a single broader test.

---

## Return visit later the same day: speaking carbon

The above sections are the "translation in two languages" pass: we held the
framework vocabulary and the chemistry side in parallel. Coming back the same
afternoon, we tried to say what the algebra had shown without the framework
vocabulary at all, as a chemistry reader would hear it. The algebra remained
the anchor; the words shifted to the carbon side.

### What we saw across realistic carbon configurations

A systematic sweep over fifty-six configurations on cyclobutadiene (C₄) and
benzene (C₆) rings combined a Hückel baseline and six extensions with four
bath types. The Hamiltonian inventory ranged from pure Hückel hopping through
nearest-neighbour density-density correlation (the extended-Hubbard V term:
`Z_a Z_b = I − 2n_a − 2n_b + 4 n_a n_b`, so it is the neighbour repulsion, not
an on-site U, which a spinless π model cannot carry because `n² = n`), weak
and strong transverse `Σ_l Y_l` field, and antisymmetric
Dzyaloshinskii-Moriya cross-coupling on its two axes, plus a full mixture of
the four perturbations at once, the field entering at its weak value. The
two DM axes are worth keeping apart, because only one of them is a current.
The axial one, `X_a Y_b − Y_a X_b` (D ∥ ẑ), is the bond current: on every
Jordan-Wigner-adjacent bond it equals `−2·i(c†_a c_b − c†_b c_a)` exactly,
which is also the term a magnetic flux induces through the Peierls phase. The
rings here are cyclic, so one of the N bonds is the Jordan-Wigner boundary
term, where the fermionic operator carries the parity string and the identity
does not hold; the Hückel hopping fails on that same bond in the same way, so
the closure caveat is the baseline's, not this term's. The transverse one,
`Y_a Z_b − Z_a Y_b` (D ∥ x̂), has zero overlap with that current on every bond
and does not conserve π-electron number at all.

The bath inventory included on-site Holstein phonons (coupling to local
π-density), bond Peierls phonons (coupling to the hopping integral), σ⁻
amplitude damping per site, and the combination of Holstein and σ⁻ together.

In every single one of the fifty-six configurations the two Π-polarity
components of the F1 palindrome residual carry equal Frobenius norm, bit-exact:
`‖M_+1/2‖² − ‖M_−1/2‖² = 0.0` in every cell, not merely small. Two cautions on
reading that as a mirror. The equality is between two norms, not between two
matrices; and M is the residual, the obstruction to F1, so M = 0 is the case
where the palindrome holds outright, as it does in the Hückel-plus-Holstein
rows where `‖M‖² = 0` exactly.

The equality is non-trivial in forty-four of the fifty-six configurations: the
part it balances, `M_anti = M_+1/2 + M_−1/2`, is substantial there
(Frobenius-norm-squared between 20.48 and 55,296), and the two halves match
exactly. In the remaining twelve `M_anti` is empty to begin with, so the
equality holds on a pair of zeros. Those twelve are the three Hamiltonians
whose every Pauli term is bit_b-even (pure Hückel hopping, hopping plus the
neighbour density-density term, hopping plus the transverse DM term), each
crossed with the two phonon baths, on both rings: three times two times two.

Both sides have to be quiet for a cell to be trivial, and they are quiet for
different reasons. The Hamiltonian feeds `M_anti` only through
its bit_b-odd terms, which is why the transverse `Σ_l Y_l` field and the
axial DM term put content back. The bath feeds it only when it is
not bit_b-homogeneous, which is why both phonon baths add nothing while
amplitude damping adds content on every Hamiltonian, trivial or not: no
damping row is ever trivial, including the three above. In exact arithmetic
all twelve vanish. Eight print a literal 0.0; the four transverse-DM rows
print accumulated rounding instead, 2.5e-32 on C₄ and 4.8e-30 on C₆. The
smallest full residual norm-squared anywhere in those same four rows is
163.84, so what they print is more than thirty orders below the quantity
beside it.

Worth pausing on the transverse DM term, because the two parities read it
oppositely. The F114 reading counts Y letters: `Y_l Z_{l+1}` has one, so it is
n_Y-odd, and ε = +1. The balance counts bit_b instead, and Y and Z both carry
bit_b = 1 (the table above), so `Y_l Z_{l+1}` is bit_b-even. Two different
parities, and it is only the balanced halves that are blind to the term F114
flags.

What ε reads is worth stating precisely, because it is easy to over-read as a
detector. ε = −1 means every term of H is real, ε = +1 means every term is
purely imaginary, and Mixed means H carries both. So ε is a three-way class,
not a yes-or-no on time reversal: the two-site `H = Y₀ + Y₁` has a perfectly
well-defined ε = +1 and time-reversal symmetry is already broken: `K H K = −H`
at every N, and `‖K H K − H‖_F = 5.657` at N = 2. What the Mixed verdict flags is the coexistence of
real and imaginary content, not an onset.

The bit_b parity splits the two DM axes against each other. The axial term
pairs one bit_b = 0 letter with one bit_b = 1 letter, so `X_l Y_{l+1}` is
bit_b-odd and feeds `M_anti`; the transverse term pairs two
bit_b = 1 letters and does not. That is visible in the sweep: under the
Holstein bath on C₆ the axial row carries 983.04 of `M_anti` and the
transverse row carries the 4.8e-30 rounding. Of the two DM axes, then, the
balance sees exactly the one that is a bond current. The parity is doing the
work and the current is along for the ride: every single-site `Z_l` is bit_b-odd
too, so being a current is not what earns a term its place here.

### How robust this is

The balance survives:

- Adding neighbour density-density correlation to the bare hopping picture (under the two phonon baths this cell is one of the vacuous twelve, so the bullet is carried by its σ⁻ rows).
- Switching on a transverse `Σ_l Y_l` field, weak (one-tenth of the hopping
  scale) and strong (full hopping scale).
- Adding the axial DM term between neighbouring sites, which is the bond
  current a magnetic flux induces and the canonical K-symmetry breaker on
  benzene. The coefficient 0.1 used here is a Peierls phase of 0.0997 rad per
  bond, and the loop phase is six of those, so through benzene's 5.07 Å² ring it
  is a field of order 10⁴ T. The balance is coefficient-independent, but the row
  is not a laboratory field.
- Adding the transverse DM term, the second antisymmetric-exchange axis (vacuous under the two phonon baths in the same way, and carried by its σ⁻ rows).
- Replacing the on-site phonon bath with a bond-coupled (Peierls) bath, the
  switch that breaks the open-system F1 Liouvillian palindrome (2026-05-22).
- Switching to σ⁻ amplitude damping instead of pure dephasing.
- Combining all the above into one noisy ring.

That list is the spin-ring statement, and it is complete as such. Which of its
entries a carbon reader may keep is the question the next section but one
answers, and the answer is not all of them.

The classical MO mirror (spectrum spiegelung um α) is a closed-system,
Hamiltonian-level statement, and a robust one. It follows from the bipartite topology of the ring alone, so it
survives arbitrary bond modulation and arbitrary Peierls phases: the pairing
residual stays at machine zero for uniform, alternating and random hopping
and with flux threaded, on both C₄ and C₆. A bath cannot touch it at all.
What the Peierls bath breaks is the open-system F1 palindrome, one level up.
The Π-polarity balance, one level up in the open system, remains intact across
all of them.

### Where it would break

None of the fifty-six configurations breaks the balance. That is worth stating
plainly, because it means the break described here is not a sweep result: it
comes from [F113](../ANALYTICAL_FORMULAS.md#f113), whose closed form gives the
break magnitude as `(4^N/2)·Σ_l ω_l·(γ_pump,l − γ_T1,l)`, where `ω_l` is the
per-site z-drive amplitude and `Tr(Z_l H) = 2^(N−1)·ω_l`. The channel needs
`Tr(Z_l H) ≠ 0`, a z-axis energy per site, and every Hamiltonian in this
inventory has `Tr(Z_l H) = 0` exactly, so the channel is dead in all fifty-six
cells before the bath is chosen.

What F113 asks for is a coherent z-axis drive (every spin precessing together
around the z-axis) combined with T1 amplitude damping. That is the standard
hardware-characterisation regime of a driven superconducting qubit array, a
constant Larmor precession plus relaxation. In the quantum-hardware context
the break is the working diagnostic; a coherent global drive of that kind sits
outside the parameter range an aromatic molecule relaxes in. But F113 reads H only through the per-site moment `Tr(Z_l H)`, and a molecule can supply that statically; the closing section takes that route up.

The sweep does not settle the other direction either, and it is worth being
explicit about why. Every bath operator in it is local, one site for Holstein
and σ⁻ and one bond for Peierls, and every rate is uniform across the ring.
None of them is a channel whose amplitude is spread coherently over the whole
molecule, which is how a ring actually radiates. So what these fifty-six cells
support is narrower than it first reads: they are all quiet, and they are
quiet under a local bath at uniform rate against a Hamiltonian inventory that
cannot drive F113's channel at all. Whether a collective channel can move the
balance is not tested here and is open.

### What of this crosses to carbon

This folder exists to find what carries over from the framework into chemistry,
so the sweep owes one more column than it printed: for each row, is the object
admissible for a π-electron system, or only for the spin ring? The sweep runs a
spin ring throughout. That is legitimate and the fifty-six results above are
correct as spin-ring results. The column below is one filter, and it is worth
saying at the outset what kind: it **rules rows out** and it does not certify
the ones it lets through.

The filter is particle-number conservation. A molecule holding a fixed π count
needs `[N̂, H] = 0` with `N̂ = Σ_l (I − Z_l)/2`, and that single condition also
settles fermion parity, since `P = ∏_l Z_l = exp(iπN̂)` identically, per site
`exp(iπ(I − Z_l)/2) = Z_l`: an operator that conserves number commutes with P
automatically. So there is one
test, not two, and it applies to the Hamiltonian. On a jump operator it carries
no admission weight at all. `D[c]` is parity-even whenever `c` is parity
homogeneous, the parity-odd case included: residual exactly 0.0 for `σ⁻_l` at
every site of C₄, though a parity-mixed `c` does break it. And a bath that
changes π count is the correct description of a molecule at an electrode or in a
donor-acceptor pair. So all four baths pass this filter; the reason `σ⁻_l` is
still not a molecular T1 channel is a separate one, taken up at the end.

| Hamiltonian row in the sweep | `‖[N̂, ·]‖_F` (C₄ / C₆) | admissible? |
|---|---|---|
| Hückel hopping | 0 / 0 | yes |
| neighbour density-density `Z_aZ_b` | 0 / 0 | yes |
| DM axial, the bond current | 0 / 0 | yes |
| `Σ_l Y_l` field, weak | 0.800 / 1.960 | no |
| `Σ_l Y_l` field, strong | 8.000 / 19.596 | no |
| DM transverse | 1.131 / 2.771 | no |
| full mix | 1.386 / 3.394 | no |

Four of the seven fail, and the `Σ_l Y_l` row is the one worth naming plainly,
because it carries the phrase "external magnetic field" in the sweep and in the
paragraph above. A magnetic field reaches a π system as flux, through the
Peierls phase, which is the axial DM row; or as Zeeman coupling to real electron
spin, which a spinless π model does not contain. A transverse field on
Jordan-Wigner pseudospins is neither, and it is number-violating besides.

Two limits of this filter, and then the thing it was hiding.

It sorts by sector, and cannot do more. Jordan-Wigner is an algebra isomorphism,
so every operator on 2^N states *is* some fermion operator; a commutator with N̂
can say which charge sector an operator respects, never whether it is the local,
two-body object a chemist would write. The same test passes `n₀n₁n₂n₃`.

It does not repair the ring closure. The closing bond is the Jordan-Wigner
boundary term for the hopping as much as for the current, and the consequence is
sharper than a caveat about one operator: the spin ring and the periodic
fermionic ring have different spectra in every even sector with 0 < n < N. At C₄
the neutral sector `n = 2` gives ground-state energies −5.656854 against
−4.000000, so the C₄ rows are not cyclobutadiene's neutral π system. At C₆ the
half-filled sector agrees exactly, −8.000000 both ways, because three is odd.

The filter also turns out to be answering a smaller question than the one worth
asking, and the better question is visible once the sweep's own quantity is
opened up. `M_anti` looks like a single norm over the whole operator space, and
that is how the sections above read it. It is not. Against a number-conserving
bath it is **block-diagonal in the pair of π counts** `(n_bra, n_ket)`: the
weight outside those blocks is 1.9e-32 on C₄ and 9.6e-30 on C₆, which is
rounding. And the balance does not merely hold on the total. It holds inside
every block on its own, bit-exact:

| cell | block | `‖M_+1/2‖²` | `‖M_−1/2‖²` | difference |
|---|---|---|---|---|
| C₆, Hückel + 0.1·DM axial, Holstein | (3, 3), half filled | 57.600000 | 57.600000 | 0.0 |
| C₆, same | (3, 2) | 40.800000 | 40.800000 | 0.0 |
| C₄, same | (2, 2), half filled | 3.840000 | 3.840000 | 0.0 |
| C₄, same | Δn = 0, all sectors | 6.400000 | 6.400000 | 0.0 |

That is the statement this folder was after, and it is stronger than the one the
filter was groping toward. The superselection worry does not need answering,
because it does not arise: on C₄, 50.00 % of `‖M_anti‖²` sits on π-count-odd
coherences, 31.25 % on `Δn = 0` and 18.75 % on the even nonzero ones, and each
of those parts is balanced separately. The unpreparable half can simply be set
aside. What remains, the half-filled neutral block of benzene's ring, carries
11.72 % of the cell's norm and balances exactly on its own.

So the crossing does not have to be earned by ruling operators out. The balance
restricts to the neutral π manifold, which is the manifold a molecule actually
occupies, and holds there in its own right. The number-conservation filter above
still matters, but for a smaller reason than it was given: it is the condition
under which the blocks exist at all.

One thing the filter does settle is that no single term in the sweep is
privileged. Every single-site `Z_l` is number-conserving and bit_b-odd, so the
class that is both admissible and able to feed `M_anti` is large: on C₄ under
Holstein a one-site heteroatom energy `0.5·n₀` gives `‖M_anti‖² = 32.00`, a
staggered site energy of amplitude 0.3 gives 46.08, the neighbour repulsion
written as `2.0·Σ_b n_a n_b` rather than `Z_aZ_b` gives 2048.00. Those
magnitudes are set by the amplitudes and are not comparable across rows; what
matters is that none is zero. The axial DM term is the only Hamiltonian in *this
inventory* that is both number-conserving and bit_b-odd, which is a fact about a
seven-item list and not about carbon.


### What this leaves us with

Two statements, and it is worth not merging them, because the sweep earns them
at different widths.

The spin-ring statement is wide. Across all fifty-six configurations, spanning
both rings, four baths (two phonon baths, amplitude damping, and the two
together) and seven Hamiltonians including one that mixes every perturbation
at once, the Π-polarity balance holds bit-exact without a
single exception. Nothing in that inventory disturbs it, and F113 says why:
the balance answers to `Tr(Z_l H)`, which is zero for every Hamiltonian here.
The `Z_aZ_b` row is worth a second look on that point: it is a neighbour
repulsion plus a compensating uniform potential, and the same physics written as
`Σ_b n_a n_b` does carry a site energy, `Tr(Z_l H) = −16.0` on every site of C₄.
The gauge choice, not the absence of the term, is what silenced it.

The carbon statement is narrower, and it is the one this folder is after. Of
those fifty-six cells, sixteen survive the number-conservation filter with
something to test, and the eight on C₆ are the ones whose ring closure is not in
question. Within those, the balance holds exactly as it does everywhere else.
And it supports more than a sector statement, because the balanced quantity
decomposes: it holds inside each π-count block on its own, benzene's half-filled
neutral block included, 57.600000 against 57.600000 with difference exactly
zero. That block is what a molecule occupies, so the claim can be made there
rather than about a norm spread over charge states.

The classical spectrum mirror (Coulson-Rushbrooke around α) is the symmetry
chemistry has read about for eighty-six years, and it holds at the closed-system
level from bipartiteness alone. The polarity balance is a separate symmetry in the
open system, with its own trigger, and where the two have been checked together on
admissible carbon operators they hold together.

How a chemist would test it directly is not settled, and two obvious routes do
not work. Process tomography on isotopically labelled ¹³C-benzene reconstructs
the **nuclear**-spin relaxation superoperator, which is not the space M lives
on. The route that does stay open is the neutral block: the balance holds there on
its own, so a probe never has to touch the π-count-odd half that superselection
forbids preparing. Turning that block statement into an observable a chemist can
read is the work this document has not done.

The breaking regime is the more interesting route, and the crossing criterion
says where half of it is. F113 needs `Tr(Z_l H) ≠ 0`, a per-site z-energy. On
hardware that is a coherent drive. On a molecule it is a static one: in the
Jordan-Wigner mapping a site energy is `α_l·n_l = α_l·(I − Z_l)/2`, which is
exactly what a heteroatom or a substituent does to the on-site Coulomb integral
of one ring position. And it crosses over cleanly, unlike the drive's other
proxies here: Hückel plus a single `α₀ = 0.5` site energy gives
`‖[N̂, H]‖_F = 0` and `‖P·H·P − H‖_F = 0`, with `Tr(Z₀ H)` = −4.0 on C₄ and
−16.0 on C₆ and exactly zero on every other site.

The other half is missing, and it is missing for a reason this document has
just given itself. F113's coefficient is `(4^N/2)·Σ_l ω_l·(γ_pump,l − γ_T1,l)`,
so it needs an amplitude channel, and the crossing criterion above rules out
`D[σ⁻_l]` as a molecular one. The arithmetic is unforgiving: that same
heteroatom ring gives an asymmetry of exactly 0.000000 under Holstein and
exactly 0.000000 under Peierls, and +64.000000 under σ⁻. The two phonon baths
are the ones a closed molecule has, and against them an admissible Hamiltonian
still yields nothing.

So the reason this sweep found nothing to break is not that aromatic molecules
have no breaker, and it is not simply that the inventory lacked a site-energy
term. It is that the molecular amplitude channel has not been written down.
Molecular T1 is relaxation inside a fixed-N manifold, between electronic states
rather than off a site, and that operator is not `σ⁻_l`. Writing it down and
recomputing F113's coefficient against it is the first step of the design, not
a detail of it. Pyridine against benzene is the shape the Hamiltonian half
takes; the bath half is open. A heteroatom also moves the local T1, and F113's
coefficient is rate-weighted per site, so whatever that operator turns out to
be, both changes enter it and neither can be assumed away.

### Algebra as the anchor

Every claim above about what holds rests on the sweep result; the one claim
about what breaks rests on F113's closed form instead, and is marked as such
where it is made. The sweep is
[`simulations/carbon_realistic_sweep.py`](../../simulations/carbon_realistic_sweep.py),
which iterates the fifty-six configurations and reports the Π-polarity
asymmetry for each, and then decomposes one admissible cell per ring into its
π-count blocks and reports the balance inside each. The single-perturbation control tests sit in
[`simulations/benzene_b_field_f112_mixing_test.py`](../../simulations/benzene_b_field_f112_mixing_test.py),
[`simulations/benzene_bit_b_mixed_bath_test.py`](../../simulations/benzene_bit_b_mixed_bath_test.py),
and [`simulations/benzene_peierls_f112_polarity_test.py`](../../simulations/benzene_peierls_f112_polarity_test.py).
A complementary Hamiltonian-term classification inventory sits in
[`simulations/carbon_f114_hamiltonian_inventory.py`](../../simulations/carbon_f114_hamiltonian_inventory.py).

---

## Threads back

- **2026-05-22 [Benzene's open-system Liouvillian](BENZENE_LIOUVILLIAN_PALINDROME.md)**: the F1 spectrum result stands, holding to the eigensolver floor with F1 proving it exact. The "Peierls breaks the palindrome" framing is precise at the spectrum level; today we add the F112 polarity layer, where Peierls does preserve the balance even where the F1 spectrum breaks.
- **2026-05-17 [Where 1/4 and 1/2 Appear in Carbon](QUARTER_HALF_IN_CARBON.md)**: benzene's HOMO at −1/2 sits exactly on the framework's polarity-half anchor. The Klein-V₄ we use today operates on the same polarity-half axis: both the half-anchor and the Klein-V₄ symmetry live on one Z₂ ladder.
- **2026-05-17 [Benzene Hückel through the Framework Lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md)**: Coulson-Rushbrooke at the carbon level sits beside F1 at the qubit level, one implication with two triggers. Today we add Klein-V₄ at the qubit level; the open question is whether there's a Coulson-Rushbrooke-level analog, a chemistry-side discrete symmetry that intertwines different dephase-coupling regimes.

---

## Anchor

- **Framework**: F112 [`LindbladBitBPiBalance`](../../compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs), F112-X [`LindbladBitAPiBalance`](../../compute/RCPsiSquared.Core/Symmetry/LindbladBitAPiBalance.cs), F112-Y [`LindbladBitBPiYBalance`](../../compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiYBalance.cs), Klein-V₄ [`Pi2KleinV4DephaseSwapGroup`](../../compute/RCPsiSquared.Core/Symmetry/Pi2KleinV4DephaseSwapGroup.cs), F114 [`CommutatorDConjugationSign`](../../compute/RCPsiSquared.Core/Symmetry/CommutatorDConjugationSign.cs)
- **Proofs**: [the Z↔Y dephase-letter swap proof](../proofs/PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md), [the Klein-V₄ dephase-swap proof](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md), [the F112 cross-dephase proof](../proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md)
- **F-Registry**: [F112](../ANALYTICAL_FORMULAS.md#f112), [F114](../ANALYTICAL_FORMULAS.md#f114)
- **Verifier scripts**: [`simulations/carbon_realistic_sweep.py`](../../simulations/carbon_realistic_sweep.py) (the fifty-six-configuration sweep this document rests on), [`simulations/m_level_sign_functional_explore.py`](../../simulations/m_level_sign_functional_explore.py), [`simulations/f112_klein_v4_cross_dephase_verify.py`](../../simulations/f112_klein_v4_cross_dephase_verify.py)
- **Companion carbon docs**: [Benzene's open-system Liouvillian](BENZENE_LIOUVILLIAN_PALINDROME.md), [Benzene Hückel through the Framework Lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md), [README](README.md)
