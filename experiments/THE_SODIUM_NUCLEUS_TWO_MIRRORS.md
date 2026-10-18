# The sodium nucleus: two mirrors in the quadrupolar relaxation of a spin 3/2

**Date:** 2026-10-17
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 5.5)
**Verifier:** [`sodium_two_mirrors.py`](../simulations/sodium_two_mirrors.py) (`python simulations/sodium_two_mirrors.py`; exact, sympy 1.14, every parameter a symbol where it can be; 105 checks)
**Rates, not ours:** A. M. Souza, A. Gavini-Viana, I. S. Oliveira, R. S. Sarthour, R. Auccaise, J. Teles, E. R. deAzevedo, T. J. Bonagamba, *Nuclear spin 3/2 electric quadrupole relaxation as a quantum computation*, [arXiv:1105.2539](https://arxiv.org/abs/1105.2539), Eqs. (15) to (27), which cite R. Auccaise et al., J. Magn. Reson. 192, 17 (2008); the spectral densities follow J. R. C. van der Maarel, J. Chem. Phys. 94, 4765 (1991), and G. Jaccard, S. Wimperis, G. Bodenhausen, J. Chem. Phys. 85, 6282 (1986).
**A word fence.**
- *Order* p = m − m′ is the coherence order of |m⟩⟨m′|.
- J₀, J₁, J₂ are the spectral densities of the electric field gradient at 0, ω₀ and 2ω₀, as they enter the secular Redfield generator.
- *Dark* means untouched by the J₀ dissipator: an operator |m⟩⟨m′| that commutes with T₂₀. This is the operator sense of DARK in MirrorWorld's `EndCount` (commuting with every jump), restricted to the one jump T₂₀. It is not a state annihilated by a jump, which is the sense of [Singlet Under a Shared Field](SINGLET_UNDER_A_SHARED_FIELD.md). *Bright* is the rest.
- A *mirror about −c* is the map λ ↦ −λ̄ − 2c on the eigenvalues of one order block, and c is called its centre. The spectrum of the order p block is symmetric about Im λ = −p·ω, not about 0, so the repository's F1 form λ ↦ −λ − 2σ would send it to the block of order −p. The order −p block has the complex-conjugate characteristic polynomial of the order p block (checked), so the map used here is F1's map followed by that conjugation, read inside one block.

**Conventions.** Spin operators I_z, I_± for I = 3/2; T₂₀ = (3I_z² − I(I + 1))/√6, T₂,±₁ = ∓(I_zI_± + I_±I_z)/2, T₂,±₂ = I_±²/2; Lindblad dissipator Σ_m J_|m|·(T₂ₘρT₂ₘ† − ½{T₂ₘ†T₂ₘ, ρ}); Hamiltonian ω_Q·T₂₀ + ω·I_z, so ω_Q is the coefficient of T₂₀ in this normalization, not the conventional quadrupolar splitting. With these conventions Souza et al.'s rates are reproduced with their constant set to C = 3; that is a normalization of the J here, not a value of their physical constant. Dynamic frequency shifts (the imaginary parts of the spectral densities) and the finite polarization of the equilibrium state are not in the model.

## What this is about

Sodium in tissue is a spin 3/2: four levels, not two. It relaxes through the fluctuating electric field gradient around it, not through a magnetic field, and sodium MRI uses that relaxation to tell tissue compartments apart, most directly through the triple-quantum coherence between m = +3/2 and m = −3/2. The repository's mirror results live at d = 2, and for d > 2 they hold only partially ([F121](../docs/ANALYTICAL_FORMULAS.md), under equidistant dephasing); [Qubit Necessity](../docs/QUBIT_NECESSITY.md) §8a leaves the non-dephasing dissipators at d > 2 open. This page reads one physical one, the quadrupolar dissipator of a single spin 3/2.

When all three spectral densities are present it has no single mirror. It has two. The slow part of the gradient, J₀, sees one property of a level: outer (|m| = 3/2) or inner (|m| = 1/2). Everything that agrees on that property, the populations, the central transition and the triple quantum, is dark to it, and this dark half of the spectrum is palindromic about −3(J₁ + J₂), with the central transition and the triple quantum sitting exactly at its centre. The bright half, the satellite and double-quantum coherences, is palindromic about a second centre 3J₀ further out. Both mirrors are exact, for every static quadrupolar splitting and every offset. The dark centre is the mean of the two longitudinal rates 6J₁ and 6J₂, and it is also the triple-quantum rate; that identity was already a line away from the published rates, and the two mirrors give it a place.

## Abstract

One spin 3/2 under the secular Redfield generator of quadrupolar relaxation at infinite temperature (jumps T₂,ₘ, m = −2 … 2, at rates J₀, J₁, J₁, J₂, J₂), with a static splitting ω_Q·T₂₀ and an offset ω·I_z. The generator is block diagonal in the coherence order. Exactly, with every parameter a symbol: the populations decay at 0, 6J₁, 6J₂, 6(J₁ + J₂); the central transition and the triple quantum at 3(J₁ + J₂), shifted only by their frequencies; at ω_Q = 0 the satellites at 3(J₀ + J₁) and 3(J₀ + J₁) + 6J₂, the double quantum at 3(J₀ + J₂) and 3(J₀ + J₂) + 6J₁ (Souza et al.'s rates with their constant set to 3); a splitting draws each pair together (after the frequency shift p·ω), and past 6ω_Q² = 9J₂² (satellites) or 6ω_Q² = 9J₁² (double quantum) turns it into a pair on the same real part, conjugate after the frequency shift p·ω. The dark set pairs under λ ↦ −λ̄ − 2c₁ with c₁ = 3(J₁ + J₂), the bright set under the same map with c₂ = 3(J₀ + J₁ + J₂) = c₁ + 3J₀; with ω_Q and ω present, the satellite and double-quantum polynomials keep real coefficients after the frequency shift p·ω and the root sum −2c₂. The two centres coincide exactly when J₀ = 0. The whole spectrum has one centre exactly when J₀ = 0, or when ω_Q = ω = 0 and one of J₁, J₂ vanishes (proved below; for non-negative densities a second, independent route gives it again). The reason is that at I = 3/2, and at no other spin, the five Hermitian components of T₂,ₘ are Dirac matrices: the dissipator is Pauli-string dephasing, and one map, X ↦ Γ₀·Xᵀ, carries both mirrors. At every spin I ≥ 1 other than 3/2, J₀ alone breaks the mirror (proved); without a static splitting, at every offset, no rates other than all zero give a global centre at any spin from I = 1 to 6 other than 3/2 (proved spin by spin). A static splitting, on resonance, changes that at spin 1 only: there it creates a centre on two tuned rays with non-negative rates, both with J₀ = 0 (proved by elimination; the two rays found again by hand).

## The sweep: what the repository held

Swept 2026-10-17 by one survey agent, stores by name.

- `docs/ANALYTICAL_FORMULAS.md`: [F121](../docs/ANALYTICAL_FORMULAS.md), the qudit partial palindrome under full-Cartan (equidistant) dephasing; at d = 4, N = 1 it counts 8 of 16 modes paired about one centre. No entry on quadrupolar relaxation, spectral densities J₀/J₁/J₂, Redfield, ²³Na, triple quantum or central transition.
- `docs/proofs/`: [PROOF_QUDIT_PARTIAL_PALINDROME](../docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md), the equidistant ladder, the ceiling, the product cap and its clock-shift mirror, checked at d = 4, N = 1, 2. Nothing on spin 3/2 as an angular momentum or on non-equidistant levels.
- `docs/`: [Qubit Necessity](../docs/QUBIT_NECESSITY.md), the trunk d² − 2d = 0, and §8a: 0 of 236 sampled qutrit dissipators paired fully, non-dephasing dissipators at d > 2 left open.
- The typed layer: `QuditPartialPalindromeCeiling`, `QuditProductMirrorCap`, `QuditMirrorProtectionScalingClaim` (Tier 1) and the live `QuditPartialPalindromeWitness` (`inspect --root qudit`), all on the equidistant model. MirrorWorld: `Formulas.F121_*`, d = 2 and 3 only.
- `docs/GLOSSARY.md`: nothing on spin > 1/2 physics; its table of the spent senses of *dark* (a long-lived mode, the record class, the lens letters, an unlit site) does not hold the operator sense used here, which the word fence names.
- `experiments/`, the OpenArcs registry, `docs/CAUGHT_ERRORS.md`, `fw.Confirmations`: nothing on spin > 1/2 physics or quadrupolar relaxation; CAUGHT_ERRORS holds the qudit counting corrections (the product cap holds for d ≤ 5, so at d = 4).

- The Clifford side, searched for Dirac matrices, SO(5), spin 3/2, ququarts, anticommuting strings, depolarizing and colourings:
  - [Qubit Necessity](../docs/QUBIT_NECESSITY.md)'s ququart row (4 immune operators against 12) assumes a dephasing with four distinct eigenvalues. T₂₀ has two, each twice, so the count here is 8 against 8, the qubit's balance on one bit of two; the row and this page agree.
  - [Why Depolarizing Noise Breaks the Palindrome](DEPOLARIZING_PALINDROME.md) §3: under dephasing along X, Y and Z on one qubit the pairing exists exactly when one of γ_X, γ_Y, γ_Z vanishes. That is the shape of this page's theorem at ω_Q = ω = 0, with three channels Γ₀, {Γ₁, Γ₂}, {Γ₃, Γ₄} in place of X, Y, Z. F158's registry entry breaks the depolarizing case by the product argument used below.
  - [The palindrome as a colouring](THE_PALINDROME_AS_A_COLOURING.md) defines colourings and anticommuting sums, the grade-one vectors of the Clifford algebra the strings generate, which the Γ's are.
  - [F158](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md) and [F118](../docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md) apply as they stand once the Γ's are read as strings.
  - MirrorWorld's `EndCount` and `PauliString` hold nothing on this nucleus, and nothing else in the repository holds Dirac matrices or SO(5).

- The moment route used beyond spin 3/2, searched for odd power sums, resultants and subresultants: the power-sum lemma of [PROOF_F87_WINDOWED_MONOMIAL_CONVERSE](../docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) §1 (a spectrum is symmetric about 0 exactly when its odd power sums vanish); F117's and F120's trace faces in `docs/ANALYTICAL_FORMULAS.md`; MirrorWorld's `Hardness.OddPowerSums` and the Sextant's OBSTRUCTED verdict, an odd centred trace that does not vanish; and [The Mirror under a Shared Field](THE_MIRROR_UNDER_A_SHARED_FIELD.md), whose cross-check reads the third and fifth moments through a resultant and its first subresultant ([`shared_field_odd_moments.py`](../simulations/shared_field_odd_moments.py)). The section below runs that route on this generator; nothing in the repository held it for a spin I. For the static splitting at spins other than 3/2 no store held the question at all; the nearest is F121's remark that for its d = 3 example the Hamiltonian helps at no centre, a property of that H and not a theorem. The two-way reading of the power-sum lemma (Newton's identities) is the lemma's own, and it is what makes the moments decisive at spin 1.

Checked for adjacency: F121 and this page are about the same Hilbert space and different dissipators; F121's single-centre count (8 of 16) and this page's two-centre pairing (16 of 16) do not contradict each other. F158's end count is a second route to this page's global-centre theorem, and F118's palindromizer is the shape of its one map. [The Mirror under a Shared Field](THE_MIRROR_UNDER_A_SHARED_FIELD.md) found the same shape for qubit chains under correlated Z noise: the global mirror breaks, mirrors survive block by block about centres that differ by what the watching charges. [The Shared Field in a Crystal](THE_SHARED_FIELD_IN_A_CRYSTAL.md) is the companion page on multiple-quantum NMR.

**Prior art.** The rates are Souza et al.'s Eqs. (15) to (27), and the reading of the four levels as two qubits, |3/2⟩, |1/2⟩, |−1/2⟩, |−3/2⟩ as |00⟩, |01⟩, |10⟩, |11⟩, is theirs too, with a quantum-circuit model of the same relaxation; in their labelling the outer/inner property below is the parity Z⊗Z. They also state that the triple quantum and the central transition are the only coherences free of J₀, a decoherence-free subspace of the J₀ channel spanned by Bell states: this page's dark coherences are theirs. That the two share the slow rate J₁ + J₂, that the longitudinal relaxation is biexponential, and that the satellites of a resolved splitting decay at the bright centre's rate, is the standard picture behind triple-quantum-filtered sodium NMR; both centres are known rates. What this page adds is reading them as the centres of two exact palindromes, and the global-centre rule below; the literature was not searched for either, and the mean-rate identity follows in one line from the cited rates.

## The prediction and how it fared

**Written before computing:** with J₀ alone the spectrum is an exact palindrome, since T₂₀ is a Z on one bit (held at ω_Q = ω = 0: two levels of rate, 0 and 3J₀, eight modes each; a splitting or an offset breaks it, as the theorem below shows); with J₁, J₂ present the palindrome breaks in general, and a guess that extreme narrowing (J₀ = J₁ = J₂) restores it. The guess failed: at extreme narrowing the two centres sit at 6J and 9J and the global mirror breaks. The prediction was also too coarse in the other direction: with J₀ = 0 the palindrome survives J₁ and J₂ at every splitting. What survives is the block structure below.

**Written before computing, for other spins.** From the J₀ cell prices alone: J₀ alone is no palindrome at I = 1 (5 zero-price cells against 4), I = 5/2 and I = 7/2 (held, and proved below for every I ≥ 1 other than 3/2). Guessed: 3/2 is the only spin where J₀ alone is a palindrome, because only there the classes of |m| make one bit of two equal halves (held, and proved below for every I ≥ 1; the reason found is sharper than the guess, the five quadrupole components are Dirac matrices only at 3/2). Guessed, as an even chance: with J₀ = 0 a global centre might survive at other spins (failed without a splitting: there no rates other than all zero give one, proved below for every spin from 1 to 6 other than 3/2; with a splitting, on resonance, held at spin 1 on two tuned rays, both with J₀ = 0, and failed at spins 2 to 6). Guessed: at I = 5/2 the cells blind to J₀ do not share one rate (ill-posed: there they are not an invariant set).

## The one bit the slow field reads

T₂₀ = (3I_z² − 15/4)/√6 is diagonal, (3, −3, −3, 3)/√6 on m = 3/2, 1/2, −1/2, −3/2. It sees one thing about a level, whether it is outer or inner. Label each level by two bits, a = outer/inner and b = the sign of m; then T₂₀ = (3/√6)·Z_a, one Z on bit a (in Souza et al.'s labelling the same operator is Z⊗Z). Its dissipator charges |m⟩⟨m′| at J₀·½·(t_m − t_m′)², with t the diagonal entries: zero when m and m′ agree on bit a, ½·(6/√6)²·J₀ = 3J₀ when they disagree. That is the repository's Hamming price, on one bit of two.

| coherence | order p | outer/inner agree? | J₀ cost |
|---|---|---|---|
| populations | 0 | yes | 0 |
| ±1/2 central transition | 1 | yes (both inner) | 0 |
| ±3/2 triple quantum | 3 | yes (both outer) | 0 |
| satellites, ±3/2 ↔ ±1/2 | 1 | no | 3J₀ |
| double quantum | 2 | no | 3J₀ |

## The two mirrors

The verifier computes the characteristic polynomial of every order block with J₀, J₁, J₂, ω_Q, ω as symbols; the blocks of order −p have the conjugate polynomials of those of order p, so the statements below for p ≥ 0 carry over.

- **Dark, centre c₁ = 3(J₁ + J₂).** The populations have roots 0, −6J₁, −6J₂, −6(J₁ + J₂) at every ω_Q and ω, two pairs about −c₁. The central transition is the factor x + c₁ + iω of the order-1 polynomial and the triple quantum the whole order-3 polynomial x + c₁ + 3iω: both sit at Re λ = −c₁, each its own mirror image.
- **Bright, centre c₂ = 3(J₀ + J₁ + J₂).** The remaining order-1 quadratic and the order-2 quadratic, written in y = x + p·iω, have real coefficients and the linear coefficient 6(J₀ + J₁ + J₂). Their roots are either a real pair about −c₂ or a complex pair on Re = −c₂ (the static splitting enters only the constant term, as 6ω_Q²), and in both cases λ ↦ −λ̄ − 2c₂ maps the pair to itself. Under the conjugated map a pair is its own mirror only if both roots sit on Re λ = −c, or share their frequency with real parts summing to −2c; the real coefficients and the root sum of the shifted quadratic give exactly that. What the two mirrors carry beyond the single block is that the blocks of one sector share one centre: the populations' four roots centred on −c₁ with the central transition and the triple quantum on it, both bright quadratics centred on −c₂.
- **The distance between the centres** is c₂ − c₁ = 3J₀, the J₀ cost of a bright cell. The two centres coincide exactly when J₀ = 0, so then the whole spectrum has one centre, for every ω_Q and ω.

## When the whole spectrum has one centre

**Theorem.** For real J₀, J₁, J₂, ω_Q, ω, the sixteen eigenvalues, counted with multiplicity, have one common centre −s, under λ ↦ −λ − 2s or equally under λ ↦ −λ̄ − 2s, exactly when J₀ = 0, or when ω_Q = ω = 0 and J₁·J₂ = 0.

The two maps give the same answer because the spectrum is closed under conjugation (the order −p block carries the conjugate polynomial of the order p block). The theorem needs no sign on the J; physical spectral densities are non-negative anyway.

*The centre is fixed by the trace.* The Hamiltonian part of the generator is traceless, so Tr L = −24J₀ − 48J₁ − 48J₂ at every ω_Q and ω. A pairing has its centre at minus the mean eigenvalue, so 2s = 3J₀ + 6J₁ + 6J₂.

*Proof.* Let P(x) be the characteristic polynomial of the generator, of degree 16, the product of the seven order blocks' polynomials. The spectrum is symmetric about −s, with multiplicities, exactly when P(−x − 2s) = P(x) holds identically in x. The difference P(−x − 2s) − P(x) has fourteen non-zero coefficients, each a polynomial in J₀, J₁, J₂, ω_Q, ω with real coefficients. (The coefficient of x¹⁶ cancels because the degree is even, that of x¹⁵ by the choice of s, that of x¹⁴ vanishes identically; the first survivor, at x¹³, is 144·J₀·ω_Q², which alone forces J₀ = 0 or ω_Q = 0.) The fourteen coefficients generate the same ideal as the two polynomials

  J₀·ω_Q²,  J₀·(9J₁² + ω²)·(9J₂² + 4ω²),

which are the ideal's reduced Gröbner basis in lexicographic order, and equal ideals have equal zero sets. Over the reals, 9J₁² + ω² vanishes only at J₁ = ω = 0, and 9J₂² + 4ω² only at J₂ = ω = 0. So the coefficients all vanish exactly when J₀ = 0, or when ω_Q = 0, ω = 0 and J₁ = 0 or J₂ = 0. ∎

*Why it looks this way.* At ω_Q = ω = 0 the sixteen rates r = −λ are explicit. The generator preserves the trace, so λ = 0 is always an eigenvalue, and the mirror must send it to rate 2s. The difference 2s − r is, rate by rate, 3·(J₀ + 2J₁ + 2J₂), 3·(J₀ + 2J₂), 3·(J₀ + 2J₁), 3J₀ (the populations), 3·(J₀ + J₁ + J₂) (the dark coherences), 3·(J₁ + 2J₂), 3J₁ (the satellites), 3·(2J₁ + J₂), 3J₂ (the double quantum). For non-negative densities each of these is non-zero unless some density vanishes, so with all three positive no rate sits at 2s, 0 has no partner, and there is no centre. The partner of 0 is the largest population rate when J₀ = 0, the upper satellite when J₁ = 0, and the upper double quantum when J₂ = 0. For J₁ = 0 the reflection r ↦ 3J₀ + 6J₂ − r sends the multiset {0², (6J₂)², (3J₂)⁴, (3J₀)², (3J₀ + 6J₂)², (3J₀ + 3J₂)⁴} onto itself (exponents are multiplicities): 0 ↔ 3J₀ + 6J₂, 6J₂ ↔ 3J₀, 3J₂ ↔ 3J₀ + 3J₂. The populations pair with the satellites, the dark coherences with the double quanta. J₂ = 0 is the same with J₁ and J₂ exchanged, and the roles swap: the populations pair with the double quanta, the dark coherences with the satellites. J₀ = 0 is the coincidence of the two mirrors (the section above).

The offset enters only because it gives every block of order p ≠ 0 the frequency p·ω, while the partner of 0 must be real. The secular dissipator commutes with rotations about z, so the frame rotating at ω removes the offset; in that frame the theorem reads: one centre exactly when J₀ = 0, or ω_Q = 0 and J₁·J₂ = 0. In the laboratory frame a rate at 2s is not enough, because the partners must also match in number: at J₁ = J₂ = 0, J₀ = 1, ω_Q = √6, ω = 6 a satellite eigenvalue is real and sits exactly at −2s = −3, but 0 is a fourfold eigenvalue there and −3 only a twofold one, so the pairing still fails. The static splitting cannot be rotated away. In the rotating frame at J₁ = 0 < J₂ it draws the upper satellite inward below 2s (while 6ω_Q² < 9J₂²; beyond, the pair is complex at the rate 3(J₀ + J₂)) and turns the double-quantum pair complex; at J₁ = J₂ = 0 the satellites are complex as well. In each case the real eigenvalue 0 has no real partner at rate 2s.

The verifier gates the trace, the identity through the Gröbner basis, and the three reflections at a vanishing density. It also checks sixteen exact points, each by two routes: P(−x − 2s) − P(x) zero or non-zero at that point, and F158's end count (in the section after next). The gates can fail: a Gröbner basis that differs from the one above, or a single point that disagrees, fails the run. The points cover both outcomes; among them are six breaking cases with J₀ > 0 and J₁·J₂ = 0 away from ω_Q = ω = 0, the counting point above one of them, and J₀ = 0 with both frequencies on.

## Why spin 3/2: five Dirac matrices

The five Hermitian components of the quadrupole operators are

  Γ₀ = T₂₀,  Γ₁ = T₂,₁ + T₂,₁†,  Γ₂ = i(T₂,₁ − T₂,₁†),  Γ₃ = T₂,₂ + T₂,₂†,  Γ₄ = i(T₂,₂ − T₂,₂†).

At I = 3/2 they anticommute pairwise and square to 3/2, 3, 3, 3 and 3 times the identity: five Dirac matrices. In Souza et al.'s labelling of the levels as two qubits, |3/2⟩, |1/2⟩, |−1/2⟩, |−3/2⟩ = |00⟩, |01⟩, |10⟩, |11⟩, they are, up to scalars, the Pauli strings ZZ, ZX, ZY, XI and YI. (In the bits a, b of the section above T₂₀ is Z_a; Souza et al.'s two bits are another pair, whose parity is a.) The algebra is known: C. Wu, J.-P. Hu and S.-C. Zhang, *Exact SO(5) symmetry in the spin-3/2 fermionic system*, Phys. Rev. Lett. 91, 186402 (2003), [arXiv:cond-mat/0302165](https://arxiv.org/abs/cond-mat/0302165), build the SO(5) of spin-3/2 fermions on five Dirac matrices whose bilinears carry spin 2, so their five Γ's span the quadrupole matrices. What follows is what the algebra does to the relaxation.

**The dissipator is Pauli-string dephasing.** Each pair T₂,±ₘ is the pair of its two Hermitian components at half the rate, so the generator's dissipator is Σₐ γₐ·(Γ̂ₐρΓ̂ₐ − ρ), with Γ̂ₐ the Γₐ scaled to square 1 and γ = 3J₀/2, 3J₁/2, 3J₁/2, 3J₂/2, 3J₂/2. At ω_Q = ω = 0 the sixteen products 1, Γₐ and ΓₐΓ_b (a < b) are its eigenoperators, and each decays at twice the summed rate of the jumps it anticommutes with. That one rule is the whole rate table of the Abstract:

| operators | coherence | rate |
|---|---|---|
| 1 | population | 0 |
| Γ₀ | population | 6(J₁ + J₂) |
| Γ₁Γ₂ | population | 6J₁ |
| Γ₃Γ₄ | population | 6J₂ |
| Γ₁Γ₃, Γ₁Γ₄, Γ₂Γ₃, Γ₂Γ₄ | central transition and triple quantum | 3(J₁ + J₂) |
| Γ₁, Γ₂ | satellites | 3J₀ + 3J₁ + 6J₂ |
| Γ₀Γ₁, Γ₀Γ₂ | satellites | 3(J₀ + J₁) |
| Γ₃, Γ₄ | double quantum | 3J₀ + 6J₁ + 3J₂ |
| Γ₀Γ₃, Γ₀Γ₄ | double quantum | 3(J₀ + J₂) |

The dark operators are the eight that commute with Γ₀, the bright ones the eight that anticommute with it.

**One map for both mirrors.** With Γ₀ scaled to diag(1, −1, −1, 1), the linear map S(X) = Γ₀·Xᵀ satisfies, exactly and at every ω_Q and ω,

  S·L·S⁻¹ = −L + 2·L_{J₀} − 6(J₁ + J₂),

and so does X ↦ Xᵀ·Γ₀. Here L_{J₀} is the J₀ part of the dissipator, L_{J₀} = −3J₀·P with P the projector on the bright cells; P commutes with L and with S. So S maps each sector to itself, with λ ↦ −λ − 2c₁ on the dark sector and λ ↦ −λ − 2c₂ on the bright one. The two mirrors of this page are one map:
- The transpose reverses the Hamiltonian part, because H is real, and keeps every jump term, because every Γ is symmetric or antisymmetric.
- Multiplying by Γ₀ keeps the Hamiltonian part, because Γ₀ commutes with H.
- It reflects the J₁ and J₂ dissipators, whose jumps anticommute with Γ₀.
- It keeps the J₀ dissipator, whose jump is Γ₀ itself, and that kept part moves the bright centre out by 3J₀.

At J₀ = 0 nothing is kept and S is an F1 palindromizer. It is built the way [F118](../docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md) builds the qubit palindromizers, Π_Z(ρ) = ρᵀ·X^⊗N and Π_Y(ρ) = X^⊗N·ρᵀ, with Γ₀ in place of X^⊗N.

**A second route, for non-negative densities.** The normalized jumps are Hermitian and square to 1. So for J ≥ 0, a jump of rate zero being simply absent, [F158](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md) applies as it stands, with σ = Σγ = s. The spectrum pairs about −s exactly when dim ker L = dim ker(L + 2s). The near kernel holds the operators that commute with H and with every jump present; the far kernel holds those that commute with H and anticommute with every jump present. One fact about H keeps the count short. H is diagonal, so an operator commutes with it exactly when it lives on the cells |m⟩⟨m′| whose energy gap E_m − E_m′ vanishes, where E_m is the diagonal of H.
- **All three densities positive.** The five Γ's generate every 4 × 4 matrix, so the near kernel holds only the scalars. Nothing anticommutes with all five, because their product is a multiple of 1, so the far kernel is empty. The count is one against none: no centre.
- **J₀ = 0.** Γ₀ anticommutes with every jump present, and it commutes with H, both being diagonal, at every ω_Q and ω. It is invertible, and by F158 an invertible element of the far kernel forces the two kernels to have equal dimension, which gives the centre. In the two-qubit labelling Γ₀ is the string ZZ: a Pauli string in the far kernel commuting with H, what `EndCount` calls a colouring after [the palindrome as a colouring](THE_PALINDROME_AS_A_COLOURING.md), whose colourings are the one-letter-per-site case. If J₁ = J₂ = 0 as well, L = −i[H, ·], and its spectrum, the energy gaps, is symmetric about 0 = s.
- **J₁ = 0 < J₀, J₂.** The operators that anticommute with Γ₀, Γ₃ and Γ₄ span Γ₁ and Γ₂, which is the span of T₂,₁ and T₂,₁†. T₂,₁ lives on the two cells |3/2⟩⟨1/2| and |−1/2⟩⟨−3/2|, with gaps √6·ω_Q + ω and −√6·ω_Q + ω, and T₂,₁† lives on their transposes. So the far kernel is two-dimensional at ω_Q = ω = 0 and empty otherwise. The near kernel is spanned by 1 and Γ₁Γ₂ at every ω_Q and ω. A centre exists exactly at ω_Q = ω = 0. The case J₂ = 0 < J₀, J₁ is the same, with Γ₃, Γ₄, T₂,₂ and the gaps ±√6·ω_Q + 2ω.
- **J₁ = J₂ = 0 < J₀.** The near kernel holds the dark operators with zero gap, the far kernel the bright ones.
  - The dark gaps are 0 (four populations), ±ω (central transition) and ±3ω (triple quantum).
  - The bright gaps are ±(√6·ω_Q + ω) and ±(√6·ω_Q − ω) for the satellites, ±(√6·ω_Q + 2ω) and ±(√6·ω_Q − 2ω) for the double quantum.
  - At ω = 0 the counts are 8 against 8 if ω_Q = 0, and 8 against 0 otherwise.
  - At ω ≠ 0 the near kernel has dimension 4. At most one of the four conditions √6·ω_Q = ±ω, ±2ω can hold, so the far kernel has dimension at most 2, and there is no centre. The counting point above, ω = √6·ω_Q, is a case with exactly 2.

This is the theorem again, for J ≥ 0. Negative densities lie outside F158's class, and the Gröbner route alone covers them. The verifier checks the rank equality at all sixteen points beside the polynomial identity: two exact routes, one answer.

**Only at 3/2.** For I ≥ 1, Γ₀² ∝ (3m² − I(I + 1))² is a multiple of 1 only if 3m² − I(I + 1) takes two opposite values. At I = 3/2 the values are ±3. At I = 1 they are 1 and −2, and every I ≥ 2 has at least three values of m². So among the spins with a quadrupole moment, 3/2 is the only one whose quadrupole operators form a Clifford algebra.

## Beyond spin 3/2

**J₀ alone breaks the mirror at every I ≥ 1 other than 3/2 (proved).** The J₀ dissipator charges |m⟩⟨m′| the price J₀·½(t_m − t_m′)², where t_m = (3m² − I(I + 1))/√6. A mirror must send the price 0 to the top price, so the two must occur equally often.
- The cells of price 0 are those with m² = m′². Each class of m² holds two levels, except m = 0 at integer I, so these cells number 4I + 2 at half-integer I and 4I + 1 at integer I.
- The top price joins only the two extreme classes, m² = I² and the smallest m². It falls on 2·2·2 = 8 cells at half-integer I and on 2·2·1 = 4 at integer I.

The counts agree only at I = 3/2 (8 against 8). At I = 1, 2, 5/2, 3 and 7/2 they are 5 against 4, 9 against 4, 12 against 8, 13 against 4 and 16 against 8.

**No global centre without a splitting, for any rates other than all zero (proved for every spin from 1 to 6 in half-integer steps except 3/2).** At each of these spins, at ω_Q = 0 and every offset ω, no real rate triple other than J₀ = J₁ = J₂ = 0 gives the whole spectrum one centre, negative densities included. The route needs no eigenvalue, and it is the repository's moment route (the sweep above):
- A centre would sit at −s, fixed by the trace as at 3/2. Pairing λ with −λ − 2s makes every odd centred moment μ_k = Tr (L + s)^k vanish, the easy half of the power-sum lemma.
- A diagonal similarity makes the generator rational: it rescales I₊ and I₋ to rational entries, commutes with the diagonal H, and keeps the spectrum. At ω_Q = ω = 0, μ₃, μ₅ and μ₇ are then exact homogeneous polynomials in J₀, J₁, J₂, of degrees 3, 5 and 7.
- At every spin listed, μ₃ and μ₅ have finitely many common real zeros: 5 projective points at I = 1, 3 at I = 2 to 7/2, 11 at I = 4 to 6. μ₇ is non-zero at each of them, so at ω_Q = ω = 0 there is no centre.
- The offset cannot create a centre. At ω_Q = 0 the generator is the dissipator (the generator at ω = 0) plus the offset term −iω·[I_z, ·]. The dissipator is self-adjoint, because every jump T₂,ₘ comes with its adjoint, ±T₂,₋ₘ, at the same rate, and it commutes with [I_z, ·], because [I_z, T₂,ₘ] = m·T₂,ₘ. So the block of order p has the real spectrum of the zero-frequency block, shifted by −i·p·ω. For ω ≠ 0 a pairing λ ↦ −λ − 2s must then send block p to block −p, whose spectrum at zero frequency is the same (ρ ↦ ρ† carries one block onto the other and conjugates a real spectrum to itself); so every block would be palindromic about −s by itself, and the zero-frequency spectrum with it.

At I = 3/2 the same route returns μ₃ = 0, μ₅ = 14580·J₀J₁²J₂² and μ₇ = 229635·J₀J₁²J₂²·(J₀² + 2J₁² + 2J₂²), which vanish exactly where the theorem above puts a centre: the route does not say no where the theorem says yes. The test runs one way only: all odd moments vanishing would prove a centre, three vanishing would not. The verifier checks the two jump identities at every spin, recomputes the moments, finds the common zeros through the resultant and the first subresultant, reads the signs and μ₇ by interval arithmetic on exact isolating intervals, and checks the counts given here.

Up to I = 7/2 every common real zero of μ₃ and μ₅ has a negative density, so for physical rates the first two odd moments already decide. At I = 4 to 6, five or six of the common zeros have all three densities positive (of them with J₀ > J₁ > J₂, the order a monotone spectral density gives: one at I = 4 to 5, two at 11/2 and 6), and only μ₇ excludes them. This is a proof spin by spin, not for every I. At I = 5/2 the dark/bright split itself is gone: the cells with |m| = |m′|, blind to J₀, are not an invariant set, because J₁ and J₂ carry them into the rest.

**A static splitting can make a centre, at spin 1 only.** Proved on resonance (ω = 0) for every spin from 1 to 6 in half-integer steps except 3/2. The splitting ω_Q·T₂₀ cannot be rotated away; it makes the generator non-normal, and its spectrum can be complex. The moment route still runs:
- The splitting is diagonal and commutes with the rational similarity. The odd moments become exact polynomials in J₀, J₁, J₂ and q = ω_Q/√6, homogeneous of degree k, with only even powers of q, because every trace is real.
- μ₃ is linear in u = q². Solving it for u and eliminating, the common zeros of μ₃, μ₅, μ₇ and μ₉ are found exactly: resultants in J₁, their common factors in J₀, a Gröbner basis for each non-linear factor, the case where μ₃'s coefficient of u vanishes, the line J₂ = 0 and its corner (1, 0, 0).
- At I = 1 the generator has nine eigenvalues. By Newton's identities four vanishing odd moments are then the whole condition, so here the route decides both ways. At I ≥ 2 they are necessary only, which is enough for a no.

At I = 2 to 6 the four moments have no common real zero other than all rates zero, where L = −i[H, ·] is trivially symmetric: a splitting gives no centre there either. At I = 1, besides that point, they vanish together on exactly six rays (J₀, J₁, J₂, ω_Q²), rays under J ↦ tJ, ω_Q ↦ tω_Q: on (4, 0, 1, −1/6), where ω_Q would be imaginary, and on five real ones, at each of which the characteristic polynomial satisfies (−1)ⁿP(−x − 2s) = P(x), n = 9, exactly. Three of the five need a negative density. With non-negative rates exactly two remain, and both have J₀ = 0:
- J₁ = 4J₂/5 at ω_Q² = J₂²/10: the spectrum is {0, −1², −(6/5)³, −(7/5)², −12/5}·J₂, palindromic about −6/5·J₂.
- J₁ = 8J₂ at ω_Q² = 21J₂²/2: {0, −5², −6³, −7², −12}·J₂, palindromic about −6·J₂.

*How it happens*, the two rays found again by hand. At spin 1, in rates r = −λ:
- the populations are 0, 3J₁/2 and (J₁ + 4J₂)/2 (the last is the deuterium longitudinal rate);
- the double quantum is (J₁ + 2J₂)/2, twice;
- the single-quantum pair is the two roots of a quadratic, each twice (orders ±1), whose sum is (3J₀ + 3J₁ + 2J₂)/2 and whose product carries 3ω_Q²/2, the only place the splitting enters;
- and 2s = 4(J₀/2 + J₁ + J₂)/3.

Put one single-quantum rate on the centre s and let the other mirror the double quantum: the sum then forces J₀ = 0. The populations must pair among themselves, which fixes J₁ = 4J₂/5 (3J₁/2 on the centre) or J₁ = 8J₂ ((J₁ + 4J₂)/2 on the centre), and the product fixes ω_Q². The splitting is tuned until a single-quantum mode, which it pushes, lands on the mirror image of the double quantum, which it cannot move. At both rays every eigenvalue is simple in its block and the generator is semisimple: no exceptional point. The pairing crosses the coherence orders, −1 from the single quantum with −7/5 from the double quantum at the first ray. This route reads one pairing pattern; that no other pattern adds a ray is the elimination's part.

**The offset, and which frame.** ω is the offset from the frame rotating at the Larmor frequency, and the secular dissipator commutes with rotations about z, so ω = 0 is that frame, on resonance. At both rays every block spectrum is real, so an offset only adds the frequencies −p·ω: the nine decay rates stay palindromic about s in every frame. The complex map λ ↦ −λ − 2s, which also flips the frequencies, holds on resonance only: for ω ≠ 0 the orders stay apart, and the single-quantum lines alone, {1, 6/5}·J₂ or {6, 7}·J₂, are not symmetric about s. Off the rays and away from all rates zero, at spin 1 and at I ≥ 2, the coefficients of (−1)ⁿP(−x − 2s) − P(x) are polynomials in ω, and one of them is non-zero at ω = 0, so at each rate triple and splitting at most finitely many offsets could give a centre; which, if any, is not decided here.

**What the rays are.** With a static splitting the medium is ordered, and J₀, J₁, J₂ are three densities, those of the field gradient's m = 0, ±1 and ±2 components at the frequencies 0, ω₀ and 2ω₀; they need not follow one curve J(ω). A rotation about the field axis leaves the m = 0 component unchanged, so motion that only turns the gradient about the field has J₀ = 0: the rays need a component that does not fluctuate, not an exotic density. Each ray is still fine-tuned, three conditions on four parameters, and any further slow motion gives J₀ > 0 and leaves it. Both sit just below the exceptional point of the single-quantum pair, 6ω_Q² = J₁² (6ω_Q²/J₁² = 15/16 and 63/64): the splitting draws the pair's two rates together, from {7/10, 3/2}·J₂ to {1, 6/5}·J₂ at the first ray, and stops short of merging them. Below that point the two modes share one frequency, so the splitting at the rays is smaller than the linewidth: no resolved quadrupolar doublet is seen. At ω_Q = 0 the dissipator alone pairs 3 of the 9 eigenvalues about −s at both rays; with the splitting all 9 pair.

**What it says beyond the nucleus.** At spin 3/2 the splitting only removes centres. Inside [F158](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md)'s class that is forced where the far end is empty, as it is at 3/2 with all three densities positive: both ends must also commute with H (Lemma 1), so a Hamiltonian can shrink the far end but not open it, while the near end keeps 1. The spin-1 jumps are not Hermitian operators squaring to 1, so they lie outside that class, and there a Hamiltonian can supply a centre the dissipator lacks. Inside the class no palindrome exists at all at odd total dimension with a jump present (F158 §(f5)); spin 1 has d = 3. [F121](../docs/ANALYTICAL_FORMULAS.md) notes for its d = 3 example, two qutrits under equidistant dephasing with the SU(3) Heisenberg coupling, that the Hamiltonian helps at no centre, and marks those sampled counts as a property of that H, not a theorem; one spin 1 under quadrupolar relaxation is a d = 3 case where a Hamiltonian does help, at tuned points.

## What a sodium measurement can read

Under this generator the triple-quantum rate, equal to the slow transverse rate, is the mean of the two biexponential longitudinal rates: c₁ = (6J₁ + 6J₂)/2, at every static splitting, since none of the three rates depends on ω_Q. In Souza et al.'s units (rates C·J) the longitudinal rates are 2CJ₁ and 2CJ₂ (the third population rate, 2C(J₁ + J₂), belongs to the quadrupolar-order mode Γ₀ of the section on Dirac matrices above, which I_z recovery does not enter), and the slow transverse and triple-quantum rate is C(J₁ + J₂). It is a property of the model, not a new measurement: a relaxation study that reports these rates of one sample can check it, and a deviation beyond the measurement error and the fit bias says that something outside the model is present (more than one sodium environment in the voxel, dynamic frequency shifts, motion that the secular Redfield form does not capture).

## What this does not show

No sodium data were compared. The model is a single isolated spin 3/2 at infinite temperature in the secular Redfield picture; multiple compartments, exchange, finite polarization and non-secular terms are outside it. Nothing here improves a sodium image. The mean-rate identity follows directly from rates published before this page; only its place as a mirror centre is ours.

## Open

- The identity against measured sodium relaxation rates (gels, cartilage, brain), whether it is already in use as a consistency check, and which deviation a two-compartment voxel produces.
- Beyond spin 3/2: a proof that no rates other than all zero give a global centre at every spin (proved only spin by spin to I = 6, without a splitting); with a splitting, which offsets ω ≠ 0 could still give a centre off the spin-1 rays (finitely many at each point away from all rates zero); and whether a weaker, block-wise mirror survives there.
- An ensemble of spins 3/2. With an independent gradient on each spin the Clifford reading makes it a Pauli-string chain on twice as many qubits, inside [F158](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md)'s class. A shared gradient gives collective jumps, sums of strings that do not square to 1, which is the setting of [The Mirror under a Shared Field](THE_MIRROR_UNDER_A_SHARED_FIELD.md) instead.
- A C# witness: MirrorWorld's `EndCount` at N = 2 could carry the theorem now, with the jumps ZZ, ZX, ZY, XI, YI and H = ω_Q·(3/√6)·ZZ + ω·(ZI + ½·IZ), since T₂₀ = (3/√6)·ZZ and I_z = ZI + ½·IZ in this labelling. `EndCount` takes integer inputs, so ω_Q is to be chosen as a multiple of √6/3, or the ZZ coefficient taken as the parameter.
- The general d > 2 question of Qubit Necessity §8a: which other physical dissipators split a palindrome into centres instead of breaking it.

## Links

- [F121](../docs/ANALYTICAL_FORMULAS.md) and [PROOF_QUDIT_PARTIAL_PALINDROME](../docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md), the equidistant d > 2 case.
- [Qubit Necessity](../docs/QUBIT_NECESSITY.md), the d² − 2d = 0 trunk and the open d > 2 dissipators.
- [The Mirror under a Shared Field](THE_MIRROR_UNDER_A_SHARED_FIELD.md), blockwise centres for qubit chains.
- [The Shared Field in a Crystal](THE_SHARED_FIELD_IN_A_CRYSTAL.md), the multiple-quantum companion page.
- [F158](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md), the two-end count, and [F118](../docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md), the palindromizer as a transpose and a multiplication.
- [Why Depolarizing Noise Breaks the Palindrome](DEPOLARIZING_PALINDROME.md) and [the palindrome as a colouring](THE_PALINDROME_AS_A_COLOURING.md), the qubit cases of the same count.
