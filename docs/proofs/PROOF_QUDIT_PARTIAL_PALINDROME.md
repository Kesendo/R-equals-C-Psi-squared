# PROOF: the partial palindrome at d>2 is the symmetric overlap of the disagreement count

**Status:** Tier 1 derived (one-line combinatorial identities, exact; verified machine-exact for the dissipator spectrum at d = 3, N = 2 and against brute enumeration on a (d, N) grid).
**Date:** 2026-06-11
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Builds on:**
- [the qubit necessity](../QUBIT_NECESSITY.md): the per-site split d : (d²−d) and the polynomial trunk d² − 2d = 0 that closes only at d = 2 (typed as `QubitNecessityPi2Inheritance`).
- [Absorption Theorem](PROOF_ABSORPTION_THEOREM.md): a qubit Pauli cell has dissipator eigenvalue −2γn_XY, a basis coherence has n_XY=Hamming(i,j), and a full eigenmode has Re(λ)=−2γ⟨n_XY⟩.
- [on the one diagonal](../../reflections/ON_THE_ONE_DIAGONAL.md): the dissipator IS the recentred disagreement-count diagonal; this proof is that diagonal read one dimension up.

## Abstract

The palindromic mirror is exact only for qubits (d = 2); [the qubit necessity](../QUBIT_NECESSITY.md) proves this from the per-site balance d = d² − d, i.e. d² − 2d = 0, which closes only at d = 2. For d > 2 the spectrum is not random: N = 2 qutrits were observed to pair 36–52 of 81 eigenvalues, a residual structure no principle captured (OQ-002). This proof captures it. Under full-Cartan dephasing (the complete diagonal Cartan subalgebra as jump operators), the d levels are mutually **equidistant**, so the dissipator eigenvalue of a coherence |i⟩⟨j| is exactly −2γ·Hamming(i,j), the **same cell ladder as the qubit**. What differs is the multiplicity per rung: the number of coherences at Hamming distance k is

  **c_k = d^N · C(N, k) · (d−1)^k.**

The palindrome reflects rung k against rung N−k. For d = 2 the factor (d−1)^k = 1 and c_k = 2^N·C(N, k) is symmetric in k, so every rung finds its mirror: 100%. For d > 2 the (d−1)^k tilts the distribution toward large k, and only the symmetric overlap pairs. The dissipator's paired ceiling is

  **paired(d, N) = Σ_k d^N · C(N, k) · (d−1)^{min(k, N−k)},**

which equals d^{2N} (full) iff d = 2. For d = 3, N = 2 it is 54/81; the qutrit fraction erodes with N (66.7%, 66.7%, 51.9% for N = 1, 2, 3); d = 4 gives 50%, 50%, 31.2%. The (d−1) tilt base is exactly the per-site decaying-to-immune ratio (d²−d) : d = (d−1) : 1 from [the qubit necessity](../QUBIT_NECESSITY.md), raised to the number of disagreeing sites. The full interacting Liouvillian is richer and is left open in §4.

## §1 The equidistant ladder

Take as jump operators the diagonal Cartan generators (for qutrits λ₃ = diag(1, −1, 0) and λ₈ = diag(1, 1, −2)/√3; generally a basis of the diagonal subalgebra). Each is diagonal, so a coherence |i⟩⟨j| is an eigenmode of the single-site dissipator with rate

  Σ_M [M_{ii}M̄_{jj} − ½(|M_{ii}|² + |M_{jj}|²)] = −½ Σ_M (M_{ii} − M_{jj})².

For the qutrit Cartan {λ₃, λ₈} this evaluates to −2 for every i ≠ j: (0,1) gives ½(4 + 0); (0,2) gives ½(1 + 3); (1,2) gives ½(1 + 3); all equal 2. The three levels are mutually equidistant. Hence per site the rate is −2γ if i ≠ j and 0 if i = j, and across N sites

  rate(|i⟩⟨j|) = −2γ · #{l : i_l ≠ j_l} = −2γ · Hamming(i, j),

the identical dissipator-cell ladder to the qubit's −2γn_XY grading ([Absorption Theorem](PROOF_ABSORPTION_THEOREM.md)). The diagonal generators are immune (rate 0); the d² − d off-diagonal generators decay. A full Liouvillian eigenmode is governed by the corresponding Hilbert-Schmidt expectation, not generally by one integer Hamming label. This is the per-site d : (d²−d) split of [the qubit necessity](../QUBIT_NECESSITY.md), now read as a rung of the rate ladder.

## §2 The multiplicity tilt

Count the coherences at Hamming distance k. Choose the k disagreeing sites: C(N, k). On each, j_l ≠ i_l is one of d−1 values: (d−1)^k. The ket i ranges over all d^N states. Hence

  **c_k = d^N · C(N, k) · (d−1)^k,**   Σ_k c_k = d^N · Σ_k C(N,k)(d−1)^k = d^N · d^N = d^{2N}.

The factor (d−1)^k is the only place the dimension enters non-trivially. Its base d−1 is the per-site decaying : immune ratio (d² − d)/d, raised to the number of disagreements. At d = 2 it is 1^k = 1 and c_k = 2^N·C(N, k), the symmetric binomial.

## §3 The ceiling, and why only d = 2 is full

The dissipator spectrum is real, with rung k at rate −2γk and multiplicity c_k. The palindrome reflects λ ↦ −2·(Nγ) − λ about the center −Nγ, i.e. rung k against rung N−k. Modes at rung k find a partner iff there is a mode at rung N−k; the number that pair is, summing each two-rung pair and the self-mirrored middle (N even),

  paired(d, N) = Σ_{k < N/2} 2·min(c_k, c_{N−k}) + [N even]·c_{N/2}.

Because C(N, k) = C(N, N−k) and d − 1 ≥ 1, min(c_k, c_{N−k}) = d^N·C(N, k)·(d−1)^{min(k, N−k)}, so

  **paired(d, N) = Σ_k d^N · C(N, k) · (d−1)^{min(k, N−k)}.**

This is d^{2N} (everything pairs) iff (d−1)^{min(k,N−k)} = (d−1)^k for all k, i.e. iff d − 1 = 1, i.e. **d = 2**. For d > 2 the high rungs (k > N/2) carry more modes than their low-rung partners; the excess

  Σ_{k > N/2} d^N·C(N, k)·[(d−1)^k − (d−1)^{N−k}]

is unpaired. For d = 3, N = 2: c = [9, 36, 36], rung 0 (×9) pairs into rung 2 (×36) leaving 27, rung 1 (×36) self-mirrors; paired = 54, excess = 27. The qutrit fraction erodes with N (66.7%, 66.7%, 51.9%); d = 4 gives 50%, 50%, 31.2%. The d² − 2d = 0 uniqueness of [the qubit necessity](../QUBIT_NECESSITY.md) reappears here as the unique fully-paired column of an N-indexed family.

## §4 The sampled interacting spectrum

The ceiling above is the **dissipator's** palindrome, taken about the physical center −Nγ (the k ↔ N−k reflection, where the qubit palindrome is exact). For the symmetric SU(3) Heisenberg at N = 2, adding H lowers the count from 54 to 48 about −Nγ = −2γ and from 72 to 60 about −3γ (where the two large dissipator rungs are equinumerous), and there is no centre at which it helps. A pair λ ↔ μ can only form about the centre −(Re λ + Re μ)/2, so the centres worth checking are the midpoints of the real-part levels, {0, −2γ, −4γ} for the dissipator and {0, −2γ, −3γ, −4γ} for the full L: eight in all, 0, −γ, −1.5γ, −2γ, −2.5γ, −3γ, −3.5γ, −4γ. At the five integer ones the dissipator pairs 9, 18, 54, 72, 36 and the full L pairs 6, 12, 48, 60, 27; at the three half-integer ones both pair nothing; and at every other centre both counts are zero. So the full count never exceeds the dissipator's, at J = 0.05, 1 and 10 alike. That is a property of this H and not of every H: `H=cI` has `L_H=0` and leaves the dissipator result unchanged. The interacting count must therefore be stated for the chosen H and center.

The [Absorption Theorem](PROOF_ABSORPTION_THEOREM.md) gives the universal
right-eigenoperator Hilbert-Schmidt reading
`Re(λ)=−2γ(v†Qv)/(v†v)` for every Hermitian H. In the SU(3) Heisenberg case,
symmetry additionally quantizes those right-HS expectations to
`{0,1,1.5,2}` with multiplicities `{6,36,12,27}`; the `−3γ` rung has
`<Q>=1.5`. A generic H spreads the right-HS means between these special
rungs, but does not break the Absorption Theorem. The producer's separate
biorthogonal quantity `w†Qv` must not be substituted for this Rayleigh ratio.

So the interacting partial palindrome has **no H-independent closed form**: the paired count floats with the Hamiltonian (60 for the SU(3) Heisenberg at the tested nonzero couplings `J=0.05,1,10`; at `J=0` it returns to the dissipator counts 54/72; 48 and 52 for less symmetric couplings; near zero for a generic H). The dissipator's closed form (54 about −Nγ, full iff d = 2) is the only invariant skeleton; the interacting count is a property of each H's symmetry. Verified in [`simulations/qutrit_interacting_palindrome.py`](../../simulations/qutrit_interacting_palindrome.py) (self-validating, including the centre scan). The counts are read where they do not depend on the matching tolerance: they are the same at every tolerance from 10⁻¹⁰ to 10⁻⁶, while the smallest physical splitting in the spectrum, the detuning of the −3γ rung from 4J, is 4J − √(16J² − γ²) ≈ γ²/(8J) = 3.1·10⁻⁵ at J = 10, so a tolerance of 10⁻⁴ would merge levels the physics keeps apart. The SU(3) Heisenberg's specific 60 is decoded by representation theory in §8.

## §5 What is ours and what is the home

The d : (d²−d) per-site split, the binomial rung count, and the d = 2 necessity are catalogued (the necessity is typed as `QubitNecessityPi2Inheritance`; the qubit rate law is the Absorption Theorem). What this proof banks is the **synthesis**: the equidistant-ladder reading that carries the qubit rate law verbatim to d > 2, the multiplicity tilt (d−1)^k as the per-site ratio raised to the disagreement count, and the closed-form ceiling whose unique full column recovers d² − 2d = 0 as an N-family. The verification anchor is [`simulations/qutrit_partial_palindrome.py`](../../simulations/qutrit_partial_palindrome.py) (self-validating); the typed claim is `QuditPartialPalindromeCeiling` (`compute/RCPsiSquared.Core/Symmetry/`, parent `QubitNecessityPi2Inheritance`).

## §6 The operator realization: the product cap, Π's formula one dimension up, and the wreath family

The sections above count the pairs; this section asks **which operator performs the pairing**, and the answer brings the trunk polynomial in a third time.

**Theorem (the product cap).** Let W = ⊗_l q_l be any per-site mirror, with arbitrary site-dependent q_l, one-sided or two-sided, linear or antilinear, that intertwines the dissipator palindrome W L_D = (−L_D − 2Nγ)·W. Then W pairs at most

  **P(d, N) = max_m (2d)^(N−2m) · (d³ − d²)^m of the d^{2N} coherences,**

which is **(2d)^N for every N when d ≤ 5** and larger from d = 6 on.

*Proof.* Per site, split the letters into dark = {|x⟩⟨x|} (d of them, rate 0) and lit = {|i⟩⟨j|, i ≠ j} (d² − d, rate −2γ), and write q_l in blocks between the two classes. Give a block the lit grade c = (lit out) + (lit in) ∈ {0, 1, 2}: dark → dark has c = 0, dark ↔ lit has c = 1, lit → lit has c = 2. The rate is additive over sites, so the intertwining asks, on every nonzero product of blocks, that the out-rung plus the in-rung equal N, i.e. Σ_l c_l = N. Fix a nonzero block choice at every other site and vary site l: two nonzero blocks of q_l with different grades would give two nonzero products with different sums, and one of them breaks the identity. So every nonzero block of q_l carries one grade c_l. The ranks are then bounded per grade, r(0) ≤ d, r(1) ≤ 2·min(d, d² − d) = 2d, r(2) ≤ d² − d, and rank(W) = Π_l rank(q_l). Σ_l c_l = N with c_l ∈ {0, 1, 2} forces as many sites at grade 0 as at grade 2, say m of each, so rank(W) ≤ (2d)^(N−2m)·(d·(d² − d))^m, and the maximum over m is P(d, N). Every pattern is attained: an integer factor of one grade reaches its bound with residual exactly zero, and a factor that mixes grades breaks the identity. ∎

Where the optimum sits is one comparison. Two sites at grades (0, 2) give d·(d² − d) = d³ − d²; two sites at (1, 1) give (2d)² = 4d². The first wins iff d²(d − 1) > 4d², i.e. iff **d ≥ 6**, with a tie at d = 5. So for d ≤ 5, the qutrit included, the strict per-site dark ↔ lit swap is an optimum (the only one below d = 5; at d = 5 the (0, 2) patterns tie with it) and the cap is (2d)^N at every N; from d = 6 on, pairing a dark-only site with a lit-only site beats the swap. The smallest case is d = 6, N = 2: q₀ = P_dark and q₁ = P_lit give W = P_dark ⊗ P_lit, supported entirely on the self-complementary rung h = 1, with rank

  **d(d² − d) = 180 > 144 = (2d)^N,**

and it obeys the identity because source and image both decay at −2γ.

The cap is the full space iff P(d, N) = d^{2N}, i.e. iff **d² − 2d = 0**, at every d: for d ≤ 5 that is (2d)^N = d^{2N}, and from d = 6 on every grade rank is below d², so nothing is full. This is [the qubit necessity](../QUBIT_NECESSITY.md) trunk as an operator bound, its third appearance after the per-site split and the ceiling column. And the cap is strictly below §3's combinatorial ceiling for every d ≥ 3, N ≥ 2 (at d = 3, N = 2: cap 36 < ceiling 54). Divide both by d^N. For d ≤ 5, P/d^N = 2^N = Σ_k C(N, k), while the ceiling is Σ_k C(N, k)·(d − 1)^{min(k, N−k)}, strictly larger once d ≥ 3 and some 0 < k < N exists, i.e. N ≥ 2. For d ≥ 6 with m = ⌊N/2⌋, P/d^N = 2^{N−2m}·(d − 1)^m, while the ceiling holds C(N, m)·(d − 1)^m plus the positive k = 0 term at even N (and C(N, m) ≥ 1), and 2·C(N, m)·(d − 1)^m ≥ 2N·(d − 1)^m > 2·(d − 1)^m at odd N ≥ 3. So: since a global, non-product partial isometry reaches the ceiling exactly (greedy rung matching, exact intertwining on its support; verified), **the gap ceiling − P(d, N) is precisely the non-product part of the partial palindrome**, 18 at d = 3, N = 2. This is the inverse of the [golden-router story](PROOF_CEILING_GOLDEN_ROUTER.md): there the suspected non-locality dissolved because window sums gave per-site maps room the per-term test could not see; here the dissipator is strictly local, the rung count is rigid, and the non-locality is provable.

**The operator attaining the cap for d ≤ 5 is the qubit palindromizer's own formula.** Define, exactly as in [F118](PROOF_PI_FACTORS_AS_R_TIMES_D.md) (Π_Z(ρ) = ρᵀ·X^⊗N),

  **Π_d(ρ) = ρᵀ · Shift^⊗N,**  Shift the cyclic clock shift |x⟩ ↦ |x+1⟩.

The full Π_d is a permutation and therefore has rank d^{2N}. Per site Π_d sends the letter (i, j) to (j, i−1): dark (x, x) ↦ (x, x−1), lit-aligned (a, a−1) ↦ (a−1, a−1) dark. On the **shift-aligned subspace**, per-site span{(x,x)} ∪ {(a, a−1)} of dimension 2d per site, hence (2d)^N in all, Π_d is closed and the intertwining residual is **exactly zero** (verified at d = 3, N = 1..3 and d = 4, N = 1..2), so Π_d P_aligned is a grade-1 mirror at every site and attains the cap whenever d ≤ 5; on the complement this realization fails at O(γ). There are two chiralities Π_d^± (the two shift directions, aligned subspaces {(a, a∓1)}); at d = 2 the two off-diagonals coincide, the chiralities merge, the aligned subspace is everything, and the mirror is full: **the qubit magic of §3 is, in operator language, the degeneracy of the two shift chiralities.**

**The mirror group becomes a wreath family.** With D the transpose (F118's reflection), the closure obeys, verified exactly at d = 2, 3, 4, 5:

  **ord(Π_d) = 2d,  |⟨Π_d, D⟩| = 2d²,  ⟨Π_d, D⟩ ≅ Z_d ≀ Z₂**

(the elements are the index swaps-or-not composed with independent two-sided shifts (a, b) ∈ Z_d × Z_d; D-conjugation exchanges the two shift factors rather than inverting them, which is what makes the extension a wreath product and not a generalized dihedral group). At d = 2 this is exactly D₄: **the [F118 mirror group](PROOF_PI_FACTORS_AS_R_TIMES_D.md) is the d = 2 column of a d-indexed wreath family**, and for d ≥ 3 the reflection D swaps the Π_d^± chiralities instead of normalizing one mirror, the group-level face of the partiality.

The verification anchor for this section is [`simulations/qudit_product_mirror_cap.py`](../../simulations/qudit_product_mirror_cap.py) (self-validating: cap and trunk arithmetic; the product lemma with integer grade factors, exact residuals and a mixed-grade negative control; P(d, N) against the enumeration of every grade pattern and the d = 6 projector; exact alignment residuals; the global ceiling-reacher; the group law; the d = 2 degeneracy). The typed claim is `QuditProductMirrorCap` (`compute/RCPsiSquared.Core/Symmetry/`, parents `QuditPartialPalindromeCeiling` + `QubitNecessityPi2Inheritance`), whose battery recomputes the lemma, the threshold and the d = 6 control at construction. The first question this section leaves open, whether a translation-invariant non-product mirror can exceed the cap, is answered in §7 for the cases computed, and the SU(3)-Heisenberg interacting count 60 of §4 is resolved as a representation-theory exercise in §8.

## §7 No intermediate at the computed cases: the non-product part is translation-invariant

§6 leaves open whether a structure between the strict per-site product and a fully global operator, a translation-invariant non-product mirror, can beat the product cap. At every case computed it can, and it goes all the way: **translation invariance alone recovers the entire ceiling.**

The palindrome intertwiner W·L_D = (−L_D − 2Nγ)·W forces W to be block-anti-diagonal in the Hamming grading (it sends rung h to rung N − h), so no intertwiner has rank above the §3 ceiling Σ_h min(c_h, c_{N−h}). Impose translation invariance, [W, T] = 0 with T the cyclic site shift: the allowed entries collapse into T-orbits, one free coefficient each. Give every orbit an integer coefficient and take the rank over the prime field GF(p). Reduction mod p can only lower a rank, so rank_p ≤ rank over ℚ ≤ the maximal translation-invariant rank ≤ the ceiling, and a rank mod p equal to the ceiling settles all three at once. At (d, N) = (3, 2), (4, 2), (3, 3):

  **translation-invariant rank = the ceiling** (54, 128, 378), strictly above the product cap (36, 64, 216).

So at these three cases the gap of §6 is **non-product, not non-translation-invariant**: the whole non-local part of the partial palindrome lives in the translation-invariant sector, and what the strict per-site product cannot do, a homogeneous entangling mirror does completely. The difference sits on the self-complementary rung (h = N/2, for N = 2 the 2d(d² − d)-dimensional h = 1 block, 36 at d = 3): a (1, 1) product reaches at most 2d² of those coherences there (d·d in each direction), a (0, 2) product at most d(d² − d) (one direction only), while the translation-invariant intertwiner reaches all of them. The hierarchy at these cases is two-tiered, product and then translation-invariant = ceiling, with no intermediate layer and the d² − 2d = 0 trunk as the gate between them. Whether translation invariance reaches the ceiling at every d and N is open; the argument above is a rank computation per case, not a general construction.

This completes the inversion of the [golden-router story](PROOF_CEILING_GOLDEN_ROUTER.md). There the per-term test *suspected* non-locality and the window sums (a translation-invariant reading) dissolved it; here the strict product is *provably* capped, and again the next structure up, translation invariance, recovers everything at every case computed. Both seams say one thing: the apparent non-locality is the home of a **homogeneous** structure, and only the strict per-site product is too rigid to see it.

The verification anchor is [`simulations/qudit_ti_intermediate.py`](../../simulations/qudit_ti_intermediate.py) (self-validating: the translation-invariant rank mod p equals the ceiling at (3, 2), (4, 2), (3, 3) and is full at d = 2, the product cap is strictly below at d = 3, 4, and a mutation that replaces a translation-invariant rank by the product cap must fail its gate). The typed claim `QuditProductMirrorCap` recomputes the rank mod p at (2, 2), (3, 2), (4, 2) in its battery and at (3, 3) in its tests.

## §8 The SU(3)-Heisenberg count, decoded

§4 found the interacting count H-dependent and gave 60/81 for the SU(3) Heisenberg without saying why 60. Representation theory says why. SU(3) does not survive the dephasing (the single-site jumps break the global symmetry, so ⟨Q⟩ is not an SU(3) class function), but it organizes the skeleton, and the count falls out of three layers.

**The operator space as an SU(3) representation.** The two-qutrit Hilbert space is 3⊗3 = 6 ⊕ 3̄ (symmetric sextet, antisymmetric antitriplet), so the Liouville space End(ℂ³⊗ℂ³) = (6⊕3̄)⊗(6̄⊕3) decomposes into the irreps 1 (×2), 8 (×4), 27, 10, 1̄0̄ (dimensions 2 + 32 + 27 + 10 + 10 = 81, verified by the quadratic Casimir on coherence space). The SU(3) Heisenberg H = Σ_a λ_a⊗λ_a is a Casimir difference, constant on each Hilbert sector: E(6) = +4/3, E(3̄) = −8/3, a single gap Δ = E(6) − E(3̄) = 4J. So L_H has exactly three eigenvalues, splitting the coherences by which sectors the ket and bra live in: the **intra-sector** (L_H = 0, 45 coherences, 6⊗6̄ ⊕ 3̄⊗3) and the **inter-sector** (L_H = ±iΔ, 36 coherences, 6↔3̄).

**The dephasing reads the energy split.** Because Re(λ) = −2γ⟨Q⟩ (§4) and ⟨Q⟩ respects the L_H blocks (it commutes with the H-energy grading even while breaking SU(3)), the Hamming average distributes as

  intra: ⟨Q⟩ ∈ {0:6, 1:18, 2:21},  inter: ⟨Q⟩ ∈ {1:18, 3/2:12, 2:6},

reconstructing §4's {0:6, 1:36, 3/2:12, 2:27}. Here `⟨Q⟩` is computed from each right eigenvector by the Hilbert-Schmidt Rayleigh ratio `v†Qv/v†v`; the script checks that its imaginary part is below numerical tolerance and separately prints the complex difference from the inapplicable biorthogonal diagonal `w†Qv` (maximum 0.025008 in this run). The half-integer ⟨Q⟩ = 3/2 (the −3γ rung) is **exactly the inter-sector**: every one of its 12 modes has right-HS inter-sector projector weight 1 within numerical precision. The half-integer Hamming average is the signature of crossing the symmetric/antisymmetric divide; intra-sector coherences carry integer ⟨Q⟩ only.

**The 60.** The palindrome about the center −3γ reflects each mode (⟨Q⟩, Im λ) to (3 − ⟨Q⟩, −Im λ): the Hamming rung must complement and the full frequency, not merely its sign, must conjugate. A global one-to-one tolerance matching of the full complex eigenvalues gives

  intra Q = 1 ↔ Q = 2: 2·min(18, 21) = 36;
  inter Q = 1 ↔ Q = 2 at |Im λ| = 4J: 12;
  inter Q = 3/2 self across ±Δ: 2·min(6, 6) = 12;

total **36 + 12 + 12 = 60**, with 21 unpaired (the 6 immortal ⟨Q⟩ = 0 populations, 3 surplus intra Q = 2, and 12 surplus inter Q = 1). Of the last 12, six sit on the detuned |Im λ| ≈ 3.99875 branch and six are the multiplicity surplus at |Im λ| = 4J. So the SU(3)-Heisenberg's 60 is (SU(3) energy sectors) × (Hamming rungs) × (full frequency conjugation), and the standout is structural: the −3γ rung is the sym↔antisym seam. This is the representation-theory account of one H's count, not a universal law (it leans on H having exactly two Hilbert sectors with the 6/3̄ structure), consistent with §4's H-dependence.

The verification anchor is [`simulations/su3_heisenberg_rep_theory.py`](../../simulations/su3_heisenberg_rep_theory.py) (self-validating: the Casimir multiplicities, the L_H sectors, the right-HS intra/inter ⟨Q⟩ split with the 3/2-is-inter-sector identity, the convention-separating complex-difference control, and the 60 from the census).
