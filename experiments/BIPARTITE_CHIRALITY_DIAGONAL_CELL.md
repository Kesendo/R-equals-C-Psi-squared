# Bipartite-Chirality of the F87 Diagonal Cell (k=3 and k=4)

**Status:** The criterion *soft ⟺ H's hopping graph is bipartite in the dephasing basis* is
verified bit-exact at k=3 (N=4 for all three dephase letters, N=5 for Z) and at k=4 (N=4, all
three letters), with 0 mismatches throughout. The direction **bipartite ⟹ soft is derived**
(the chiral K, modulo the F80 one-sidedness M = −2i(H⊗I), itself bit-exact). The converse
splits by support: **at full support (k=N) it closes** , a Mixed+Mixed pair has only two flip
generators, so it is always bipartite, hence soft (modulo M), which settles F111's blocked
"Mixed+Mixed = soft" , while the **windowed regime (k<N)** is reduced to a validated first-order statement (the degenerate
D̂-block set-asymmetry, c=0 ⟺ bipartite, bit-exact; moment route ruled out), the set-level proof
still open.
**Date:** 2026-05-30
**Regenerate:**
- [`simulations/f87_42_8_bipartite_fullcell.py`](../simulations/f87_42_8_bipartite_fullcell.py) `[N] [letters]` , k=3 criterion over the whole diagonal cell (default N=4, all letters; pass `5 Z` for the N=5 Z check)
- [`simulations/f87_k4_bipartite_bridge.py`](../simulations/f87_k4_bipartite_bridge.py) , k=4 criterion + the F111 template cross-check (N=4)
- [`simulations/f87_bipartite_chiral_witness.py`](../simulations/f87_bipartite_chiral_witness.py) , the three derivation links and the optimal λ↔−λ−2σ pairing residual
- [`simulations/f87_flip_generators.py`](../simulations/f87_flip_generators.py) , the flip-generator count |S|, full support (|S|≤2) vs windows (|S|≥3), and the GF(2) φ⟺bipartite check
**Anchors:** [PROOF_F103 §7](../docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md),
[PROOF_F111](../docs/proofs/PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md),
[ChiralKClaim](../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs),
[F87DiagonalCellBipartiteWitness](../compute/RCPsiSquared.Diagnostics/F87/F87DiagonalCellBipartiteWitness.cs).

## The question

F87 sorts a Pauli pair under single-letter dephasing into truly / soft / hard by whether the
Liouvillian spectrum pairs each λ with −λ−2σ. Inside the *diagonal Klein cell* (the cell whose
Klein index matches the dephase letter: Z → (0,1), X → (1,0), Y → (1,1)), PROOF_F103 §7 found a
classical reading of soft vs hard. Read H in the dephasing letter's eigenbasis, where it becomes
a real hopping matrix; let G_H be its graph (basis states as nodes, nonzero off-diagonal entries
as edges). The pair is soft exactly when G_H is **bipartite**, i.e. when H admits a chiral
sublattice K (diagonal in that basis, KHK = −H), the same K = diag((−1)^sublattice) of AZ class
BDI that [ChiralKClaim](../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs) carries.

Two questions this experiment answers:

1. Does that criterion hold beyond the k=3 body count where it was found?
2. [F111](../docs/proofs/PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md) has its own k=4 rule (hard ⟺
   at least one term is a pure-D template), with a derivation that stalled on the converse,
   "Mixed+Mixed = soft." Is F111's rule the same criterion seen at k=4, and does the bipartite
   reading say anything about that blocker?

## Method

For each dephase letter, enumerate the diagonal cell's k-body Pauli templates, form the
y-parity-homogeneous unordered pairs (the F87 dissipator-resonance window), and for each pair:

- compute the **actual** class with the canonical classifier (`PauliPairTrichotomy.Classify` in
  C#, `classify_pauli_pair` in the Python framework, both the same greedy λ↔−λ−2σ pairing);
- independently compute the **predicted** class from the bipartite criterion: build H via the
  sliding-window k-body builder, rotate into the dephase-letter eigenbasis (Hadamard for X, the
  Y→Z rotation for Y, identity for Z), and 2-colour G_H (a zero diagonal plus no odd cycle).

Then count mismatches between *soft* and *bipartite*. At k=4 the bridge script additionally
counts hard pairs carrying **no** pure-D template, the direct F111 cross-check.

## Results

### k=3 (the F103 cell), N=4 and N=5

```
N=4   Z-deph:  hard=50  soft=26   mismatches(soft⟺bipartite)=0
N=4   X-deph:  hard=50  soft=26   mismatches=0
N=4   Y-deph:  hard=50  soft=26   mismatches=0
N=5   Z-deph:  hard=50  soft=26   mismatches=0
```

The cell counts are N-independent (the pair set is alphabet-only); the criterion holds bit-exact
at both N.

### k=4 (the F111 cell), N=4

```
Z-deph:  hard=228  soft=828   mismatches=0   hard-without-template=0
X-deph:  hard=228  soft=828   mismatches=0   hard-without-template=0
Y-deph:  hard=228  soft=828   mismatches=0   hard-without-template=0
```

The 228 hard matches F111's frozen count. `hard-without-template=0` says every hard pair carries
a pure-D template, and the 0 mismatch says soft ⟺ bipartite here too.

### The derivation links (witness, N=4, Z-deph)

```
SOFT  XXZ+ZXX :  M=-2i(H⊗I) True | KHK=-H True | (K⊗I)(-i{H,.}-D)(K⊗I)=-L True
                 marginals: Re-about-(-σ) 3.9e-15,  Im-about-0 1.7e-14
                 optimal |λ↔-λ-2σ| residual = 2.7e-14   (paired)
HARD  XXZ+XZX :  M=-2i(H⊗I) True | chiral K does NOT exist (odd cycle)
                 marginals: Re 5.5e-2,  Im 3.4e-14        optimal residual = 1.6e-1
HARD  ZZZ+XXZ :  M=-2i(H⊗I) True | chiral K does NOT exist (template lifts diagonal)
                 marginals: Re 8.0e-2,  Im 4.8e-14        optimal residual = 1.0e-1
```

The frequency marginal {Im λ} stays mirror-symmetric even when hard (Lindbladian spectrum is
conjugation-closed); the decay marginal {Re λ} is what breaks. The break is genuine: the best
possible pairing leaves a residual of order 10⁻¹, not a near-miss.

### Why no odd cycle at full support: the flip-generator count

In the dephasing basis each Pauli term acts as a single bit-flip mask (its off-diagonal X/Y
positions), so H's hopping graph is the Cayley graph of the distinct edge masks S = {a⊕b}, and it
is bipartite iff a **linear** φ over 𝔽₂ sends every mask to 1 (that φ *is* the chiral K). At full
support (k=N) each term places once, so a Mixed+Mixed pair has only **|S| ≤ 2** masks, and two
nonzero 𝔽₂ vectors always admit such a φ. Measured
([`f87_flip_generators.py`](../simulations/f87_flip_generators.py)):

```
FULL     k=N=3:   max|S|=2   all bipartite   soft=42   hard=0    (φ-exists ⟺ bipartite: 0 mismatches)
FULL     k=N=4:   max|S|=2   all bipartite   soft=828  hard=0    (0 mismatches)
WINDOWS  k=3<N=4: max|S|=4   not all bip.    soft=26   hard=16   (0 mismatches)
```

A third mask, hence a possible odd cycle, appears only when a term is placed at more than one
**window**, i.e. when k < N. The contrast is sharpest at body count k=3: full support at N=3 (all
soft) versus windowed at N=4 (16 hard) , the same body count, opposite behaviour. N=3's long-noted
specialness here is simply that k=3 is already full support, so there are no windows and no odd
cycles.

## What it means

**1. The criterion is k-universal across the two body counts tested.** soft ⟺ bipartite, 0
mismatches, k=3 and k=4, all three dephase letters. The mechanism PROOF_F103 §7 derived is not a
k=3 accident.

**2. At k=4 only one of the two k=3 obstructions fires.** At k=3, hardness has two sources: a
pure-D template lifting the diagonal (rule a), or an odd hopping cycle from opposite
position-parity (rule b). At k=4 (full support, k=N), `hard-without-template=0` says only the
diagonal-lift survives, there are no odd-cycle hard pairs. So at k=4 the bipartite criterion
reduces to *non-bipartite ⟺ has a template*, which is exactly F111's rule, now with a mechanism:
the template is diagonal, it lifts H's diagonal, and a diagonal entry cannot sign-flip under a
diagonal K.

**3. F111's blocked converse closes at full support (modulo M).** "Mixed+Mixed = soft" decomposes
into three links:

- **(i) Mixed+Mixed ⟹ zero diagonal.** Derived: a non-template term has at least one off-diagonal
  letter, so it is off-diagonal; two of them give a zero-diagonal H.
- **(ii) zero diagonal ⟹ bipartite.** **Derived at full support** by the flip-generator count
  above: k=N gives |S| ≤ 2 masks, and two nonzero 𝔽₂ vectors always admit the linear φ. (Genuinely
  full-support-specific: at k=3 *with windows* (N=4) it is false , XXZ+XZX is zero-diagonal with a
  3-cycle, hard.)
- **(iii) bipartite ⟹ soft.** Derived (the chiral K), modulo the F80 one-sidedness M = −2i(H⊗I),
  bit-exact.

F111's cell is k = N = 4, full support, so all three links hold: **"Mixed+Mixed = soft" is derived
there, modulo M** , the converse that stalled F111 across three derivation paths, now closed, with
the chiral K exhibited constructively as the separating functional φ. What stays open is the
*windowed* converse (k < N), where a third mask can close an odd cycle.

## The windowed converse as a first-order perturbation (the clock's tick)

The windowed converse (non-bipartite ⟹ hard) resisted a direct proof , ruling out every operator
that could pair the spectrum. The clock reframes it. A mode winds as e^(λt) = e^(−αt)·e^(iωt), and
the tick is γ. At **γ = 0** (the Takt stopped) L = −i[H,·] is purely imaginary, hence symmetric
about 0 = −σ: **soft**. So the break is tick-driven, and it is **first order**: the pairing
residual is γ-linear, residual = c·γ as γ → 0 (e.g. the soft pair stays at machine zero for all γ;
the hard pairs grow linearly).

So the converse becomes a first-order perturbation question. At γ = 0 the eigenmodes are the
H-eigenbasis coherences |E_a⟩⟨E_b| at −i(E_a−E_b); the dephasing D̂ = Σ_l(Z_l·Z_l − I) shifts them,
and because H's frequencies are degenerate the shifts are the **eigenvalues of D̂ restricted to each
degenerate-frequency block**, not its diagonal (the naive diagonal formula even gives the soft pair
a spurious asymmetry). The degenerate first-order spectrum reproduces the actual break **bit-exact**,
under a stable optimal-transport measure (the min-sum-assignment `.max()` is unstable on
near-degenerate spectra; that instability, not the physics, produced the earlier ragged 2.0/1.6/3.2):

```
SOFT  XXZ+ZXX:  c(degenerate 1st order) = 0.0000 = c(measured)   ✓
FLUX  IXY+XIY:  c = 0.1191 = 0.1191   ✓
REAL  XXZ+XZX:  c = 0.2559 = 0.2559   ✓     (μ-vs-L mismatch ~10⁻⁴/γ, so the shifts are exact)
```

So the converse is a validated first-order statement: **c = the set-asymmetry of the degenerate
first-order D̂-block spectrum about −σ, and c = 0 ⟺ bipartite.**
([`f87_first_order_degenerate.py`](../simulations/f87_first_order_degenerate.py),
[`f87_break_gamma_scaling.py`](../simulations/f87_break_gamma_scaling.py).)

**The low moments are blind (the early dead-end Phase B reopened).** In the computational basis D̂ is diagonal,
D̂|i⟩⟨j| = −2·Hamming(i⊕j)·|i⟩⟨j| (the Absorption Theorem), which made a moment proof tempting:
soft ⟺ Tr((L+σ)^{odd}) = 0. But for the genuine-cycle pairs the scout tested, the first-order odd moments **vanish for soft and hard alike**
(Tr((L+σ)³) = Tr((L+σ)⁵) = 0 to first order, verified for both;
[`f87_moment_condition.py`](../simulations/f87_moment_condition.py); a lifted diagonal is different,
it breaks already at first order, p₃ = 9216·γ for IIZ+IZI). For cycles the low moments do not see the
break; it is *moment-invisible at low order*. The hardness is a **set-pairing asymmetry** (the optimal-transport
distance between Spec(L) and Spec(−L−2σ)), finer than any low moment. At low order this read as the moment route being dead, forcing a set-level (combinatorial) statement about the degenerate D̂-block spectrum. But it is only blind at low order: Phase B (below) shows the first NONVANISHING odd moment, at m\* = 2ℓ + deg, is a positive monomial, so the moment route closes the converse after all, built on that same combinatorial structure. That set-level statement is independently derived in PROOF_F103 §7.5 (2026-06-04): the odd cycle obstructs the chiral functional that would supply the gain channel's reflection-floor mode; the first-order reduction it rested on is itself closed in §7.6 (degenerate PT + analyticity), so the converse is fully derived modulo standard perturbation theory.

## Honest status

- **Derived:** bipartite ⟹ soft (modulo M = −2i(H⊗I), verified bit-exact); and, at full support
  (k=N), Mixed+Mixed ⟹ bipartite (the |S| ≤ 2 flip-generator / linear-φ argument). Together these
  close F111's "Mixed+Mixed = soft" modulo M.
- **Derived modulo the first-order premise (2026-06-04, PROOF_F103 §7.5):** the windowed converse
  non-bipartite ⟹ hard. The ω=0 first-order block is the gain channel Σ_l Z_l(·)Z_l restricted to
  H's commutant; its +N population-Perron mode is always present (Σ_l Z_l² = N·I), and soft ⟺ the
  −N reflection-floor is also attained ⟺ an anti-diagonal commutant element A = FK exists (F = X^⊗N,
  using FHF = −H) ⟺ bipartite ⟺ no K3 triangle. ω=0 is decisive (the +N mode pairs only at ω=0).
  Separately the operator-search is dissolved: any palindromizer forces spec(L) = spec(−L−2σ), so no
  non-chiral similarity escapes (verified 236 pairs, N=4 all letters + N=5 Z). The low-order moment route stays
  ruled out (Phase B, below, closes the converse at the first nonvanishing moment). The first-order reduction is itself closed
  (§7.6, 2026-06-04): L₀ = −i[H,·] is normal so degenerate PT is exact (O(γ) shifts = M_ω eigenvalues to
  1e-9, N=4/5), and analyticity (char-poly coefficients polynomial in γ) turns the nonzero first-order
  break into a break for all but finitely many γ, so non-bipartite ⟹ hard at the generic physical γ.
  Adversarially stress-tested (5 holes ruled out, N=4/5). So the converse is fully derived modulo
  standard perturbation theory.
- **Typed split (2026-06-08):** the genericity result (non-bipartite ⟹ hard for all but finitely many
  γ) is Derived, carried by F87DiagonalCellBipartiteWitnessSet; the remaining all-γ closure including the
  physical operating point is isolated as the Tier1Candidate WindowedConverseAllGammaClaim (upgraded to
  "proven modulo R-deg + R-sign" in the Phase B bullet below).
- **Proven modulo R-deg + R-sign (Phase B, 2026-06-09,
  [PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md](../docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md)):** recentre
  the Liouvillian as M = A + γQ; two involutions 𝓕 = F⊗F and R = I⊗F force every surviving odd power-sum word
  to have #A_L, #A_R, #Q all odd, giving the threshold #A ≥ 2ℓ (ℓ the unsigned odd-girth) and a second,
  independent re-proof of bipartite ⟹ soft (no odd word survives). The first nonvanishing odd power-sum closes
  to a **positive monomial** c·γ^deg (c > 0, deg ∈ {1, 3}), which has no positive real root, so non-bipartite is
  **hard at every operating point γ > 0**, upgrading "all but finitely many γ" to "all γ > 0". This spine (the
  sign table, the threshold, the soft re-proof, the deg-1 positivity closed form P_{3,1} = 6·4^N·Σ_l c_l²
  over H's single-site-Z Pauli coefficients, manifestly ≥ 0)
  is Tier1Derived (`WindowedConverseThresholdClaim`); the monomial-and-positive lemma stays Tier1Candidate
  (`WindowedConverseAllGammaClaim`), gated on the two residuals R-deg (the genuine-cycle degree lift
  m\* = 2ℓ + 3) and R-sign (the genuine-cycle coefficient P_{m\*,3} > 0), both verified bit-exact cell-wide at
  N = 4 and at N = 5 / N = 6 representatives.
- **Retraction (the "flux moment-invisible exception" was a discarded-Im H bug).** An earlier reading treated
  the complex-H flux pair (odd #Y, X…Y bond, Gaussian-integer H) as a genuine all-γ exception, on the grounds
  that its signed third moment cancels, (H³)_{ii} = 0. That cancellation is real but irrelevant, and the
  earlier read mis-handled the imaginary part of H. The two-reflection argument is governed by
  *unsigned* index-trajectory existence: Q is diagonal, so a trace's support is the unsigned three-walk count
  (|H|³)_{ii} = 6 > 0, which is immune to the signed XX+YY cancellation. The flux pairs therefore obey the
  identical monomial law (the IXY+XIY flux pair has p_9 = 589824·γ³, pure γ³, hard ∀γ > 0), with no
  exception. The sign table is checked bit-exact on this complex-H pair, so the converse covers flux pairs the
  chiral-K route never had to face.
- **Not yet tested:** windowed k=4 (N > 4) and k=5. The windowed converse is proven modulo R-deg + R-sign;
  full support is settled.

## Links

- Criterion and proof: [PROOF_F103 §7](../docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md)
- k=4 rule it unifies with: [PROOF_F111](../docs/proofs/PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md)
- The chiral symmetry behind it: [ChiralKClaim](../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs) (AZ class BDI, KHK = −H)
- The F80 one-sidedness it rests on: [PROOF_F80](../docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md) (M = ±2i times a Hamiltonian object for Π²-odd cell terms)
- Typed witness: [F87DiagonalCellBipartiteWitness](../compute/RCPsiSquared.Diagnostics/F87/F87DiagonalCellBipartiteWitness.cs) (Tier1Derived) + [BipartiteChirality](../compute/RCPsiSquared.Diagnostics/F87/BipartiteChirality.cs) primitive
