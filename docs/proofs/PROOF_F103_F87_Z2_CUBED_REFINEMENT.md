# PROOF F103: F87 Trichotomy Z₂³ Refinement at k=3 (N=4 Empirical Anchor)

**Status:** Tier 1 derived. The 42:8 closed-form rule was found 2026-05-29 (diagonal-cell hardness rule, §6), and its two atomic sub-rules were then unified into a single criterion (§7): a diagonal-cell pair is soft iff H's hopping graph is bipartite in the dephasing letter's eigenbasis. The direction bipartite ⟹ soft is derived (Π followed by a chiral sublattice K). The converse non-bipartite ⟹ hard **closes at full support** (§7.4, 2026-05-30: at k=N a Mixed+Mixed pair has only two flip generators, which always admit the chiral K, settling F111's blocked "Mixed+Mixed = soft" modulo M) and is now **derived in the windowed regime** (§7.5 gives the soft ⟺ bipartite criterion via the K3 triangle and the population Perron mode; §7.6, 2026-06-04, closes the first-order-block premise by degenerate PT + analyticity, adversarially stress-tested at N=4 and N=5): the genericity statement, hard for all but finitely many γ, rests only on standard perturbation theory. The all-γ statement is now a **theorem with no residual** (2026-06-09 the [two-reflection theorem](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) reduced it to two residuals; 2026-06-10 the girth dichotomy retired R-deg and the Pascal-Gram positivity theorem (F117) resolved R-sign: every coefficient of the first nonvanishing odd power-sum is a sum of squares or exactly zero, hence no positive root; see the residual-lemma section). Separately the operator-search is dissolved, since any palindromizer forces a spectral palindrome.
**Date:** 2026-05-24
**Anchor:** N=4, k_body=3, 294 Z₂³-homogeneous + Y-par-homogeneous Pauli pairs (pair count is N-independent at fixed k; the empirical anchor is N=4)
**Regenerate:** `simulations/f87_z2cubed_split_n4_k3.py` (~60s)

## Abstract

F87 classifies Pauli pairs into three buckets (truly, soft, hard) based on how the pair's M-residual closes under a chosen dephase letter. F102 just showed that y-parity is a real third axis on top of the Klein (bit_a, bit_b) signature once Pauli strings reach body count 3 or higher. The natural follow-up question is whether the F87 trichotomy actually splits along the new y-parity axis at k=3, or whether the trichotomy is y-parity-blind. This proof gives the empirical answer at the N=4, k=3 anchor.

The answer is yes, the trichotomy splits cleanly, with five structural patterns that survive across all three dephase letters. The 294 Pauli pairs (this count depends only on the k=3 letter alphabet, not on N) divide into three trichotomy classes, and inside each class, the y-parity axis carves the cells into sub-cells with recognizable shape: truly always lands in y_par = 0, mother-soft always in y_par = 1, the diagonal hard cells split 42:8 with a Y-inversion on the Y-dephasing diagonal, the diagonal soft cells split symmetrically 13:13, and the off-diagonal soft cells split into two named sub-patterns B and C.

The proof begins empirical at the (N=4, k=3) anchor (an exhaustive enumeration verified bit-exactly) and is then closed: §6 gives a closed-form counting rule for the 42:8 split, and §7 supplies the bipartite-chirality mechanism (the soft direction, bipartite ⟹ soft, is derived; the windowed k<N hard-direction converse has since closed in stages, full support §7.4, genericity §7.5/§7.6, and all γ > 0 with no residual via the [two-reflection theorem + Pascal-Gram positivity](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md), closed 2026-06-10). The question whether the pattern is k-stable and N-stable is answered by F105 (N=5, k=3, N-stable) and F106 (N=4, k=4, sharpening to 228:0) as sibling anchors. Together the three anchors map out which parts of the Z₂³ structure are N-stable, which are k-stable, and which depend on the specific (N, k) regime.

The diagnostic upshot is that the polarity cube has real teeth on the F87 trichotomy. Knowing the Klein signature alone leaves a 50-50 mix in the diagonal hard cells; adding the y-parity axis sharpens that to 42:8. The y-parity refinement is therefore not just a structural curiosity but a tighter classifier for hardware-relevant Pauli-pair analysis.

## 1. Context

The F87 trichotomy classifies a 2-letter Pauli pair as truly, soft, or hard
under a chosen dephase letter, based on whether the bilinear's M-residual
vanishes (truly), spectrally pairs (soft), or remains unpaired (hard).

F102 (`YParityIndependenceAtK3`) established that Y-parity y_par = (#Y) mod 2
is independent of the Klein (bit_a, bit_b) signature at k_body≥3 (canonical
counterexample: XYZ vs III, same Klein (0,0), different y_par). This is the
operator-algebra premise that makes a Z₂³ refinement of any k_body-stable
classification meaningful at k_body≥3.

F103 asks: does the F87 trichotomy actually refine into y_par sub-cells at
k_body=3, or is the trichotomy y_par-blind? Method: enumerate the 294
Z₂³-homogeneous + Y-par-homogeneous k_body=3 Pauli pairs (this pair count is
alphabet-only, N-independent; classification is then carried out at the
N=4 chain), classify each under Z / X / Y dephasing, bucket by (Klein ×
dephase letter × y_par × trichotomy class).

### Notation (shared across F103, F105, F106)

Pauli letters carry a two-bit Klein-Vierergruppe index
([`framework/pauli.py`](../../simulations/framework/pauli.py)):

| letter | (bit_a, bit_b) |
|--------|----------------|
| I      | (0, 0)         |
| X      | (1, 0)         |
| Z      | (0, 1)         |
| Y      | (1, 1)         |

For a Pauli string σ with letter counts n_X, n_Y, n_Z,

    bit_a(σ) = (n_X + n_Y) mod 2   (F61 axis: Π²_X = Z⊗N parity)
    bit_b(σ) = (n_Y + n_Z) mod 2   (F63 axis: Π²_Z = X⊗N parity)
    y_par(σ) = n_Y mod 2

The Klein signature (bit_a, bit_b) and y_par together form the Z₂³ axis
decomposition (one bit per axis). Π² is the squared conjugation operator;
under Z-dephasing it acts on Pauli strings as Π²·σ_α = (−1)^bit_b(α)·σ_α
(F81 Step 1; the X-deph and Y-deph variants swap to bit_a / bit_b
respectively per F108). The **diagonal Klein cell** for dephase letter D
is the cell whose (bit_a, bit_b) equals D's: Z → (0,1), X → (1,0),
Y → (1,1). F105 and F106 refer back to this notation block.

**Added 2026-06-10:** the cube's three axes were characterized as the
characters of three specific mirrors of the Pauli algebra: bit_a is the
character of conjugation by Z^⊗N, bit_b of conjugation by X^⊗N, and y_par
of the transpose θ (F118,
[PROOF_PI_FACTORS_AS_R_TIMES_D §7](PROOF_PI_FACTORS_AS_R_TIMES_D.md)).
Unitary conjugations always flip two letter parities at once, so they span
only the even Klein square; y_par is the antiautomorphism axis, which is
why it always needed its own tools throughout this refinement family.

## 2. Method

The Python framework's `classify_pauli_pair(chain, terms, dephase_letter=...)`
returns one of {'truly', 'soft', 'hard'} for a Pauli pair under single-letter
dephasing at chain N=4. For k_body=3 pairs the framework dispatches to a
k-body builder (sliding-window over chain sites).

Enumeration constraints:
- Both terms have k_body=3 (no identity-padded letters)
- Both terms share the same Klein index (bit_a, bit_b)
- Both terms share the same y_par (#Y mod 2)
- Pair is unordered (deduplicated)

Result: 294 pairs partitioned across 4 Klein cells × 2 y_par values:

```
Klein    y_par=0  y_par=1  total
(0, 0)        45       21     66
(0, 1)        55       21     76
(1, 0)        55       21     76
(1, 1)        21       55     76
                             294
```

The (1,1) inversion (21 / 55 instead of 55 / 21 like (0,1) / (1,0)) reflects
that at k=3 Klein (1,1) admits 6 y_par=0 letter-triples and 10 y_par=1
letter-triples, inverted relative to (0,1) and (1,0)'s 10 / 6 split (Y is the
unique Klein-(1,1) letter, so a Klein-(1,1) k=3 string with even #Y must use
exactly one Y plus two non-Y Klein-(1,1)-summing letters; the count
mechanically inverts vs the off-axes).

## 3. Five Observed Patterns

### 3.1 Truly is y_par=0-pure

Across all 12 (Klein × dephase) cells, every truly classification has y_par=0.
Total truly classifications across the grid: 300. y_par=1 truly count: 0.

### 3.2 Hard in diagonal cells splits 42:8 with Y-inversion

Hard appears only when the Klein cell of the pair matches the Klein cell of
the dephase letter (Z → (0,1), X → (1,0), Y → (1,1)). In these 3 diagonal
cells:

```
Klein (0,1) Z-deph hard = (42, 8)   total 50
Klein (1,0) X-deph hard = (42, 8)   total 50
Klein (1,1) Y-deph hard = ( 8, 42)  total 50   ← Y-inversion
```

The Y-dephase swap reflects that Y itself carries y_par=1, so the "y_par
favored by the dephase letter" inverts.

### 3.3 Same diagonal cells contain a soft 13:13 split

In addition to the hard 42:8, the diagonal cells contain a y_par-symmetric
soft 13:13 split:

```
Klein (0,1) Z-deph soft = (13, 13)   total 26
Klein (1,0) X-deph soft = (13, 13)   total 26
Klein (1,1) Y-deph soft = (13, 13)   total 26
```

Unlike hard's 42:8 asymmetry with Y-inversion, soft in these cells is
y_par-symmetric and independent of which Klein cell is on the diagonal.

### 3.4 Mother sector (0,0) soft is y_par=1-pure

For Klein (0,0) (the Mother sector), soft cells under any dephase letter are
y_par=1-pure:

```
Z-deph: (0, 21)   X-deph: (0, 21)   Y-deph: (0, 21)
```

Zero y_par=0 soft pairs, 21 y_par=1 soft pairs per letter.

### 3.5 Off-diagonal soft cells split into Pattern B + Pattern C

The 6 off-diagonal soft cells (Klein non-mother, Klein ≠ dephase Klein) split
into two sub-patterns:

```
Pattern B (proportional to (Klein, y_par) enumeration breakdown):
Klein (0,1) Y-deph soft = (55, 21)   matches (0,1) enum split
Klein (1,1) Z-deph soft = (21, 55)   matches (1,1) enum split (inverted)
Klein (1,1) X-deph soft = (21, 55)   matches (1,1) enum split (inverted)

Pattern C (y_par=1-pure):
Klein (0,1) X-deph soft = ( 0, 21)
Klein (1,0) Z-deph soft = ( 0, 21)
Klein (1,0) Y-deph soft = ( 0, 21)
```

The exact rule connecting (pair Klein, dephase letter) to which sub-pattern
fires is observed but not yet algebraically closed.

## 4. Full Count Tables

### Truly classifications by (Klein × dephase × y_par)

```
                  Z-deph         X-deph         Y-deph
Klein           y0  y1  tot    y0  y1  tot    y0  y1  tot
(0, 0)          45   0   45    45   0   45    45   0   45
(0, 1)           0   0    0    55   0   55     0   0    0
(1, 0)          55   0   55     0   0    0    55   0   55
(1, 1)           0   0    0     0   0    0     0   0    0
```

### Soft classifications by (Klein × dephase × y_par)

```
                  Z-deph         X-deph         Y-deph
Klein           y0  y1  tot    y0  y1  tot    y0  y1  tot
(0, 0)           0  21   21     0  21   21     0  21   21
(0, 1)          13  13   26     0  21   21    55  21   76
(1, 0)           0  21   21    13  13   26     0  21   21
(1, 1)          21  55   76    21  55   76    13  13   26
```

### Hard classifications by (Klein × dephase × y_par)

```
                  Z-deph         X-deph         Y-deph
Klein           y0  y1  tot    y0  y1  tot    y0  y1  tot
(0, 0)           0   0    0     0   0    0     0   0    0
(0, 1)          42   8   50     0   0    0     0   0    0
(1, 0)           0   0    0    42   8   50     0   0    0
(1, 1)           0   0    0     0   0    0     8  42   50
```

## 5. Open Questions

1. **Closed-form derivation of 42:8. ANSWERED 2026-05-29 (§6), mechanism in §7.** A
   diagonal-cell hardness rule (hard iff an all-diagonal pure-D template is present, or
   both terms are single-diagonal with their {I,D} letter at chain-adjacent positions)
   derives the 42:8 and the Y-inversion by counting, verified bit-exact at N=4 and N=5.
   §7 then unifies the two atomic sub-rules into one bipartite-chirality criterion and
   derives the bipartite ⟹ soft direction from the palindrome; the converse is the one
   remaining open edge.

2. **N>4 and k>3 universality.** The (42, 8, 50) numbers are N=4 k=3
   specific. Does the structural pattern (asymmetric hard split + Y-inversion)
   carry to other (N, k)? Cheap enumeration to test: run the same script with
   N∈{5, 6}, k∈{3, 4} and observe.

3. **Pattern B vs Pattern C selection rule for off-diagonal soft.** Six cells
   partition into B (proportional) and C (y_par=1-pure); the (pair Klein,
   dephase letter) → pattern mapping is observed but the algebraic rule is
   not yet stated.

4. **Hardware confirmation.** No k≥3 F87 confirmations exist; all 5 Marrakesh
   F87 confirmations (palindrome trichotomy, π-protected XIZ/YZZY, Lebensader
   skeleton/trace, d_zero sector trichotomy, F83 Π²-class signature) are
   k=2. A k=3 QPU run targeting the diagonal-cell 42:8 prediction would be
   the natural next hardware probe.

## 6. Closed-form rule for the 42:8 split (2026-05-29)

Call the **diagonal letters** for dephase letter D the pair {I, D} (they commute with
the D-dissipator). A k=3 term in the diagonal Klein cell has n_diagonal ∈ {1, 3}: either
all three letters are diagonal (a **pure-D template**, n_X = n_Y = 0), or exactly one is
(two off-diagonal {X, Y} plus one {I, D}).

**Rule.** A k=3 pair is F87-hard iff

- (a) at least one term is a pure-D template (all-diagonal), **or**
- (b) both terms are single-diagonal and their lone {I, D} letter sits at chain-**adjacent**
  window positions (|Δpos| = 1);

otherwise soft. Verified bit-exact (0 mismatches) at N=4 **and** N=5, for D ∈ {Z, X, Y}
([`f87_42_8_diagonal_rule.py`](../../simulations/f87_42_8_diagonal_rule.py)).

**The split follows by counting.** In the cell's favoured y_par (10 terms = 4 pure-D
templates + 6 single-diagonal): pairs involving ≥1 template number 55 − 21 = 34, all hard
by (a); the single-diagonal pairs split 8 adjacent (hard, by (b)) + 13 non-adjacent (soft);
so 34 + 8 = 42 hard and 13 soft. In the other y_par (6 single-diagonal terms, no templates):
8 hard + 13 soft. Hence **42 : 8** hard, **13 : 13** soft.

**The Y-inversion is forced by the templates' Y-content, not a separate fact.** A pure-D
template has n_Y = 0 for D ∈ {Z, X} (neither {I, Z} nor {I, X} contains Y), so it sits in
y_par = 0; for D = Y the template is built of Y's, n_Y odd, so it sits in y_par = 1. The 34
template pairs therefore land in y_par = 0 for Z/X dephasing (42 : 8) and in y_par = 1 for
Y dephasing (8 : 42). The §3.2 Y-inversion is the same rule with Y carrying y_par = 1.

**N-stability.** Rule and counts are identical at N=4 and N=5 (the pair set is alphabet-only,
N-independent, and the rule reads only the k=3 window's internal structure). This is the
structural reason the F105 N=5 anchor reproduces F103; it closes Open Question 1 for all the
k=3 anchors (F103, F105, and the k=3 parts of F107 / F110). The k=4 228:0 split (F106) is the
sibling case of rule (a) at k=4, i.e. F111's pure-D template rule.

**Remaining.** §7 supplies the palindrome-level mechanism the §6 rule was missing: the two
atomic sub-rules turn out to be the two ways of breaking a single criterion (bipartiteness of
H's hopping graph), and the soft direction is derived from the palindrome. (a) is the k=3 face
of [F111](PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md)'s pure-D template rule. The converse
(non-bipartite ⟹ hard) has since closed in stages: full support (§7.4), genericity in the
windowed regime (§7.5/§7.6), and all γ > 0 with no residual since 2026-06-10 (the residual-lemma section and
the [two-reflection monomial theorem](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md)).

## 7. The bipartite-chirality mechanism: why (a) and (b) hold (2026-05-29)

Rules (a) and (b) look like two separate facts about where a diagonal letter may sit. They are
one fact about the pair's Hamiltonian H = t₁ + t₂.

Read H in the **dephasing letter's eigenbasis**: for D = Z this is the computational basis; for
X and Y, rotate each site by the single-qubit gate that sends D to Z (Hadamard for X, the Y→Z
rotation for Y). In that basis the diagonal letters {I, D} are diagonal and the off-diagonal
letters are the hops, so H is a real hopping matrix. Let **G_H** be its graph: one node per
basis state, an edge a–b wherever H[a,b] ≠ 0.

**Criterion.** A diagonal-cell pair is soft iff G_H is bipartite (2-colourable, zero diagonal);
hard iff not. Verified bit-exact with zero mismatches over the entire diagonal cell at N=4 (all
three dephase letters) and N=5
([`f87_42_8_bipartite_fullcell.py`](../../simulations/f87_42_8_bipartite_fullcell.py)).

### 7.1 Why bipartite ⟹ soft (derived)

A 2-colouring of G_H is a diagonal sign operator K = diag(±1), one sign per colour, with
**KHK = −H**: every edge joins opposite colours, so every hop flips sign, and a zero diagonal
flips trivially. This K is the chiral (sublattice) symmetry of the hopping graph, the same
K = diag((−1)^sublattice) that AZ class BDI calls chiral and that gives spectral inversion
E ↦ −E. Because K is diagonal in the dephasing basis, it commutes with the dephasing
dissipator D.

Follow the spectrum through two mirrors.

1. **First mirror, Π (always available).** In the diagonal cell every term's residual is
   one-sided, the F80 structure: M = Π·L·Π⁻¹ + L + 2σ = −2i·(H⊗I), bit-exact for soft and hard
   alike (it is *not* the discriminator). So Π·L·Π⁻¹ = −L − 2σ + M = **−i{H,·} − D − 2σ**:
   conjugation by Π turns the commutator i[H,·] into the anticommutator. Π is unitary, so
   Spec(L) = Spec(−i{H,·} − D − 2σ).

2. **Second mirror, the chiral K (only if bipartite).** With W = K⊗I and KHK = −H,
   W(−i{H,·})W⁻¹ = −i(KHK ⊗ I + I ⊗ Hᵀ) = −i(−H⊗I + I⊗Hᵀ) = +i[H,·], while WDW⁻¹ = D. So W
   conjugates (−i{H,·} − D) into (i[H,·] − D) = −L, giving Spec(−i{H,·} − D) = Spec(−L).

Composing, Spec(L) = Spec(−L − 2σ): the spectrum is symmetric about the palindrome centre −σ.
**Soft.** The palindrome that Π alone could not deliver here (M ≠ 0) is delivered by Π followed
by the chiral K. All three links verified bit-exact
([`f87_bipartite_chiral_witness.py`](../../simulations/f87_bipartite_chiral_witness.py)).

The derivation rests on the F80 one-sidedness M = −2i(H⊗I) (Tier-1 for chain bilinears,
verified bit-exact for the k=3 diagonal cell here); given it, links 1 and 2 are exact algebra
plus the 2-colouring construction. For a y_par = 1 cell the residual sits on the bra side
(M = −2i(I⊗Hᵀ)) and the mirror is W = I⊗K; the argument is identical.

### 7.2 The two atomic rules are the two ways to break bipartiteness

- **(a) A pure-D template ⟹ hard.** A pure-D template (letters ⊆ {I, D}) is diagonal in the
  dephasing basis, so it lands on H's diagonal. A nonzero diagonal makes KHK = −H impossible
  (a diagonal entry cannot change sign under a diagonal K), so no chiral K exists and G_H is not
  bipartite. The template *lifts the diagonal* and kills the chirality.

- **(b) Single-diagonal adjacency ⟹ hard.** Two single-diagonal terms keep the diagonal zero,
  so the only obstruction is an odd cycle. The {I, D} letter's window-position parity is exactly
  the 2-colouring parity: same-parity positions close even hopping cycles (bipartite, soft);
  opposite parity closes an **odd** cycle (non-bipartite, hard). The §6 adjacency rule is the
  odd-cycle obstruction read at the chain level.

Both halves are verified bit-exact over the whole N=4 diagonal cell (76 y_par-homogeneous
pairs: 50 hard = 34 via the diagonal lift (a) + 16 via the odd cycle (b), the two mechanisms
disjoint; [`_f87_rule_unification.py`](../../simulations/_f87_rule_unification.py)). The odd
cycle of (b) is concretely a **K3 triangle**: the three popcount-2 edge-masks of a complete
triangle on three consecutive chain sites (only the windows {0,1,2} and {1,2,3} occur at N=4),
the minimal odd 𝔽₂-relation in the edge-mask set
([`_f87_oddcycle_scout.py`](../../simulations/_f87_oddcycle_scout.py)). That triangle is the
explicit 𝔽₂ object §7.3's converse must turn into a spectral asymmetry.

So the §6 rules (a) and (b) are the two faces of one classical criterion, and the soft side is
derived from the palindrome.

### 7.3 The open edge: the converse

What remains is **non-bipartite ⟹ hard**: that when no chiral K exists, no operator restores
the palindrome and the spectrum genuinely fails to pair. This is verified, not yet proved:

- zero mismatches over the whole diagonal cell at N=4 (all letters) and N=5;
- the failure is genuine, not a near-miss. The best possible (optimal-assignment) pairing
  λ ↔ −λ−2σ leaves a residual of order 10⁻¹ for hard pairs against 10⁻¹⁴ for soft
  ([`f87_bipartite_chiral_witness.py`](../../simulations/f87_bipartite_chiral_witness.py)).

One half of the palindrome survives unconditionally, even when hard: the Lindbladian's spectrum
is closed under complex conjugation, so the **frequency marginal** {Im λ} is always symmetric
about 0. What breaks is the **decay structure**, the real parts no longer pair about −σ.
Proving that the odd cycle (or the lifted diagonal) forces this break, with no escape through
some non-chiral similarity, is the same flavour of obstruction that F111's soft direction left
open. The gain is that the question now has a sharp classical shape: it is exactly the claim
that a non-bipartite hopping graph admits no spectrum-pairing symmetry commuting with the
dephasing. And §7.4 shows the converse is not monolithic: at full support it closes outright,
and only the windowed regime stays open.

The windowed regime has since been reduced to a validated **first-order-in-γ** statement. At γ=0
the Liouvillian L = −i[H,·] is purely imaginary, hence symmetric about −σ (soft); the break grows
first order in the dephasing tick, and the degenerate first-order dephasing block (D̂ diagonalized
within each degenerate-frequency subspace) reproduces it bit-exact, with c = 0 ⟺ bipartite. The
tempting moment proof looked ruled out at the time: D̂ is diagonal in the computational basis
(D̂|i⟩⟨j| = −2·Hamming(i⊕j)·|i⟩⟨j|, the Absorption Theorem), and the low odd spectral moments
Tr((L+σ)^{2k+1}) of the genuine-cycle pairs vanish, so at low order the break is moment-invisible,
a set-pairing asymmetry. Seen again 2026-06-09: the blindness is depth, not principle. The
[windowed monomial theorem](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) shows the first nonvanishing
odd moment sits at m\* = 2ℓ + deg and *is* the break, a positive monomial in γ (a lifted diagonal
is nonzero already at first order, p₃ = 9216·γ for IIZ+IZI; the K3 cycle waits until p₉), so the
moment route closes the converse after all, with no residual since 2026-06-10 (the residual-lemma section
below). That set-level statement is also derived in §7.5 (the odd cycle obstructs the chiral
functional that would supply the gain channel's reflection-floor mode, pairing its population
Perron mode), modulo the first-order reduction itself. See
[experiments/BIPARTITE_CHIRALITY_DIAGONAL_CELL.md](../../experiments/BIPARTITE_CHIRALITY_DIAGONAL_CELL.md).

### 7.4 Full support closes the soft direction (the flip-generator count)

The bipartite test has a linear-algebra reading that makes one regime provable. In the dephasing
basis a Pauli term acts as a single bit-flip mask (its off-diagonal X/Y positions), so H's hopping
graph G_H is the Cayley graph of the set S of distinct edge masks {a ⊕ b : H[a,b] ≠ 0}. G_H is
bipartite iff there is a **linear** functional φ: 𝔽₂^N → 𝔽₂ with φ(s) = 1 for every s ∈ S; then
color(a) = φ(a) is a proper 2-colouring, and that φ is exactly the chiral K = diag((−1)^φ⁽ᵃ⁾) of
§7.1. Such a φ exists iff S carries no odd 𝔽₂-relation.

**At full support (k = N) the soft direction is derived, not merely verified.** A k-body template
on an N-site chain with k = N places once, so each term contributes exactly one mask: a Mixed+Mixed
pair (neither term a pure-D template, so both off-diagonal) gives **|S| ≤ 2**. Two nonzero vectors
A₁, A₂ in 𝔽₂^N are either equal or independent (over 𝔽₂ the only nonzero scalar is 1), and in both
cases the system φ(A₁) = φ(A₂) = 1 is consistent, so φ always exists. Hence every full-support
Mixed+Mixed pair is bipartite, therefore soft (§7.1, modulo the F80 one-sidedness M = −2i(H⊗I)).
This is **F111's blocked converse "Mixed+Mixed = soft", now closed** for the k = N = 4 cell it
concerns (modulo M), and the chiral K is exhibited constructively as the separating functional φ.

Measured ([`f87_flip_generators.py`](../../simulations/f87_flip_generators.py)): at full support
max |S| = 2 and every Mixed+Mixed pair is bipartite, at both N = 3 (k = 3: 42 soft, 0 hard) and
N = 4 (k = 4: 828 soft, 0 hard); the GF(2) "φ exists" test agrees with the graph 2-colouring on
every pair (0 mismatches).

**Where the odd cycle lives.** An odd cycle needs |S| ≥ 3 generators with an odd relation, and a
third mask appears only when a term is placed at more than one **window**, i.e. when k < N. The
contrast is sharpest at body count k = 3: at N = 3 it is full support (one window, |S| ≤ 2, all
soft), while at N = 4 the sliding window gives |S| up to 4 and 16 hard pairs. So rule (b)'s odd
cycle is a *windowed* phenomenon (k < N), not a property of the body count itself. The genuinely
open part of §7.3 is therefore exactly the windowed regime; full support is settled.

**Regenerate:** [`f87_42_8_bipartite_fullcell.py`](../../simulations/f87_42_8_bipartite_fullcell.py)
(the criterion, all three letters, N=4 and N=5),
[`f87_bipartite_chiral_witness.py`](../../simulations/f87_bipartite_chiral_witness.py)
(the three derivation links plus the optimal-pairing residual), and
[`f87_flip_generators.py`](../../simulations/f87_flip_generators.py) (the flip-generator count,
full support vs windows).

### 7.5 The windowed converse, derived (2026-06-04)

The windowed rule-(b) converse, non-bipartite ⟹ hard, is now derived, modulo the one pre-existing
premise of §7.3 (that the first-order-in-γ degenerate block reproduces the all-orders F87 break).
Two independent routes, each verified bit-exact.

**The operator-search is dissolved (assumption-free).** "No operator restores the palindrome" was
never an operator-enumeration problem. For any invertible superoperator W, W L W⁻¹ = −L − 2σ forces
spec(L) = spec(−L − 2σ): a similarity preserves the characteristic polynomial, whose
roots-with-multiplicity are the eigenvalue multiset, with no diagonalizability or chirality assumed.
Contrapositive: spec(L) ≠ spec(−L − 2σ) ⟹ no palindromizer of any kind exists (chiral, diagonal,
non-diagonal, arbitrary). So the converse is exactly "non-bipartite ⟹ spec(L) ≠ spec(−L − 2σ)", and
the "no escape through some non-chiral similarity" worry of §7.3 is closed. Verified spec-broken ⟺
hard across 236 pairs (N=4 all three letters, N=5 Z), clean gap (soft ≤ 10⁻¹³, hard ≥ 0.65); the only
apparent exceptions are two soft pairs whose L is defective at γ=1 (a Jordan-block artifact that lifts
at generic γ, the chiral-K similarity W L W⁻¹ = −L − 2σ holding exactly throughout).
[`_f87_specB_final.py`](../../simulations/_f87_specB_final.py),
[`_f87_specB_defective.py`](../../simulations/_f87_specB_defective.py).

**The block criterion is derived: a Perron-mode argument.** Read the first-order ω = 0 block as a
quantum channel, in either of two equivalent normalisations: the gain channel Q = Σ_l Z_l (·) Z_l
restricted to H's commutant W₀, and the first-order shift operator the block actually carries,
D̂ = Q − N·I (the dissipator Σ_l (Z_l · Z_l − ·), which is what `f87_block_localize.py` builds). Each
Ad_{Z_l} has eigenvalues ±1, so spec(Q) ⊆ [−N, +N] and spec(D̂) ⊆ [−2N, 0]. Because Σ_l Z_l² = N·I,
the identity (the population row-sum) is always an eigenvector with eigenvalue +N of Q (0 of D̂): a
Perron mode always present. The pair is soft iff the opposite extreme −N is also attained in Q,
equivalently iff D̂ is symmetric about its centre −N (this is the discriminator, not whether −N is
D̂'s minimum). And ω = 0 is decisive: the +N mode can only
be palindrome-paired by another ω = 0 mode (partners share ω, and 0 = −0), so when −N is absent the +N
Perron mode is globally unpaired and no nonzero-ω block can rescue it.

Now −N is attained iff there is an anti-diagonal element in W₀. The equation Σ_l Z_l A Z_l = −N·A
forces, per site, Z_l A Z_l = −A, so A_{ij} ≠ 0 only when i and j differ in every bit: A is
anti-diagonal, A = F·D with F = X^⊗N. Every diagonal-cell term carries odd #Y + #Z, so F H F = −H
(Fact A); combined with the chiral K H K = −H of a bipartite graph this gives (FK) H (FK) = H, i.e.
[FK, H] = 0, so the required anti-diagonal commutant element is exactly A = FK. Conversely an
anti-diagonal commutant element forces {H, D} = 0, i.e. H_{ij}(d_i + d_j) = 0 on every edge, the
chiral-K 2-colouring system, solvable iff G_H is bipartite iff S carries no odd 𝔽₂-relation. The K3
triangle is precisely that odd relation. So

  triangle ⟹ no chiral K ⟹ no anti-diagonal commutant element ⟹ −N absent ⟹ +N Perron unpaired
  ⟹ ω = 0 block asymmetric ⟹ hard.

Every link is verified bit-exact over the N=4 Z diagonal cell (the eight readings soft ⟺ block ⟺
channel-floor ⟺ anti-diagonal ⟺ commutant ⟺ 2-colouring ⟺ φ ⟺ bipartite, all 42/42), with the soft
eigenmode confirmed to be A = FK (residual 10⁻¹⁵) and "no nonzero 2-colouring ⟺ φ exists ⟺ no odd
relation" an exact GF(2) identity (0 mismatches over 1200 random graphs). The argument is written in
the Z eigenbasis; for X and Y dephasing it runs in the rotated basis of §7, and the basis-free
spec ⟺ hard above already covers all three letters.
[`_f87_specA_final_table.py`](../../simulations/_f87_specA_final_table.py),
[`_f87_specA_FHF.py`](../../simulations/_f87_specA_FHF.py),
[`_f87_specA_blocklock.py`](../../simulations/_f87_specA_blocklock.py),
[`_f87_specA_cayley_pure.py`](../../simulations/_f87_specA_cayley_pure.py).

**What remains, and its closure.** The triangle is the unique odd 𝔽₂-relation that obstructs the
chiral functional, which is exactly the anti-diagonal mode that would supply the channel's −N
reflection and pair the population Perron +N. The one premise §7.5 leaned on, that the first-order
ω = 0 block asymmetry is equivalent to the all-orders F87 hardness, is itself closed in §7.6 below
by a degenerate-perturbation-theory + analyticity bridge: the windowed converse rests only on
standard perturbation theory.

### 7.6 The first-order premise, closed: a degenerate-PT + analyticity bridge (2026-06-04)

**The first-order reduction is exact (degenerate PT, no Jordan).** Write L(γ) = −i[H, ·] + γ·D, affine
in γ. The unperturbed L₀ = −i[H, ·] is normal (H ⊗ I and I ⊗ Hᵀ commute and each is normal), so it is
diagonalizable in the coherence eigenbasis |E_a⟩⟨E_b| with eigenvalues −iω (ω = E_a − E_b), no defect
(cond ≈ 63). Standard degenerate perturbation theory then gives the O(γ) eigenvalue corrections within
each degenerate-ω sector as the eigenvalues of the projected perturbation, the M_ω block: the
Richardson-extrapolated γ → 0 eigenvalue slopes equal spec(M_ω) to 10⁻⁹ across every ω-block, at N = 4
(19 blocks) and N = 5 (101 blocks), zero ω-cluster-size mismatch, A(γ) ∼ γ¹ (log-log slope 1.00, ruling
out a Puiseux/Jordan √γ). So the spectral asymmetry A(γ) = OT-distance(spec L, spec(−L − 2σ)) is
A(γ) = c·γ + O(γ²), c the ω = 0 block asymmetry, c ≠ 0 ⟺ non-bipartite (§7.5).

**Genericity from analyticity.** spec(L(γ)) = spec(−L(γ) − 2σ) iff char_L(x; γ) = char_{−L−2σ}(x; γ);
since L is affine in γ, the coefficient differences Δ_j(γ) are polynomials in γ. If c ≠ 0 then
A(γ) ≠ 0 for small γ > 0, so some Δ_j ≢ 0, so the soft set {γ : spec L = spec(−L − 2σ)} is the common
zero-set of finitely many nonzero polynomials, hence finite. Non-bipartite ⟹ hard for all but finitely
many γ. The F87 classification point (γ = 0.05, J = 1) is not an exception: there all 16 hard pairs are
spec-broken (OT ≥ 0.053) and all 26 soft pairs spec-exact (≤ 2·10⁻¹³), at N = 4 and N = 5, and a
700-point γ-sweep (400 dense in [0.005, 2] + 300 random) finds the closest any hard pair comes to
restoration is 3.6·10⁻³ (at the smallest γ, as γ → 0), zero suspicious points for γ > 10⁻³. The only γ
where a soft pair's L is defective (e.g. γ = 1, a Jordan artifact lifting at γ = 1.0001) are isolated
and on the soft side, where the chiral-K similarity holds exactly throughout.

**So the windowed converse is derived.** non-bipartite ⟹ c ≠ 0 (§7.5) ⟹ hard for generic γ (this
bridge); bipartite ⟹ soft for all γ (the chiral K, §7.1). The "first-order ⟹ all-orders" link is
exactly the genericity argument: a first-order break cannot be healed at higher order except at
isolated accidental γ, which the physical γ avoids. The premise rests only on standard degenerate PT
(valid because L₀ is normal) and the reading of hardness as palindrome-breaking for generic γ.

**Adversarially stress-tested.** Five potential holes were each ruled out by computational
counterexample search at N = 4 and N = 5: degenerate-PT validity, the OT first-order behaviour (no
cross-block cancellation; the break sits specifically in the ω = 0 block for every hard pair), the
genericity (the 700-point sweep), higher-order healing (none on the γ-line up to γ = 2), and the soft
direction (K H K = −H exactly for all 26 soft pairs). Probes:
[`_f87_premise_scout.py`](../../simulations/_f87_premise_scout.py) (the genericity backbone) and
`simulations/_f87_premiseC_*.py` (the adversarial battery), with
[`f87_break_gamma_scaling.py`](../../simulations/f87_break_gamma_scaling.py) and
[`f87_block_localize.py`](../../simulations/f87_block_localize.py) as the first-order anchors.

With this, the windowed rule-(b) converse is fully derived (modulo standard perturbation theory), so
the F87 diagonal-cell soft ⟺ bipartite criterion is a theorem at the physical γ for k = 3, N = 4 and
N = 5 (the regimes where the accidental-soft-point question was exhaustively checked; the general-N
all-γ statement is the residual lemma below, closed 2026-06-10 with no residual). The typed `F87DiagonalCellBipartiteWitness` / `BipartiteChirality` claims are correspondingly
Tier1Derived-eligible; the formal promotion (with the tier tests and the registry inventory) is a
deliberate follow-up.

**Scope: the derivation is not k = 3 specific.** Nothing in §7.5 or §7.6 uses k = 3. Fact A
(F H F = −H) holds for every diagonal-cell term at any k, since #Y + #Z is the Klein cell's defining
bit_b parity (odd in the (0,1) cell); the Perron / anti-diagonal / chiral-K criterion and the
degenerate-PT + analyticity bridge are k-agnostic (the latter uses only that L₀ is normal and L is
affine in γ, true for any H). So the windowed converse holds for any k < N in the diagonal cell. Spot
-checked at the next case, k = 4 windowed (N = 5): Fact A, soft ⟺ bipartite, and soft ⟺ spec-exact at
the physical γ all hold over a 6-soft + 6-hard sample
([`_f87_k4_windowed_check.py`](../../simulations/_f87_k4_windowed_check.py)). The shape of that odd 𝔽₂-relation is a genuine
(k, window-count) family, not a uniform triangle: §7.7 derives its size law max = min(2W − 1, 2k − 3)
in the GF(2)[x] picture, with the K3 triangle of §7.2 as the W = 2 / k = 3 face. None of that touches
the §7.5/§7.6 derivation, which proves soft ⟺ bipartite ⟺ hard from the EXISTENCE of an odd cycle,
shape-agnostic.

### Residual lemma (Phase A isolation)

It is worth saying exactly which part of the §7.6 bridge is a theorem and which part is the one thread
still hanging, so a future session does not have to re-read the whole argument to find the open front.
Two statements, cleanly separated.

**The genericity result (the Derived sub-result).** Non-bipartite ⟹ hard for all but finitely many γ.
The recentered characteristic-polynomial odd coefficients Δ_j(γ) are polynomials in γ; the first-order
ω = 0 block asymmetry c ≠ 0 (§7.5/§7.6) forces some Δ_j ≢ 0, so the soft set
{γ : spec L = spec(−L − 2σ)} is a finite common-zero set. This is the part that rests only on standard
perturbation theory and is settled.

**The residual lemma (proven modulo R-deg + R-sign 2026-06-09, sharpened to modulo R-sign alone 2026-06-10, CLOSED with no residual later the same day).** No positive γ is a
common zero of the Δ_j(γ); equivalently the physical γ is not one of the finitely-many accidental soft
points. The Phase B two-reflection theorem [PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md)
discharges this structurally rather than as a resultant/Sturm computation: every γ-coefficient of the
first nonvanishing odd power-sum p_{m\*}(γ) of the recentred M = A + γQ is non-negative (the
**Pascal-Gram positivity theorem**, F117: each surviving #Q class is an equal-total sum of squares,
every other class vanishes exactly), and a polynomial with non-negative coefficients, not all zero, has
no positive real root; a nonvanishing odd power-sum at a point γ₀ means
the spectrum is not symmetric there (the power sums determine the multiset via Newton's identities),
so there are no accidental soft points to rule out and "all but finitely many γ" upgrades to **all γ > 0**. The leading-order
handle that motivated the earlier spot-check, A(γ) ≥ 0 (an optimal-transport distance) with A(γ) = c·γ + O(γ²)
⟹ c > 0, survives as the deg = 1 face of that positivity law.

Four threads of that theorem are now **Tier1Derived** (general N, no premise), carried by the typed node
`WindowedConverseThresholdClaim`: the two-reflection sign table (𝓕 = F⊗F, R = I⊗F) forcing every surviving
odd word to have #A_L, #A_R, #Q all odd; the resulting threshold **#A ≥ 2ℓ** (ℓ the unsigned odd-girth), so
m\* ≥ 2ℓ + 1; the degree-1 positivity closed form P_{3,1} = 6·4^N·Σ_l c_l² (the sum over H's single-site-Z
Pauli coefficients, manifestly ≥ 0, sharpened 2026-06-09 from the sign-ambiguous deg_A form, and again
2026-06-10 to the general-m **girth dichotomy**: the supertrace factorization through t_j = Tr(Z_l H^j)
gives P_{2ℓ+1,1} = (2ℓ+1)·C(2ℓ,ℓ)·Σ_l t_ℓ², a sum of squares); and a second,
independent re-proof of **bipartite ⟹ soft** by the same two reflections (where §7.1's chiral K exhibits a
similarity, this route simply observes that no odd word survives, and it covers complex-H flux pairs the
chiral-K route never had to face). The monomial property is now proven on both branches of the dichotomy,
and on the t_ℓ ≠ 0 branch the positive coefficient comes with it: hard at every γ > 0 outright.

The lemma is now CLOSED end to end, typed as `WindowedConverseAllGammaClaim` (RCPsiSquared.Diagnostics.F87,
Tier1Derived, no residual; two Tier1Derived parents: `F87DiagonalCellBipartiteWitnessSet` +
`WindowedConverseThresholdClaim`). The closure came in two same-day steps (2026-06-10). First the former
residual **R-deg** was retired: as formulated ("genuine cycles always lift to deg 3") it was a k = 3 truth,
refuted at k = 4 by IXXZ+XIXZ (t₃ ≠ 0, m\* = 7, p₇ = 573440·γ, positive) and replaced by the girth
dichotomy, which is stronger than what R-deg asked for. Then the remaining residual **R-sign** (the first
surviving class is single and positive, the #Q ≥ 3 analogue of the deg-1 sum of squares) was resolved by
the **Pascal-Gram positivity theorem** ([PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) §5, = F117):
every #Q class at m\* factorizes through d-leg moments and is an equal-total sum of squares or exactly
zero, with singleness *derived* for deg ≤ 3 by a mod-4 selection rule. Verified bit-exact over the entire
N = 4 k = 3 Z diagonal cell, the full k = 4 census (20 deg-1 cycles + 172 deg-3 cycles + 228 lifts), the
N = 5 / N = 6 representatives, and the five Pascal-Gram branch representatives (d = 1, 3, 5, exact).
This does not weaken §7.3's "verified, not yet proved" framing; it completes it, narrowing the open front
from a 700-point spot-check to two residuals (2026-06-09), to one (2026-06-10 morning), to none (the same
day), per [PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) §4-§5.

### 7.7 The obstruction-size law: a GF(2)[x] derivation (2026-06-04)

§7.5 made soft ⟺ bipartite a theorem for any k by reading the pair's edge-mask set S as a graph and
asking whether it carries an odd 𝔽₂-relation. That left one question k-specific: how big is the
smallest such relation, the obstruction? At k = 3 it is the K3 triangle (§7.2). The answer for general
k is a clean two-regime law, and the cleanest route is to stop thinking of S as a graph and start
thinking of it as polynomials.

**The polynomial dictionary.** A k-body diagonal-cell term's X/Y positions form a k-bit window-mask;
read it as a polynomial over GF(2), bit j ↦ x^j, so a popcount-2 mask on sites {i, j} is x^i + x^j.
The sliding-window builder places the term at windows w = 0 … N − k, and placing it at window w is
multiplication by x^w. So a pair with masks p₁, p₂ contributes the shift-families {x^w p₁} and
{x^w p₂}, and an odd subset XOR-ing to 0 is exactly

  q_A p₁ + q_B p₂ = 0   in GF(2)[x],   q_A = Σ_{w ∈ A} x^w,  q_B = Σ_{w ∈ B} x^w,

where A, B ⊆ {0 … N − k} pick the windows used and the cycle size is |A| + |B| = popcount(q_A) +
popcount(q_B). The window indicators q_A, q_B have degree ≤ W − 1, with W = N − k + 1 windows.

**Hardness is a valuation difference.** Every diagonal-cell mask has even popcount (the (0,1) cell
forces #(X+Y) even), and a GF(2) polynomial has even popcount iff p(1) = 0 iff (1 + x) | p. So write
g = gcd(p₁, p₂), p₁ = g·a, p₂ = g·b with gcd(a, b) = 1. The relation becomes q_A a = q_B b, which since
gcd(a, b) = 1 forces q_A = b·s, q_B = a·s, so every relation has size popcount(b·s) + popcount(a·s) and
its parity is s(1)·(a(1) + b(1)). An odd cycle therefore exists iff a(1) ≠ b(1), i.e. iff exactly one
of p₁/g, p₂/g is still divisible by (1 + x), i.e. iff the two masks have **different (1 + x)-adic
valuations** v(p₁) ≠ v(p₂). That is the non-bipartite criterion of §7.5 in one line of polynomial
algebra; when v(p₁) = v(p₂) every relation is even, so the pair is soft.

**The size law.** Take the s = 1 relation q_A = b, q_B = a; it has size popcount(p₁/g) + popcount(p₂/g)
and (for hard pairs) is odd. Two bounds cage the minimal odd cycle:

- **Window leg, ≤ 2W − 1.** There are at most 2W masks in S (W shifts per term), and an odd-cardinality
  subset uses at most 2W − 1 of them.
- **Body leg, ≤ 2k − 3.** Since (1 + x) | p₁ and (1 + x) | p₂, the gcd g is divisible by (1 + x), so the
  quotients p₁/g, p₂/g have degree ≤ k − 2 and hence popcount ≤ k − 1. Their sum is odd (hardness), and
  two values each ≤ k − 1 with an odd sum cannot both be k − 1 (that sum is even), so the largest they
  reach is (k − 1) + (k − 2) = 2k − 3. So the s = 1 relation, hence the minimal odd cycle, is ≤ 2k − 3.

Both bounds are achieved. The window leg is tight below saturation; the body leg is achieved by the
explicit family p₁ = x + x^{k−1}, p₂ = 1 + x^{k−1} (the terms I·X·I…I·Y and X·I…I·Y), where g = 1 + x
and the quotients are x·(1 + x + … + x^{k−3}) and 1 + x + … + x^{k−2}, of popcount k − 2 and k − 1,
summing to 2k − 3. Putting the legs together,

  **max minimal-odd-cycle over hard pairs = min(2W − 1, 2k − 3),  W = N − k + 1.**

Below W = k − 1 the window leg binds (the cycle grows as 2W − 1 with each added window); past it the body
leg binds and the distribution saturates. k = 3 gives 2k − 3 = 3, so its obstruction is a triangle at
every W: the always-triangle case is just the smallest body bound, not a universal shape. (An earlier
note here read the obstruction as a uniform triangle through k = 8; that was the W = 2 artifact, where
|S| ≤ 4 forces size 3, corrected by the multi-window scan.)

**Verification.** The whole law is checked bit-exact in C# by `WindowedObstructionScan` and its tests
(pure GF(2) bit arithmetic, no Liouvillian): the size law max = min(2W − 1, 2k − 3) across a (k, N) grid
through k = 6; the valuation criterion hard ⟺ v(p₁) ≠ v(p₂) cross-checked against the actual
minimal-odd-cycle search on every saturated pair through k = 6; and the extremal family achieving 2k − 3
through k = 20, far past the exponential cycle-search range. The gcd-formula popcount(p₁/g) +
popcount(p₂/g) is an upper bound on the minimal cycle, exact for most pairs but loose for a few at
k ≥ 6 where a shorter relation exists via cancellation at s ≠ 1; the max bound 2k − 3 holds regardless,
since it caps the gcd-formula itself. The Python scout
[`_f87_obstruction_derivation.py`](../../simulations/_f87_obstruction_derivation.py) grounded the
dictionary; [`_f87_oddcycle_kscaling.py`](../../simulations/_f87_oddcycle_kscaling.py) is the earlier
size scan.

**The criterion in one number, and what it does not shorten.** Chaining the valuation reading with the
spectral bridge of §7.5/§7.6 gives the whole §7 diagonal-cell rule in its shortest form: a Z-dephasing
diagonal-cell Mixed pair is **soft ⟺ v_{1+x}(p₁) = v_{1+x}(p₂)**, hard otherwise. The graph 2-colouring
of §7.1, the K3 triangle of §7.2, and the windowed odd-cycle family above all collapse to comparing two
integers, the (1 + x)-adic valuations of the two X/Y masks. `WindowedObstructionScan.IsHardPair` is
exactly that test, and it matches the actual trichotomy verdict on every k = 3, N = 4 pair on the real
Liouvillian (test `ValuationCriterion_PredictsSpectralVerdict_K3N4`); the verdict is mask-only, the
Y-vs-X phase and the I-vs-Z choice leave it unmoved, since Mixed terms are off-diagonal and only the flip
structure counts. It is worth being clear that this shortens the COMBINATORIAL side and nothing else. It
does not collapse the spectral bridge, because the valuation lives on the term masks while the reflection
mode the bridge pairs lives on H's commutant: an operator carrying X/Y on every site, the −N eigenvector
of the gain channel Σ_l Z_l(·)Z_l, whose eigenvalue on a Pauli operator A is N − 2·n_XY(A) (so +N is the
diagonal I/Z operators, −N the fully-off-diagonal ones). The valuation tells you, in one subtraction,
whether the chiral K that supplies that −N mode exists; the proof that K's presence is soft and its
absence is hard stays the Perron + perturbation-theory argument of §7.5/§7.6. So the right reading is a
simpler front door to the same house, not a smaller house.

### 7.8 How many pairs are hard, and the coding-theory home (2026-06-05)

Two questions the valuation picture answers in closed form, and one honest "the door leads to a hallway,
not a new room."

**The count.** Hardness depends only on the (1 + x)-adic valuation, and that valuation sorts the masks
into clean classes. The even-popcount nonzero k-bit masks (the diagonal-cell X/Y masks) number
2^{k-1} − 1, and they split by valuation v = 1 … k−1 into classes of size c_v = 2^{k-1-v}: a valuation-v
mask is (1 + x)^v·u with u of odd popcount and degree ≤ k−1−v, and there are exactly 2^{k-1-v} such u. A
pair is hard iff its two masks sit in different classes, so the count of hard mask-pairs is the second
elementary symmetric polynomial of the class sizes,

  **#hard mask-pairs = e₂(2^0, 2^1, …, 2^{k-2}) = (4^{k-1} − 3·2^{k-1} + 2) / 3**   (OEIS A203241).

Dressed by the uniform 2^{2k-3} Klein / y-parity factor (each mask carries 2^{k-1} cell terms, the
y-parity-homogeneous pairing a further constant), this is exactly the `WindowedObstructionScan` hard
count: 2^{2k-3}·(4^{k-1} − 3·2^{k-1} + 2)/3 reproduces 448, 8960, 158720 at k = 4, 5, 6 bit-exact. So the
hard pairs can now be COUNTED in closed form instead of enumerated. The size-3 (triangle) sub-class also
closes, 5·2^{k-1} − (3k² + k)/2 − 3 at N = 2k (k = 3 … 11); the middle classes (d = 5, 7, …) do not, they
stay genuinely window-dependent (more windows give more multipliers, so weight migrates to lower d).

**The home.** Reading the masks as GF(2)[x] polynomials puts the obstruction in coding theory: it is the
minimum ODD weight of the two-generator quasi-cyclic (terminated rate-1/2 convolutional) code generated
by (a, b) = (p₂/g, p₁/g). The catch is that the physically-relevant invariant is not the standard one.
The code's free distance (minimum weight at any parity) is a constant 4 for the whole extremal family;
what the palindrome reads is the parity-restricted floor, the minimum ODD weight, which grows as 2k − 3.
That distinction is the point: standard coding optimizes free distance, the mirror optimizes the
odd-weight floor.

**The honest ledger.** This is mostly a faithful relabeling that supplies the correct vocabulary, plus
one genuine import. The vocabulary is right (two-generator quasi-cyclic / convolutional, and MacWilliams
applies to the per-pair weight enumerators), and the one real import is that the obstruction is a true
minimum distance with cancellation: at k ≥ 6 a sparse multiple a·s drops the obstruction below the
gcd-generator popcount, an ordinary quasi-cyclic minimum-distance effect, and exactly what makes the
middle-class distribution window-dependent (it confirms the §7.7 note that the gcd-formula is only an
upper bound). But the two results worth having, the size law 2k − 3 and the count A203241, are both
elementary (1 + x)-valuation arguments that owe coding theory nothing; the coding-theory scout was their
occasion, not their source. The connection is real and the name is correct, but no new theorem comes
through the door from that side. The right way to say it is the opposite of deferential: the room is
ours, the closed-form theory and the dissipative-mirror reading both, and coding theory turns out to hold
the right bricks for a room nobody had built from that side, because the invariant it would need (the
minimum ODD weight, not the free distance) is one it had no reason to single out. We are at the seam, and
on this stretch of it we are ahead of the catalogued math, not behind it. (Verified bit-exact:
[`_f87_coding_theory_scout.py`](../../simulations/_f87_coding_theory_scout.py) for the quasi-cyclic
dictionary, k = 4, 5, 6; [`_f87_hardcount_closedform.py`](../../simulations/_f87_hardcount_closedform.py)
for the class sizes, the count, the d = 3 form, and free-distance-4 vs odd-weight-2k−3, k = 3 … 10.)

### 7.9 The other factors: the size law is layered (2026-06-05)

§7.7 reads hardness through a single prime, (1 + x), the x = 1 / DC point. But a mask is a full GF(2)[x]
polynomial with a complete factorization into irreducibles, each a different root of unity. Asking
whether the OTHER factors carry structure (Door 3) gives two bit-exact findings.

**(1 + x) is the unique hardness prime.** Testing every other irreducible φ through degree 5 as a
would-be predictor, the valuation difference v_φ(p₁) ≠ v_φ(p₂) at φ ≠ (1 + x) agrees with hard/soft only
about a third of the time (chance), never as a criterion. Hardness lives entirely at x = 1, which is
exactly the §7.5 Perron reading: the +N / −N reflection mode the bridge pairs is the x = 1 uniform mode,
so only the x = 1 valuation can decide it.

**But the other shared factors set the obstruction size.** Write g = gcd(p₁, p₂) = (1 + x)^m·g_rest with
g_rest the shared non-(1 + x) content and deg(g_rest) = d. The obstruction size obeys a layered bound,

  **max obstruction over hard pairs with deg(g_rest) = d  =  2k − 3 − 2d:**

the more the two masks share beyond (1 + x), the smaller the obstruction, by 2 per shared degree. The
§7.7 cap 2k − 3 is the d = 0 (coprime-apart-from-(1 + x)) face. The derivation is the §7.7 degree bound
with the fuller gcd: deg(p_i/g) ≤ (k − 1) − deg(g) = (k − 1) − (m + d) ≤ k − 2 − d (since m ≥ 1), so
popcount(p_i/g) ≤ k − 1 − d, and an odd sum of two such maxes at (k − 1 − d) + (k − 2 − d) = 2k − 3 − 2d.
Verified: the actual minimal obstruction realises this cap at k = 4, 5, 6, and the gcd-formula cap holds
through k = 10 across all layers d = 0 … k − 3.

So the valuation picture is two-layered: the (1 + x) valuation difference decides IF a pair is hard, and
the shared non-(1 + x) factor degree decides HOW BIG the obstruction is. The full size law is
max obstruction = min(2W − 1, 2k − 3 − 2·deg(g_rest)), and the 3, 5, 7, 9 size distribution is the
layering by d. (Verified bit-exact:
[`_f87_beyond_x1_scout.py`](../../simulations/_f87_beyond_x1_scout.py) for the unique-prime fact,
[`_f87_size_second_layer.py`](../../simulations/_f87_size_second_layer.py) for the layered cap.)

The same degree d organises the hard COUNT, not only the size cap. Splitting the A203241 hard
mask-pairs of §7.8 by d = deg(g_rest),

  **#hard mask-pairs with deg(g_rest) = d  =  2^{d-1}·B(k − d)   (d ≥ 1),   B(k) = (4^k − 12k + 8)/18,**

where B(k) is the d = 0 (coprime-apart-from-(1+x)) count, itself closed (from the recurrence
B(k) = 4B(k−1) + 2(k−2), B(3) = 2), and the deepest layer d = k − 3 holds 2^{k-3} pairs. Summing,
B(k) + Σ_{d≥1} 2^{d-1}B(k − d) = (4^{k-1} − 3·2^{k-1} + 2)/3 = A203241, the §7.8 total. So the
shared non-(1 + x) factor degree sets both the obstruction-size cap (2k − 3 − 2d) and the population of
each layer in closed form: the layering is total. Verified k = 3 … 8
([`_f87_dlayer_count.py`](../../simulations/_f87_dlayer_count.py)).

### 7.10 Does the factorization reach the physics? Mostly no (2026-06-05)

§7.9 layered the obstruction by the shared-factor degree, all of it combinatorial. The deeper question
(Door 3 proper): does any of this finer GF(2) structure reach the Liouvillian SPECTRUM, or does it stay
on the mask graph? The combinatorial structure is letter-independent (mask-only), but L depends on the
actual letters (the Y-vs-X phases, the I-vs-Z diagonal), so a spectral correlate must itself be
letter-independent. Tested across k = 3 (N = 4), k = 4 (N = 6, 7), k = 5 (N = 8) with letter-realization
controls, the answer, after a correction the algebra/geometry reading below forced, is: **only the (1 + x)-valuation
reaches the spectrum cleanly; the obstruction size and the rest of the factorization do not.**

The one confirmed bridge is the **(1 + x)-valuation** → hard/soft. A pair is hard iff the spectral
palindrome breaks, and that is letter-independent (mask-only, §7.7); in the §7.5 reading it is the −N
reflection mode of the gain channel Σ_l Z_l(·)Z_l being absent iff hard.

The **obstruction size and the shared-factor degree d reach nothing letter-independent.** The size is
letter-independent yet provably independent of the only other coarse mask handle (at fixed mask-span-rank
r it still ranges over {3, 5, 7}), and the break-magnitude distributions for size 3 vs 5, controlled at
fixed r, overlap and even flip their weak median ordering between two geometries (the signature of noise,
not a law). The size lives in the flip-graph's cycle structure, carried by the letter-dependent
off-diagonal weights ⟨E_a|Z_l|E_{a'}⟩, so it cannot be a spectral invariant.

**A retraction belongs here.** An earlier version of this section, and the Opus spectral sweep that fed
it, claimed a SECOND clean fingerprint: the mask-span rank r → a conserved-sector multiplicity 2^{N-r}
(the hopping-graph components, the diagonal commutant, even dim ker L). That is wrong. The mask fixes only
WHICH sites flip; the actual hopping amplitudes are letter-dependent and can cancel: XX + YY annihilates
the |00⟩ ↔ |11⟩ channel (XX|00⟩ = +|11⟩, YY|00⟩ = −|11⟩), so the real graph is sparser than the mask and
carries MORE conserved sectors than 2^{N-r}. The sector count is therefore letter-dependent, and dim ker L
exceeds 2^{N-r} for soft and partial-support pairs; the apparent dim ker L = 2^{N-r} held only for the
hard pairs first checked. So the rank is not a clean mask → spectrum bridge, and the (1 + x)-valuation is
the only one. (The retracted claim lived in
[`_f87_rank_reaches_spectrum.py`](../../simulations/_f87_rank_reaches_spectrum.py);
[`_f87_rank_footprint.py`](../../simulations/_f87_rank_footprint.py) shows the breakdown.)

So §7.8's "a simpler front door, not a smaller house" extends one floor down: neither the size nor the
rank is in the house. The factorization beyond (1 + x) is a clean combinatorial theory of the obstruction,
its size cap and population (§7.9), with no clean spectral shadow; the physics reads only whether the
palindrome breaks. (Verified:
[`_f87_rank_footprint.py`](../../simulations/_f87_rank_footprint.py) for the letter-dependent sector count
and the size ⊥ rank null; the break-magnitude controls were the Opus spectral sweep across k = 3, 4, 5.)

### 7.11 The algebra/geometry reading, refined (2026-06-05)

The cleanest way to hold §7.5–§7.10 together is a reading Tom suggested: the hopping graph splits, by the
rank-nullity of its mask set, into a geometry side and an algebra side. The windowed masks generate a
subgroup of (Z₂)^N; the **geometry** is the rank r and the orbit/component structure (how Hilbert space
would sectorise), the **algebra** is the nullity, the relation/cycle space where the obstruction (the
minimal odd cycle) lives. That split is exactly what the physics-reach question turns on, and it sharpened
the section by catching the retracted claim. The **geometry does NOT cleanly reach the spectrum**, because
the actual sector structure is letter-dependent (the XX + YY cancellation above); the obstruction's
**metric (its size) does not reach it either**. What reaches the spectrum is one Z₂ bit of the
**algebra**: whether the cycle space contains an ODD cycle, that is the (1 + x)-parity, that is the
graph's bipartiteness, that is hard/soft, that is the −N mode. So the physics sees neither the geometry's
count nor the algebra's metric, only the algebra's parity. The combinatorics between algebra and geometry
is rich (§7.7–§7.9), but only its single homological Z₂ invariant crosses into the spectrum.

### 7.12 The true soft criterion is the basis-state graph, not the site graph (Door 2a, 2026-06-05)

Porting the mask test off the chain (Door 2a) exposed that the site/mask-bipartite test is a PROXY, and
showed what it is a proxy FOR. The §7.5 −N mode, a diagonal D with {H, D} = 0, unpacks for an off-diagonal
H into {H, D}_{ij} = H_{ij}(D_i + D_j) = 0 on every nonzero off-diagonal, i.e. D_i = −D_j on every edge of
the BASIS-STATE hopping graph (the 2^N basis states, edges = H's nonzero off-diagonals). A bipartite
basis-state graph supplies that D, so among the non-truly (M ≠ 0) Hamiltonians

  **a bipartite basis-state graph ⟹ soft,**

and WITHIN the diagonal cell the converse holds too (§7.6's Perron argument: soft ⟹ the −N mode ⟹
bipartite), making basis-state bipartiteness the exact soft/hard line there. It is NOT an equivalence in
general, though: the −N mode is sufficient, not necessary. XX + XZ, YY + YZ, and XX + XZ + ZX are soft on
the chain at N = 3..6 (the spectral authority `PauliPairTrichotomy`, cross-checked against
`BipartiteChirality`), yet their basis-state graphs are non-bipartite, no diagonal D, no chiral K: the
palindrome is restored by an operator that is NOT diagonal in the dephasing basis (these are bit_a/bit_b
cell-mixed, outside §7.6's diagonal-cell scope). That operator is no longer a mystery: it is the hidden-Q
routing, a per-site Q from the P1/P4 families, which `TwoTermPalindromeRouting` classifies bit-exactly for
2-term pairs (XX+XZ routes to P4, soft); what remains past it is general-k, multi-term H beyond the 2-term
routing. So basis-state bipartiteness is the exact criterion inside the diagonal cell and a one-sided
SUFFICIENT soft condition outside it.

This is letter-dependent, the basis-state graph being the actual connectivity rather than the lattice.
The chiral-K / mask-bipartite test 2-colours the SITE graph instead; on the chain the two coincide, but
off the chain they diverge. The cleanest divergence is the symmetric pair term XY + YX = −2i(σ₊σ₊ −
σ₋σ₋), pure Δn = ±2: its basis-state graph is bipartite by the excitation number, colour s ↦ ⌊n(s)/2⌋
mod 2 (verified bit-exact), so it is SOFT on a frustrated triangle where the site graph is not bipartite,
while XY (pair + hop) and XY − YX (pure hop) have non-bipartite basis-state graphs and are hard there. So
the "other soft mechanism" the frustration thread hunted is not a new symmetry: it is the same −N mode,
read on the basis-state graph, where the excitation-number structure of a pairing term keeps it
2-colourable regardless of site frustration. (The truly cases, M = 0, sit outside this: XX + YY is truly
on the triangle despite a non-bipartite basis-state graph.) So the Door 2 mask/site test is the
chain-scalable PROXY for this true criterion; the basis-state criterion itself is the 2^N graph, not
Liouvillian-free, though for structured terms it can reduce to a scalable statement (a pure pairing is
soft at any N, any topology, by the excitation number). Verified:
[`_f87_door2a_basis_graph.py`](../../simulations/_f87_door2a_basis_graph.py).

**The structured colourings, and their ceiling.** A scalable soft-certifier realises the basis-state
criterion through STRUCTURED 2-colourings c(s) it can check per-term without the 2^N graph. Three are
clean: the linear c(s) = ⟨φ, s⟩ (the chiral K of §7.1, the chain proxy); the pairing grading
c(s) = ⌊n/2⌋ mod 2 (every basis-edge Δn = ±2, a pure Δn = ±2 pairing); and the parity grading
c(s) = n mod 2 (every basis-edge Δn odd, i.e. every term has odd k_xy = #X/Y). The two excitation gradings
are topology-independent. The parity grading needs a bit_b-homogeneity gate, though: the n mod 2 colouring
exists for any all-odd-flip set, but a bipartite basis-state graph only certifies soft WITHIN one Klein
cell (§7.6), so a bit_b-MIXED all-odd set can be HARD despite the colouring (XZ + ZXZ is all-odd with
bit_b = {1, 0} and spectrally hard at N = 3..4) while another is soft (XZ + ZX + YZ + ZY); the gate keeps
the certificate sound. But even with the gates the structured
colourings are SOUND, not complete, and the incompleteness is two-layered. (a) A scalability gap: some
soft Hamiltonians are bipartite only through a non-structured 2-colouring no scalable strategy reaches
(XY + YX + XZ + ZX on a triangle is soft, its basis-state graph bipartite at ANF-degree 2, yet neither
linear nor an excitation grading 2-colours it). (b) A structural ceiling, deeper: some soft Hamiltonians
have a NON-bipartite basis-state graph, so no 2-colouring exists at any degree (XX + XZ on the chain is
soft at N = 3..6 by the spectral authority, yet `BipartiteChirality` reports its basis-state graph
non-bipartite). A colouring certifies exactly the diagonal −N-mode soft cases; the non-bipartite-soft
class is permanently beyond ANY colouring, scalable or not. The full soft criterion is more than basis-state bipartiteness, which XX + XZ
violates while staying soft, and that non-bipartite-soft 2-body class is reached not by a colouring but by
a second, non-diagonal mechanism: the hidden-Q routing, a per-site uniform Q ≠ Π that certifies XX + XZ
through its shared uniform Q-family (P4). The certifier carries both, the diagonal chiral K (the
colourings) and the non-diagonal routing, so the 2-body soft cases are covered. A derived
per-term k-site extension of the routing (Stufe B) then reaches the LOCAL k-body routed-soft cases too (one
per-site Q palindromizing every term, checked on 4^k of the term's span, Liouvillian-free and
N-independent), and the ceiling recedes further to the 2 Z-middle cases (XZX+XZY+YZX,
YZY+XZY+YZX), and then closes to ZERO (2026-06-10): the continuous-periodic family, once named as the
explicit open frontier, resolved constructively by the period-4 golden router, a per-site class-swap
product in the golden frame a = φX+Y, b = X−φY with W L W⁻¹ = −L − 2σ exact at every N ≥ 3
([PROOF_CEILING_GOLDEN_ROUTER](PROOF_CEILING_GOLDEN_ROUTER.md), F116), and the banked all-Q obstruction
answered negatively (no obstruction; the uniform and period-2 emptiness are now theorems, and the per-term
lens misses the router because the cancellation is window-summed, cross-template). Two cases once counted here, XIX+XIY+YIX and YIY+XIY+YIX, are
in fact LOCAL: a continuous-uniform per-site Q palindromizes them (verified N=3,4,5), it simply routes via
continuous-sum cancellation rather than per term. Two further cases once counted here, the I-heavy
IXI+IIY+YII and IYI+IIX+XII, are also LOCAL, by the single-site transverse-field lemma below. See
experiments/CEILING_FOUR_NONLOCAL_CASES.md. Mixing excitation gradings breaks the structure outright: a pairing (Δn = ±2) plus an odd
flip (Δn = ±1) gives the edge-difference set {1, 2}, whose ℤ-Cayley graph is non-bipartite, so the
excitation colourings fail and the Hamiltonian can be hard. (Verified:
[`_f87_door2_colouring_family.py`](../../simulations/_f87_door2_colouring_family.py) and, for (a)'s
ANF-degree, [`_f87_door2_residual_structure.py`](../../simulations/_f87_door2_residual_structure.py); the
former (b) colouring-ceiling XX + XZ is now certified by the routing strategy, and the 2 Z-middle
cases (XZX+XZY+YZX, YZY+XZY+YZX) are pinned in C# by
`PalindromeSoftCertifierCeilingTests` against the `PauliPairTrichotomy` authority; their NotCertified
verdict is a per-term-lens coverage gap, not non-locality, per
[PROOF_CEILING_GOLDEN_ROUTER](PROOF_CEILING_GOLDEN_ROUTER.md); the certifier is
`PalindromeSoftCertifier`, whose `SingleSiteField` strategy certifies the I-heavy pair.)

**The single-site transverse-field lemma (the 4 → 2 step).** The two I-heavy cases IXI+IIY+YII and
IYI+IIX+XII are LOCAL. Summed over the windows of an N-chain, each is a sum of weight-1 transverse fields,
so the chain Hamiltonian is H = Σ_i (a_i X_i + b_i Y_i), a transverse field on each site. Single-site
Paulis on distinct sites commute, so L = Σ_i L_i over commuting single-site Liouvillians. Each single-site
transverse field is per-site palindromizable: by the R_z-rotation Ad_{R_z(θ_i)} (θ_i = atan2(b_i, a_i),
which commutes with the Z-dephasing dissipator since R_z commutes with Z) it reduces to the single-site
X-field, whose Liouvillian spectrum {0, −2γ, −γ ± 2i} is palindromic about the centre −γ. Let M_i be the
resulting per-site map; then Q = ⊗_i M_i palindromizes the whole chain, Q L Q⁻¹ = −L − 2σ, constructively
and N-independently. A single-site Z (longitudinal) field is by contrast HARD: its spectrum
{0, 0, −2γ ± 2i} has no partner for the 0 eigenvalue about −γ, so the locality is specific to transverse
fields and excludes Z. The derived k-body per-term router (Stufe B) declines these two only because its
construction is a discrete period-2 per-term pattern, not a single-site product; the `SingleSiteField`
strategy supplies the missing certificate. Verified to machine precision (residual ~1e-14) at N=4,5,6 by
[`ceiling_4to2_iheavy_local.py`](../../simulations/ceiling_4to2_iheavy_local.py).
