# PROOF F103: F87 Trichotomy Z₂³ Refinement at k=3 (N=4 Empirical Anchor)

**Status:** Tier 1 derived. The 42:8 closed-form rule was found 2026-05-29 (diagonal-cell hardness rule, §6), and its two atomic sub-rules were then unified into a single criterion (§7): a diagonal-cell pair is soft iff H's hopping graph is bipartite in the dephasing letter's eigenbasis. The direction bipartite ⟹ soft is derived (Π followed by a chiral sublattice K). The converse non-bipartite ⟹ hard **closes at full support** (§7.4, 2026-05-30: at k=N a Mixed+Mixed pair has only two flip generators, which always admit the chiral K, settling F111's blocked "Mixed+Mixed = soft" modulo M) and **stays verified-not-derived in the windowed regime** (k<N, all three letters, N=4 and N=5).
**Date:** 2026-05-24
**Anchor:** N=4, k_body=3, 294 Z₂³-homogeneous + Y-par-homogeneous Pauli pairs (pair count is N-independent at fixed k; the empirical anchor is N=4)
**Regenerate:** `simulations/f87_z2cubed_split_n4_k3.py` (~60s)

## Abstract

F87 classifies Pauli pairs into three buckets (truly, soft, hard) based on how the pair's M-residual closes under a chosen dephase letter. F102 just showed that y-parity is a real third axis on top of the Klein (bit_a, bit_b) signature once Pauli strings reach body count 3 or higher. The natural follow-up question is whether the F87 trichotomy actually splits along the new y-parity axis at k=3, or whether the trichotomy is y-parity-blind. This proof gives the empirical answer at the N=4, k=3 anchor.

The answer is yes, the trichotomy splits cleanly, with five structural patterns that survive across all three dephase letters. The 294 Pauli pairs (this count depends only on the k=3 letter alphabet, not on N) divide into three trichotomy classes, and inside each class, the y-parity axis carves the cells into sub-cells with recognizable shape: truly always lands in y_par = 0, mother-soft always in y_par = 1, the diagonal hard cells split 42:8 with a Y-inversion on the Y-dephasing diagonal, the diagonal soft cells split symmetrically 13:13, and the off-diagonal soft cells split into two named sub-patterns B and C.

The proof begins empirical at the (N=4, k=3) anchor (an exhaustive enumeration verified bit-exactly) and is then closed: §6 gives a closed-form counting rule for the 42:8 split, and §7 supplies the bipartite-chirality mechanism (the soft direction, bipartite ⟹ soft, is derived; the windowed k<N hard-direction converse is the one remaining edge, §7.3). The question whether the pattern is k-stable and N-stable is answered by F105 (N=5, k=3, N-stable) and F106 (N=4, k=4, sharpening to 228:0) as sibling anchors. Together the three anchors map out which parts of the Z₂³ structure are N-stable, which are k-stable, and which depend on the specific (N, k) regime.

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
(non-bipartite ⟹ hard) is the one open edge; see §7.3.

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
tempting moment proof is ruled out: D̂ is diagonal in the computational basis
(D̂|i⟩⟨j| = −2·Hamming(i⊕j)·|i⟩⟨j|, the Absorption Theorem), but the odd spectral moments
Tr((L+σ)^{2k+1}) vanish for soft and hard alike, so the break is moment-invisible , a set-pairing
asymmetry, not a moment identity. The remaining gap is therefore set-level: that the degenerate
first-order block spectrum fails to pair iff an odd cycle is present. See
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
