# The Windowed F87 Converse, Per Block: Where the Odd Cycle Breaks the Pairing

**Status:** A sharpening of the then-open windowed converse, not a proof. The per-block decomposition
is analytic (the transpose relation M(−ω) = M(ω)^T) and verified bit-exact at N=4; the localization
(which block breaks) is read off the same N=4 anchor. At the time of writing the windowed converse
non-bipartite ⟹ hard stayed open exactly as [PROOF_F103 §7.3](../docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md)
left it. *Superseded 2026-06-09:* the converse is now proven modulo R-deg + R-sign by the
[two-reflection monomial theorem](../docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md); this
document remains the per-block localization record.
**Date:** 2026-05-30
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Script:** [`simulations/f87_block_localize.py`](../simulations/f87_block_localize.py)
**Builds on:** [PROOF_F103 §7.3](../docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) (the
first-order reduction, [`f87_first_order_degenerate.py`](../simulations/f87_first_order_degenerate.py)),
[BIPARTITE_CHIRALITY_DIAGONAL_CELL](BIPARTITE_CHIRALITY_DIAGONAL_CELL.md) (the criterion + the
clock's tick), the chiral K ([PROOF_K_PARTNERSHIP](../docs/proofs/PROOF_K_PARTNERSHIP.md)), and
the F80 one-sidedness M = −2i·(H⊗I) ([PROOF_F80](../docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md)).

---

## The question

The F87 trichotomy sorts a Pauli pair, under single-letter dephasing, into truly / soft / hard by
whether the Liouvillian spectrum pairs each λ with −λ − 2σ. §7.3 reduced the open windowed
converse (non-bipartite ⟹ hard, body count k < N) to a first-order-in-γ statement: at γ = 0 the
Liouvillian L = −i[H, ·] is purely imaginary and trivially symmetric about −σ; the break grows
first order in the dephasing tick, and the **degenerate first-order dephasing block** , the
dissipator D̂ diagonalized within each degenerate-frequency subspace of the H-eigenbasis coherences
|E_a⟩⟨E_b| , reproduces the actual break bit-exact, with c = 0 ⟺ bipartite. §7.3's remaining gap is
set-level: *"that the degenerate first-order block spectrum fails to pair iff an odd cycle is
present."* This note asks a smaller question: not yet *why* it fails, but **where**.

## The decomposition: the pairing is per block

The first-order modes are μ = −iω + γ·s, where ω = E_a − E_b is the frequency and s is an
eigenvalue of the per-frequency dephasing block

```
M_ω[(a',b'),(a,b)] = Σ_l ⟨E_a'|Z_l|E_a⟩ ⟨E_b|Z_l|E_b'⟩  −  N·δ.
```

The palindrome sends μ to −μ − 2σ = −i(−ω) + γ·(−s − 2N) (σ = Nγ): a mode in the **−ω block** with
shift −s − 2N. Now the backbone. Because H is real in the dephasing basis and each Z_l is
symmetric, the minus-frequency block is the transpose of the plus-frequency one,

```
M(−ω) = M(ω)^T,
```

so the two share a spectrum (verified to machine precision, max ±ω mismatch ≈ 2·10⁻¹⁴). The
partner shifts −s − 2N therefore have to land back in the *same* block's spectrum, and the
palindrome collapses to a condition **on each block separately**:

> **soft ⟺ for every frequency ω, the block spectrum {s ∈ eig M_ω} is symmetric about −σ**
> (equivalently, the pure Z-overlap part Z_ω := M_ω + N·I is symmetric about 0).

A hard pair cannot pair its full spectrum because at least one of its frequency blocks is itself
asymmetric about −σ. This is finer than "the spectrum fails to pair": it says the failure is
block-local.

## What the blocks show (N = 4)

The size-weighted sum of the per-block asymmetries about −σ reproduces the measured break c
exactly ([`f87_block_localize.py`](../simulations/f87_block_localize.py)):

```
SOFT  XXZ+ZXX  (bipartite)      c = 0.0000   no block breaks
FLUX  IXY+XIY  (non-bipartite)  c = 0.1191   ω=0 (size 40)  +  (±2,     each 0.369)
REAL  XXZ+XZX  (non-bipartite)  c = 0.2559   ω=0 (size 40)  +  (±1.236, each 0.667)  +  (±3.236, each 0.298)
```

Three facts fall out:

1. **The ω = 0 block always breaks.** This is the diagonal-frequency sector: the populations
   |E_a⟩⟨E_a| together with any same-energy coherences. In both hard cases it is the robust
   breaker (FLUX 0.32, REAL 0.48), present whenever the pair is hard.
2. **The ω ≠ 0 blocks break in equal ± pairs.** This is exactly the transpose relation
   M(−ω) = M(ω)^T: the +ω and −ω blocks share a spectrum, so they share an asymmetry.
3. **A bipartite pair breaks nothing**, block by block, consistent with c = 0.

## Why ω = 0 is the cleanest place to attack the converse

The chiral K = diag((−1)^ℓ) inverts the single-particle spectrum, E ↦ −E, so it maps the ω block
to the **−ω** block: for ω ≠ 0 it relates two *different* blocks. The exception is ω = 0, where K
maps the block to itself ((a,a) ↦ (N+1−a, N+1−a), and same-energy coherences accordingly). So on
the ω = 0 block alone, the chiral K , when it exists , acts *internally*, and the natural reading
is that it chiralizes Z₀ (forces the symmetry about 0); when no K exists (an odd cycle), ω = 0 is
precisely the block that the table shows breaking. The set-level question , Z_ω symmetric about 0
⟺ bipartite , is thus isolated, for ω = 0, to a single block on which the would-be pairing operator
acts within rather than across. That is the cleanest entry for a proof, and it is exactly where the
empirics put the robust break.

## Honest status

- **Sharpening, not proof.** At the time of writing the windowed converse non-bipartite ⟹ hard
  remained open, as in §7.3 (since 2026-06-09: proven modulo R-deg + R-sign, see the Status note).
  What is added: the pairing is *per block* (analytic, via the transpose relation), and the break
  *localizes* (ω = 0 always; ω ≠ 0 in equal ± pairs). The open set-level claim becomes per-block:
  Z_ω symmetric about 0 ⟺ bipartite, with ω = 0 the natural first target.
- **The transpose relation is N-general; the localization table is N = 4** (the §7.3 /
  BIPARTITE_CHIRALITY anchor). The per-block reduction does not depend on N.
- **First-order (small-γ) statement.** It reads the converse off the degenerate first-order
  blocks; §7.3's result (first-order reproduces the break bit-exact) is what licenses that.
- **Windowed only.** The full-support converse (k = N) is separately closed in §7.4 via the
  flip-generator |S| ≤ 2 argument; this note is only about the windowed regime k < N, where an
  odd cycle can appear.

## Anchor and open work

- Script: [`f87_block_localize.py`](../simulations/f87_block_localize.py) (per-block asymmetry,
  the size-weighted sum = c, and the transpose check M(−ω) = M(ω)^T).
- The open per-block crux: Z_ω symmetric about 0 ⟺ bipartite, for each ω.
- The ω = 0 sub-case (K acts within the block) as the proof entry; then the ω ≠ 0 blocks, where K
  relates a block to its transpose partner, are the remaining set-level step.
