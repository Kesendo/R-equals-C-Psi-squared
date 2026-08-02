# Schur round (round 4): the read-off corrected, the law placed as an affine sp(12) reflection

Scratch, local-only (`_rb_schur{1..5}.py`). Run from repo root, `PYTHONIOENCODING=utf-8`. Exact over ℤ.
This round CORRECTS round 3's Lever-A framing (which was broken at the read-off) and places the
two-row law as a verified affine-Weyl statement.

## The correction (load-bearing): round 3's e_n read-off was never valid
- `content5(e_n)` (round 3, `_rb_trans2.py` `__main__`) prints **0 for all n** — because `e_n(w)` is
  fully S₅-symmetric (gate A1), any W(C₅) antisymmetrization of it vanishes. The README claim
  "L(e_n) = c_k(n)" was asserted but its own code shows 0.
- A type-A (S₅-alternant) read-off of `e_n` at any tested (δ, M²-shift) also gives 0 or no constant
  ratio (`_rb_schur1.py`, `_rb_schur2.py`). **`e_n(w)` alone does not carry `c_k(n)` by any tested
  read-off.** The Lever-A "Schur palindrome over the 36-weight alphabet" was thus not on solid ground.

## What IS true and gated this round
1. **The true m-slice of X is S₅-antisymmetric** (`_rb_schur3.py`, 0/12000 across m = 24, 20, 16):
   `X_m = A_δ · S_m` (5-var alternant × symmetric cofactor) is genuine — the alternant part of round
   3's factorization is real. Slice built by exact convolution `X = P1 ⊛ P2` restricted to
   coord-1 = m (packing: coord c in bits [16c, 16c+16), bias 32768).
2. **But `F_k(m)` is NOT the 5-var dominant coefficient** `[X]_{(m,tail_k)}` (`_rb_schur4.py`: ratio
   non-constant, supports disjoint). The read-off flips the **6th** sign (the m/coord-1 sign) too —
   it is the full W(C₆) antisymmetrization, not a frozen-coord-1 5-var content. This is why the
   5-variable `e_n` reduction cannot capture it.

## The correct framing (verified, `_rb_schur5.py`)
`n_λ` is the **symplectic sp(12) = C₆ character coefficient** of the virtual character
`S := X / A_{ρ,C6}`, ρ = (6,5,4,3,2,1), read off by the committed
`n_raw(λ) = Σ_{ε∈{±1}⁶} sgn(ε) [X]_{ε∘2(λ+ρ)}` (= 2 n_λ, f133 gate G5).

- **The two-row law is exactly the affine reflection `μ₁ ↦ 22 − μ₁`** in μ = λ+ρ (μ₁ = j+6, the
  partner 16−j = 22−μ₁), holding μ₂..μ₆ fixed. Gate: **0/36** mismatches. In raw t₁-exponent units
  (m = 2μ₁) this is the anti-period m ↦ 44−m of round 2b (44 = 2·22).
- **This reflection is AFFINE, not finite:** the reflection center μ₁ = 11 (m = 22) is NOT the center
  of X's coord-1 support (which is 0). The finite Weyl reflection about 0 is the PROVEN total-negation
  evenness (`Xc(−ν) = Xc(ν)`, round 2b); the open part is the off-center affine node at 11.
- **Exact domain, reconfirmed on the correct read-off** (`_rb_schur5.py` part B): the reflection
  `μ₁ ↦ 22−μ₁` on three-row (j,k,l) holding (k,l) HOLDS at **l = 0 (0/36), l = 1 (0/25), l = 3
  (0/9)** and FAILS at **l = 2 (8/16)**. This is the ρ₄-tail resonance boundary the earlier rounds
  saw — now verified via the symplectic read-off, not the 6-var slice.

## Residual, precisely named (the honest hand-off)
Prove that the virtual sp(12) character `S = X / A_ρ` is invariant, on two-row weights, under the
**affine Weyl reflection s₀ of C₆ at level 21** (μ₁ ↦ 22−μ₁). Equivalent forms already banked:
anti-period `F_k(m+44) = −F_k(m)` in-wall (round 2b); palindrome `c_k(n) = c_k(14−n)` (round 3).
The finite half (reflection about 0 = total negation) is PROVED; the affine half is the open core.

**Why this is the right target, and the literature/rep-theory pointer:** X = Δ̂ · ∏₃₁ (deleted-sheet
spinor product; the deleted sheet 1⁶ is the TOP spin weight of the B₆ spinor system). Affine-Weyl
quasi-periodicity with a fixed off-center node and a sign is the signature of a **theta-function /
affine character** factor. The precise statement to prove or cite: removing the top spin weight from
the B₆ spinor sine-product `SP` induces, on the sp(12)-dual two-row coefficients, the level-21 affine
C₆ reflection s₀. This is the "affine C₆/B₆ wall reflection" round 2 conjectured, now pinned to
s₀ at level 21 with its exact resonance domain (l ∈ {0,1,3}). No published identity was found (round
1's King/Rains–Warnaar gaps); the object is the deleted-top-weight spinor discriminant, treated by
neither.

## Net
Round 3's e_n/type-A Schur reformulation is **withdrawn** (read-off invalid). The law is instead the
affine sp(12) Weyl reflection s₀ (μ₁ ↦ 22−μ₁, level 21) of X/A_ρ on two-row weights — verified 0/36,
domain l ∈ {0,1,3} pinned, finite half proved, affine half open with a precise theta/affine-character
target. Scratch: `_rb_schur1.py` (W(C₅) read-off = 0), `_rb_schur2.py` (type-A search, no hit),
`_rb_schur3.py` (true slice S₅-antisym 0/12000), `_rb_schur4.py` (F_k ≠ dominant coeff),
`_rb_schur5.py` (sp(12) affine reflection 0/36 + domain l∈{0,1,3}).
