# The seed exists because N is odd: a nullity reduction of the census input

*Existence side of the codim-1 containment corollary. Complements the exclusion shell census of
[F89_MULTI_SECTOR_MONODROMY.md](F89_MULTI_SECTOR_MONODROMY.md) and the count-scan census table of
[F89_PATH_K_DIABOLIC.md](F89_PATH_K_DIABOLIC.md). 2026-07-04.*

## The one input the corollary cannot derive

The containment corollary of [PROOF_CODIM1_BY_ADDITIVITY](../docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md)
transports and folds a seed across the whole diamond, but it takes the seed's *existence* as given:
**a real defective exceptional point must exist on the (1,2) block at each odd N.** Until now that
input was supplied by census, the PT-break count scan, checked to N = 11 (the table in
[F89_PATH_K_DIABOLIC.md](F89_PATH_K_DIABOLIC.md); the counts 4/6/7/9 at N = 5/7/9/11). A census is a
lower bound over a window, not a law, and it grows more expensive every N (the N = 11 run took 2 h 31 m).
This note replaces the census question with an exact identity, and reduces the remaining work to two
named lemmas. It does *not* close the existence theorem; it makes precise what is left, so no future
session re-runs the scan to "confirm".

## The pencil, and two spectral endpoints

The (1,2) block is the affine pencil L(q) = A + q·C on coherences |a⟩⟨b| (a a two-excitation ket, b a
one-excitation bra) of the XY chain under uniform Z-dephasing:

- **A** = −2·n_diff, the dephasing diagonal, with only two values: rung −2 (n_diff = 1, the bra site
  sits on a ket site) and rung −6 (n_diff = 3, the three sites all distinct).
- **C** = −i(H₂ ⊗ I − I ⊗ H₁), the coherent XY hop; C is anti-Hermitian, C = i·K with K self-adjoint in
  the orbit metric, so K has real spectrum.
- At odd N the spectrum of L(q) is self-conjugate for every real q (the unmirrorable central site;
  `DiabolicReflectionParityWitness`), so the characteristic polynomial is real in Λ and the eigenvalues
  are real or in complex-conjugate pairs.

Count the real eigenvalues at the two ends of the q-axis (exact SVD nullities, not eigenvalue thresholds,
which mis-read the tail at large q):

- **r(∞) = nullity(C)**, the asymptotically-real modes: at large q, λ ∼ q·(i·κ), κ ∈ spec K, so a mode
  stays real iff κ = 0, iff it lies in ker C.
- **r(0⁺) = nullity(P₋₂ C P₋₂) + nullity(P₋₆ C P₋₆)**, the modes with a real (that is, zero) first-order
  shift on each degenerate dephasing level. This is measured just above q = 0, past the immediate
  imaginary lift-off of the levels whose first-order shift is nonzero.

## The identity, and why it forces the seed

The verifier `simulations/seed_existence_nullity_check.py` computes both endpoints exactly, N = 3..13:

| N | dim | n₂ = nullity(P₋₂CP₋₂) | n₆ = nullity(P₋₆CP₋₆) | r(∞) = nullity(C) | r(0⁺) − r(∞) |
|---|-----|-----|-----|-----|-----|
| 3 | 9 | 2 | 3 | 3 | **2** |
| 5 | 50 | 4 | 6 | 6 | **4** |
| 7 | 147 | 6 | 9 | 9 | **6** |
| 9 | 324 | 8 | 12 | 12 | **8** |
| 11 | 605 | 10 | 21 | 21 | **10** |
| 13 | 1014 | 12 | 18 | 18 | **12** |

**r(0⁺) − r(∞) = N − 1, exactly, for every odd N.** Since N − 1 > 0, the real-eigenvalue count strictly
drops between small and large q, and the drop happens entirely at finite q > 0 (r(0⁺) is read past the
q = 0 lift-off). A drop in the real count is a real-to-complex transition; the following lemma turns it
into a seed.

**Lemma (defective at a simple zero; Tier 1 derived).** At a simple real zero q\* of the discriminant
disc_Λ χ, exactly two eigenvalue branches meet in a square-root branch point, so L(q\*) carries a 2×2
Jordan block (geometric 1 < algebraic 2): a defective, order-2 exceptional point, at a real double
eigenvalue. This is the classical Kato fact (a single order-k defective Puiseux cycle makes the
discriminant vanish to order k − 1; a semisimple or real-real crossing gives an *even* order, hence no
sign change); it makes "a real-to-complex transition is defective" a theorem at simple zeros, and pins
the "defective, not semisimple" character the census could only read numerically per locus.

So the existence question reduces to the surplus N − 1, and the surplus splits into three pieces, two of
them now proved.

## Piece 1: n₂ = N − 1 is a path count (Tier 1 derived)

Parametrize the −2 states by ordered pairs (p, q), p ≠ q: |p,q⟫ = |{p,q}⟩⟨{p}|, with p the shared
ket-and-bra site and q the ket-only site. The only XY hops that stay on the −2 rung give
P₋₂ C P₋₂ = i·K_red with two edge types: a **q-hop** (the ket-only site moves, |p,q⟫ ↔ |p,q±1⟫, weight
−1) that cannot cross the shared site p, and a **swap** (|p,q⟫ ↔ |q,p⟫, weight +1) only when |p − q| = 1.
Because the q-hop cannot cross p, each fixed-p fiber splits into two segments, and the N − 1 swaps glue
them into

  **Path_k = R_k ∪ L_{k+1},   k = 0, …, N − 2,   a simple path on (N − 1 − k) + (k + 1) = N vertices.**

The graph of K_red is therefore a disjoint union of **N − 1 simple paths, each on N vertices**. Every
path is a tree, so its signed adjacency gauges to the unsigned path P_N, whose nullity is 1 iff N is odd.
Hence nullity(P₋₂ C P₋₂) = (N − 1)·[N odd] = N − 1 for odd N (and 0 for even N). Each fiber-path has
exactly N vertices, and a path carries a zero mode iff its length is odd: **the seed-forcing kernel is
nonempty precisely because N is odd**, the combinatorial face of the same unmirrorable central site that
gives the odd-N self-conjugacy. (For even N the kernel is empty, consistent with the even-N absence of
real diabolics.)

## Piece 2: r(∞) is the free-fermion fusion-resonance count (Tier 1 derived)

ker C is the space of intertwiners {ρ : H₂ρ = ρH₁}. Equivalently, under Jordan-Wigner the (1,2)
coherence |{j,k}⟩⟨{i}| is the cubic operator c†_j c†_k c_i, and the coherent part C = −i[H, ·] is
diagonal in the single-particle eigenmode basis: [H, c†_a c†_b c_c] = (λ_a + λ_b − λ_c)·c†_a c†_b c_c,
with λ_k = 2 cos(kπ/(N + 1)) the one-magnon energies. The one-magnon spectrum is simple, and the
two-magnon energies are the Slater sums λ_a + λ_b (the same free-fermion additivity E_DE = ε_j + ε_k
read for the frequencies in [F89_PATH_K_DIABOLIC.md](F89_PATH_K_DIABOLIC.md)). Therefore

  **r(∞) = nullity(C) = #{ (c, {a,b}) : a < b, λ_a + λ_b = λ_c }, the fusion-resonance count.**

The asymptotically-real Liouvillian modes are exactly the free-fermion fusion resonances, where a
two-magnon energy equals a one-magnon energy. The count is number-theoretic in the cosines (extra
resonances at special angles make it jump: 3, 6, 9, 12, 21, 18 at N = 3..13, not a clean 3(N−1)/2).

## The narrative: the wild resonances cancel, the odd-N kernel survives

The fusion-resonance count Σ appears in *both* endpoints: r(∞) = Σ and r(0⁺) = n₂ + n₆ with n₆ = Σ
(Piece 3, below). It cancels, and the surplus is carried entirely by the odd-N path kernel:

  **r(0⁺) − r(∞) = n₂ = N − 1.**

The number-theoretically irregular resonances contribute equally to the small-q and large-q real counts
and force nothing; the seed is forced by the odd-N structure of the −2 rung alone.

## Piece 3: the remaining lemma (open), and the spectral-inheritance lead

The cancellation needs **(N1′): n₆ = nullity(P₋₆ C P₋₆) equals the same fusion-resonance count Σ.**
Verified N = 3..15. It is not a projected intertwiner (the −2/−6 split lives in the dephasing/position
basis, conjugate to the eigenmode basis where C is diagonal, and the projection of a resonant eigenmode
is not a −6 zero mode). The lead is a strong, non-generic structural fact (verified N = 3..7):

  **spec(P₋₂ C P₋₂) ⊆ spec(C)  and  spec(P₋₆ C P₋₆) ⊆ spec(C), every eigenvalue exactly.**

A compression whose spectrum is an exact sub-multiset of the original is a supersymmetry-like spectral
inheritance, not mere Cauchy interlacing. Given it, and with n₂ = N − 1 known and the conservation
Σ_μ [mult_μ(C₂) + mult_μ(C₆) − mult_μ(C)] = 0, (N1′) is equivalent to the multiplicity statement
**mult₀(C₆) = mult₀(C)**: the coupling adds the surplus N − 1 to the kernel entirely through the −2
channel. The open task is to prove the spectral inheritance (an intertwining relation among C₂, C₆ and
the coupling) and read off the zero multiplicity.

## Reproduce

```bash
python simulations/seed_existence_nullity_check.py
```

Self-validating: it asserts (F1) the surplus N − 1 exactly (N = 3..13), (N2) the (N − 1)-paths-of-N
decomposition, (FF) nullity(C) = the fusion-resonance count = n₆, and (SI) the spectral inheritance.

## Status

The existence input is **reduced, not closed**. Proven: the defective-at-a-simple-zero lemma, n₂ = N − 1
(the path count), and r(∞) = the fusion-resonance count (free fermions). Open, and now precisely named:

1. **(N1′)** n₆ = the fusion-resonance count, equivalently mult₀(C₆) = mult₀(C) under the verified
   spectral inheritance; the free-fermion / third-quantization structure of the whole quadratic
   Liouvillian is the natural tool.
2. **The codim-2 β-exotic:** a count-dropping transition is defective unless it is the non-generic
   order-3 point where the local 2×2 has a nilpotent linear term; ruling this out for all odd N is a
   codimension-2 genericity statement (the fixed gap 4 and the structurally nonzero inter-rung hop argue
   against it, but it is unproven).

When (N1′) and the β-exotic close, the census input becomes a law for all odd N, and the containment
diamond membership follows at every N with no further scan. The two proved lemmas are natural candidates
to promote to typed Claims with live witnesses (the `DiabolicReflectionParityWitness` template) once the
theorem is complete.
