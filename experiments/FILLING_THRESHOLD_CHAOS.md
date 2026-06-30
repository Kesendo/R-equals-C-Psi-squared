# Dissipative Quantum Chaos is a Filling Threshold, not an Integrability One

**Date:** June 30, 2026
**Tier:** 1 (computational, live-witnessed)
**Arc:** `f89_galois_open_doors` door C, decisive follow-up (the `diabolic_over_higher_n` neighbour)
**Witness:** `inspect --root fillcsr` (`FillingThresholdWitness`)
**Engine:** `FillingThresholdCsr` (Diagnostics) on `WeightCoherenceBlock.Build(n, wKet, wBra, q, Δ, field)` (Core)

## What this is about

Door C of the F89 Galois arc asked a sharp question and got a clean **null**: the part of the
(SE,DE) relaxation spectrum that carries the non-solvable Galois group S_d (the H_B-mixed half,
S_8/18/32/53 for path 3..6) does **not** read as dissipative quantum chaos. Its complex spacing
ratio (CSR, Sá-Ribeiro-Prosen) is Poisson-like, not GinUE: algebraic chaos over the coupling q and
spectral chaos at fixed q are different things here (`inspect --root galoischaos`,
[RANDOM_MATRIX_THEORY.md](RANDOM_MATRIX_THEORY.md) Result 3).

The two follow-up stages then showed the null is **robust**: breaking the Liouvillian's free-fermion
additivity (XXZ anisotropy Δ) or the Hamiltonian's integrability (a random Z-field) — with or without
interactions — never drove the (SE,DE) CSR to GinUE (the `IntegrabilityBreakingCsr` Door-C sweep). The
reading that survived: the null is **structural / kinematic**. The (SE,DE) block is a *dilute*
2-excitation sector (one excitation in the ket, two in the bra); a dilute sector cannot thermalize,
so it stays non-chaotic however hard you break integrability.

This experiment is the decisive test of that reading. If the dilute block is non-chaotic *because*
it is dilute, then a **dense** coherence sector of the **same** Liouvillian — one with extensive
excitation content — should reach GinUE under the same knobs. That is exactly what we find. Chaos
here is a **filling threshold**, not an integrability one.

## Terms used here

- **Coherence block (wKet, wBra)** — the Liouville-space sector spanned by |a⟩⟨b| with popcount(a)=wKet,
  popcount(b)=wBra. The Z-dephasing XXZ Liouvillian L = −i[H,·] + D is closed on it (the XX+YY hopping
  conserves each leg's weight; the Δ·ZZ, the dephasing, and a Z-field are diagonal). (SE,DE) = (1,2) is
  the Door-C block; **dilute** = small total weight, **dense** = wKet, wBra near N/2 (extensive filling).
- **CSR** — for each eigenvalue λ, z = (NN − λ)/(NNN − λ) with NN/NNN its nearest / next-nearest
  neighbour. ⟨|z|⟩ is radial rigidity, ⟨cos θ⟩ angular repulsion. 2D-Poisson: ⟨|z|⟩≈0.66, ⟨cos θ⟩≈0.
  GinUE (dissipative quantum chaos, class A): ⟨|z|⟩≈0.74, ⟨cos θ⟩≈−0.24 (finite-size references are
  computed live, never hardcoded).
- **Class A** — the symmetry class whose reference is GinUE. We use **unequal** weight (p, p+1): the F1
  palindrome Π maps the (p,p+1) block to the *conjugate* (p+1,p) block, not to itself, so no residual
  antiunitary survives — the GinUE target is the right one (not AI⁺/AII⁺). Confirmed live: under a
  random field the block spectrum's conjugation-match fraction is ≈ 0.

## The setup

The general block builder `WeightCoherenceBlock.Build(n, wKet, wBra, q, Δ, field)` (Core) is the same
verbatim physics as the Door-C `XxzCoherenceBlock` extended to any (wKet, wBra) and carrying a per-site
longitudinal field Σ_k w_k Z_k (the diagonal frequency −i·q·(fe(ket) − fe(bra)), fe(c)=Σ_k w_k·z_k).
`FillingThresholdCsr.DisorderSweep` draws w_k ~ U[−W, W] per realization, pools the per-spectrum z's of
the **off-real** bulk (the valid CSR domain once conjugation symmetry is broken), and bootstraps a 95%
CI. References are finite-size-matched to the measured per-spectrum z-count. Canonical operating point:
γ = 1, q = J/γ = 1, interacting Δ = 1, ergodic disorder window W = 0.75.

Methodology is inherited verbatim from the Door-C harness `IntegrabilityBreakingCsr` (it reuses that
class's `Reduce` and finite-size references): pool per-spectrum z's, never raw eigenvalues across spectra
(that superimposes independent point processes and fakes Poisson); bootstrap the CI; compare against
finite-size-matched Poisson/GinUE references (not the asymptotic values, which carry the wrong edge bias).

## The result: the dilute block stays Poisson, the dense block reaches toward GinUE

Walking the filling ladder at fixed N (the unequal blocks (1,2) → (3,4)/(4,5), references size-matched
to the per-spectrum z-count):

| N | block | dim | ⟨\|z\|⟩ | ⟨cos θ⟩ | GinUE ⟨cos θ⟩ ref | % of GinUE angle |
|---|-------|-----|--------|---------|-------------------|------------------|
| 6 | dilute (1,2) | 90   | 0.670 | −0.040 | −0.171 | ~23% |
| 6 | **dense (3,4)** | 300  | **0.719** | **−0.089** | −0.206 | **43%** |
| 7 | dilute (1,2) | 147  | 0.674 | −0.042 | −0.185 | ~23% |
| 7 | **dense (3,4)** | 1225 | **0.712** | **−0.129** | −0.229 | **56%** |
| 8 | dilute (1,2) | 224  | 0.680 | −0.035 | −0.243 | ~14% |
| 8 | **dense (4,5)** | 3920 | **0.718** | **−0.162** | −0.243 | **67%** |

(N=8 with a 95% bootstrap CI: dilute ⟨|z|⟩=0.680 [0.677,0.684], dense ⟨|z|⟩=0.718 [0.715,0.721], over ~14k/16k pooled z's; references finite-size-matched at the 3920-point per-spectrum size.)

Read across the two filling regimes:

- **The dilute (1,2)=(SE,DE) block stays Poisson.** ⟨cos θ⟩ ≈ −0.04 at every N, ⟨|z|⟩ ≈ 0.67 — no
  angular repulsion, no N-trend. The Door-C null, reproduced through the general builder. A 2-excitation
  sector does not thermalize, however hard the disorder + interactions break integrability.
- **The dense (p,p+1) block near half-filling is chaotic.** Its radial statistic ⟨|z|⟩ ≈ 0.71–0.72
  sits **at the GinUE value already**, and its angular repulsion ⟨cos θ⟩ is **negative and climbs toward
  GinUE with the block size**: −0.089 → −0.129 → −0.162 at N = 6/7/8 (≈ 43% → 56% → 67% of the
  size-matched GinUE angle). The dilute block stays flat at ≈ 0 (~14–23%) across the same N.

The radial statistic reaching GinUE while the angular one converges more slowly from above is the
expected finite-size signature of genuine class-A (GinUE) statistics: ⟨|z|⟩ saturates first, ⟨cos θ⟩
(a finer correlation) trails and catches up with size. The dense block is on the GinUE side and walking
in; the dilute block is parked on the Poisson side.

(The isospectral pairs (2,3) ≅ (3,4) at N=6 read identically, as they must — particle-hole / the global
spin-flip QP relate them; see [F89d cross-fold](F89_PATH_K_DIABOLIC.md).)

## Why it lands: extensive filling is the missing ingredient

Door-C established that breaking the Galois (over q) or Hamiltonian (Bethe) integrability does not move
the dilute CSR. This experiment isolates what does: **filling**. The same Liouvillian, the same disorder,
the same interactions — only the excitation content changes — and the CSR crosses from Poisson to GinUE.
So:

- Galois chaos (the non-solvable S_d over the coupling field) and spectral chaos (GinUE at fixed q) are
  **distinct objects** and merge **only at extensive filling**. The dilute (SE,DE) block where the S_d
  lives is the wrong place to look for GinUE — not because the algebra is too simple, but because the
  sector is too dilute to thermalize.
- The MBL caveat holds: the strong-disorder corner (W = 2) relaxes the dense block back toward Poisson
  (localization), so the GinUE window is the intermediate-W ergodic band — exactly where many-body
  dissipative chaos is expected. Interactions (Δ = 1) deepen the repulsion over the free-fermion case
  (Δ = 0), as a non-integrable many-body sector should.

This closes Door C's last open edge with a structural verdict, not the originally-anticipated
integrability gradient: the fixed-q chaos line is drawn by **filling**, and the dilute block's persistent
Poisson statistics are the kinematic shadow of a sector that holds too few excitations to scramble.

## Reproduce

```bash
dotnet run --project compute/RCPsiSquared.Cli -c Release -- inspect --root fillcsr
```

recomputes the dilute-vs-dense contrast live (N=6 and N=7) with finite-size-matched references and the
class-A guard, and prints the CONFIRMED verdict. The full filling ladder and the N=6→7→8 size scaling are
the reconnaissance tests in
`compute/RCPsiSquared.Diagnostics.Tests/Foundation/FillingThresholdCsrTests.cs`
(`Reconnaissance_FillingLadder_N6`, `Reconnaissance_DenseSizeScaling`, `Reconnaissance_DenseN8_Full`),
tagged `[Trait("Category", "SLOW_FILLCSR")]` (the N=8 EVDs run minutes); run them with
`--filter "Category=SLOW_FILLCSR"`.

## See also

- [RANDOM_MATRIX_THEORY.md](RANDOM_MATRIX_THEORY.md) — Result 3 (the galoischaos null, the dilute control)
  and Result 5 (this finding).
- [F89_TOPOLOGY_ORBIT_CLOSURE.md](F89_TOPOLOGY_ORBIT_CLOSURE.md) — the live-Galois section (S_d over q).
- `inspect --root galoischaos` (the Δ=0 / dilute control) and `inspect --root fillcsr` (this result);
  the Door-C staged sweep + methodology live in the `IntegrabilityBreakingCsr` harness.
