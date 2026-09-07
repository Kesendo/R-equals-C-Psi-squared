# Finite-Size Filling-Associated Crossover in Dissipative Spectral Statistics

**Date:** June 30, 2026
**Tier:** 1 (computational, live-witnessed)
**Arc:** `f89_galois_open_doors` door C, decisive follow-up (the `diabolic_over_higher_n` neighbour)
**Witness:** `inspect --root fillcsr` (`FillingThresholdWitness`)
**Engine:** `FillingThresholdCsr` (Diagnostics) on `WeightCoherenceBlock.Build(n, wKet, wBra, q, Δ, field)` (Core)

This is filling-associated crossover evidence at N=6..8: movement toward GinUE,
not a causal or thermodynamic threshold theorem. The file/API names retain the research label.

## In plain words (the ground truth first)

The physical object is a **chain of N spins** (qubits) with nearest-neighbour XXZ interactions (the standard
anisotropic spin-exchange model), each spin also continuously **watched by its environment** (local
Z-dephasing, the decoherence a real device feels).
An open quantum system like this is generated not by a Hamiltonian alone but by a **Liouvillian** L, the
superoperator that plays the Hamiltonian's role once dissipation is present: L = −i[H,·] + D, where H is
the spin Hamiltonian and D is the dephasing. Its eigenvalues live in the complex plane (real part = a
relaxation rate, imaginary part = an oscillation frequency), and that cloud of eigenvalues is the
**spectrum** we study.

"**Dissipative quantum chaos**" is a statement about the *local spacing statistics* of that complex cloud.
Two reference fingerprints bracket it: a structureless **2D-Poisson** cloud (the signature of an integrable
/ non-chaotic system, eigenvalues don't "feel" each other) and **GinUE** (the Ginibre Unitary Ensemble, the
random-matrix reference for a fully chaotic dissipative system, whose eigenvalues repel). The diagnostic
that tells them apart is the complex spacing ratio (**CSR**), summarised below. The question of this
experiment is which fingerprint a given piece of the spectrum wears, and what controls it.

Hamiltonian integrability and Liouvillian spacing statistics are different questions.
Poisson-like spacings are compatible with integrable structure, but do not prove it.
Here the controlled comparison changes the coherence sector at fixed interacting disorder.

"**Filling**" is how many excitations (flipped spins) the piece of the spectrum under study carries,
relative to N: a **dilute** sector holds a handful of excitations, a **dense** sector holds ~N/2
(extensive, i.e. growing with the system size).
The result, in one line: the sampled dense sectors show stronger GinUE-like repulsion than the dilute sectors.

The rest of this document is the precise version. (Acronyms expanded on first use; EVD = eigenvalue
decomposition, the numerical step that produces the spectrum.)

## What this is about

Door C of the F89 Galois arc asked a sharp question and got a clean **null**: the part of the
(SE,DE) relaxation spectrum that carries the non-solvable Galois group S_d (the H_B-mixed half,
S_8/18/32/53 for path 3..6) does **not** read as dissipative quantum chaos. Its complex spacing
ratio (CSR, Sá-Ribeiro-Prosen) is Poisson-like, not GinUE: algebraic chaos over the coupling q and
spectral chaos at fixed q are different things here (`inspect --root galoischaos`,
[RANDOM_MATRIX_THEORY.md](RANDOM_MATRIX_THEORY.md) Result 3).

The Door-C sweeps distinguish three Hamiltonian cases:

- Nonzero Delta breaks free-fermion additivity, but uniform XXZ remains Bethe-integrable.
- A random longitudinal Z disorder at Delta=0 remains quadratic (Anderson/free fermions).
- A generic random field plus Delta!=0 is the interacting disordered nonintegrable test.

The open-chain Bethe solution and the disordered interaction/Anderson boundary are independent
Hamiltonian inputs ([open XXZ](https://arxiv.org/abs/0707.1995),
[random-field XXZ](https://arxiv.org/abs/2403.09608)). Reflection, conjugation or cross-fold symmetry can
break without each knob breaking Hamiltonian integrability. The finite executed CSR evidence is a
dilute non-GinUE reading over the sampled knobs, followed by the **same** Liouvillian's dilute-vs-dense
comparison at canonical Delta=1 plus disorder. It supports a filling dependence in that regime,
not a universal inability of dilute sectors to thermalize.

## Terms used here

- **Coherence block (wKet, wBra)**: the Liouville-space sector spanned by |a⟩⟨b| with popcount(a)=wKet,
  popcount(b)=wBra. The Z-dephasing XXZ Liouvillian L = −i[H,·] + D is closed on it (the XX+YY hopping
  conserves each leg's weight; the Δ·ZZ, the dephasing, and a Z-field are diagonal). (SE,DE) = (1,2) is
  the Door-C block; **dilute** = small total weight, **dense** = wKet, wBra near N/2 (extensive filling).
- **CSR**: for each eigenvalue λ, z = (NN − λ)/(NNN − λ) with NN/NNN its nearest / next-nearest
  neighbour. ⟨|z|⟩ is radial rigidity, ⟨cos θ⟩ angular repulsion. 2D-Poisson: ⟨|z|⟩≈0.66, ⟨cos θ⟩≈0.
  GinUE (dissipative quantum chaos, class A): ⟨|z|⟩≈0.74, ⟨cos θ⟩≈−0.24 (finite-size references are
  computed live, never hardcoded).
- **Class A**: the symmetry class whose reference is GinUE. We use **unequal** weight (p, p+1): the F1
  palindrome Π maps the (p,p+1) block to the *conjugate* (p+1,p) block, not to itself, so no residual
  antiunitary survives; the GinUE target is the right one (not AI⁺/AII⁺). Confirmed live: under a
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

- **The dilute (1,2)=(SE,DE) block stays Poisson.** ⟨cos θ⟩ ≈ −0.04 at every N, ⟨|z|⟩ ≈ 0.67: no
  strong angular repulsion or clear N-trend in this table. The Door-C null is reproduced through the
  general builder at the executed interacting-disorder operating point; it is not a thermalization proof.
- **The sampled dense (p,p+1) blocks move toward GinUE.** Their radial statistic ⟨|z|⟩ ≈ 0.71–0.72
  is **near the GinUE reference**, and their angular repulsion ⟨cos θ⟩ is **negative and climbs toward
  GinUE with the block size**: −0.089 → −0.129 → −0.162 at N = 6/7/8 (≈ 43% → 56% → 67% of the
  size-matched GinUE angle). The dilute block stays flat at ≈ 0 (~14–23%) across the same N.

The radial statistic is closer to the GinUE reference than the angular statistic.
Their differing finite-size trends support a crossover reading; they do not establish
asymptotic convergence or a thermodynamic phase boundary.

(The isospectral pairs (2,3) ≅ (3,4) at N=6 read identically, as they must: particle-hole / the global
spin-flip QP relate them; see [F89d cross-fold](F89_PATH_K_DIABOLIC.md).)

## What the comparison establishes

At the same Liouvillian, disorder and interactions, changing excitation content changes the observed
CSR: the sampled dense sectors approach the GinUE reference more closely than the dilute ones.
Galois structure over the coupling field and spectral statistics at fixed coupling remain distinct objects.
These finite sizes do not establish a universal filling threshold or a cause of thermalization.

The strong-disorder corner (W = 2) moves the dense block back toward Poisson-like statistics; this is
compatible with localization but is not an MBL proof. In the sampled disorder window, interactions
(Δ = 1) deepen repulsion relative to the quadratic free-fermion case (Δ = 0). The distinction is
between those executed operating points, not between every disordered and every clean Hamiltonian.

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

- [RANDOM_MATRIX_THEORY.md](RANDOM_MATRIX_THEORY.md): Result 3 (the galoischaos null, the dilute control)
  and Result 5 (this finding).
- [F89_PATH_K_GALOIS.md](F89_PATH_K_GALOIS.md): the live-Galois section (S_d over q).
- `inspect --root galoischaos` (the Δ=0 / dilute control) and `inspect --root fillcsr` (this result);
  the Door-C staged sweep + methodology live in the `IntegrabilityBreakingCsr` harness.
