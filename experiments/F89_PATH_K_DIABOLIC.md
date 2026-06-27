# The Path-k Diabolic at the N=4→N=5 Transition: integrability keeps the crossings, the self-fold placed them on the real axis

**Status:** Tier 2 (computed, gate-validated at path-3; one honest open item below). The instrument reproduces the known path-3 diabolic structure exactly; the path-4 result (diabolics exist but sit at complex q, none at physical real q) is decisive and confirmed by a fine real-axis scan plus a tight zoom. All 15 in-region coalescences classify consistently across three signals (monodromy loop, gap-scaling exponent, EpCharacter): 11 diabolics, 4 defective EPs. The earlier character-vs-loop disagreement on two near-axis pairs was resolved (loop contamination by a neighbour EP at the dense path-4 spectrum, fixed by reading the small intrinsic loop radius). The complex-q set is not claimed complete. The interpretive synthesis (integrability = existence, self-fold = real-axis placement) is a grounded reading, not yet a theorem.

**Date:** 2026-06-27
**Authors:** Thomas Wicht, Claude (Opus 4.8)

## What this is about

At N=4 (the path-3 chain) the watched (SE,DE) coherence block has a single silent **diabolic** point at a physical coupling q = J/γ = q_EP ≈ 0.659: two relaxation modes cross with their eigenvectors staying independent (semisimple), at λ_EP = −4 + 2iJ. The [zeros_connecting_structure](../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs) arc traced it from below ([`gmscan --trace`](../compute/RCPsiSquared.Diagnostics/Foundation/GaloisMonodromyWitness.cs)): the +2.349 σ_T twin pair slides onto the fold and merges there. That picture is N=4-only, because the within-block self-fold (the half-filling self-complement DE = bar(DE), the rung balance 2 = N−2) exists only at N=4. The forward question: **does an analogue survive to N=5 (path-4), or does removing the self-fold kill the diabolic?**

The honest prior (the reviewed investigation plan's R-1/R-2): a diabolic is a *node* of the spectral curve, codim-3-complex = 6 real conditions, overdetermined by 4 in a 2-real-DOF complex-q scan, so **generically absent unless a symmetry auto-satisfies it**. At N=4 the self-fold (plus the reality of q_EP) did exactly that. Remove it at N≥5 and the generic expectation is: no diabolic. The answer turned out richer than either "yes" or "no", and it reconciles the codimension argument with integrability.

## The instrument, validated at path-3

[`pkmono --diabolic`](../compute/RCPsiSquared.Cli/Commands/PathKMonodromyScanCommand.cs) (engine: [`PathKMonodromyScout.FindDiabolics`](../compute/RCPsiSquared.Diagnostics/Foundation/PathKMonodromyScout.cs)) hunts the residual block's coalescences: it scans the complex-q plane for full-block min-gap minima (no per-q AT stripping, which is valid only at q=2; the residual SET is tracked from q=2 by continuity, since F_AT·F_residual are distinct factors), then classifies each candidate by three independent signals:

- the **monodromy loop** (a small q-loop: identity ⟹ no √-branch ⟹ diabolic candidate; transposition ⟹ defective EP),
- the **gap-scaling exponent** (gap ∝ |q−q\*|^p: p ≈ 1 linear ⟹ two sheets crossing ⟹ diabolic; p ≈ ½ ⟹ defective √-branch),
- the **EpCharacter** verdict (Riesz projector: geo = alg, departure ≈ 0 ⟹ semisimple/diabolic; geo < alg ⟹ defective Jordan block).

At path-3 the three agree and reproduce the **known** diabolic structure of the squared discriminant factor (3q⁴+q²−1)² exactly: the real pair q ≈ ±0.659 (λ = −4 ± 1.318i) **and** the imaginary pair q ≈ ±0.876i (λ = −2.248 / −5.752), all identity-loop, exponent 1.00; the genuine EPs are correctly transposition-loop, exponent ½. This is the [`PathKDiabolicTests`](../compute/RCPsiSquared.Diagnostics.Tests/Foundation/PathKDiabolicTests.cs) gate (the path-k path reproduces `OcticRootsAt`; the classifier re-finds q_EP as diabolic and a generic eigenvalue as not).

## The path-4 (N=5) result

Diabolics **exist** at path-4, and there are several. But the decisive observation is *where in q*: **none sits at real (physical) q.** The clean diabolics (all three signals agreeing) are conjugate pairs at complex q:

| q | merged λ | signals |
|---|---|---|
| ± 0.879 i | −6.662 | identity loop, exp 1.00 (imaginary-q, analytic continuation, like path-3's ±0.876i) |
| ± 1.695 i | −0.383 | identity loop, exp 1.00 (imaginary-q) |
| 0.6407 ± 0.180 i | −4.077 ± 1.115 i | identity loop, exp 0.98 (clean) |
| 0.7654 ± 0.024 i | −4.371 ∓ 2.056 i | identity loop, exp 0.98 (clean) |
| 1.9447 ± 1.217 i | −2.455 ∓ 3.473 i | identity loop, exp 1.02 (clean) |

A fine real-axis strip scan (`--re 0.3,2.5 --im -0.25,0.25 --cell 0.01`) finds **no diabolic at im(q) = 0**; the only real-q coalescence is a **defective** EP (q ≈ 1.0776, λ ≈ −3.792). The closest approach to the physical axis is a conjugate pair at q = 0.6118 ± 0.012i, and a tight zoom (`--cell 0.002`) confirms it stays **split** at ± 0.012i rather than merging onto the axis. So at every physical (real) q the N=5 chain has only loud defective EPs; the silent diabolics have retreated into complex q.

## The synthesis: two causes that coincided at N=4

From below, the (SE,DE) block eigenvalue is λ = −2γ⟨n_XY⟩ + i(E_SE − E_DE), and by free-fermion integrability E_DE = ε_j + ε_k is a **sum** of single-particle energies ε_k = 2J cos(kπ/(N+1)) ([the AT](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) for the rate, the Slater additivity for the frequency, documented for path-4 in [F89_TOPOLOGY_ORBIT_CLOSURE.md](F89_TOPOLOGY_ORBIT_CLOSURE.md)). An additive spectrum has abundant accidental coincidences, and they are **semisimple** because the Hamiltonian and the diagonal dissipator do not mix the occupation-number sectors: two crossing modes come from different sectors, so their eigenvectors stay independent. These are the level **crossings** of an integrable (Poisson-like, already read by the galoischaos witness) spectrum, the opposite of the level repulsion (avoided crossings = defective EPs) of a chaotic one. Turning on XXZ anisotropy Δ ≠ 0 breaks the sectors and flips diabolic → defective, exactly [DIABOLIC_BY_INTEGRABILITY](../hypotheses/DIABOLIC_BY_INTEGRABILITY.md)'s gate.

So **integrability is the symmetry R-2 named**: it auto-satisfies the node conditions and makes diabolics abundant, at all N, throughout the complex-q plane. What it does *not* do is put them at real q. That was the **self-fold's** job, and only at N=4: the half-filling DE = bar(DE) antiunitary fold makes one crossing-q self-conjugate, i.e. real (q_EP = 0.659, the real root of 3q⁴+q²−1). Without the self-fold (N≥5) the crossings stay in conjugate pairs straddling the real axis, never landing on it; q = 0.6118 ± 0.012i is the ghost of where one would be.

This reconciles the surprise with the reviewed prior: **R-1/R-2 (codimension and the self-fold) govern real-axis placement, integrability governs existence.** Both were right about different things. The N=4 "diabolic on the rung-2 line −4 at physical q" was the coincidence of the two; at N=5 they separate, the diabolics scatter in λ and leave the physical axis, and a physicist tuning J/γ meets loud defective EPs where the four-site chain had a silent crossing.

## Resolved during the investigation

- **The two near-axis pairs (q ≈ 0.5915 ± 0.091i, 0.6118 ± 0.012i) ARE true diabolics.** Their EpCharacter and gap-exponent (≈ 1) read diabolic, but the fixed 0.02 monodromy loop read a transposition. A loop-radius sweep settled it: the loop is the **identity at small radius** (r ≤ 0.008) and only becomes a transposition once a neighbouring defective EP enters the annulus (r ≥ 0.012), while a genuinely defective EP (the control at q ≈ 0.9938) is a transposition down to r = 0.002. So it was loop contamination by a dense-spectrum neighbour, not a true defect. The fix reads the candidate's **intrinsic** (small-radius) monodromy; all 11 diabolics now classify loop-identity, consistent with character and exponent. Locked by the `Path4_NearAxisDiabolic_IsLoopIdentity` regression test.

## Open item (do not over-read past this)

- **Completeness is not claimed.** The scanned region was bounded; the full count of complex-q diabolics (the analogue of path-3's four) is not established. An exact discriminant factorization of F_18 (Route B) would settle existence and count definitively but is likely infeasible at degree 18.

## Reproduce

```
# path-3 validation (the known 4 diabolics: 2 real + 2 imaginary)
dotnet run --project compute/RCPsiSquared.Cli -c Release -- pkmono --diabolic --k 3 --re -1.2,1.2 --im -1.2,1.2 --cell 0.04
# path-4: diabolics exist at complex q
dotnet run --project compute/RCPsiSquared.Cli -c Release -- pkmono --diabolic --k 4 --re 0.2,3 --im -1.5,1.5 --cell 0.05
# the decisive real-q check (none at im(q)=0) + the tight zoom
dotnet run --project compute/RCPsiSquared.Cli -c Release -- pkmono --diabolic --k 4 --re 0.3,2.5 --im -0.25,0.25 --cell 0.01
dotnet run --project compute/RCPsiSquared.Cli -c Release -- pkmono --diabolic --k 4 --re 0.595,0.63 --im -0.05,0.05 --cell 0.002
```

Gate tests: `dotnet test compute/RCPsiSquared.Diagnostics.Tests -c Release --filter "FullyQualifiedName~PathKDiabolicTests"`.

## See also

- [zeros_connecting_structure arc](../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs) (this is its forward edge, now answered)
- [DIABOLIC_BY_INTEGRABILITY](../hypotheses/DIABOLIC_BY_INTEGRABILITY.md) (the N=4 mechanism this generalizes)
- [F89_BRANCH_LOCUS_PALINDROME](F89_BRANCH_LOCUS_PALINDROME.md) (the palindrome that folds the locus; path-3's 4 diabolics)
- [F89_TOPOLOGY_ORBIT_CLOSURE](F89_TOPOLOGY_ORBIT_CLOSURE.md) (the path-k residual degrees and Slater additivity)
- the reviewed investigation plan (local, `docs/superpowers/plans/`; R-1/R-2 and the gate-first design)
