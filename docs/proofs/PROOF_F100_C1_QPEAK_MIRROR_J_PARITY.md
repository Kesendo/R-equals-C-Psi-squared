# PROOF F100: F71 c₁/Q_peak bond-mirror deviation is exactly odd in the F71-anti-palindromic J (observable-side twin of F92)

**Status:** Tier 1 derived (algebraic R-equivariance argument on the PROOF_C1 apparatus + numerical empirical witness at N=3, 4, 5, residuals < 1e-9)
**Date:** 2026-05-20; 2026-05-21 (κ-obstruction section added)
**Authors:** Thomas Wicht, Claude (Anthropic)
**Typed claim:** [`C1QPeakMirrorJParity.cs`](../../compute/RCPsiSquared.Core/F71/C1QPeakMirrorJParity.cs)

---

## Statement

For an N-qubit XY chain with uniform Z-dephasing and a bond-coupling profile J = (J_0, ..., J_{N−2}), the F71 bond-mirror deviation of the closure-breaking coefficient c₁ (and of the F86c per-bond Q_peak observable),

    D(b) := c₁(b) − c₁(N−2−b),

is an **exactly odd function** of the F71-anti-palindromic component of J. It is the observable-side twin of F92 (the spectrum-side statement).

For palindromic J the deviation vanishes identically: D(b) = 0 for every bond b, however non-uniform the palindromic part of J is. F71 never required uniform J; it requires palindromic J. The non-uniform breakdown is therefore graceful (leading-order linear in the asymmetry), not a hard violation.

The identical conclusion holds for the F86c per-bond Q_peak observable K_b(Q, t): ΔQ_peak(b) is odd in J_anti and zero for palindromic J.

## The J_sym / J_anti structure

The F71 chain-mirror R maps site i ↔ N−1−i, hence bond b ↔ N−2−b. Write the mirrored profile

    F71(J)_b := J_{N−2−b}.

Decompose J into its F71-palindromic and F71-anti-palindromic components,

    J = J_sym + J_anti,
    J_sym  := (J + F71(J)) / 2     (F71-palindromic:      F71(J_sym)  =  J_sym),
    J_anti := (J − F71(J)) / 2     (F71-anti-palindromic: F71(J_anti) = −J_anti).

This is the same orthogonal split that F92 uses on the spectrum side. There the J_anti component is the *invisible* direction: the diagonal-block eigenvalue multiset depends only on J_sym. F100 is the complementary face: the c₁/Q_peak bond-mirror deviation lives **entirely** in J_anti and is exactly odd in it.

c₁(b; J) is the EQ-018 closure-breaking coefficient at bond b for base profile J (see [PROOF_C1_MIRROR_SYMMETRY](PROOF_C1_MIRROR_SYMMETRY.md)). The asymmetry of a bond pair is the difference B_b = J_b − J_{N−2−b} = 2(J_anti)_b; F100 says D(b) is leading-order linear in B_b.

## Algebraic proof

### Step 1. F71-conjugation parity of c₁

The superoperator R_sup(·) = R(·)R satisfies

    R_sup · L(J) · R_sup = L(F71(J)),

because R relabels the coupling at bond b to bond N−2−b. The PROOF_C1 apparatus (per-site purity is invariant under R up to coherence sign-flips that square away, plus the partial-trace-under-reflection lemma) carries through unchanged for a non-uniform base profile, yielding R-equivariance of the entire c₁ pipeline:

    c₁(b; J) = c₁(N−2−b; F71(J)).

For palindromic J (F71(J) = J) this is exactly PROOF_C1's c₁(b) = c₁(N−2−b).

### Step 2. Oddness via the J_sym / J_anti split

With D(b; J) := c₁(b; J) − c₁(N−2−b; J), apply Step 1 with b → N−2−b: c₁(N−2−b; J) = c₁(b; F71(J)). Hence

    D(b; J) = c₁(b; J) − c₁(b; F71(J)).

Evaluating D at the mirrored profile:

    D(b; F71(J)) = c₁(b; F71(J)) − c₁(b; F71(F71(J)))
                 = c₁(b; F71(J)) − c₁(b; J)
                 = −D(b; J).

Since F71(J) = J_sym − J_anti, in (J_sym, J_anti) coordinates:

    D(b; J_sym, −J_anti) = −D(b; J_sym, J_anti).

D is exactly odd in J_anti at fixed J_sym, to all orders. ∎

### Consequences

- **Palindromic survival:** J_anti = 0 ⟹ D = −D ⟹ D = 0. The F71 c₁/Q_peak bond-mirror holds for **every** palindromic J, however non-uniform J_sym is. F71 never required uniform J; it requires palindromic J. Uniform is merely the simplest palindromic profile.
- **Graceful breakdown:** the Taylor series of D in J_anti has odd powers only, so D is leading-order linear in the asymmetry parameter B_b = J_b − J_{N−2−b} = 2(J_anti)_b. Graceful, not a hard violation.
- **J_sym-dependence (Tier 2 empirical):** the leading coefficient κ_b is the c₁-gradient evaluated at J_sym and generically depends on J_sym. The parity argument fixes the oddness, NOT the coefficient. κ_b admits no closed form; see the section "The leading coefficient κ" below.
- **Q_peak:** the identical R-conjugation argument applies to F86c's per-bond observable K_b(Q, t), built R-equivariantly: K_b(Q, t; J) = K_{N−2−b}(Q, t; F71(J)), so ΔQ_peak(b) is odd in J_anti, zero for palindromic J. c₁ is numerically verified; Q_peak follows by the identical argument.

## Empirical witness

Witness script: [`simulations/_f71_nonuniform_j_verification.py`](../../simulations/_f71_nonuniform_j_verification.py). c₁ is extracted via the α-rescaling pipeline on per-site purity, probe states ψ_1+vac and ψ_2+vac. The base profile J = J_sym + s·J_anti_dir is swept over s ∈ {0, ±0.04, ±0.08, ±0.12} with a linear-ramp J_anti direction, across 4 palindromic J_sym profiles (uniform 0.8 / 1.0 / 1.2 plus one non-uniform palindromic "valley"). The full 4^N Liouvillian is used, no truncation.

| N | palindromic survival max\|D(s=0)\| | oddness max\|D(+s)+D(−s)\| | typ \|D\| at max s |
|---|---|---|---|
| 3 | 1.5e−11 | 7.8e−11 | 0.54 |
| 4 | 6.1e−11 | 5.5e−10 | 0.99 |
| 5 | 4.0e−10 | 1.0e−9 | 4.05 |

Even-power coefficients (constant and quadratic, from a cubic fit of D vs s) are below ~3e−8, confirming only odd powers survive. The leading coefficient κ shows 76% / 62% / 143% relative spread across the 4 J_sym profiles at N=3 / 4 / 5, confirming the J_sym-dependence. N=3, 4, 5 verified; N=6 not required (three independent N suffice).

## Connection to F92

F92 and F100 are the two faces of the same J_sym / J_anti split, read on different observables of the same Liouvillian:

- **F92 (spectrum side):** the F71-refined diagonal-block eigenvalue multiset depends only on J_sym. J_anti is invisible to the diagonal-block spectrum; the breaking it induces lives entirely in the F71-cross-block eigenvectors. F92 is the *invariance* statement: the spectrum does not see J_anti.
- **F100 (observable side):** the c₁/Q_peak bond-mirror deviation D depends only on J_anti, and depends on it as an exactly odd function. F100 is the *deviation* statement: the bond-mirror observable sees J_anti, and sees only J_anti.

Together they account for the full split. The spectrum-side twin keeps J_anti out of the eigenvalues; the observable-side twin localises the entire bond-mirror deviation in J_anti and proves it odd. On the palindromic orbit J_anti = 0 the deviation vanishes and F100 reduces to PROOF_C1's c₁(b) = c₁(N−2−b), now seen to hold for every palindromic J, not only uniform J.

## The leading coefficient κ: no closed form

F100 fixes the *parity* of D(b) exactly: its Taylor series in J_anti has odd powers only (Tier 1 derived, above). It does not fix the *magnitude*. Writing the leading term

    D(b; s) = κ_b · s + O(s³)

with s the amplitude along a fixed anti-palindromic direction, κ_b is the leading coefficient, equivalently the c₁-gradient ∂c₁(b)/∂(J_anti) at the palindromic base J_sym. This section records why κ_b admits no closed form (analysis 2026-05-21).

**The obstruction, by inheritance from c₁.** κ_b is a derivative of c₁. The closure-breaking coefficient c₁ is bilinear in the initial state with a kernel K; EQ-018 was the dedicated search for K's closed form, and it closed only the endpoint bonds (via F73), leaving the interior open: "interior kernel values remain non-closed", "their closed form remains open". Operationally the c₁ the F100 pipeline uses is the LSQ α-rescaling fit on per-site purity ([PROOF_C1_MIRROR_SYMMETRY](PROOF_C1_MIRROR_SYMMETRY.md)), which EQ-018 found is "rational, not bilinear" and probe-specific, not a universal kernel entry. Either way, c₁ has no closed form at an interior bond: the bilinear kernel's interior entries are non-closed, and the operational LSQ c₁ is not a closed form at all. A closed form for κ_b in (N, b, J_sym) would require a closed-form, differentiable expression for c₁(b; J); none is available. κ_b can be no more closed than c₁, and c₁ is a fit. The argument is self-standing: it does not route through F86's g_eff, the exceptional point, or the OBC sine-sum structure.

**The empirical signature confirms it.** The witness [`_f71_nonuniform_j_verification.py`](../../simulations/_f71_nonuniform_j_verification.py) measured κ_b across N=3, 4, 5 and four palindromic J_sym profiles. κ shows no bond-position law, no J_sym law, no N-scaling, sign changes between J_sym magnitudes, and a relative J_sym-spread that reaches 143% at N=5 (76% / 62% / 143% at N=3 / 4 / 5). That is the signature of a structureless residue, not of a function awaiting its formula.

**Same family as F86's g_eff; F89's closure does not rescue it.** κ_b is the same kind of object as F86's coupling g_eff: a non-primitive, bond-position-dependent residue downstream of a projection or fit (cf. [PROOF_F86B_OBSTRUCTION](PROOF_F86B_OBSTRUCTION.md)). F89's path-polynomial closed form D_k (Tier 1 derived) does not transfer to it: F89 closes D_k *via* the orbit-polynomial reduction, which exploits Bloch-orbit symmetry and discards exactly the per-bond cross matrix elements a per-bond coefficient needs (verified negative, [`f86_geff_via_f90_bridge_probe.py`](../../simulations/f86_geff_via_f90_bridge_probe.py)). The route that closes an orbit-symmetric quantity is structurally the route that strips a per-bond one; this is the reduced-model route, obstruction L4 of the F86b catalogue, where every finite reduction proved insufficient for the per-bond signal.

**What κ does have: the open positive direction.** "No closed form" is not "nothing". κ_b is accessible to first-order Lindblad perturbation theory: it can be written as a formal expression in the eigenvalues, slow-mode eigenvectors, and first-order J-derivatives of L(J_sym), the slow-mode eigenvector-mixing series. That expression would *explain* κ's J_sym-dependence (through the eigenvector-mixing matrix elements) and would show whether the large N=5 spread reflects a near-degeneracy in the slow-mode sector, lifting κ from a Tier-2 empirical residue to a Tier-1 structural object expressed in the primitive basis. It is still not a closed form: the palindromic-J slow modes the series is built from have themselves no closed form for general N (open-boundary XY with arbitrary palindromic couplings is not exactly solvable). The honest open direction for κ is that formal spectral expression, not a closed form.

**Methodological frame.** This is F100's instance of the lesson that what survives a closed-form effort is the symmetry, not the number ([`ON_THE_Q_AXIS_AND_THE_PTF_LESSON`](../../reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md)): F100's parity, D odd in J_anti, is the Tier-1-derived symmetry that survived; κ_b is the number that, by construction, is the residue. κ_b stays Tier 2 empirical.

## Anchors

- Typed claim: [`C1QPeakMirrorJParity.cs`](../../compute/RCPsiSquared.Core/F71/C1QPeakMirrorJParity.cs) (Tier1Derived; static helpers `PalindromicComponent`, `AntiPalindromicComponent`, `PalindromicDeviation`, `IsPalindromic`).
- Source proof for the uniform-J base case: [PROOF_C1_MIRROR_SYMMETRY](PROOF_C1_MIRROR_SYMMETRY.md).
- Spectrum-side twin: [PROOF_F92_BOND_ANTI_PALINDROMIC_J](PROOF_F92_BOND_ANTI_PALINDROMIC_J.md).
- F86c per-bond Q_peak mirror (uniform J): [PROOF_F86C_F71_MIRROR](PROOF_F86C_F71_MIRROR.md).
- Empirical witness: [`simulations/_f71_nonuniform_j_verification.py`](../../simulations/_f71_nonuniform_j_verification.py).
- F-entry: [F100 in ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md).
- Inventory: [`docs/SYMMETRY_FAMILY_INVENTORY.md`](../SYMMETRY_FAMILY_INVENTORY.md).
