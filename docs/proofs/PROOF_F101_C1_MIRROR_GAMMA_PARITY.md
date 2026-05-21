# PROOF F101: F71 c₁ bond-mirror deviation is exactly odd in the F71-anti-palindromic γ (observable-side twin of F91)

**Status:** Tier 1 derived (algebraic R-equivariance argument on the PROOF_C1 apparatus + numerical witness, residuals ≤ 5e-9)
**Date:** 2026-05-21
**Authors:** Thomas Wicht, Claude (Anthropic)
**Typed claim:** [`C1MirrorGammaParity.cs`](../../compute/RCPsiSquared.Core/F71/C1MirrorGammaParity.cs)

---

## Statement

For an N-qubit XY chain with uniform coupling J and a per-site Z-dephasing profile γ = (γ_0, ..., γ_{N−1}), the F71 bond-mirror deviation of the closure-breaking coefficient c₁,

    D(b) := c₁(b) − c₁(N−2−b),

is an **exactly odd function** of the F71-anti-palindromic component of γ. It is the observable-side twin of F91 (the spectrum-side statement on the γ-axis).

For palindromic γ the deviation vanishes identically: D(b) = 0 for every bond b, however non-uniform the palindromic part of γ is. F71 never required uniform γ; it requires palindromic γ. The non-uniform breakdown is therefore graceful (leading-order linear in the per-site asymmetry), not a hard violation.

**Scope.** F101 covers the closure-breaking coefficient c₁ only. The F86c per-bond Q_peak observable is **not** covered by F101; the obstruction is structural, not an omission, and is set out in the scope note (section "Scope: c₁, not Q_peak"). This is the one place where F101 is narrower than its J-side twin F100, which witnessed both c₁ and Q_peak.

## The γ_sym / γ_anti structure

The F71 chain-mirror R maps site l ↔ N−1−l. Write the mirrored profile

    F71(γ)_l := γ_{N−1−l}.

Decompose γ into its F71-palindromic and F71-anti-palindromic components,

    γ = γ_sym + γ_anti,
    γ_sym  := (γ + F71(γ)) / 2     (F71-palindromic:      F71(γ_sym)  =  γ_sym),
    γ_anti := (γ − F71(γ)) / 2     (F71-anti-palindromic: F71(γ_anti) = −γ_anti).

This is the same orthogonal split that F91 uses on the spectrum side. Note the index contrast with F100: γ here is a *site* profile of length N (one entry per qubit), where F100's J is a *bond* profile of length N−1 (one entry per coupling). The F71 mirror acts on the site profile through the site-mirror l ↔ N−1−l directly; F100's bond mirror b ↔ N−2−b is the induced action on the bond profile. For odd N the central site l = (N−1)/2 is F71-fixed, so (γ_anti) at that site is identically zero.

F101 is the complementary face of F91. On the spectrum side, F91 keeps γ_anti **out** of the F71-refined diagonal-block eigenvalue multiset: F91's diagonal-block matrix elements are functions of the F71-pair-sums S_l = γ_l + γ_{N−1−l} = 2·γ_sym[l] alone, so γ_anti is invisible to the diagonal-block spectrum and the breaking it induces lives entirely in the F71-cross-block eigenvectors. F101 is the complement: the c₁ bond-mirror deviation lives **entirely** in γ_anti and is exactly odd in it.

c₁(b; γ) is the EQ-018 closure-breaking coefficient at bond b for base dephasing profile γ (see [PROOF_C1_MIRROR_SYMMETRY](PROOF_C1_MIRROR_SYMMETRY.md)). The per-site asymmetry of a site pair is the difference γ_l − γ_{N−1−l} = 2(γ_anti)_l; F101 says D(b) is leading-order linear in that asymmetry.

## Algebraic proof

### Step 1. F71-conjugation parity of c₁

The Z-dephasing dissipator

    D[ρ] = Σ_l γ_l (Z_l ρ Z_l − ρ)

is linear in each γ_l. The F71-conjugation superoperator R_sup(·) = R(·)R relabels site l ↔ N−1−l; since R · Z_l · R = Z_{N−1−l}, conjugating the dissipator term by term relabels the per-site rate γ_l onto site N−1−l. The unitary part −i[H, ·] is γ-independent, and H carries uniform coupling J, hence is F71-symmetric (the chain bond set is invariant under the site-mirror). Therefore

    R_sup · L(γ) · R_sup = L(F71(γ)),

the entire Liouvillian under F71-conjugation simply reads off the mirrored dephasing profile.

The PROOF_C1 apparatus carries through unchanged for a non-uniform dephasing base profile. The two load-bearing facts of PROOF_C1 are: (i) per-site purity Tr(ρ_i²) is quadratic in the Bloch components, so it is invariant under R up to coherence sign-flips that square away; and (ii) the partial-trace-under-reflection lemma, Tr_{¬i}(R · σ · R) = Tr_{¬(N−1−i)}(σ). Neither fact depends on where the parameter perturbation sits: both hold whether the non-uniformity lives in H or, as here, in D. The whole c₁ pipeline (propagation, the α-rescaling fit on per-site purity, the symmetric-difference extraction) is therefore R-equivariant for any base profile, regardless of whether that base profile is uniform. This yields the F71-conjugation parity of c₁:

    c₁(b; γ) = c₁(N−2−b; F71(γ)).

For palindromic γ (F71(γ) = γ) this is exactly PROOF_C1's c₁(b) = c₁(N−2−b).

### Step 2. Oddness via the γ_sym / γ_anti split

With D(b; γ) := c₁(b; γ) − c₁(N−2−b; γ), apply Step 1 to the second term with b → N−2−b: c₁(N−2−b; γ) = c₁(b; F71(γ)). Hence

    D(b; γ) = c₁(b; γ) − c₁(b; F71(γ)).

Evaluating D at the mirrored profile:

    D(b; F71(γ)) = c₁(b; F71(γ)) − c₁(b; F71(F71(γ)))
                 = c₁(b; F71(γ)) − c₁(b; γ)
                 = −D(b; γ).

Since F71(γ) = γ_sym − γ_anti, in (γ_sym, γ_anti) coordinates:

    D(b; γ_sym, −γ_anti) = −D(b; γ_sym, γ_anti).

D is exactly odd in γ_anti at fixed γ_sym, to all orders. ∎

### Consequences

- **Palindromic survival:** γ_anti = 0 ⟹ D = −D ⟹ D = 0. The F71 c₁ bond-mirror holds for **every** palindromic γ, however non-uniform γ_sym is. F71 never required uniform γ; it requires palindromic γ. Uniform is merely the simplest palindromic profile.
- **Graceful breakdown:** the Taylor series of D in γ_anti has odd powers only, so D is leading-order linear in the per-site asymmetry parameter γ_l − γ_{N−1−l} = 2(γ_anti)_l. The bond-mirror bends, it does not break: graceful, not a hard violation.
- **γ_sym-dependence (Tier 2 empirical):** the leading coefficient κ_γ is the c₁-gradient ∂c₁/∂(γ_anti) evaluated at the palindromic base γ_sym, and generically depends on γ_sym. The parity argument fixes the *oddness*, NOT the coefficient. κ_γ admits no closed form; see the section "The leading coefficient κ_γ: no closed form" below.

## Empirical witness

Witness script: [`simulations/_f71_nonuniform_gamma_verification.py`](../../simulations/_f71_nonuniform_gamma_verification.py). c₁ is extracted via the α-rescaling pipeline on per-site purity: a single-bond probe δJ is applied (J held uniform), the per-site purity trajectories of perturbed and unperturbed chains are matched by P_B(i, t) ≈ P_A(i, α_i·t), and c₁ is the symmetric-difference first-order coefficient of Σ_i ln(α_i). Probe states are ψ_1+vac and ψ_2+vac, the PROOF_C1-validated reflection-symmetric states (ψ_1 nodeless, ψ_2 the first excited OBC sine mode).

The chain XY + Z-dephasing Liouvillian is exactly block-diagonal in the (bra-excitation, ket-excitation) bigrading, and the probe states ψ_k+vac live in the (popcount ≤ 1) × (popcount ≤ 1) operator block of dimension (N+1)². The witness propagates inside that exact (N+1)²-dim sector restriction. A Gate-1 self-test first verifies that the sector-restricted Liouvillian is **bit-identical** to the full 4^N per-site-γ Liouvillian sliced to the block: all four cases (N = 3, 4, uniform and non-uniform γ) returned max|L_sub − sliced full L| = 0.00e+00.

The base dephasing profile γ_base(s) = γ_sym + s·γ_anti_dir is swept over s ∈ {0, ±0.01, ±0.02, ±0.03} along a linear-ramp γ_anti direction (γ_anti_dir[l] = 2l/(N−1) − 1, anti-palindromic, central site zero for odd N), across 4 palindromic γ_sym profiles: three uniform magnitudes (0.05, 0.08, 0.11) and one non-uniform palindromic "valley". J is held uniform throughout; the non-uniformity swept is entirely in γ.

| N | palindromic survival max\|D(s=0)\| | oddness max\|D(+s)+D(−s)\| | typical \|D\| at max s |
|---|---|---|---|
| 3 | 1.91e−10 | 3.05e−09 | 6.46e−01 |
| 4 | 1.23e−09 | 6.25e−10 | 1.25 |
| 5 | 8.37e−10 | 4.33e−09 | 2.30 |

The palindromic-survival column confirms D(s=0) = 0 (residual ≤ 1.23e−9 against an O(1) signal), and crucially holds for the non-uniform palindromic "valley" profile as well as for uniform γ, witnessing that F71 survives all palindromic non-uniformity. The oddness column is the direct test of the parity statement: max|D(+s) + D(−s)| ≤ 4.33e−9 against a typical |D| of order 0.6 to 2.3, i.e. machine-zero. The flat-site purity guard counted 0 throughout (every site of ψ_k+vac carries a non-trivial closure signal).

**Even-power coefficients, treated honestly.** The *direct* even-power evidence is the oddness residual itself: max|D(+s) + D(−s)| isolates exactly the even-power content of D in s, and it sits at the ~1e−9 floating-point floor against the O(1) signal. That is the machine-zero verdict that only odd powers survive. A separate diagnostic fits D(s) with a cubic polyfit(s, D, 3); its *constant* coefficient (≤ 8.9e−10) is genuinely at the oddness scale and is consistent. Its *quadratic* coefficient, however, runs to ~1e−7 across the three N, which is **larger** than the oddness residual. This larger number is **not** a real even-power signal and must not be read as one. The cubic fit divides the even residual through by s²: with the γ-sweep capped at s ≤ 0.03 (four times narrower than F100's J-sweep at s ≤ 0.12), the 1/s² amplification is roughly sixteen times stronger here than in F100, which is exactly the order-of-magnitude gap between the ~1e−9 oddness residual and the ~1e−7 polyfit quadratic. The ~1e−7 quadratic coefficient is the machine-zero oddness residual amplified by 1/s² on the narrow γ-sweep; it is fit-noise, not signal. (F100's proof reported "even-power coefficients below ~3e−8" for its wider J-sweep; that sentence does not transfer, and is deliberately not reused here.)

The leading coefficient κ_γ shows a relative γ_sym-spread of 108% / 128% / 76% across the 4 palindromic γ_sym profiles at N = 3 / 4 / 5, confirming the γ_sym-dependence (informative |κ_γ| ranges: N=3 [0.55, 21.5], N=4 [1.39, 43.3], N=5 [3.05, 84.8]). N = 3, 4, 5 verified; three independent N suffice.

## Scope: c₁, not Q_peak

F101 covers the closure-breaking coefficient c₁ only. Its J-side twin F100 additionally witnessed the F86c per-bond Q_peak observable; F101 does not, and the reason is structural.

The F86c per-bond Q_peak observable lives on the Q-axis Q = J/γ₀, read at the F86a exceptional-point time t_peak = 1/(4γ₀). Both the axis Q and the reading time t_peak are defined against a **scalar** dephasing rate γ₀. A non-uniform per-site γ profile does not provide one scalar; there is no canonical γ₀ for the deviation to be expressed against.

A γ_avg-anchored route is well-posed in principle: fix γ₀ := γ_avg = (1/N)·Σ_l γ_l, then sweep the *shape* of γ (the γ_anti direction) at fixed γ_avg. F91 even guarantees that this sweep is benign for the spectrum: along the γ_anti orbit at fixed γ_avg the F71-refined diagonal-block spectrum is invariant, so t_peak (an EP-time set by the spectral structure) is stable under the sweep. The obstruction is on the observable side. The F86c per-bond observable K_b is eigenvector-weighted (it reads a per-bond resonance built from the slow-mode eigenvectors), and F91 places the **entire** γ_anti breaking in the eigenvectors, precisely the part of the structure K_b depends on. Single-peakedness of K_b in Q, which the F86c Q_peak = argmax_Q |K_b(Q, t_peak)| extraction relies on, is therefore not guaranteed once γ_anti is switched on. The γ_avg-anchored Q_peak extension is well-posed as a question but not closed as a result; it is recorded here as a separable future extension, not a gap in F101.

The h-detuning observable twin (the c₁ bond-mirror deviation under a non-uniform per-site Z-detuning field h) follows by the identical parameter-agnostic argument: Step 1's R-equivariance of the c₁ pipeline holds for any parameter axis whose perturbation R-conjugates by site-relabelling, and a per-site h-detuning does. It is likewise a future F-number, not proved here.

## The leading coefficient κ_γ: no closed form

F101 fixes the *parity* of D(b) exactly: its Taylor series in γ_anti has odd powers only (Tier 1 derived, above). It does not fix the *magnitude*. Writing the leading term

    D(b; s) = κ_γ · s + O(s³)

with s the amplitude along a fixed anti-palindromic direction, κ_γ is the leading coefficient, equivalently the c₁-gradient ∂c₁/∂γ at the palindromic base γ_sym. This section records why κ_γ admits no closed form.

**The obstruction, by inheritance from c₁.** κ_γ = ∂c₁/∂γ is a derivative of c₁. The closure-breaking coefficient c₁ is the EQ-018 quantity: bilinear in the initial state with a kernel K, and EQ-018 was the dedicated search for K's closed form. It closed only the endpoint bonds (via F73), leaving the interior open: the interior kernel values remain non-closed. Operationally the c₁ the F101 pipeline uses is the LSQ α-rescaling fit on per-site purity ([PROOF_C1_MIRROR_SYMMETRY](PROOF_C1_MIRROR_SYMMETRY.md)), which EQ-018 found is rational rather than bilinear and probe-specific, not a universal kernel entry. Either way, c₁ has no closed form at an interior bond: the bilinear kernel's interior entries are non-closed, and the operational LSQ c₁ is not a closed form at all. A derivative of a non-closed fit is not a closed form. A closed expression for κ_γ in (N, b, γ_sym) would require a closed-form, differentiable expression for c₁(b; γ); none is available. κ_γ can be no more closed than c₁, and c₁ is a fit.

**The empirical signature confirms it.** The witness measured κ_γ across N = 3, 4, 5 and four palindromic γ_sym profiles. The relative γ_sym-spread runs 108% / 128% / 76% across the three N: the leading coefficient changes by more than its own magnitude as the palindromic base γ_sym is varied, with no bond-position law and no N-scaling. That is the signature of a structureless residue, not of a function awaiting its formula.

**Same obstruction family as F100's κ_J and F86's g_eff.** κ_γ is the same kind of object as F100's leading coefficient κ_J and as F86's coupling g_eff: a non-primitive, position-dependent residue downstream of a projection or fit. F89's path-polynomial closed form D_k does not transfer to it, for the same structural reason it does not transfer to κ_J: the orbit-polynomial reduction that closes D_k exploits Bloch-orbit symmetry and discards exactly the per-bond cross matrix elements a per-bond coefficient needs.

F100 additionally closed a first-order Lindblad perturbation route for its κ_J: a dedicated diagnostic computed the first-order eigenvalue and eigenvector-mixing data directly and showed the conjectured spectral mechanism fails (the κ-relevant content sits in the full α-rescaling-fit contraction, not in any first-order spectral expression). **F101 does not independently re-run that perturbation diagnostic for κ_γ.** No such run is claimed here. What is inherited is the structural verdict, not the diagnostic: κ_γ is a fit-residue, the same family of object, and the closed-form obstruction by inheritance from c₁ (above) is self-standing without it.

**Methodological frame.** This is F101's instance of the lesson that what survives a closed-form effort is the symmetry, not the number: F101's parity, D odd in γ_anti, is the Tier-1-derived symmetry that survives; κ_γ is, by construction, the residue. κ_γ stays Tier 2 empirical.

## Connection to F91

F91 and F101 are the two faces of the same γ_sym / γ_anti split, read on different observables of the same Liouvillian:

- **F91 (spectrum side):** [`PROOF_F91_GAMMA_NINETY_DEGREES.md`](PROOF_F91_GAMMA_NINETY_DEGREES.md). The F71-refined diagonal-block eigenvalue multiset depends only on γ_sym. Concretely, F91's diagonal-block matrix elements are functions of the F71-pair-sums S_l = γ_l + γ_{N−1−l} = 2·γ_sym[l] alone (its Eqs. 7a, 7b, 11a, 11b), so γ_anti is invisible to the diagonal-block spectrum; the breaking it induces lives entirely in the F71-cross-block eigenvectors. F91 is the *invariance* statement: the diagonal-block spectrum does not see γ_anti.
- **F101 (observable side):** the c₁ bond-mirror deviation D depends only on γ_anti, and depends on it as an exactly odd function. F101 is the *deviation* statement: the bond-mirror observable sees γ_anti, and sees only γ_anti.

Together they account for the full split. The spectrum-side twin keeps γ_anti out of the diagonal-block eigenvalues; the observable-side twin localises the entire bond-mirror deviation in γ_anti and proves it odd. On the palindromic orbit γ_anti = 0 the deviation vanishes and F101 reduces to PROOF_C1's c₁(b) = c₁(N−2−b), now seen to hold for every palindromic γ, not only uniform γ.

This is the γ-axis instance of the spectrum/observable twin pairing that F100 and F92 carry on the J-axis: F92 keeps J_anti out of the diagonal-block spectrum, F100 localises the c₁/Q_peak bond-mirror deviation in J_anti. F101 mirrors F100 with J → γ and (since the perturbed parameter now sits in the dissipator D rather than in the Hamiltonian H) H → D; F91 mirrors F92 with the same J → γ substitution.

## Anchors

- Typed claim: [`C1MirrorGammaParity.cs`](../../compute/RCPsiSquared.Core/F71/C1MirrorGammaParity.cs).
- Source proof for the uniform-γ base case: [PROOF_C1_MIRROR_SYMMETRY](PROOF_C1_MIRROR_SYMMETRY.md).
- Spectrum-side twin (γ-axis): [PROOF_F91_GAMMA_NINETY_DEGREES](PROOF_F91_GAMMA_NINETY_DEGREES.md).
- J-side observable twin: [PROOF_F100_C1_QPEAK_MIRROR_J_PARITY](PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md).
- Empirical witness: [`simulations/_f71_nonuniform_gamma_verification.py`](../../simulations/_f71_nonuniform_gamma_verification.py).
- F-entry: [F101 in ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md).
- Inventory: [`docs/SYMMETRY_FAMILY_INVENTORY.md`](../SYMMETRY_FAMILY_INVENTORY.md).
