# PROOF: the K₁ chiral mirror rate law is a site-wise trajectory identity

**Status:** Tier 1 derived (four-step operator argument, exact for every N, site, and time; verified to machine precision at N = 5 and N = 7, worst deviation 8.9·10⁻¹⁶).
**Date:** 2026-06-10
**Authors:** Thomas Wicht, Claude (Fable 5)
**Builds on:**
- [review/EQ014_FINDINGS.md](../../review/EQ014_FINDINGS.md): the EQ-014 audit that retracted the PTF closure law Σ ln α_i = 0 and isolated the mirror law Σ f_i(ψ_k) = Σ f_i(ψ_{N+1−k}) as the machine-exact survivor.
- [reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md](../../reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md): named the survivor "PTF's surviving Tier-1 law" and identified K₁ as the symmetry behind it.
- [hypotheses/PERSPECTIVAL_TIME_FIELD.md](../../hypotheses/PERSPECTIVAL_TIME_FIELD.md): the PTF time-rescaling framing whose per-site fits α_i define the rates f_i.
- [PROOF_K_PARTNERSHIP.md](PROOF_K_PARTNERSHIP.md) + `compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs`: the same sublattice chirality K₁ on the eigenvalue side (spectrum inversion E_{N+1−k} = −E_k; bipartite ⟹ soft). This proof is its eigenvector/dynamics sibling.

## Abstract

Let H be the uniform XY chain on N sites (J = 1) under uniform Z-dephasing γ₀, V = δJ·½(X_bX_{b+1} + Y_bY_{b+1}) a single J-bond defect, and define the **odd-sublattice Z product**

  **K₁ = Π_{l odd} Z_l**  (sites 0-indexed; K₁² = I).

K₁ anticommutes with every XY bond term (exactly one endpoint of each bond is odd), so K₁HK₁ = −H and K₁VK₁ = −V, while K₁ commutes with every dephasing operator Z_l. Combining this conjugation with complex conjugation (H and V are real in the computational basis, and the Z-dephasing action is real) yields the **trajectory identity**: for any real initial state ψ,

  **P_i(t; H+V, K₁ψ) = P_i(t; H+V, ψ)**  for every site i and every time t, exactly,

where P_i(t) = Tr(ρ_i(t)²) is the site-i purity under the full Lindblad evolution. The single-excitation modes ψ_k(l) ∝ sin(πk(l+1)/(N+1)) satisfy K₁ψ_k = ψ_{N+1−k} exactly, and for the PTF pair states φ_k = (|vac⟩ + |ψ_k⟩)/√2 any relative sign a sublattice product may introduce is absorbed by the U(1) phase e^{iπN̂}, which commutes with H, V, and the dissipator. Hence

  **P_i(t; φ_k) = P_i(t; φ_{N+1−k})**  for every i and t.

The PTF per-site rescaling fits α_i, and therefore the first-order rates f_i = (α_i − 1)/(δJ/J), are functionals of these purity trajectories, so they are **site-wise** equal between k and N+1−k. The published mirror law

  Σ_i f_i(ψ_k) = Σ_i f_i(ψ_{N+1−k})  (EQ-014, machine-exact at N = 5, 7, 8)

is the corollary obtained by summing over sites. The law that survived the EQ-014 retraction is not a statement about fitted rates at all; it is a site-wise identity of the trajectories the fits are read from.

## §1 Step 1: the algebra (exact)

The XY chain is bipartite: every bond (l, l+1) joins an even site to an odd site. The bond term X_lX_{l+1} + Y_lY_{l+1} contains exactly one Pauli letter (X or Y) on an odd site, and Z anticommutes with X and Y on that site while the even-site letter is untouched, so K₁ flips the sign of every bond term:

  K₁ H K₁ = −H,  K₁ V K₁ = −V.

(The defect bond (b, b+1) is itself a bond, so the same one-odd-endpoint count applies; no uniformity is used.) K₁ is a product of Z's, so [K₁, Z_l] = 0 for every l: the dephasing dissipator D(ρ) = Σ_l γ_l(Z_lρZ_l − ρ) is K₁-invariant. Verified machine-exact at N = 5 (block 1 of the verifier).

## §2 Step 2: unitary conjugation

Let ρ(t) solve the Lindblad equation with Hamiltonian H+V and dissipator D. Then ρ̃(t) = K₁ρ(t)K₁ solves the Lindblad equation with Hamiltonian K₁(H+V)K₁ = −(H+V) and the same dissipator (D commutes with conjugation by K₁ termwise). K₁ acts on each site as I or Z, so each single-site reduced state transforms unitarily (ρ_i → Z ρ_i Z or ρ_i) and site purities are K₁-invariant:

  P_i(t; −(H+V), K₁ψ) = P_i(t; H+V, ψ).

## §3 Step 3: complex conjugation flips the sign back

H and V are real matrices in the computational basis (XX+YY hopping has real entries), and the Z-dephasing dissipator has a real action (Z real, rates real). Conjugating the Lindblad equation elementwise maps a solution for Hamiltonian −H' to a solution for +H' on the conjugated state: the commutator term −i[−H', ρ] conjugates to −i[H', ρ*] and D(ρ)* = D(ρ*). For a real initial state ρ(0)* = ρ(0), so the two evolutions produce conjugate trajectories, and purities are conjugation-invariant:

  P_i(t; −H', ψ) = P_i(t; H', ψ)  for real ψ.

The pair states φ_k are real (sine amplitudes plus a real vacuum component), and K₁φ_k is again real, so Step 3 applies on top of Step 2.

## §4 Step 4: mode mapping and the U(1) sign absorption

Chaining Steps 2 and 3 gives the trajectory identity P_i(t; H+V, K₁ψ) = P_i(t; H+V, ψ). It remains to identify K₁φ_k. The sine identity

  sin(π(N+1−k)l/(N+1)) = −(−1)^l sin(πkl/(N+1))  (l = 1, …, N the 1-indexed site)

says ψ_{N+1−k} carries the phase pattern +1 on odd 1-indexed sites (= even 0-indexed) relative to ψ_k. The 0-indexed odd-sublattice product applies Z exactly where the excitation phase must flip, so

  **K₁ψ_k = ψ_{N+1−k} exactly, with no sign, for every k**

(this is the gauge the typed `ChiralKClaim` records; block 1 verifies it at machine precision for all k at N = 5). Since K₁|vac⟩ = +|vac⟩, the pair states map cleanly: K₁φ_k = φ_{N+1−k}.

The complementary even-sublattice product K₁ᶜ = Π_{l even} Z_l anticommutes with H and V just as well but picks the minus sign, K₁ᶜψ_k = −ψ_{N+1−k}, giving K₁ᶜφ_k = (|vac⟩ − |ψ_{N+1−k}⟩)/√2. That relative sign is absorbed by the U(1) phase

  e^{iπN̂} = Π_l Z_l  (N̂ the excitation number operator),

which acts as +1 on the vacuum and −1 on the single-excitation sector, commutes with H and V (each bond term flips two sites) and with every Z_l, and acts per site as Z, so it changes neither the dynamics class nor any site purity. Either sublattice product therefore proves the same identity:

  **P_i(t; φ_k) = P_i(t; φ_{N+1−k})  for every site i and time t, exactly.**

The PTF fits α_i (and the rates f_i, the closure sums Σ ln α_i, and every other functional of the per-site purity curves) inherit the equality site-wise; summing over i gives the EQ-014 Σ-mirror law as a corollary.

## §5 Scope

The argument needs exactly four ingredients and nothing else:

1. H and V real in the computational basis;
2. K₁(H+V)K₁ = −(H+V) (one chirality for the unperturbed chain AND the defect);
3. a K₁-invariant dissipator ([K₁, Z_l] = 0 holds for any Z-string conjugator and any site rates γ_l, uniform or not);
4. a real initial state, modulo a U(1) phase that commutes with the dynamics.

It does NOT need: uniformity of J, the specific XY form (any K₁-odd real bipartite hopping works, including site-dependent couplings and multiple defects, as long as every term joins the two sublattices), single-excitation states (any real ψ with K₁ψ identified gives a trajectory identity between ψ and K₁ψ), or first-order perturbation theory (the identity is exact at every δJ, which is why the fitted f_i mirror at finite δJ = 0.1 came out machine-exact in EQ-014 rather than merely first-order-exact).

What breaks it: a Z-field or any K₁-even Hamiltonian component (Step 1 fails), complex hopping phases (Step 3 fails), or a dissipator that is not a Z-string algebra element (T1 amplitude damping breaks Step 2's dissipator invariance).

## §6 History

The law was found 2026-04 in EQ-014 as a Σ-law: the audit that retracted PTF's closure claim (Σ ln α_i = 0 is not a first-order theorem; Σ f_i is state-dependent and O(1)) noticed that the sums were nonetheless **pairwise equal** between ψ_k and ψ_{N+1−k}, machine-exact at N = 5, 7, 8, and `ON_THE_Q_AXIS_AND_THE_PTF_LESSON` named K₁ as the enforcing symmetry without writing the derivation down. It lived untyped for five weeks: in those two documents and a comment in `compute/RCPsiSquared.Diagnostics/Ptf/PerturbationMatrixElements.cs`.

Derived 2026-06-10 as a site-wise trajectory identity, using the involution-plus-sign-table idiom of the windowed-converse wave (PROOF_F87_WINDOWED_MONOMIAL_CONVERSE: pick the involutions, read the sign table, let reality or a U(1) phase close the loop). The derivation strengthened the law twice in passing: from Σ-level to site-wise, and from fitted-rate-level to trajectory-level. The eigenvalue side of the same K₁ (spectrum inversion, bipartite ⟹ soft) was already typed as `ChiralKClaim`; this proof adds the eigenvector/dynamics side, typed as `ChiralMirrorTrajectoryClaim`.

## §7 Verification

- [`simulations/ptf_chiral_mirror_trajectory.py`](../../simulations/ptf_chiral_mirror_trajectory.py): self-validating, four blocks. Algebra at N = 5 (machine-exact); trajectory identity at N = 5 (γ₀ = 0.05, δJ = 0.1, dense expm, worst site-purity deviation 8.9·10⁻¹⁶ over k ∈ {1, 2} and t ≤ 8); the U(1) minus-sign branch via the even-sublattice product (deviation 0 against φ₄, 5.6·10⁻¹⁶ against φ₂); N = 7 sparse spot check (worst 7.8·10⁻¹⁶).
- `compute/RCPsiSquared.Diagnostics.Tests/Ptf/ChiralMirrorTrajectoryClaimTests.cs`: C# algebra tests at N = 4, the trajectory identity through the Core spectral propagator at t ∈ {0.5, 2.0}, and registry wiring (typed parent `ChiralKClaim`).
- The original Σ-law numbers: `review/EQ014_FINDINGS.md` (N = 5, 7, 8; ψ-state table).
