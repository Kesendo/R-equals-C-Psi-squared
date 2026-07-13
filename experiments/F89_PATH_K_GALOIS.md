# F89 path-k Galois: the algebraic certificate that the H_B-mixed factors are unwritable in radicals

*Extracted 2026-07-13 from [F89_TOPOLOGY_ORBIT_CLOSURE.md](F89_TOPOLOGY_ORBIT_CLOSURE.md) as its own
concern (the F89 corpus de-monolith, step 2). That note derives WHERE the octic comes from (the path-3
(SE, DE) S_2-sym factorisation F_a · F_b · F_8 and its path-k siblings, the AT-lock mechanism, the
closed-form amplitudes); this one carries the algebraic verdict on the residual H_B-mixed factors:
**Gal(F_8/Q(i)(q)) = S_8** via the specialization + Dedekind + Jordan certificate, the path-3..6
generalisation S_8/S_18/S_32/S_53 with the live isolate-before-DDF witness, the discriminant's diabolic
double zero at q ≈ 0.659, and the amplitude q-obstruction that follows. The geometric route to the same
S_8 (eigenvalue braids) is the sibling [F89_MONODROMY_MIRROR.md](F89_MONODROMY_MIRROR.md); this note is
the from-above algebraic half of that pair. Notation (q = J/γ, F_a/F_b/F_8, AT-locked, overlap /
no-overlap, diabolic vs defective) is the parent's §Notation.*

**Status:** Tier 1 derived for the Galois verdicts (Gal(F_8) = S_8; path-4/5/6 = S_18/S_32/S_53, all
non-solvable), for the diabolic-degeneracy location q_EP = √((−1+√13)/6) with λ_EP = −4γ + 2iJ, and for
the live isolate-before-DDF pipeline; Tier 2 empirical for the octic-mode amplitude q-dependence.
**Date:** 2026-07-13 (sections moved from the parent, faithful modulo small review-driven precision fixes; their content dates from 2026-06-21..25)
**Authors:** Thomas Wicht, Claude (Fable 5)
**Scripts:**
- [`f89_path3_octic_galois.py`](../simulations/f89_path3_octic_galois.py): the gate-first q0 = 2 Frobenius certificate for Gal(F_8) = S_8 (+ multi-q0 and base-field robustness).
- [`f89_pathk_galois.py`](../simulations/f89_pathk_galois.py): the path-3..6 engine of record (`full-d`): isolate F_d, certify Gal(F_d) = S_d for d = 8/18/32/53.
- [`f89_path3_ep_locator.py`](../simulations/f89_path3_ep_locator.py): the diabolic-degeneracy locator (3q⁴+q²−1 = 0).
- [`f89_jordan_definitive.py`](../simulations/f89_jordan_definitive.py): the artifact-free diabolic-vs-defective test at q_EP (g1=2, dep=0).
- [`f89_path3_octic_amplitude_q_scan.py`](../simulations/f89_path3_octic_amplitude_q_scan.py): the octic-mode amplitude q-scan (no polynomial fit).

**Typed claims / live witnesses:** [`F89Path3OcticGaloisClaim`](../compute/RCPsiSquared.Core/Symmetry/F89Path3OcticGaloisClaim.cs) (live at `inspect --root f89galois`), [`F89PathKHbMixedDegreesClaim`](../compute/RCPsiSquared.Core/Symmetry/F89PathKHbMixedDegreesClaim.cs), [`F89Path3OcticEpClaim`](../compute/RCPsiSquared.Core/Symmetry/F89Path3OcticEpClaim.cs) (`inspect --root epcharacter`, `inspect --root f89octic`).
**Related:** [F89_TOPOLOGY_ORBIT_CLOSURE.md](F89_TOPOLOGY_ORBIT_CLOSURE.md) (the parent: where F_8 is born), [F89_MONODROMY_MIRROR.md](F89_MONODROMY_MIRROR.md) (the same S_8 from eigenvalue braids), [F89_BRANCH_LOCUS_PALINDROME.md](F89_BRANCH_LOCUS_PALINDROME.md) (the branch locus is a palindrome), [F89_TOPOLOGY_CONTROLS_GALOIS_WRITABILITY.md](F89_TOPOLOGY_CONTROLS_GALOIS_WRITABILITY.md) (the four-topology sequel: graph symmetry caps writability), [F89_PATH_K_DIABOLIC.md](F89_PATH_K_DIABOLIC.md) (diabolic points across path-k), [`hypotheses/DIABOLIC_BY_INTEGRABILITY.md`](../hypotheses/DIABOLIC_BY_INTEGRABILITY.md) (why the degeneracy is diabolic), [`reflections/ON_WHAT_CANNOT_CLOSE.md`](../reflections/ON_WHAT_CANNOT_CLOSE.md) (the plain-words reading).

---

## Path-3 octic non-solvability: Gal(F_8) = S_8 (Tier 1 derived)

**The Galois group of the octic over Q(i)(q) is the full symmetric group S_8** ([`f89_path3_octic_galois.py`](../simulations/f89_path3_octic_galois.py), gate-first). The foundation is the degree-52 discriminant:

    disc(F_8) = const · q²⁴ · (3q⁴ + q² − 1)² · P_10(q²)   (const a nonzero normalization-dependent rational; only the factor structure is load-bearing)

where P_10(q²) is degree 10 in q² (degree 20 in q, even powers only) and is NOT a perfect square in Q. The square factor (3q⁴+q²−1)² locates the **diabolic degeneracy** at q² = (−1+√13)/6 ≈ 0.434, q ≈ 0.659. The non-square disc gives **Gal(F_8) ⊄ A_8**, but that alone does not pin the group (S_4 is solvable yet ⊄ A_4).

**The method (specialization + Dedekind + Jordan).** For a good q0 ∈ Q (disc(q0) ≠ 0), the specialized group G_{q0} = Gal(F_8(·,q0)/Q(i)) is a SUBGROUP of the generic G = Gal(F_8/Q(i)(q)) (specialization can only shrink), so proving G_{q0} = S_8 at one good q0 forces G = S_8. Factoring F_8(·,q0) modulo a split prime p ≡ 1 (mod 4) (Z[i]/𝔭 = F_p, i ↦ r with r²≡−1) gives, when squarefree, a Frobenius cycle type. **Certificate at q0 = 2** (F_8(·,2) monic over Z[i]): irreducible over Q(i) (⇒ transitive); the split prime 𝔭|5 (F_5, i↦2) factors it to cycle type **(5,2,1)**, whose square is a 5-cycle (⇒ primitive, since a 5-orbit fits no degree-8 block system, and no proper primitive degree-8 group has order divisible by 5 ⇒ ⊇A_8 by Jordan, 5 ≤ 8−3) and which is itself odd (⇒ ⊄A_8). Hence G_{q0=2} = S_8 ⇒ **Gal(F_8/Q(i)(q)) = S_8**. Confirmed at q0 ∈ {2, 3, ½, 3/2} over Q(i) (14–16 distinct cycle types each) and **robust to enlarging the base** to Q(i,√5) (still irreducible + 5-cycle at q0 ∈ {2,3}).

**Consequence (radical-non-solvability).** S_8 is non-solvable, so the eight roots λ_k(q) admit **no expression by radicals** over Q(i)(q) (Abel-Ruffini). This does NOT exclude non-radical special-function expressions (theta / hypergeometric, the higher-degree analogues of the quintic's Bring radical), which exist for any algebraic function.

**What this means (the content is negative).** S_8 is the *generic* Galois group of an irreducible degree-8 polynomial (van der Waerden 1936; Bhargava, *Annals* 201, 2025), it is not exotic. The point is that free-fermion integrability **spends itself entirely on the factorisation**: the AT-locked F_a/F_b quadratics carry the single-particle frequencies −1±√5 in radicals, and the diabolic point sits on the *solvable* quartic factor (3q⁴+q²−1); the residual octic carries no further algebraic structure. The closed-form program for path-3 terminates exactly at the AT-protected half. (Contrast: the SIC-POVM spectral polynomials, Appleby-Yadsan-Appleby-Zauner 2012, gave a *solvable* Galois group, opposite polarity at the same seam.) Scope: this is the group of the path-3 (SE,DE) S_2-sym octic *factor*, not of "the Liouvillian spectrum"; it is a similarity-invariant of that invariant sub-block. Method reference: K. Conrad, "Recognizing Galois groups S_n and A_n". (Not to be confused with differential Galois theory / "Liouvillian solutions", a different object.)

This S_8 is the small-symmetry (chain) row of a wider four-topology law: on a more symmetric wiring the graph automorphism group caps the (SE,DE) factor degrees, so the relaxation is radically writable exactly when that cap is ≤ 4 (only the complete graph, every N) and grows into S_n only where the symmetry is small (chain, ring). See [F89_TOPOLOGY_CONTROLS_GALOIS_WRITABILITY.md](F89_TOPOLOGY_CONTROLS_GALOIS_WRITABILITY.md).

## Path-3 octic-mode amplitude q-dependence: no closed-form fit (Tier 2)

For each of the 8 octic-derived modes, sigs(N) follows const(q)/[N²(N−1)] (degree-0 polynomial in N). The constant **does NOT admit a polynomial fit ≤ degree 5 in q** ([`f89_path3_octic_amplitude_q_scan.py`](../simulations/f89_path3_octic_amplitude_q_scan.py)):

| q | Σ_8 octic sigs · N²(N−1) | Dominant mode (largest sigs · N²(N−1)) |
|---|---|---|
| 0.5 | 1.68 | mode at rate ≈ 3.6γ, sigs ≈ 1.34 (near the q=0.659 diabolic degeneracy) |
| 0.75 | 1.50 | mode at rate ≈ 3.78γ, sigs ≈ 0.99 (past the degeneracy) |
| 1.0 | 2.20 | mode at rate ≈ 3.35γ, sigs ≈ 1.91 |
| 1.5 | 2.73 | mode at rate ≈ 3.35γ, sigs ≈ 2.44 |
| 2.0 | 2.92 | mode at rate ≈ 3.35γ, sigs ≈ 2.47 |
| 3.0 | 3.00 | mode at rate ≈ 3.35γ, sigs ≈ 2.22 |
| 5.0 | 2.76 | mode at rate ≈ 4.65γ, sigs ≈ 1.94 (rate-crossing through dominant) |

The Σ has no monotone behavior; it rises from 1.68 (q=0.5) to ≈3.0 (q=2.5−3) then declines at q→∞. Mode-by-mode tracking is fragile due to rate crossings; pair-summing by Hamming-complement (Γ_a + Γ_b = 8γ at fixed |ω|/J) shows the dominant pair_1 sum monotonically rising q=0.75 → q=2 then declining. **Diabolic-degeneracy locus at q ≈ 0.659**: pair_1 mode (sigs=1.34 at the nearest sampled q=0.5), the two octic eigenvalues coalesce there (a genuine double root), consistent with the (3q⁴+q²−1)² discriminant-factor zero. The crossing is diabolic (semisimple): the eigenvectors stay independent (g1=2, dep=0), NOT a defective EP (see the location subsection below). This connects path-3's (SE, DE) octic structure to the F86 rate-channel phenomenology.

**Status**: Tier 2 empirical for the amplitudes (no polynomial-in-q fit), now *explained* by the Tier-1 result that the octic Galois group is non-solvable (**Gal(F_8) = S_8**, see above): the eigenvalues themselves have no radical closure in q, so neither do amplitudes built from them. The closed-form analytical layer ends at the AT-locked F_a/F_b quadratics (4 of 12 S_2-sym eigenvalues) and the F_a amplitudes (F_b is invisible to S(t)). Path-3 is "half-solved": exactly the AT-protected half admits radical closure, and the other half provably cannot.

## Path-3 octic diabolic-degeneracy location (Tier 1 derived)

The (3q⁴+q²−1)² perfect-square factor of disc(F_8) locates a **diabolic degeneracy at q = √((−1+√13)/6) ≈ 0.658983** (verified bit-exact: 3q⁴+q²−1 = O(10⁻¹⁶) at this q in [`f89_path3_ep_locator.py`](../simulations/f89_path3_ep_locator.py)). The factor is a *perfect square*: the EP-condition (3q⁴+q²−1) appears to even multiplicity 2, a *double* zero of the discriminant in q, so the two eigenvalues cross linearly/analytically (a defective √-branch EP would force a *simple* zero, so the perfect square is consistent with an analytic crossing but is corroborating, not decisive, on its own); the crossing is semisimple (diabolic), established decisively by the scalar-λI restriction of the octic onto the coalescing span and confirmed artifact-free (g1=2, dep=0, ‖P‖≈3.88 finite; `f89_jordan_definitive.py`, `inspect --root epcharacter`). Numerical sweep around q_EP identifies the merging pair: two octic eigenvalues with rates approaching 4γ and 4γ (from above and below the spectral midpoint of rate 2γ and rate 6γ) and frequencies converging to 2J. Together:

    λ_EP ≈ −4γ + 2iJ

The merged-eigenvalue rate Re(λ_EP) = −4γ sits exactly at the AT-spectral midpoint of the (SE, DE) sector (between rate 2γ for overlap and rate 6γ for no-overlap). This is the structural signature of a 2-level Liouvillian coalescence at the midpoint of two AT-quantized rate channels spanning a 4γ gap. The number g_eff = 2/q_EP ≈ 3.034 is the EP-location relation Q_EP = 2/g_eff of the SEPARATE F86a 2-level rate-channel reduction (it fixes the eigenVALUE location); it is NOT a genuine coupling within the octic: the octic's own 2×2 restriction onto the coalescing span is scalar λ·I, so the eigenvectors stay independent (diabolic), not hybridized.

**Preliminary observation (not Tier-locked)**: the numerical value Re(λ_EP) = −4γ coincides with F86's reported t_peak = 1/(4γ₀) for the (n, n+1)-block 2-level reduction. Whether this reflects a genuine shared 2-level rate-midpoint/coalescence structure across F89 and F86, or is an algebraic coincidence of the AT-rate-midpoint, requires F86 itself to be on a settled foundation before drawing inheritance conclusions. F86 is a collection of partial results that has not been closed; cross-framework bridges should be revisited when F86 is restarted with a clean slate. The F89-side observation stands on its own: **path-3's (SE, DE) octic has a diabolic degeneracy at q ≈ 0.659 with merged eigenvalue −4γ + 2iJ**, derived purely from F89 internal structure. *Why* it is diabolic and not the generic defective (free-fermion integrability, the twin scalar restriction, and the decisive XXZ Δ ≠ 0 gate that defects it) is the hypothesis [`DIABOLIC_BY_INTEGRABILITY.md`](../hypotheses/DIABOLIC_BY_INTEGRABILITY.md); the same diabolic-on-the-line is the spectral face of the [branch-locus palindrome](F89_BRANCH_LOCUS_PALINDROME.md).

**Status**: Tier 1 derived for the q_EP location (analytical from disc factorisation) and the merged-eigenvalue Re(λ_EP) = −4γ identification (numerical sweep). The empirical match with F86's t_peak number is logged as a preliminary cross-framework observation; promotion to a structural bridge claim is deferred until F86 is on stable footing.

## Computing the path-3..6 Galois groups live: isolate-before-DDF via the rate-confined invariant subspace (Tier 1 derived)

The S_8 certificate above generalises to all four paths, and the whole pipeline is now a **live C# witness** (`inspect --root f89galois`): for each k ∈ {3, 4, 5, 6} it recomputes the H_B-mixed factor F_d from the block at inspect time (F_d is never imported) and reads off Gal(F_d/Q(i)(q)) = S_d (d = 8/18/32/53, non-solvable). Per path: build the (SE, DE) S_2-sym block at q0 = 2 over Z[i] (integer-mirror basis, ×2-cleared to Gaussian integers), take its **division-free Berkowitz** characteristic polynomial C, isolate F_d = C / AT, reduce F_d modulo split primes 𝔭 (p ≡ 1 mod 4), and apply the same generalised-Jordan certificate (a d-cycle ⟹ transitive; a prime cycle in (d/2, d−3] ⟹ ⊇A_d; an odd cycle ⟹ S_d). This runs live even for path-6 (75×75 block, degree-53 DDF); the full render recomputes all four in ~8 s. The Python engine of record is [`f89_pathk_galois.py`](../simulations/f89_pathk_galois.py) (`full-d`); the committed F_d(λ,2) literals are kept only as a test cross-check, and the full-D tests assert that the live-reconstructed F_d equals them.

**Isolate-before-DDF (the degree-pollution principle).** You cannot read F_d's cycle type off the *full* block charpoly C, even after taking its squarefree part. C's roots are the AT eigenvalues **plus** F_d's, so a Frobenius cycle type read from C permutes d + (#AT roots) points: a subdirect product of the wrong degree, which breaks the "transitive on exactly d points" premise that Jordan's theorem requires. F_d must first be isolated as a standalone irreducible of degree exactly d. (This is why DDF-ing the whole charpoly is wrong, independently of the AT degeneracy: the failure is structural degree pollution, not "the discriminant vanishes".) Isolation is validated by the **exact triple** over Q(i): exact division with remainder 0, deg F_d = d, and gcd(AT, F_d) = 1, which together force F_d to be exactly the degree-d non-AT factor.

**Reconstructing the AT factor: the naive Slater rule fails; the real structure is rate-confined invariant subspaces.** Isolating F_d needs AT(λ) = ∏(λ − λ_AT) as an exact Z[i] polynomial. The natural guess, that the AT spectrum is the single-particle F_a Bloch frequencies (rate 2γ) plus the 2-particle DE-Slater sums E_j + E_k (rate 6γ), is correct for path-3/4 but **fails from path-5 on** (the no-overlap F_b sector is not a clean Slater-sum multiset). The robust construction reads the AT factor off the **rate-confined invariant subspaces** of the parent's §"AT-lock mechanism" ([F89_TOPOLOGY_ORBIT_CLOSURE.md](F89_TOPOLOGY_ORBIT_CLOSURE.md)) directly: writing the cleared block as M = D + iK (D the ×2 rates ∈ {−4, −12}, K the real-symmetric hopping), for each rate sector take W = the largest M-invariant subspace inside it and set AT = ∏ over sectors of charpoly(M restricted to W). The largest invariant subspace is found by the iterative shrink W ← {v : Kv ∈ span W}, which applies K once per step and so avoids the coefficient blow-up of forming Kᵐ (essential at path-6, dim 75). Validated exactly (sympy and the committed C# oracle) for path-4/5/6: AT·F_d reproduces the full charpoly, with AT degrees 8/13/22 (split by rate as 2+6, 3+10, 3+19 over the overlap/no-overlap sectors) and F_d degrees 18/32/53.

**Tier honesty.** The per-class degree table {8, 18, 32, 53} is Tier-1-*candidate* (switch-enumerated by N_block, no closed form); the Galois verdict S_d on each is Tier-1-*derived* via the q0 = 2 specialization certificate (specialization can only shrink the group). path-3 additionally carries the all-q ⊄A_8 discriminant above; path-4/5/6 rest on the q0 certificate alone, but equally firmly (the certificate gives the full S_d, not just non-solvability). Typed as `F89PathKHbMixedDegreesClaim`, breadcrumbed to the live witness.

## Open / deferred

- **F89 → F86 bridge: deferred.** Path-3's octic diabolic degeneracy at q ≈ 0.659 has merged eigenvalue Re(λ_EP) = −4γ which numerically matches F86's reported t_peak = 1/(4γ₀); whether this is a structural inheritance or an algebraic AT-rate-midpoint coincidence requires F86 itself to be on settled foundations first. (Note: F86a's own real-axis EP was retracted 2026-06-21, the `LocalGlobalEpLink` retraction note in [`ANALYTICAL_FORMULAS.md`](../docs/ANALYTICAL_FORMULAS.md), so this is a coincidence between two now-clarified objects.) Cross-framework bridge claims to be revisited after F86 is restarted with a clean slate.

## Reproduce

```bash
dotnet run --project compute/RCPsiSquared.Cli -- inspect --root f89galois     # live: isolate F_d, certify Gal(F_d)=S_d for path-3..6 (~8 s)
dotnet run --project compute/RCPsiSquared.Cli -- inspect --root epcharacter   # the diabolic-vs-defective compass at q_EP
dotnet run --project compute/RCPsiSquared.Cli -- inspect --root f89octic      # the ported octic block + the live DIABOLIC verdict
python simulations/f89_path3_octic_galois.py                                  # the octic q0=2 Frobenius certificate, gate-first
python simulations/f89_pathk_galois.py                                        # the path-k engine of record (full-d)
python simulations/f89_path3_ep_locator.py                                    # the diabolic locator
python simulations/f89_jordan_definitive.py                                   # g1=2, dep=0 at q_EP (artifact-free)
python simulations/f89_path3_octic_amplitude_q_scan.py                        # the amplitude q-scan
```

---

*The chain's integrability spends itself entirely on the factorisation: the AT-protected half closes in radicals, and the residual factors are as unwritable as algebra allows, S_d at every path length we computed (path-3..6). What survives of structure past that wall is not a formula but a location: one diabolic double zero on the solvable quartic, sitting at the exact midpoint of the two absorption rates.*
