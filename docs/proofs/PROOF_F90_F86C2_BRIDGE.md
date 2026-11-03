# PROOF F90: F86 c=2 ↔ F89 Bridge Identity

**Status:** Tier 1 derived (algebraic identity by construction, operator-exact at N=5..8; the per-bond readouts agree with F86's own scan to 10⁻¹⁰ at all 22 bonds of N=5..8)
**Date:** 2026-05-11
**Authors:** Thomas Wicht, Claude (Anthropic)
**Probe:** [`simulations/f89_to_f86_kbond_via_eigendecomp.py`](../../simulations/f89_to_f86_kbond_via_eigendecomp.py) (the F89 block builders); envelope readout [`simulations/f86_envelope_second_route.py`](../../simulations/f86_envelope_second_route.py)
**Typed claim:** [`compute/RCPsiSquared.Core/Symmetry/F90F86C2BridgeIdentity.cs`](../../compute/RCPsiSquared.Core/Symmetry/F90F86C2BridgeIdentity.cs)

---

## Abstract

F86 lives on the resonance side of the project: Q_peak, EP location, per-bond response curves. F89 lives on the path-polynomial side: closed-form D_k polynomials computed from Chebyshev orbits over the chain dispersion. The two families look like they describe different things; one is dynamic (time-dependent K observables), the other is algebraic (polynomial closed forms in N). This proof says they are, in a precise sense, the same thing.

The bridge is an exact identity. The F86 c = 2 per-bond response observable K_b(Q, t) on the (1, 2) coherence block of an N-qubit chain equals the per-bond Hellmann-Feynman derivative of the F89 path-(N − 1) sub-block dynamics at bond b, with everything else (probe state, kernel, dephasing rates) algebraically matching, and a single Hamiltonian-convention factor of 2 relating the Q-axis scales. The two computations land on the same operator. They are different ways of writing down the same dynamics.

The identity is Tier 1 derived by construction: both sides reduce algebraically to the same SE-DE sub-block Liouvillian, and the Hellmann-Feynman derivative on the F89 side matches the per-bond observable construction on the F86 side. The operators agree exactly at N = 5, 6, 7, 8 (integer-coefficient matrices), and the per-bond readouts agree with F86's own scan to 10⁻¹⁰ at all 22 bonds; the numerical agreement was the empirical entry point before the algebraic identity was written down.

The diagnostic upshot is that closed-form work on either side transfers to the other. F89 closed D_k in closed form via the path-polynomial pipeline; that gives the F86 c = 2 K_b observable a closed-form route through the F90 bridge. The closed-form attempt for g_eff on the F86b obstruction side could then be checked against the bridge identity (L6 of the obstruction proof); F89 closure was confirmed in 2026-05-15, but the predicted transfer to g_eff was tested and refuted, surfacing that g_eff is genuinely a fit-residue rather than an F89-derivable quantity. F90 is the explicit identity that makes those cross-family experiments well-defined.

---

## Statement

For all N ≥ 3 and bond b ∈ {0, ..., N−2}, the F86 c=2 K_b(Q, t) observable on the (n=1, n+1=2) coherence block of an N-qubit XY chain with Z-dephasing equals the per-bond Hellmann-Feynman derivative of F89 path-(N−1) (SE, DE) sub-block dynamics applied at bond b, modulo the Hamiltonian convention factor:

    K_b^{F86 c=2}(Q_F86, t) = K_b^{F89 path-(N−1) (SE,DE)}(Q_F89 = Q_F86 / 2, t)

with all other ingredients (probe, S_kernel, dephasing rates, Liouvillian construction) algebraically identical.

---

## Convention difference (the only structural offset)

| Framework | Hamiltonian | Single-particle SE Bloch | Q_peak typical |
|---|---|---|---|
| **F86** (per `BlockLDecomposition.cs:13-14`, C# typed) | `H_b = (J/2)·(X_b X_{b+1} + Y_b Y_{b+1})` | `ε_k = 2J·cos(πk/(N+1))` (hopping J) | Endpoint ≈ 2.5 |
| **F89** (per `simulations/f89_pathk_lib.py` `build_block_H` (Python) + [`compute/RCPsiSquared.Core/F89PathK/F89BlockHamiltonian.cs`](../../compute/RCPsiSquared.Core/F89PathK/F89BlockHamiltonian.cs) (C#)) | `H = J·(XX + YY)` | `ε_k = 4J·cos(πk/(N+1))` (hopping 2J) | Endpoint ≈ 1.27 |

The conversion factor itself is typed in C# as `F90F86C2BridgeIdentity.JConventionFactor = 2.0` with helpers `F86JToF89J(double)` and `F89JToF86J(double)`. The F89-convention Hamiltonian/Liouvillian *builder* now has a C# counterpart in [`compute/RCPsiSquared.Core/F89PathK/`](../../compute/RCPsiSquared.Core/F89PathK/) (commit `b9e5fe9`, simplified in `298f51c`): `F89BlockHamiltonian.BuildBlockH(nBlock, J)` returns `2·PauliHamiltonian.XYChain(nBlock, J)` (the factor-of-2 reflecting `JConventionFactor`); `F89BlockLiouvillian.BuildBlockL(nBlock, J, γ)` returns the column-major super-operator. F86's `BlockLDecomposition` and F89PathK's `BuildBlockL` produce convention-equivalent operators under the bridge identity.

F89's effective hopping amplitude (per unit J) is 2× F86's (F86 carries the inline 1/2). Hence F89 needs *half* the J for the same physical hopping: J_F89 = J_F86 / 2 (equivalently F86's J = 2·F89's J), and Q_F89 = Q_F86 / 2. This is a one-time relabeling, not a deeper structural difference.

---

## Why the bridge holds (algebraic argument)

F86's `BlockLDecomposition.Build` constructs a Liouvillian on the (n=1, n+1=2) coherence-block flat basis of dimension N·C(N,2):

1. **Diagonal D[i, i] = −2γ₀·HD(p, q)** where (p, q) is the basis pair at flat index i. For c=2: HD(p, q) ∈ {1 (overlap), 3 (no-overlap)}, giving rates 2γ (overlap) and 6γ (no-overlap). This is **identical** to F89's (SE, DE) dephasing diagonal.

2. **Per-bond M_h_per_bond[b]** entries: for state p with adjacent pair flip at bond b that swaps two opposite bits, the bond contributes ±i to off-diagonal entries linking p ↔ p_flipped at the same q (and q ↔ q_flipped at the same p). This is **identical** to F89's per-bond hopping action on (SE, DE) basis pairs (i, jk): the SE-side hop modifies the SE state via single-particle adjacent-bond swap; the DE-side hop modifies the DE pair via the same single-particle adjacent-bond swap on either of the two DE sites.

3. **Total uniform-J L = D + Σ_b J_b·M_h_per_bond[b]** is the same operator that F89 path-(N−1) builds for the (SE, DE) sub-block, modulo the J/(J/2) convention factor.

4. **K_b(Q, t) = 2·Re ⟨ρ(t)| S_kernel | ∂ρ/∂J_b⟩** where ∂ρ/∂J_b = Duhamel-integral with kernel `M_h_per_bond[b]`. This is the **per-bond Hellmann-Feynman** derivative of F89's eigendecomposition: applying perturbation theory (Hellmann-Feynman + first-order correction) to the F89 spectrum gives precisely F86's K_b expression.

5. **Probe ρ_0 (`DickeBlockProbe`)** = `1/(2·√(N·C(N,2)))` per basis pair, identical to F89's `compute_rho_block_0` (SE, DE) Term-1 which has uniform `pre/2 = 1/(2·√(N²(N-1)/2)) = 1/√(2N²(N-1))` per basis pair (algebraically identical: `1/(2·√(N·C(N,2))) = 1/(2·√(N·N(N-1)/2)) = 1/√(2N²(N-1)) = pre/2`).

6. **S_kernel (`SpatialSumKernel`)** = `Σ_site 2·|w_site⟩⟨w_site|` where w_site picks basis pairs (p, q) differing at exactly one site (the site index) with p_site = 0, q_site = 1. For (SE, DE) basis pairs (i, jk), this means w_site picks pairs where site ∈ jk AND i = the OTHER element of {j, k}: exactly F89's `per_site_reduction_within_block_se_de` matrix (per `simulations/f89_path3_at_locked_amplitude_symbolic.py:per_site_reduction_within_block_se_de`; C# port available in [`F89PathK/F89BlockSiteReduction.cs`](../../compute/RCPsiSquared.Core/F89PathK/F89BlockSiteReduction.cs)).

All ingredients identical. The only difference is the J convention factor, hence Q_F89 = Q_F86 / 2.

This is not merely an argument: the two constructions are **bit-identical as operators**. On the shared (j, pq) = (SE, DE) basis, `‖L_F86(J) − L_F89(J/2)‖_F = 0` exactly (integer-coefficient matrices, not a small residual) at N = 5, 6, 7, 8. The bridge is an exact operator identity, not a numerical coincidence; everything below is a readout of this one shared Liouvillian.

---

## Numerical check (all 22 bonds of N=5..8 agree to 10⁻¹⁰)

[`simulations/f86_envelope_second_route.py`](../../simulations/f86_envelope_second_route.py) reads c = 2 through F89's full (SE, DE) eigendecomposition (the block builders of [`simulations/f89_to_f86_kbond_via_eigendecomp.py`](../../simulations/f89_to_f86_kbond_via_eigendecomp.py), F89 J = F86 J / 2), reproduces F86's K_b(Q, t) per bond via the per-bond Hellmann-Feynman + Duhamel integral, and reads the envelope max_t |K_b(Q, t)| the way F86's `ResonanceScan` defines Q_peak. On the same Q grid (153 points on [0.2, 4], parabolic vertex, linear half-maximum) the two routes agree on Q_peak and on HWHM_left/Q_peak to 10⁻¹⁰ at every bond of N = 5, 6, 7, 8 (to 10⁻¹¹ in Q_peak; the script's `N 2` mode reproduces the per-bond comparison), in both readings of the peak: the canonical one (the first local maximum from below) and the global one on a named window. The committed gates carry part of this: the class values at N = 5..8 (`ValidateAgainstPythonEnvelopeRouteTests` against the fixture the second route writes) and the per-orbit values of `PerF71OrbitKTableTests` (the Endpoint orbit at N = 5..10, all orbits at N = 6). None agrees bit for bit: at the canonical peaks the peak value K agrees between the routes to ~10⁻¹³ (to ~10⁻¹² at the far maxima below), and the parabolic vertex through three grid points turns that into a few 10⁻¹¹ in Q_peak.

The canonical readings, one value per mirror pair (F86c pairs bond b with N−2−b exactly), Q_peak in F86-J:

| N | b = 0 (Endpoint) | b = 1 | b = 2 | b = 3 (central) |
|---|---|---|---|---|
| 5 | 2.5171 / 0.7766 | 1.4831 / 0.7523 | | |
| 6 | 2.5806 / 0.7815 | 1.6367 / 0.7596 | 1.4283 / 0.7505 (central) | |
| 7 | 2.5601 / 0.7813 | 1.6989 / 0.7632 | 1.5011 / 0.7546 | |
| 8 | 2.5431 / 0.7808 | 1.7272 / 0.7644 | 1.5207 / 0.7563 | 1.5698 / 0.7585 |

(Each cell: Q_peak / HWHM_left/Q_peak.) Every bond has its canonical peak inside [0.2, 4]. Past a minimum after it the envelope rises again, and at N ≥ 7 some bonds reach a higher maximum on a wider window; the global reading names its window:

| N | bonds | Q_peak | HWHM_left/Q_peak | window |
|---|---|---|---|---|
| 7 | 1, 4 | 7.2554 | 0.9166 | [0.2, 20] |
| 8 | 1, 5 | 8.2247 | 0.8920 | [0.2, 20] |
| 8 | 2, 4 | 38.208 | 0.9888 | [0.2, 60] |
| 8 | 3 | 20.092 | 0.6428 | [0.2, 60] |

The first two rows agree between the routes to 10⁻¹⁰ (for N = 8 the second route's value comes from its [0.2, 60] run, whose maximum for these bonds is the same point); the two values beyond Q = 20 are read by the second route alone. The global reading is a property of the window, since the curves keep further local maxima.

**The readings are readouts of the one shared operator.** Since L_F86(J) = L_F89(J/2) bit-identically (above), any disagreement between the routes can only be in how the curve is sampled and read, not in the physics. The grid reading is itself not the continuous envelope peak: against a continuous search in Q (`--continuous N` in the same script) it sits high by 7·10⁻⁵ to 2.2·10⁻⁴ in Q_peak and off by up to 3.5·10⁻⁵ in the ratio at c = 2, N = 5..8.

---

## Implications

### For F86 open work (per [`PROOF_F86_QPEAK.md`](PROOF_F86_QPEAK.md))

1. **Item 1' (HWHM_left/Q_peak per bond class closed form):**

   Open. `F86HwhmClosedFormClaim`
   ([compute/RCPsiSquared.Core/F86/Item1Derivation/F86HwhmClosedFormClaim.cs](../../compute/RCPsiSquared.Core/F86/Item1Derivation/F86HwhmClosedFormClaim.cs))
   describes the per-bond ratio as the bare floor plus a lift, HWHM_ratio = 0.7264892771 + α_subclass · g_eff + β_subclass (g_eff the fit's own coordinate 4.39382/(Q_peak + 2)),
   with sub-classes per the `BondSubClass` enum. The floor IS analytically derived: it is the
   envelope HWHM_left/x_peak of the bare doubled-PTF 2×2 model
   ([`C2BareDoubledPtfClosedForm`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BareDoubledPtfClosedForm.cs),
   whose fixed-τ = 1 slice gives the separate value 0.671535). The (α, β) are a least-squares
   line through the canonical anchors of N = 5..8, read live by
   [`F86HwhmAlphaExtraction`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/F86HwhmAlphaExtraction.cs),
   and the anchors do not support a line: Endpoint and Flanking, the two sub-classes with four
   chain lengths, miss it by up to 1.28·10⁻³ and 4.2·10⁻⁴, 12 to 36 times the anchors' offset
   from the continuous envelope peak; Mid and CentralSelfPaired have two points each, so their
   line is exact and tests nothing. The lift above the floor is observed (+0.026 to +0.034
   Interior, +0.050 to +0.055 Endpoint at N = 5..8), and no exact law for it is held.
   Plan: [`docs/superpowers/plans/2026-05-13-f86-hwhm-closed-form-attack.md`](../superpowers/plans/2026-05-13-f86-hwhm-closed-form-attack.md).

   Open analytical step (PROOF_F90 perspective): F89's AT-locked F_a/F_b modes have closed-form eigenvectors (overlap-only / no-overlap-only support, per [`F89PathKAtLockMechanismClaim`](../../compute/RCPsiSquared.Core/Symmetry/F89PathKAtLockMechanismClaim.cs)). The bare doubled-PTF floor (envelope 0.7264892771) is the 2×2 model of the F_a/F_b pair. The lift to the class values (Interior 0.7523..0.7602, Endpoint 0.7766..0.7815 at N = 5..8) is located in the **intra-channel dispersion** of the rank-1-bridge + intra-dispersion model (see Direction (b'') below), NOT in a residual of the inter-channel coupling (per the two-dial scout below, read on the 21-point t grid). The cyclotomic Φ_{N_block+1} pattern of F89's `F89UnifiedFaClosedFormClaim.PathPolynomial(k)` table still supplies the N-scaling structure of the AT-locked / intra contribution. Deriving the lift from this rank-1-bridge + intra-dispersion structure (treated exactly, the F89 Chebyshev terrain) is the open step.

   Phase D probe (2026-05-16) refuted the multi-mode-per-cluster-pair internal mixing hypothesis:
   the 10×10 JW cluster-pair sub-block ‖xB(Q)‖_F observable shows only a 30% bump on flat
   baseline ~2.0, no Lorentzian shape and no left half-max crossing. The lift mechanism is
   therefore NOT internal to a single cluster pair. (Phase D then guessed cross-cluster-pair
   inter-channel structure, e.g. H_B-mixed octic residual interference between pairs; the
   two-dial scout in Direction (b'') below, on the 21-point t grid, puts the lift in the
   intra-channel dispersion, the inter-channel coupling only renormalizing Q_peak / g_eff.)

2. **Item 4' (c≥3 extension):** F89 path-k machinery generalises to chromaticity c ≥ 3 with HD ∈ {1, 3, ..., 2c−1} channels. The (SE, DE) 2-channel structure of c=2 generalises to (n, n+1) (c+1)-channel structure. F89 path-k for path-k ≥ 6 covers extended-N data already (N_block up to 7 verified).

3. **Item 5 (σ_0(c, N→∞) closed-form):** σ_0 is a specific F89 path-k singular value, extractable from F89's typed primitives. The 2√(2(c−1)) trajectory crossing (retracted 2026-05-08) was a finite-N artifact; the true asymptote follows from F89's cyclotomic Φ_{N_block+1} structure.

4. **Direction (b'') (full block-L derivation, NOT 4-mode):** **achieved numerically** via F89's full eigendecomposition (all 22 bonds of N=5..8 agree with F86's scan to 10⁻¹⁰). This was the explicit "active path" per `C2HwhmRatio.PendingDerivationNote`. Closed-form analytical lift via F89 AT-locked structure is the next step.

   **Localized (2026-06-11, the two-dial scout):** [`simulations/f86b2_shape_invariance_dial.py`](../../simulations/f86b2_shape_invariance_dial.py) (read on a 21-point t grid, [0.6, 1.6]·t_peak, not on the envelope; its findings are differences between dial settings on that one reading, and they have not been re-read on the envelope) splits the full block generator into the top inter-channel SVD pair (σ₀·u₀v₀†), the inter-channel tail, and the intra-channel dispersion, and dials each. Finding 1: **the ratio is a shape invariant of the tail**. Scaling the inter-channel tail to zero moves HWHM_left/Q_peak by ≤ 0.003 while Q_peak moves by 0.26-0.34 (b0: 2.84 → 2.50; b1: 1.74 → 1.48): the tail renormalizes g_eff (through Q_peak) and leaves the dimensionless lineshape alone. Finding 2: **the lift IS the intra-channel dispersion, non-perturbatively**: without it the response vanishes identically (peak-K ~ 10⁻¹⁶; the Dicke probe is orthogonal to the EP pair, and the intra modes are the only transport from probe to resonance), and at half strength the ratio overshoots (b0 ≈ 0.80). Consequence: the earlier suspicion "lift = H_B-mixed octic residual *of the inter coupling*" is refuted; the (α, β) derivation target is the **rank-1-bridge + intra-dispersion model** (two HD channels with their internal path-like sine-mode dynamics, one rank-1 coupling), treated exactly in the intra sector, which is precisely the F89 Chebyshev terrain. The response anatomy probe → intra modes → rank-1 bridge → response says where g_eff enters (only through the bridge resonance position) and where the sub-class structure enters (through the per-bond observable's overlap with the intra mode lattice); it does not make the lift linear in g_eff, and on the envelope it is not.

### For the typed-knowledge layer

F90 lives in `RCPsiSquared.Core/Symmetry/` with two F89 parents (`F89TopologyOrbitClosure` + `F89PathKAtLockMechanismClaim`). It does not declare F86 typed primitives as parents because F86 is a collection of partial results not yet closed (per `feedback_f86_is_collection_basin` memory); the bridge stands on the F89 side as a Tier-1 statement about how F89's machinery encompasses F86 c=2's K_b observable.

---

## What this is NOT

- **Not** a stand-alone closed-form derivation of HWHM_left/Q_peak per bond class. Item 1' is open: `F86HwhmClosedFormClaim` describes the lift by a line in g_eff per sub-class, fitted rather than derived and not an exact law; this bridge proof is the F89-side analytical handle on which a full Tier 1 derivation would rest.
- **Not** a claim that F86 entirely reduces to F89. F86's c≥3 strata, F86's hardware-confirmed witnesses, F86's PTF-family connections all stand on their own.
- **Not** a deprecation of F86's `C2HwhmRatio` or any F86 typed primitive. F86's empirical pipeline (full-block ResonanceScan) remains the canonical numerical witness pipeline.
- **Not** a reformulation of F89. F89 path-k stands as derived; F90 just identifies F86 c=2's K_b as a per-bond Hellmann-Feynman application of F89.

The bridge is structural: it explains WHY F86 c=2's per-bond-class HWHM constants exist (the γ₀-invariant cyclotomic structure of F89's eigendecomposition; on the canonical reading the per-bond ratio spreads from 0.7505 to 0.7815 at N = 5..8, two bands rather than one value, cf. PROOF_F86B_UNIVERSAL_SHAPE) and gives a Tier-1 numerical predictor of K_b for any N via F89's eigendecomposition.

---

## Anchors

- Probe: [`simulations/f89_to_f86_kbond_via_eigendecomp.py`](../../simulations/f89_to_f86_kbond_via_eigendecomp.py) (F89 block builders) + [`simulations/f86_envelope_second_route.py`](../../simulations/f86_envelope_second_route.py) (envelope readout)
- Typed claim: [`F90F86C2BridgeIdentity.cs`](../../compute/RCPsiSquared.Core/Symmetry/F90F86C2BridgeIdentity.cs)
- F89 parents: [`F89TopologyOrbitClosure.cs`](../../compute/RCPsiSquared.Core/Symmetry/F89TopologyOrbitClosure.cs), [`F89PathKAtLockMechanismClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/F89PathKAtLockMechanismClaim.cs)
- F86 anchor (bridged target): [`C2HwhmRatio.cs`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs), [`ResonanceScan.cs`](../../compute/RCPsiSquared.Core/Resonance/ResonanceScan.cs), [`PROOF_F86_QPEAK.md`](PROOF_F86_QPEAK.md)
- F89 deep-dive: [`F89_TOPOLOGY_ORBIT_CLOSURE.md`](../../experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md)
