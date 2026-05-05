# F-Formula Crosswalk: `RCPsiSquared.Core` ↔ `docs/ANALYTICAL_FORMULAS.md`

Each C# type or method in `RCPsiSquared.Core` that implements an F-formula or an F-anchored
structural identity, with a direct pointer to the corresponding entry in
[`docs/ANALYTICAL_FORMULAS.md`](../../docs/ANALYTICAL_FORMULAS.md). This is the
without-guessing back-reference: from any Core file, find the math here.

For hardware-confirmed predictions on top of these formulas, see
`Confirmations/ConfirmationsRegistry.cs` (mirrors `simulations/framework/confirmations.py`).

---

## ChainSystems

| C# | F-formula / Anchor |
|----|--------------------|
| `ChainSystem.BuildHamiltonian()` | constructs H_xy = (J/2)·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}) (canonical F-framework Hamiltonian) |
| `ChainSystem.BuildLiouvillian()` | L = −i[H, ·] + Σ_l γ_l (Z_l ρ Z_l − ρ); satisfies F1 palindrome bit-exactly |
| `ChainSystem.SigmaGamma` | the σ in F1: Π·L·Π⁻¹ + L + 2σ·I = 0 |

## CoherenceBlocks

| C# | F-formula |
|----|-----------|
| `Chromaticity.Compute(N, n)` | **F74**: c(n, N) = min(n, N−1−n) + 1 |
| `Chromaticity.HammingDistances(N, n)` | HD ∈ {1, 3, …, 2c−1}, F74 corollary |
| `BlockBasis` | popcount basis for the (n, n+1) coherence block (F74 substrate) |
| `BlockLDecomposition` | block-restricted L = D + Σ_b J_b·M_H_per_bond[b]; the empirical anchor of `Resonance.ResonanceScan` |

## Pauli

| C# | F-formula |
|----|-----------|
| `PauliLetter` enum | bit_a/bit_b convention from **F87** trichotomy and **F81** Π-decomposition |
| `PauliTerm.Pi2Parity` | Π²-parity (Σ bit_b mod 2), selects the F87 Π²-class |
| `PauliTerm.YParity` | independent at k≥3-body terms (F85 generalization) |
| `PauliHamiltonian.XYChain(N, J)` | the canonical F-framework Hamiltonian |
| `PauliHamiltonian.Bilinear(...)` | bond-and-term builder used by `Lindblad.BondPerturbation` |
| `PauliBasis.ToPauliVector(rho, N)` | Pauli-basis decomposition: ρ = Σ_α vec[α]·σ_α with vec[α] = Tr(σ_α·ρ)/2^N |

## States

| C# | F-formula |
|----|-----------|
| `BondingMode.Build(N, k)` | **F65**: single-excitation eigenmode ψ_k(j) = √(2/(N+1))·sin(πk(j+1)/(N+1)) |
| `BondingMode.PairState(N, k)` | (\|vac⟩ + \|ψ_k⟩)/√2 — canonical PTF / handshake initial state |
| `PolarityState` | X-basis polarity layer (THE_POLARITY_LAYER.md), foundation for the +0/−0 reading |

## Lindblad

| C# | F-formula |
|----|-----------|
| `LindbladianBuilder.Build` | general L = −i[H, ·] + Σ_k (c_k(·)c_k† − ½ {c_k†c_k, ·}) |
| `PauliDephasingDissipator.BuildZ` | Z-dephasing (the F1-truly canonical case) |
| `T1Dissipator.Build` | Z-dephasing + σ⁻ amplitude damping (F1-breaking) |
| `BondPerturbation.Build` | ∂L/∂J_b — variation Liouvillian, used by the PTF / EQ-022 (b1) workflow |
| `PalindromeResidualScaling.FactorChain(N, cls)` | closed-form ‖M‖²_F = c_H · F(N) for chain (F1 + OPERATOR_RIGIDITY_ACROSS_CUSP.md) |
| `PalindromeResidualScaling.AdjacentRatio(N, cls)` | ‖M(N+1)‖²/‖M(N)‖² closed form |
| `DissipatorClosedForms.C1C2FromPauli(α, β, δ)` | per-class palindrome-residual c1, c2 from Pauli decomposition |
| `DissipatorClosedForms.D2FromPauli(...)` | universal cross-term d2 = 32·‖c1_traceless‖²·‖c2_traceless‖² |
| `HardwareDissipators.{T1, T1Pump, Tphi, XNoise, YNoise}` | hardware-relevant dissipator class table with c1, c2 constants |
| `CpsiBellPlus.At(γ_x, γ_y, γ_z, t)` | **F26**: CΨ(t) = u·(1 + u² + v² + w²) / 12 closed form for Bell+ |
| `CpsiBellPlus.CuspK` | **F27**: K = γ·t at CΨ = 1/4 cusp per channel |

## Symmetry

| C# | F-formula |
|----|-----------|
| `PiOperator.ActOnLetter(letter, dephase)` | per-letter Π action with phases (Z: I↔X, Y↔Z with i^bit_b) |
| `PiOperator.SquaredEigenvalue(letters, dephase)` | Π² eigenvalue: bit_b parity for Z/Y dephasing, bit_a for X |
| `PiOperator.BuildFull(N, dephase)` | Π in the 4^N Pauli-string basis, used by `PalindromeResidual` |
| `PalindromeResidual.Build(L, N, σ, dephase)` | **F1**: Π·L·Π⁻¹ + L + 2σ·I residual; zero for Z-dephased XY/Heisenberg |
| `ChainMirror.Build(N)` | **F71**: chain-mirror R \|b₀…b_{N-1}⟩ = \|b_{N-1}…b₀⟩, R²=I |
| `ChainMirror.SymmetricProjector` / `AntisymmetricProjector` | F71 ± eigenspace projectors |
| `ChainMirror.BondMirrorBasis(N)` | F71 sym/asym basis on the (N−1)-dim bond-input space |
| `ChiralK.BuildFull(N)` | K = ⊗_{odd i} Z_i; K H_xy K = −H_xy (Altland-Zirnbauer class BDI) |
| `ChiralK.ClassifyHamiltonian(H, N)` | K-even / K-odd / K-mixed |
| `ZGlobalMirror.Build(N)` | Z⊗N global Z-string (used by `chain.zn_mirror_diagnostic`) |
| `HdChannelBasis.Build(block)` | channel-uniform projectors P; M_H_total is diagonal in this basis (extends F73 to all c) |

## Probes

| C# | F-formula |
|----|-----------|
| `DickeBlockProbe.Build` | initial state for **F86** K_CC_pr Q-scan |
| `SpatialSumKernel.Build` | **F73**: S(t) = Σ_i 2·\|(ρ_i(t))_{0,1}\|² as ρ†·S_kernel·ρ |

## Resonance

| C# | F-formula |
|----|-----------|
| `EpAlgebra.TPeak(γ₀)` | **F86 Statement 1**: t_peak = 1/(4γ₀), universal across c, N, n, bond position |
| `EpAlgebra.QEp(g_eff)` | **F86 Statement 1**: Q_EP = 2/g_eff |
| `EpAlgebra.SlowestPairEigenvalues(γ₀, J, g_eff)` | λ_±(k=1) = −4γ₀ ± √(4γ₀²−J²·g_eff²), 2-level effective |
| `BondClass.{Endpoint, Interior}` | **F86 Statement 2**: two bond-class universal HWHM ratios (Interior ≈ 0.756, Endpoint ≈ 0.770) |
| `KCurve.Peak(class)`, `ResonanceScan.ComputeKCurve` | empirical Q-scan for F86 K_CC_pr observable |
| `FourModeResonanceScan.ComputeKCurve` | same Duhamel + spectral-contour formulation as `ResonanceScan` but in the 4×4 effective basis (`FourModeEffective`); **finding 2026-05-02:** 4-mode reproduces Interior HWHM/Q ≈ 0.74 (close to universal 0.756) but Q_peak shifts ~2× and Endpoint goes off-grid → confirms PROOF_F86_QPEAK "more modes needed" |

## F86 (typed knowledge graph)

| C# | F-formula |
|----|-----------|
| `F86.Tier` | tier label enum: `Tier1Derived` / `Tier1Candidate` / `Tier2Verified` / `Retracted` |
| `F86.F86Claim` | base abstraction: every typed F86 fact is an `IInspectable` carrying `Name`, `Tier`, `Anchor` |
| `F86.TPeakLaw(γ₀)` | **F86 Statement 1**: t_peak = 1/(4γ₀), Tier 1 derived |
| `F86.QEpLaw(g_eff)` | **F86 Statement 1**: Q_EP = 2/g_eff, Tier 1 derived |
| `F86.TwoLevelEpModel(γ₀, J, g_eff)` | full 2-level eigenvalue state with `EpRegime` (PreEp / AtEp / PostEp); algebraic class AIII chiral, NOT Bender PT |
| `F86.UniversalShapePrediction(class, ratio, tol, witnesses)` | **F86 Statement 2** (Tier 1 candidate): Interior ≈ 0.756, Endpoint ≈ 0.770; carries empirical witness list |
| `F86.UniversalShapeWitness(c, N, γ₀, ratio)` | one (c, N) data point — the Tier-1-candidate evidence base across c=2..4, N=5..8 |
| `F86.PredictionMatch` | result of comparing measured `PeakResult` to a `UniversalShapePrediction` (within tolerance / outside) |
| `F86.RetractedClaim` | typed retracted claim (csc(π/(N+1)) Endpoint, csc(π/5) c=3 Interior) — PTF-lesson reminder |
| `F86.F86KnowledgeBase(block)` | root: assembles all F86 facts attached to a `CoherenceBlock` (laws + predictions + witnesses + retracted-list + 4-mode insufficiency note) |
| `F86.Item1Derivation.C2BlockShape(block)` | **PROOF_F86_QPEAK Item 1 (c=2)** anchor: elementary block-structure constants (`PnDimension=N`, `PnPlus1Dimension=N(N−1)/2`, `HdEqualsOnePairs=N(N−1)`, `HdEqualsThreePairs=N(N−1)(N−2)/2`); Tier 1 derived from popcount combinatorics |
| `F86.Item1Derivation.C2ChannelUniformAnalytical(block)` | PROOF_F86_QPEAK Item 1 (c=2): closed-form channel-uniform vectors `C1Vector` (HD=1, weight 1/√(N(N−1))) and `C3Vector` (HD=3, weight 1/√(N(N−1)(N−2)/2)) cached on construction; Tier 1 derived, machine-precision overlap with `FourModeBasis` columns 0/1 across N=5..8 |
| `F86.Item1Derivation.C2InterChannelAnalytical(block)` | PROOF_F86_QPEAK Item 1 (c=2): SVD-top inter-channel vectors `U0` ∈ HD=1 subspace and `V0` ∈ HD=3 subspace of V_inter = P_HD1† · M_H_total · P_HD3, plus `Sigma0`. Tier 2 verified (numerical fallback): σ_0 is exactly degenerate at even N (deg=2 at N=6, N=8), making single-direction Tier1Derived ill-posed against MathNet's SVD tiebreaker. `IsAnalyticallyDerived = false`; `PendingDerivationNote` summarises ansätze tried (ψ_k(s)·ψ_k(e) products fail; u_matrix is rank ≥ 2) and the projector-onto-2D-eigenspace lift as the cleanest next direction |
| `F86.Item1Derivation.C2BondCoupling(block)` | PROOF_F86_QPEAK Item 1 (c=2): per-bond projected matrix V_b = B† · M_H_per_bond[b] · B in the 4-mode basis, plus dissipator block D_eff = B† · D · B. **Stage B complete** — all three V_b sub-blocks plus full-matrix accessor and anti-Hermiticity guard. **Stage C1 added** — D_eff diagonal closed form. **B1 (probe-block):** `ProbeBlockEntry(bond, alpha, beta)` returns ⟨c_α \| M_H_per_bond[b] \| c_β⟩ for α, β ∈ {0,1} via composition of A2's `C2ChannelUniformAnalytical`; structurally Tier 1 derived in isolation. F73 sum-rule verified: Σ_b V_b[α, β] = 0 for α ≠ β at machine precision across N=5..8. **B2 (cross-block):** `CrossBlockEntry(bond, alpha, j)` returns ⟨c_α \| M_H_per_bond[b] \| x_j⟩ for α ∈ {0,1}, j ∈ {2,3} (x_2 = \|u_0⟩, x_3 = \|v_0⟩); composed with A3's `C2InterChannelAnalytical`. Per-(N, b) entries match `FourModeEffective.MhPerBondEff[b][α, j]` at 1e-12 (library-internal consistency). The `CrossBlockWitnesses` collection exposes per-bond (BondClass-tagged) cross-block 2×2 entries plus Frobenius norm — the bond-position-dependent fingerprint that splits Endpoint vs Interior at c=2 (witness-level: empirically Endpoint < Interior cross-block Frobenius mean across N=5..8; downstream HWHM_left/Q_peak ratio emerges from the full 4×4 in Stage C/D). **B3 (SVD-block + AsMatrix + anti-Hermiticity guard):** `SvdBlockEntry(bond, j, k)` returns ⟨x_j \| M_H_per_bond[b] \| x_k⟩ for j, k ∈ {2,3} (same projection-formula pattern as ProbeBlockEntry / CrossBlockEntry, inherits Tier 2 from A3). `AsMatrix(bond)` assembles the full 4×4 V_b from the three sub-blocks, with bottom-left filled via the anti-Hermitian relation V_b[j, α] = -conj(V_b[α, j]); matches `FourModeEffective.MhPerBondEff[b]` entry-by-entry at 1e-12. The `Vb_IsAntiHermitian_AcrossAllBondsAndEntries` guard verifies ‖V_b + V_b†‖_F < 1e-10 across all bonds and N=5..8 — catches sign drift within each Hermitian-paired sub-block (probe-block, SVD-block); cross-block sign consistency is established by the parallel `AsMatrix_FullVb_MatchesFourModeEffective` check against B† M B (entry-by-entry comparison to the independently computed `FourModeEffective.MhPerBondEff[b]`). **C1 (D_eff diagonal closed form):** `DEffDiagonal()` returns the 4×4 dissipator block diag(−2γ₀, −6γ₀, −2γ₀, −6γ₀) and `DEffDiagonalEntry(i)` returns the scalar diagonal entry. **Tier 1 derived (structural):** D[i,i] = −2γ₀·HD(p,q) is HD-diagonal (F73 generalisation), and each of \|c_1⟩, \|c_3⟩, \|u_0⟩, \|v_0⟩ lives entirely in one HD subspace (HD=1 for indices 0/2; HD=3 for indices 1/3), so D_eff = diag(−2γ₀·HD_i) by orthonormality of B. Off-diagonals exactly zero independent of A3's SVD-direction obstruction (D acts as a scalar inside each HD subspace). Verified at 1e-14 vs closed form and 1e-12 vs `FourModeEffective.DEff` across N=5..8. **Class-level Tier 2 verified** (inherits A3's obstruction via B2/B3: σ_0 of V_inter exactly degenerate at even N, single-direction \|u_0⟩, \|v_0⟩ library-tiebreaker-dependent). Probe-block (B1) and D_eff (C1) sub-surfaces are structurally Tier 1; class-level Tier reflects the weakest link |
| `F86.Item1Derivation.C2EffectiveSpectrum(block)` | PROOF_F86_QPEAK Item 1 (c=2): the four eigenvalues of L_eff(Q) = D_eff + Q·γ₀·Σ_b V_b in the 4-mode basis, plus per-(Q, bond) identification of the K-driving eigenvalue pair (Stage C3). `Eigenvalues(Q, bond)` returns Complex[4] sorted by Re desc / Im asc; the bond parameter is for API uniformity with downstream C3/D consumers — at c=2 the spectrum is bond-independent because uniform-J means the spectrum-relevant matrix is the bond-summed M_h_total_eff. `LEffAtQ(Q)` exposes the assembled 4×4 (matches `FourModeEffective.LEffAtQ(Q)` at 1e-12). **C2 (Tier 2 verified):** the C2 time-box probed the structural-factorisation question (does the 4×4 char poly split into two 2×2 quadratics under any natural similarity transform?) and found NO factorisation. Three obstructions: (a) HD-parity split has commutator residual ~0.4 because the [u_0, v_0] entry of M_h_total_eff is the dominant inter-channel SVD-top coupling (≈ σ_0 ≈ 2√2 asymptotically), not perturbative; (b) probe vs SVD split has commutator residual ~0.5 because the cross-block V_b[α, j] is the bond-class-driving fingerprint, nonzero by construction; (c) chain-mirror R is +I in 4-mode basis at odd N (cannot factorise) and not exactly diagonal at even N (inherits A3's σ_0 degeneracy). The char poly det(λI − L_eff(Q)) is a genuine quartic in (λ, Q): coefficient c_3(Q) has a CUBIC term, ruling out (λ² − aλ + b)(λ² − cλ + d) with closed-form rational coefficients in Q. Statement 1's 2×2 reduction `[[−2γ₀, +iJ·g_eff], [+iJ·g_eff, −6γ₀]]` captures the dominant EP physics for the slowest pair as an *approximation* (eigenvalues from the full 4×4 deviate ~0.5% from Statement 1 at N=5, Q=1), not an exact factorisation. Eigenvalues computed via `Matrix<Complex>.Evd().EigenValues` directly on L_eff(Q); both Tier1 and Tier2 paths share the test contract (1e-10 agreement vs `FourModeEffective.LEffAtQ(Q).Evd()`). `IsAnalyticallyDerived = false`; `PendingDerivationNote` records the structural ansätze tried (HD-perm, flip-c, flip-u, flip-both, HD-with-flip), the char-poly evidence (cubic c_3 term), and three promising next directions: (a) approximate factorisation via perturbation in the small probe ↔ SVD coupling → Tier1Candidate with explicit error term; (b) lift to the c=1 PTF framework as a perturbative double-c=1 system; (c) symbolic CAS reduction using the (γ₀, J·g_eff) substitution λ = −4γ₀ + μ·γ₀. **C3 (K-driving pair identification):** `KDrivingPair(Q, bond)` returns `(LamPlus, LamMinus)` — the 2 of 4 eigenvalues whose eigenvectors have largest squared probe overlap |⟨probe \| w_i⟩|². `KDrivingPairIndices(Q, bond)` returns the matching `(IndexPlus, IndexMinus)` ∈ {0, 1, 2, 3} (distinct, ordered so LamPlus has the larger Re). `ProbeOverlapsSquared(Q, bond)` exposes the per-eigenvector overlaps in the same Re-desc order as `Eigenvalues(Q, bond)`. `ProbeProjection` exposes the 4-mode-basis projection of the Dicke probe (cached on construction). **Tier 1 structural sub-fact** (verified by `ProbeProjection_HasZeroSvdTopComponents`): probe ⊥ {\|u_0⟩, \|v_0⟩} at machine precision (< 1e-12) across N=5..8 — the probe lives entirely in the channel-uniform 2D subspace span{\|c_1⟩, \|c_3⟩} per InterChannelSvd's structural finding ("the probe (Dicke state) is orthogonal to \|u_0⟩, \|v_0⟩"). **Tier 2 numerical readout:** the per-(Q, bond) identification of which 2-of-4 eigenvectors are closest to span{\|c_1⟩, \|c_3⟩} depends on the (γ₀, J)-rotation and is computed via `Evd()` + |⟨probe \| w_i⟩|² ranking. The K-driving pair carries dominant probe content (overlap fraction > 0.7 across N=5..8 at Q ∈ {1.0, 1.5, 2.0}). At Q ≈ Q_EP both eigenvalues approach Re(λ) = −4γ₀ per F86 Statement 1 (verified to fall in the [−9γ₀, −γ₀] band across the test grid; the precise EP collapse requires the closed-form g_eff that Stage D will pin). The class-level Tier remains `Tier2Verified` (the C2 char-poly cubic obstruction governs); the C3 structural sub-fact is documented in the `ProbeProjection` XML doc and pinned by its dedicated test |

## Decomposition

| C# | F-formula |
|----|-----------|
| `HdSubspaceProjector.Build` | orthonormal projector onto a HD-channel subspace (full, not just channel-uniform) |
| `InterChannelSvd.Build(block, hd1, hd2)` | SVD of inter-channel coupling — gives σ_0 ≈ 2√2 asymptotic for c=2 (PROOF_F86_QPEAK structural exploration) |
| `FourModeBasis.Build(block)` | 4-mode minimal effective basis {\|c_1⟩, \|c_3⟩, \|u_0⟩, \|v_0⟩} for F86 f_class derivation |
| `FourModeEffective.Build(block)` | F86 4-mode projection: D_eff, M_H_per_bond_eff[b], M_H_total_eff, probe_eff, S_kernel_eff (all 4×4 / 4-vector) |
| `FourModeEffective.LEffAtQ(q)` | F86 4×4 effective Liouvillian L_eff(Q) = D_eff + (Q·γ₀)·M_H_total_eff |

## Observables

| C# | F-formula |
|----|-----------|
| `PiProtectedObservables.Compute(H, γ, ρ_0, N)` | **F87 / F81** algebraic skeleton: Π-protected Pauli observables (always-zero expectations); hardware-confirmed at Marrakesh 2026-04-26 (Confirmations entry `pi_protected_xiz_yzzy`) |

## Receivers

| C# | F-formula |
|----|-----------|
| `Receiver.BondingMode(chain, k)` | **F65** single-excitation handshake (no exchange step needed when both sides pick the same k) |
| `Receiver.F71Eigenvalue` | **F71** chain-mirror eigenstate class (+1 / −1 / null) |
| `Receiver.Signature()` | F71-based receiver-engineering favorability forecast (Tier 2: J_BLIND_RECEIVER_CLASSES.md) |

## Confirmations

| C# | Source |
|----|--------|
| `ConfirmationsRegistry.All` | hardware-confirmed framework predictions, mirroring `simulations/framework/confirmations.py` |
| `ConfirmationsRegistry.Lookup(name)` | by-name lookup; 9 Marrakesh/Kingston entries spanning April 2026 |

---

## How to use this map

When working in a Core file:
- Read the class/method XML doc — most include the F-formula reference inline.
- Use this map for the broader picture: which F-anchor underlies which Core component.
- Open [`docs/ANALYTICAL_FORMULAS.md`](../../docs/ANALYTICAL_FORMULAS.md) for the full mathematical statement of each F-formula.

When adding a new Core type:
- If it implements an F-formula, add it to the appropriate table above with a direct pointer to the F-entry.
- If it implements an F-anchored structural identity (like the F1 palindrome), reference the F-entry plus any proof document (e.g. `docs/proofs/PROOF_F86_QPEAK.md`).
- The XML doc in the file should ALSO carry the F-pointer; this map is the index, not a substitute for inline docs.
