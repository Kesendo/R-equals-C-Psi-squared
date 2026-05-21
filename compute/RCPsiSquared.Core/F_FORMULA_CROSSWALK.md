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
| `F1.F1T1ResidualClosedForm` | **F1 T1 block**: ‖M(T1)‖²_F = 4^(N−1)·[3·Σγ² + 4·(Σγ)²]; per-site M_l with ‖M_l‖²=7, |tr(M_l)|²=16. Tier 1 derived (PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md) |
| `F1.F1T1ResidualPi2Decomposition` | **F1 T1 block Π²-split**: ‖M_anti‖²=4^(N−1)·Σγ² (F82/F84 amplitude-damping side); ‖M_sym‖²=4^(N−1)·[2·Σγ²+4·(Σγ)²]. Pythagorean orthogonal split of F1T1ResidualClosedForm. Tier 1 derived |
| `F1.F1DepolResidualClosedForm` | **F1 depol block**: ‖M(depol)‖²_F = 4^(N−1)·[(16/9)·Σγ² + 16·(Σγ)²]; per-site M_l = diag(−4/3, −4/3, −8/3, −8/3) with ‖M_l‖²=160/9, |tr(M_l)|²=64. Π²-decomposition trivial (M_anti=0); σ-shift=0; H/Z/topology-independent. Closes earlier F5 scalar diagnostic to full Frobenius. Tier 1 derived (PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md) |
| `F1.F49NonUniformCrossTermClaim` | **F49 non-uniform γ cross-term**: ‖{L_H, L_Dc}‖²_F = 4·Σ_b ‖L_H^bond_b‖²·Σ_{m∉bond_b} γ_m² + Σ_b G(bond_b, H)·(γ_{i_b}−γ_{j_b})²; L_Dc = L_D + σ·I. Per-class G-fraction (G/‖L_H^bond‖²) = 4/3 (Heisenberg) / 4 (Ising) / 0 (XY, soft XY+YX); ZZ-fraction of the bond Hamiltonian is the bond-asymmetry source. Recovers F49's 4γ²·(N−2)·‖L_H‖² at uniform γ. Closes PROOF_F1_NONUNIFORM_GAMMA "Open follow-ups" item. Tier 1 derived (PROOF_F49_NONUNIFORM_GAMMA_EXTENSION.md) |
| `F1.F1GeneralTopologyVerifiedClaim` | **F1 general-topology verification record**: the (B, D2) parameterisation of ‖M(N, G)‖²_F = c_H · F(N, G) extends bit-exactly to disconnected components (B and D2 sum across components), weighted edges (B → Σ_b J²_b), random connected Erdős-Rényi graphs, and the single-body class; verified at N=5, 6 (Python, 18 named + 60 random + disconnected + weighted + single-body), N=7 (C# F1 palindromic-pairing via `LiouvillianBlockSpectrum.ComputeSpectrumPerBlock` on chain/ring/star/K_4 + disjoint-3-chain), and N=8 (opt-in SLOW_N8 dogfood across chain/ring/star/K_4 + disjoint-4-chain Heisenberg+Z-deph with full `F1SpectrumStatistics` metric capture under `simulations/results/f1_n8_n9_metrics/`). N=9 chain wired as `SkippableFact` but blocked at the LP64 MKL P/Invoke marshalling ceiling (max joint-popcount block 15876² ≈ 4 GB > 2 GB native-array limit); `ScaleFrontierBlockedAtN = 9` field carries the diagnostic. Analytic anchor: PROOF_CROSS_TERM_FORMULA Lemma 3 + Corollary (bond-disjointness universal across any graph). Closes the last F1 OpenQuestion ("general topology beyond chain/ring/star/K_N"). Tier 2 verified (PROOF_F1_GENERAL_TOPOLOGY.md) |
| `F1.F1SpectrumStatistics` | **F1 spectrum-statistics utility**: `TopologyMetrics` record + `Compute` + JSON serialisation supporting the SLOW_N8 + SLOW_N9 dogfood capture. Five metric groups per (N, topology, H, γ): (1) wall-time profile (total + ComputeSpectrumPerBlock + effective speedup over dense), (2) palindromic-pairing precision (max / mean / median / 99p / min / outlier count), (3) spectrum-structure invariants (min/max Re/Im, dissipation gap, kernel dim, pure-imag/real counts, distinct-binned count), (4) block-decomposition cost picture (sector count, primary count, max block + sector label, top-3 sizes, total cubic cost), (5) Hamiltonian + dissipator setup (J, γ, σ shift, class name, bond list). Used by `F1GeneralTopologyN8BlockSpectrumTests` (4 N=8 systems) and `F1GeneralTopologyN9BlockSpectrumChainTests` (1 N=9 chain, currently SkippableFact). Outputs to `simulations/results/f1_n8_n9_metrics/<topology>_N<n>.json` |
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
| `F89TopologyOrbitClosure` (Symmetry/) | **F89**: S(t) for ρ_cc + uniform-J multi-bond XY depends only on the S_N-orbit of the bond set; for chain B, orbit = topology class (sorted multiset of connected-path-lengths). Tier 1 derived |
| `F89TopologyOrbitClosure.S0ClosedForm(N)` | **F89**: S(0) = (N−1)/N for ρ_cc, probe-only closed form |
| `F89TopologyOrbitClosure.ChainTopologyClass(N, bonds)` | **F89**: canonical orbit label = sorted multiset of connected-path-lengths in the bond-graph |
| `F89TopologyOrbitClosure.AreInSameChainOrbit(N, bondsA, bondsB)` | **F89**: predicate that S(t; bondsA) = S(t; bondsB) per orbit transitivity |
| `F89AdditiveIdentityClaim` (Symmetry/) | **F89 mixed-topology additive identity**: S_T(t) = Σ_i S_(k_i)(t) − (m−1)·N·S_bare(t; N); reduces 14 per-class closed forms to 6 pure-path-k forms + 1 rule (Lindbladian factorisation; verified 27/27 N=7 CSVs at 5.013e-7 precision floor). Tier 1 derived |
| `F89AdditiveIdentityClaim.BarePerSite(N, γ, t)` | per-bare-site closed form (N−1)/N²·exp(−4γt) |
| `F89AdditiveIdentityClaim.OvercountingCoefficient(m, N)` | the (m−1)·N subtraction coefficient |
| `F89AdditiveIdentityClaim.Combine(kValues, N, γ, t, sPathKFunc)` | apply identity to combine per-pure-path-k contributions into S_T(t) |
| `F89PathKVacSeParsevalClaim` (Symmetry/) | **F89 path-k (vac, SE) self-contribution**: S^(vac,SE)_block(t; k, N) = (k+1)(N−k−1)²/(N²(N−1))·exp(−4γt); pure exp(−4γt), no oscillation, via Parseval orthogonality of H_B^SE Bloch eigenstates (machine-precision verified across 15 (k, N) pairs). Tier 1 derived |
| `F89PathKVacSeParsevalClaim.Coefficient(k, N)` | rational prefactor (k+1)(N−k−1)²/(N²(N−1)) |
| `F89Path2CardanoClaim` (Symmetry/) | **F89 path-2 (SE,DE) S_2-sym closed form**: char(λ) = −(λ+2γ)(λ+6γ)·[cubic]; cubic μ³ + 10μ² + (28+32q²)μ + 24(1+4q²) = 0 in dimensionless μ=λ/γ, q=J/γ; solvable in radicals via Cardano (sympy verified). Tier 1 derived |
| `F89Path2CardanoClaim.LinearFactorRates(γ)` | the two pure-AT eigenvalues (−2γ, −6γ) |
| `F89Path2CardanoClaim.CubicCoefficients(q)` | dimensionless cubic coefficients (a, b, c) at given q = J/γ |
| `F89Path2CardanoClaim.CubicEigenvaluesNumerical(q)` | Cardano-solved cubic eigenvalues at q (1 real + 1 complex pair for typical q) |
| `F89PathKAtLockMechanismClaim` (Symmetry/) | **F89 universal AT-lock mechanism**: F_a (rate 2γ) eigvecs supported entirely on overlap basis pairs; F_b (rate 6γ) supported on no-overlap; F_b sigs ≈ 0 universally to S(t) per-site reduction; F_a count = floor(N_block/2) = number of SE-anti single-particle Bloch modes. Tier 1 derived (path-3..6 verified). |
| `F89PathKAtLockMechanismClaim.FaCount(N_block)` | F_a mode count = floor(N_block/2) |
| `F89PathKAtLockMechanismClaim.SeAntiBlochOrbit(N_block)` | The orbit {2, 4, ..., 2·floor(N_block/2)} of SE-anti Bloch indices |
| `F89PathKAtLockMechanismClaim.BlochEigenvalueY(N_block, n)` | y_n = 4·cos(πn/(N_block+1)) (in J=1 units) |
| `F89UnifiedFaClosedFormClaim` (Symmetry/) | **F89 unified F_a closed form across path-3..6**: sigs[F_a:n](N) = P_path(y_n)/[D_path · N²(N−1)] with y_n = 4cos(πn/(N_block+1)); (P,D) = {(14y+47, 9), (10y+25, 4), (13y²+82y+129, 25), (17y²+72y+80, 18)} for path-{3,4,5,6}; sum F_a rational across all paths via Newton's identities. Tier 1 derived (bit-exact). |
| `F89UnifiedFaClosedFormClaim.PathPolynomial(k)` | (coefs low-to-high, denom) per path k ∈ {3, 4, 5, 6} |
| `F89UnifiedFaClosedFormClaim.Sigma(k, n, N)` | sigs[F_a:n](N) for given path-k, Bloch index n, total qubits N |
| `F89UnifiedFaClosedFormClaim.SigmaSum(k, N)` | Σ_n sigs[F_a:n](N) over the SE-anti Bloch orbit (rational) |
| `F89Path3SeDeFactorisationClaim` (Symmetry/) | **F89 path-3 (SE,DE) S_2-sym factorisation**: char(λ) = F_a · F_b · F_8 (deg 2·2·8); F_a roots λ = −2γ + iJ(−1±√5), F_b roots λ = −6γ + iJ(−1±√5); F_8 octic irreducible over Q[i, √5] (sympy verified). Tier 1 derived. |
| `F89Path3SeDeFactorisationClaim.FaRoots(γ, J)` | the two F_a quadratic roots |
| `F89Path3SeDeFactorisationClaim.FbRoots(γ, J)` | the two F_b quadratic roots |
| `F89Path3OcticEpClaim` (Symmetry/) | **F89 path-3 octic exceptional point**: q² = (−1+√13)/6 (from disc factor (3q⁴+q²−1)²); merged eigenvalue λ_EP = −4γ + 2iJ; Re(λ_EP) sits at AT-spectral midpoint of rates 2γ (overlap) and 6γ (no-overlap). Tier 1 derived (analytical + machine-precision numerical). |
| `F89Path3OcticEpClaim.QEpSquared` | (−1+√13)/6 |
| `F89Path3OcticEpClaim.QEp` | √((−1+√13)/6) ≈ 0.658983 |
| `F89Path3OcticEpClaim.MergedEigenvalue(γ, J)` | λ_EP = −4γ + 2iJ |
| `F89Path3OcticGaloisClaim` (Symmetry/) | **F89 path-3 octic Galois group**: disc(F_8) = const · q²⁴ · (3q⁴+q²−1)² · P_10(q²) is NOT a perfect square in Q[q]; therefore Gal(F_8) ⊄ A_8 (Tier 1). Conjecture (Tier 2 open): Gal is non-solvable (likely S_8). |
| `F89Path3OcticGaloisClaim.GalNotInA8` | true (Tier 1 derived) |
| `F89Path3OcticGaloisClaim.NonSolvableConjecture_IsOpen` | true (Tier 2 open) |
| `F89PathKHbMixedDegreesClaim` (Symmetry/) | **F89 path-k H_B-mixed sub-factor degrees**: {8, 18, 32, 53} for path-{3, 4, 5, 6}; conjecturally Galois-non-solvable for degree ≥ 5 (Tier 2 conjecture in docstring). Tier 1 derived (combinatorial). |
| `F89PathKHbMixedDegreesClaim.HbMixedSubFactorDegree(k)` | S_2-sym dim minus AT-locked count per path |
| `F90F86C2BridgeIdentity` (Symmetry/) | **F90 F86 c=2 ↔ F89 bridge identity**: F86 c=2 K_b(Q,t) IS F89 path-(N−1) (SE,DE) per-bond Hellmann-Feynman, modulo F89-J = 2·F86-J convention. Verified bit-exact at 20/22 bonds across N=5..8 including orbit escapes. Tier 1 derived |
| `F90F86C2BridgeIdentity.JConventionFactor` | 2.0 (F89 J / F86 J) |
| `F90F86C2BridgeIdentity.F86JToF89J(f86J)` | f86J / 2.0 |
| `F90F86C2BridgeIdentity.F89JToF86J(f89J)` | f89J · 2.0 |
| `F71AntiPalindromicGammaSpectralInvariance` (BlockSpectrum/) | **F91 F71-anti-palindromic γ spectral invariance**: F71-refined diagonal-block spectrum of XY+Z-deph L is invariant under γ-distributions satisfying γ_l + γ_{N−1−l} = 2·γ_avg (= 90°-rotation of γ-distribution around its mean = γ-parameter side of Pi2-Z₄'s rotational axis). F71 broken as L-symmetry; invariance lives in diagonal blocks only. Bit-exact verified N=4,5,6. **Tier 1 derived (2026-05-12):** proof shows diagonal-block matrix elements depend only on F71-pair-sum multiset {S_l = γ_l + γ_{N−1−l}}; cross-block entries depend on pair-differences D_l = γ_l − γ_{N−1−l}; 90°-rotation invariance is corollary on the orbit S_l = 2γ_avg ∀l |
| `F71AntiPalindromicGammaSpectralInvariance.AntiPalindromicDeviation(γ)` | sqrt(Σ((γ_l + γ_{N-1-l}) − 2·γ_avg)²); zero iff γ is F71-anti-palindromic |
| `F71AntiPalindromicGammaSpectralInvariance.IsAntiPalindromic(γ, tol)` | predicate: γ within tol of F71-anti-palindromic |
| `F92BondAntiPalindromicJSpectralInvariance` (SymmetryFamily/) | **F92 F71-anti-palindromic J spectral invariance**: F71-refined diagonal-block spectrum of chain XY+Z-deph L invariant under J-distributions satisfying J_b + J_{N-2-b} = 2·J_avg. J-side Pi2-Z₄ twin of F91 (γ-side); same algebraic mechanism (diagonal blocks depend on F71-pair-sums T_b = J_b + J_{N-2-b}; cross-blocks on pair-differences B_b = J_b − J_{N-2-b}). Tier 1 derived (algebraic proof PROOF_F92 + bit-exact N=4,5). |
| `F92BondAntiPalindromicJSpectralInvariance.AntiPalindromicDeviation(bondJ)` | sqrt(Σ_b ((J_b + J_{N-2-b}) − 2·J_avg)²); zero iff J is F71-anti-palindromic |
| `F92BondAntiPalindromicJSpectralInvariance.IsAntiPalindromic(bondJ, tol)` | predicate: J within tol of F71-anti-palindromic |
| `F93DetuningAntiPalindromicSpectralInvariance` (SymmetryFamily/) | **F93 F71-anti-palindromic h spectral invariance**: F71-refined diagonal-block spectrum of chain XY+Z-deph+h_l Z_l L invariant under h-distributions satisfying h_l + h_{N-1-l} = 2·h_avg. Tier 1 derived (algebraic proof PROOF_F93 + bit-exact N=4,5). |
| `F93DetuningAntiPalindromicSpectralInvariance.AntiPalindromicDeviation(hPerSite)` | sqrt(Σ_l ((h_l + h_{N-1-l}) − 2·h_avg)²); zero iff h is F71-anti-palindromic |
| `F93DetuningAntiPalindromicSpectralInvariance.IsAntiPalindromic(hPerSite, tol)` | predicate: h within tol of F71-anti-palindromic |
| `C1QPeakMirrorJParity` (F71/) | **F100 c₁/Q_peak bond-mirror deviation parity**: D(b) = c₁(b) − c₁(N−2−b) is exactly odd in the F71-anti-palindromic J-component J_anti; zero for palindromic J, leading-order linear (graceful breakdown of F71 under non-uniform J). Observable-side twin of F92 (F92 = diagonal-block spectrum depends on J_sym; F100 = c₁/Q_peak deviation lives in J_anti). Tier 1 derived (algebraic R-equivariance + numerical N=3,4,5; PROOF_F100). |
| `C1MirrorGammaParity` (F71/) | **F101 c₁ bond-mirror deviation parity (γ)**: D(b) = c₁(b) − c₁(N−2−b) is exactly odd in the F71-anti-palindromic γ-component γ_anti; zero for palindromic γ, leading-order linear (graceful breakdown of F71 under non-uniform per-site dephasing). Observable-side twin of F91 (F91 = diagonal-block spectrum depends on γ_sym; F101 = c₁ deviation lives in γ_anti). c₁ only (Q_peak needs a scalar γ₀). Tier 1 derived (algebraic R-equivariance + numerical N=3,4,5; PROOF_F101). |
| `C1MirrorGammaParity.PalindromicComponent(siteGamma)` / `.AntiPalindromicComponent(siteGamma)` | the F71-palindromic / anti-palindromic split γ_sym = (γ + F71(γ))/2, γ_anti = (γ − F71(γ))/2 of a per-site γ profile |
| `C1MirrorGammaParity.PalindromicDeviation(siteGamma)` / `.IsPalindromic(siteGamma, tol)` | ‖γ_anti‖ = deviation of γ from F71-palindromic; the c₁ bond-mirror deviation D vanishes iff this is zero |
| `C1QPeakMirrorJParity.PalindromicComponent(bondJ)` / `.AntiPalindromicComponent(bondJ)` | the F71-palindromic / anti-palindromic split J_sym = (J + F71(J))/2, J_anti = (J − F71(J))/2 of a bond profile |
| `C1QPeakMirrorJParity.PalindromicDeviation(bondJ)` / `.IsPalindromic(bondJ, tol)` | ‖J_anti‖ = deviation of J from F71-palindromic; the c₁/Q_peak mirror-deviation D vanishes iff this is zero |
| `AbsorptionTheoremClaim` (Symmetry/) | **PROOF_ABSORPTION_THEOREM**: Re(λ) = −2γ₀·⟨n_XY⟩; absorption quantum 2γ₀ = a_0·γ₀; rate-quantization root of F33/F50/F55/F64-F68/F74/F89 (typed parent edge in Schicht-1 Registrations). Tier 1 derived |
| `AbsorptionTheoremClaim.Rate(nXY, γ₀)` | **AT**: α = 2γ₀·n_XY; per-Pauli-string rate |
| `AbsorptionTheoremClaim.PerCoherenceRateComputationalBasis(nDiff, γ₀)` | **AT**: α = 2γ₀·n_diff for |A⟩⟨B| basis-pair coherence (n_XY = n_diff via {I,Z}/{X,Y} per-site decomposition) |
| `AbsorptionTheoremClaim.HammingComplementPairSum(blockSize, γ₀)` | **AT + F89c**: pair-sum 2γ₀·N under column bit-flip ρ[a, b] → ρ[a, bar(b)] (n_diff ↔ N − n_diff Hamming complement); spectral maximum reached pairwise |
| `AbsorptionTheoremClaim.MaxRate(N, γ₀)` | **AT**: max rate 2γ₀·N (full XOR mode); section 4.1 of proof |
| `AbsorptionTheoremClaim.NXyFromRate(rate, γ₀)` | **AT** read backwards: n_XY = α / (2γ₀); inverse map on the 2γ₀-grid |

## Probes

| C# | F-formula |
|----|-----------|
| `DickeBlockProbe.Build` | initial state for **F86** K_CC_pr Q-scan |
| `SpatialSumKernel.Build` | **F73**: S(t) = Σ_i 2·\|(ρ_i(t))_{0,1}\|² as ρ†·S_kernel·ρ |

## Resonance

| C# | F-formula |
|----|-----------|
| `EpAlgebra.TPeak(γ₀)` | **F86a**: t_peak = 1/(4γ₀), universal across c, N, n, bond position |
| `EpAlgebra.QEp(g_eff)` | **F86a**: Q_EP = 2/g_eff |
| `EpAlgebra.SlowestPairEigenvalues(γ₀, J, g_eff)` | λ_±(k=1) = −4γ₀ ± √(4γ₀²−J²·g_eff²), 2-level effective |
| `BondClass.{Endpoint, Interior}` | **F86b**: two bond-class universal HWHM ratios (Interior ≈ 0.756, Endpoint ≈ 0.770) |
| `KCurve.Peak(class)`, `ResonanceScan.ComputeKCurve` | empirical Q-scan for F86 K_CC_pr observable |
| `FourModeResonanceScan.ComputeKCurve` | same Duhamel + spectral-contour formulation as `ResonanceScan` but in the 4×4 effective basis (`FourModeEffective`); **finding 2026-05-02:** 4-mode reproduces Interior HWHM/Q ≈ 0.74 (close to universal 0.756) but Q_peak shifts ~2× and Endpoint goes off-grid → confirms PROOF_F86_QPEAK "more modes needed" |

## F86 (typed knowledge graph)

| C# | F-formula |
|----|-----------|
| `F86.Tier` | tier label enum: `Tier1Derived` / `Tier1Candidate` / `Tier2Verified` / `Retracted` |
| `F86.F86Claim` | base abstraction: every typed F86 fact is an `IInspectable` carrying `Name`, `Tier`, `Anchor` |
| `F86.TPeakLaw(γ₀)` | **F86a**: t_peak = 1/(4γ₀), Tier 1 derived |
| `F86.QEpLaw(g_eff)` | **F86a**: Q_EP = 2/g_eff, Tier 1 derived |
| `F86.TwoLevelEpModel(γ₀, J, g_eff)` | full 2-level eigenvalue state with `EpRegime` (PreEp / AtEp / PostEp); algebraic class AIII chiral, NOT Bender PT |
| `F86.LEffMirrorAxisClaim` | **PROOF_F86B_OBSTRUCTION (the diagnosis)**: L_eff real part −4γ₀ = −2γ₀·2 is the mirror axis of the (−2γ₀, −6γ₀) channel pair; the EP is the coalescence onto it; Absorption Theorem exact at ⟨n_XY⟩=2; g_eff lives in Im(λ) as the branch's relative clock. Tier 1 derived; parent edge ← AbsorptionTheoremClaim |
| `F86.UniversalShapePrediction(class, ratio, tol, witnesses)` | **F86b** (Tier 1 candidate): Interior ≈ 0.756, Endpoint ≈ 0.770; carries empirical witness list |
| `F86.UniversalShapeWitness(c, N, γ₀, ratio)` | one (c, N) data point — the Tier-1-candidate evidence base across c=2..4, N=5..8 |
| `F86.PredictionMatch` | result of comparing measured `PeakResult` to a `UniversalShapePrediction` (within tolerance / outside) |
| `F86.RetractedClaim` | typed retracted claim (csc(π/(N+1)) Endpoint, csc(π/5) c=3 Interior) — PTF-lesson reminder |
| `F86.F86KnowledgeBase(block)` | root: assembles all F86 facts attached to a `CoherenceBlock` (laws + predictions + witnesses + retracted-list + 4-mode insufficiency note + c=2 derivation when applicable + local–global EP link + polarity inheritance link). Exposes `C2UniversalShape` for c=2 blocks (null otherwise), `LocalGlobalEpLink` (Locus 5 EP-side meta-claim, always non-null), and `PolarityInheritanceLink` (Locus 6 symmetry-side meta-claim, always non-null) |
| `F86.LocalGlobalEpLink` | **F86 ↔ FRAGILE_BRIDGE meta-claim** (Tier 2 verified): the c=2 local EP at real Q_EP and FRAGILE_BRIDGE's complex-γ-plane global EP are the same exceptional-point structure under shared AIII chiral algebra. Empirically validated by the 2026-05-06 c=2 N=5..8 Petermann-K sweep: K grows monotonically with N within each parity class (odd: ~1.79× per step; even: ~2.36× per step), odd-N peaks dominate even-N peaks by 2–4× at every step (matches A3's σ_0 R-even/R-odd-degeneracy prediction), and by N=7 K=2384 sits ~6× above FRAGILE_BRIDGE's K=403 ballpark on the real Q axis. Properties: `SharedAlgebra` (same-sign-imaginary 2×2; AIII chiral), `LocalInstanceAnchor` (F86a, PROOF_F86_QPEAK Statement 1), `GlobalInstanceAnchor` (hypotheses/FRAGILE_BRIDGE), `Witnesses` (4 pinned `PetermannSpikeWitness` rows from the c=2 N=5..8 sweep at γ₀=0.05, 121-pt Q-grid on [0.5, 4.0]), `PendingDerivationNote` (Tier1Derived promotion path: complex-γ infrastructure in `LindbladPropagator` OR closed-form K(N) at the EP). Tier2Verified honesty: shared algebra + real-axis-hit empirically validated; explicit complex-γ analytic continuation NOT carried out in code. **Locus 5 of inheritance-through-layers** (EP-side closure); companion at Locus 6 is `PolarityInheritanceLink` (symmetry-side closure) |
| `F86.PolarityInheritanceLink` | **F86 ↔ polarity-layer-pair meta-claim** (Tier 2 verified, with composition-reading + empirical-stability witness from Direction (α)): the c=2 bond-class split (Q_peak 2 ± 0.5, HWHM/Q* 1/2 + r·(1/2)) inherits structurally from the polarity-layer pair {−0.5, +0.5} at d=2 named in `Pi2KnowledgeBase` (`PolarityLayerOriginClaim` + `QubitDimensionalAnchorClaim` + `HalfAsStructuralFixedPointClaim`). The decomposition Q_peak ≈ 2 + r and HWHM/Q* ≈ 1/2 + r·(1/2) holds across c=2 N=5..8 with mean Q_peak = 2.04 ± 0.06 and Interior r_HWHM ≈ 0.502 (close to `HalfAsStructuralFixedPoint` 1/2 anchor). **Composition-reading from 2026-05-07 Direction (α) attempt:** r_Q(N, b) = `BareDoubledPtfXPeak` · Q_EP(N, b) − 2 = 4.39382 / g_eff(N, b) − 2; algebraic identity composing two existing Tier 1 facts (`C2HwhmRatio.BareDoubledPtfXPeak` = 2.196910 universal + F86a's Q_EP = 2/g_eff). Caveat: mathematically tautological if g_eff is defined via Q_peak inversion; the genuine Tier-1 content is the universality of BareDoubledPtfXPeak, so the bond-class signature must live entirely in g_eff(N, b). Encoded in `ClosedFormCompositionNote`. Direction (α)-test g_eff_E ≈ σ_0(N)·√(3/8) matches Δ ≤ 0.01 for N ≥ 6 but does not pin at tolerance 0.005. Properties: `PolarityRootAnchor` (parent Pi2 claims), `ParallelLocusReference` (Locus 5 EP-side closure pointer), `Witnesses` (4 pinned `PolarityWitness` rows from the c=2 N=5..8 `C2HwhmRatio` pipeline at γ₀=0.05), `ClosedFormCompositionNote` (composition-reading text), `EmpiricalSumQPeakAsymptote = 4.12` (the per-class sum 1/g_eff_E + 1/g_eff_I ≈ 0.937; Direction (γ') tested whether this is a structural constant and refuted it 2026-05-14: a per-class orbit-mixing artefact, not a closure), `IsAnalyticallyDerived` (false; class-Tier stays Tier2Verified because closed-form g_eff(N, b) is open), `PendingDerivationNote` (Tier1Derived promotion path: (α') full block-L per-bond Q_peak derivation joining C2HwhmRatio direction (b''); (β') Locus 5 EP-rotation × Locus 6 polarity inheritance composition; (γ') tested and refuted 2026-05-14). Tier2Verified honesty: closed-form g_eff(N, b) per bond class is open (inherits A3 / `C2InterChannelAnalytical` Tier 2 obstruction). 2026-05-07 Direction (α) finding documented in `docs/superpowers/syntheses/2026-05-07-direction-alpha-attempt.md`. **Locus 6 of inheritance-through-layers** (symmetry-side closure); companion at Locus 5 is `LocalGlobalEpLink` (EP-side closure) |
| `F86.PolarityWitness(N, QPeakInterior, QPeakEndpoint, HwhmRatioInterior, HwhmRatioEndpoint, RQpeakInterior, RQpeakEndpoint, RHwhmInterior, RHwhmEndpoint)` | one row of the pinned c=2 N=5..8 polarity-layer table held by `PolarityInheritanceLink`. R-fields are pre-computed polarity content: `RQpeak* = QPeak* − 2` (Q_peak decomposition around d=2), `RHwhm* = 2·(HwhmRatio* − 1/2)` (HWHM decomposition around 1/2 baseline). `QPeakMean` and `HwhmRatioMean` are derived properties. Frozen empirical anchor (the data is fixed at the 2026-05-07 `C2HwhmRatio` pipeline readout) |
| `F86.PetermannSpikeWitness(N, blockDim, MaxKGlobal, ArgMaxQ, MaxKInterior, MaxKEndpoint)` | one row of the pinned c=2 N=5..8 Petermann-K table held by `LocalGlobalEpLink`. `Parity` is a derived property (`N % 2 == 0 ? "even" : "odd"`). Frozen empirical anchor (not a live-recomputed witness, the data is fixed at the 2026-05-06 sweep) |
| `F86.C2UniversalShapeDerivation(block)` | **PROOF_F86_QPEAK Item 1 (c=2), Statement 2** top-level synthesis of Stages A–D (E1 integration). Wraps the Stage D2 `C2HwhmRatio` and registers the c=2 universal-shape outcome in the typed-knowledge graph. **Tier 1 candidate** (this session): empirical pipeline reproduces all 8 anchor cases (c=2 N=5..8 × {Endpoint, Interior}) at ≤ 0.005, directional Endpoint > Interior split derived empirically (gap ≈ 0.022), closed-form constant NOT yet derived. Properties: `Block`, `HwhmRatio` (Stage D2 primitive), `InteriorMean` / `EndpointMean` / `DirectionalGap` (canonical-pipeline class means), `IsClosedFormDerived` (false; flips to true on future Tier1Derived promotion via cross-block perturbation, projector-overlap lift, or char-poly factorisation), `Witnesses` (per-bond `HwhmRatioWitness` collection inherited from `C2HwhmRatio`), `PendingDerivationNote` (next-session hand-off documenting the three documented next directions). Promotion path is honest about what the algebra has not yet closed |
| `F86.Item1Derivation.C2BlockShape(block)` | **PROOF_F86_QPEAK Item 1 (c=2)** anchor: elementary block-structure constants (`PnDimension=N`, `PnPlus1Dimension=N(N−1)/2`, `HdEqualsOnePairs=N(N−1)`, `HdEqualsThreePairs=N(N−1)(N−2)/2`); Tier 1 derived from popcount combinatorics |
| `F86.Item1Derivation.C2ChannelUniformAnalytical(block)` | PROOF_F86_QPEAK Item 1 (c=2): closed-form channel-uniform vectors `C1Vector` (HD=1, weight 1/√(N(N−1))) and `C3Vector` (HD=3, weight 1/√(N(N−1)(N−2)/2)) cached on construction; Tier 1 derived, machine-precision overlap with `FourModeBasis` columns 0/1 across N=5..8 |
| `F86.Item1Derivation.C2InterChannelAnalytical(block)` | PROOF_F86_QPEAK Item 1 (c=2): SVD-top inter-channel vectors `U0` ∈ HD=1 subspace and `V0` ∈ HD=3 subspace of V_inter = P_HD1† · M_H_total · P_HD3, plus `Sigma0`. Tier 2 verified (numerical fallback): σ_0 is exactly degenerate at even N (deg=2 at N=6, N=8), making single-direction Tier1Derived ill-posed against MathNet's SVD tiebreaker. `IsAnalyticallyDerived = false`; `PendingDerivationNote` summarises ansätze tried (ψ_k(s)·ψ_k(e) products fail; u_matrix is rank ≥ 2) and the projector-onto-2D-eigenspace lift as the cleanest next direction |
| `F86.Item1Derivation.C2BondCoupling(block)` | PROOF_F86_QPEAK Item 1 (c=2): per-bond projected matrix V_b = B† · M_H_per_bond[b] · B in the 4-mode basis, plus dissipator block D_eff = B† · D · B. **Stage B complete** — all three V_b sub-blocks plus full-matrix accessor and anti-Hermiticity guard. **Stage C1 added** — D_eff diagonal closed form. **B1 (probe-block):** `ProbeBlockEntry(bond, alpha, beta)` returns ⟨c_α \| M_H_per_bond[b] \| c_β⟩ for α, β ∈ {0,1} via composition of A2's `C2ChannelUniformAnalytical`; structurally Tier 1 derived in isolation. F73 sum-rule verified: Σ_b V_b[α, β] = 0 for α ≠ β at machine precision across N=5..8. **B2 (cross-block):** `CrossBlockEntry(bond, alpha, j)` returns ⟨c_α \| M_H_per_bond[b] \| x_j⟩ for α ∈ {0,1}, j ∈ {2,3} (x_2 = \|u_0⟩, x_3 = \|v_0⟩); composed with A3's `C2InterChannelAnalytical`. Per-(N, b) entries match `FourModeEffective.MhPerBondEff[b][α, j]` at 1e-12 (library-internal consistency). The `CrossBlockWitnesses` collection exposes per-bond (BondClass-tagged) cross-block 2×2 entries plus Frobenius norm — the bond-position-dependent fingerprint that splits Endpoint vs Interior at c=2 (witness-level: empirically Endpoint < Interior cross-block Frobenius mean across N=5..8; downstream HWHM_left/Q_peak ratio emerges from the full 4×4 in Stage C/D). **B3 (SVD-block + AsMatrix + anti-Hermiticity guard):** `SvdBlockEntry(bond, j, k)` returns ⟨x_j \| M_H_per_bond[b] \| x_k⟩ for j, k ∈ {2,3} (same projection-formula pattern as ProbeBlockEntry / CrossBlockEntry, inherits Tier 2 from A3). `AsMatrix(bond)` assembles the full 4×4 V_b from the three sub-blocks, with bottom-left filled via the anti-Hermitian relation V_b[j, α] = -conj(V_b[α, j]); matches `FourModeEffective.MhPerBondEff[b]` entry-by-entry at 1e-12. The `Vb_IsAntiHermitian_AcrossAllBondsAndEntries` guard verifies ‖V_b + V_b†‖_F < 1e-10 across all bonds and N=5..8 — catches sign drift within each Hermitian-paired sub-block (probe-block, SVD-block); cross-block sign consistency is established by the parallel `AsMatrix_FullVb_MatchesFourModeEffective` check against B† M B (entry-by-entry comparison to the independently computed `FourModeEffective.MhPerBondEff[b]`). **C1 (D_eff diagonal closed form):** `DEffDiagonal()` returns the 4×4 dissipator block diag(−2γ₀, −6γ₀, −2γ₀, −6γ₀) and `DEffDiagonalEntry(i)` returns the scalar diagonal entry. **Tier 1 derived (structural):** D[i,i] = −2γ₀·HD(p,q) is HD-diagonal (F73 generalisation), and each of \|c_1⟩, \|c_3⟩, \|u_0⟩, \|v_0⟩ lives entirely in one HD subspace (HD=1 for indices 0/2; HD=3 for indices 1/3), so D_eff = diag(−2γ₀·HD_i) by orthonormality of B. Off-diagonals exactly zero independent of A3's SVD-direction obstruction (D acts as a scalar inside each HD subspace). Verified at 1e-14 vs closed form and 1e-12 vs `FourModeEffective.DEff` across N=5..8. **Class-level Tier 2 verified** (inherits A3's obstruction via B2/B3: σ_0 of V_inter exactly degenerate at even N, single-direction \|u_0⟩, \|v_0⟩ library-tiebreaker-dependent). Probe-block (B1) and D_eff (C1) sub-surfaces are structurally Tier 1; class-level Tier reflects the weakest link |
| `F86.Item1Derivation.C2KShape(block)` | PROOF_F86_QPEAK Item 1 (c=2), Stage D1: the K_CC_pr observable evaluated in the 4-mode basis as a closed-form Duhamel formula. **Tier 1 derived** (the Duhamel formula is exact in its inputs; the upstream Tier 2 numerical eigenvalues do not propagate — Tier is per-claim, analogous to C1's D_eff Tier1Derived sub-fact inside the Tier2 BondCoupling class). Formula: K_b(Q, t) = 2·Re⟨ρ(t) \| S_kernel \| ∂ρ/∂J_b⟩ with ρ(t) = R·diag(e^(λ_i·t))·R⁻¹·ρ_0 and ∂ρ/∂J_b = R·(X_b ⊙ I_jk)·R⁻¹·ρ_0 (X_b = R⁻¹·V_b·R; I_jk = (e^(λ_k·t)−e^(λ_j·t))/(λ_k−λ_j) or t·e^(λ_j·t) at degeneracy). `KAt(Q, bond, t)` → real K value at a single (Q, bond, t); `PeakOverT(Q, bond, tGrid)` → (peak \|K\|, t at peak) over a t-grid; `PeakOverDefaultT(Q, bond)` → peak \|K\| over the default 21-point [0.6, 1.6]·t_peak grid. **Verification:** `PeakOverT` matches `FourModeResonanceScan.ComputeKCurve` bit-exact (1e-10 across N=5, 7 over the default Q × t grids); the two paths share the same Duhamel structure but C2KShape repackages the inner loop as a callable per-Q primitive so D2 can call it at any Q without spinning up the full Q-scan |
| `F86.Item1Derivation.C2HwhmRatio(block)` | PROOF_F86_QPEAK Item 1 (c=2), Stage D2: HWHM_left/Q_peak ratio per bond + class-mean for c=2 universal-shape statement. **Tier 1 candidate** (directional Endpoint > Interior split derived empirically across N=5..8; closed-form constant NOT pinned this session). Tier-1-derived universal constants exposed as `BareDoubledPtfXPeak = 2.196910` (post-EP location in dimensionless x = Q/Q_EP) and `BareDoubledPtfHwhmRatio = 0.671535` (HWHM_left/Q_peak SVD-block floor); empirical Interior/Endpoint sit above this floor by ~0.08-0.10 (Direction (b) finding). `ComputePerBond(bond)` returns (Q_peak, K_max, HWHM_left) for a single bond; `HwhmLeftOverQPeakMean(BondClass)` returns the canonical-pipeline class-averaged ratio (average bond-class K(Q) curves first, then peak/HWHM; matches Python `_eq022_b1_step_e_resonance_shape.py` to ≤ 0.005 across all 8 anchor cases at γ₀ = 0.05). `Witnesses` exposes per-bond (BondClass-tagged) (Q_peak, K_max, HWHM_left, HWHM_left/Q_peak): 16 witnesses pinned at N=5..8, preserving the F71-orbit substructure (mid-chain vs flanking Interior bonds have distinct Q_peak per PROOF_F86_QPEAK "Per-F71-orbit substructure" note). Empirical pipeline uses the full-block `ResonanceScan` (not the 4-mode `C2KShape`; see F86KB's `FourModeInsufficiencyNote`: 4-mode reduced model is bond-class-blind for the HWHM ratio). `TwoLevelEpDecaySanity` = sanity-check value from the 2-level EP-rotation propagator-magnitude profile; bond-class-blind by construction (the V_b cross-block does not enter at this 2-level order), so it is NOT a derivation of the K-resonance HWHM. `IsAnalyticallyDerived = false`; `PendingDerivationNote` documents (i) the 2-level reduction + 4×4 char poly factorisation context, (ii) the 2026-05-06 Direction (b) doubled-PTF baseline (universal x_peak 2.197, HWHM/Q* 0.6715), (iii) the 2026-05-06 Direction (a') probe-block 2-level resonance attempt: STRUCTURALLY FALSIFIED, V_b probe-block is bond-class-blind for every bond at every N=5..8 (diagonal +i·c·I scalar; off-diagonal exactly zero per bond by F73 sum-rule applied per-bond). Cross-block Frobenius is unstable across N (varies 0.640, 1.318, 0.815, 0.143 due to A3's σ_0 degeneracy). SVD-block off-diagonal V_b[2, 3] IS the bond-class signature carrier (Endpoint < Interior consistently; direction OPPOSITE the empirical HWHM/Q* split). The 4-mode reduction reproduces only the ~0.673 floor, not the empirical 0.7506 lift. **Directions historical** (2026-05-06 letter-labeled taxonomy, superseded 2026-05-11 by the F90 bridge identity + F89 path-k closure): (a'') SVD-block 2-level resonance via V_b[2, 3] — retired (4-mode-blind, see C2DirectionAFalsificationProbe); (b'') full block-L derivation, not 4-mode — **SURVIVING direction**, numerically closed via F90 bridge (`F90F86C2BridgeIdentity`), analytical step is F89 AT-locked F_a/F_b + H_B-mixed octic residual per `F86HwhmClosedFormClaim`; (c'') three-block superposition — retired (see C2DirectionCFalsificationProbe); (d'') lift |u_0⟩, |v_0⟩ to projector-overlap — foundational lift (precondition for (a''); see `C2InterChannelProjector`); (e'') symbolic char-poly factorisation at Q_EP — retired (C2EffectiveSpectrum's cubic-c_3 obstruction). The active runtime `PendingDerivationNote` (see `C2HwhmRatio.BuildPendingDerivationNote`) points at `F86HwhmClosedFormClaim` and `PROOF_F90_F86C2_BRIDGE.md` directly. Historical 2026-05-06 context: `docs/superpowers/syntheses/2026-05-06-direction-a-prime-attempt.md` |
| `F86.Item1Derivation.C2EffectiveSpectrum(block)` | PROOF_F86_QPEAK Item 1 (c=2): the four eigenvalues of L_eff(Q) = D_eff + Q·γ₀·Σ_b V_b in the 4-mode basis, plus per-(Q, bond) identification of the K-driving eigenvalue pair (Stage C3). `Eigenvalues(Q, bond)` returns Complex[4] sorted by Re desc / Im asc; the bond parameter is for API uniformity with downstream C3/D consumers — at c=2 the spectrum is bond-independent because uniform-J means the spectrum-relevant matrix is the bond-summed M_h_total_eff. `LEffAtQ(Q)` exposes the assembled 4×4 (matches `FourModeEffective.LEffAtQ(Q)` at 1e-12). **C2 (Tier 2 verified):** the C2 time-box probed the structural-factorisation question (does the 4×4 char poly split into two 2×2 quadratics under any natural similarity transform?) and found NO factorisation. Three obstructions: (a) HD-parity split has commutator residual ~0.4 because the [u_0, v_0] entry of M_h_total_eff is the dominant inter-channel SVD-top coupling (≈ σ_0 ≈ 2√2 asymptotically), not perturbative; (b) probe vs SVD split has commutator residual ~0.5 because the cross-block V_b[α, j] is the bond-class-driving fingerprint, nonzero by construction; (c) chain-mirror R is +I in 4-mode basis at odd N (cannot factorise) and not exactly diagonal at even N (inherits A3's σ_0 degeneracy). The char poly det(λI − L_eff(Q)) is a genuine quartic in (λ, Q): coefficient c_3(Q) has a CUBIC term, ruling out (λ² − aλ + b)(λ² − cλ + d) with closed-form rational coefficients in Q. F86a's 2×2 reduction `[[−2γ₀, +iJ·g_eff], [+iJ·g_eff, −6γ₀]]` captures the dominant EP physics for the slowest pair as an *approximation* (eigenvalues from the full 4×4 deviate ~0.5% from Statement 1 at N=5, Q=1), not an exact factorisation. Eigenvalues computed via `Matrix<Complex>.Evd().EigenValues` directly on L_eff(Q); both Tier1 and Tier2 paths share the test contract (1e-10 agreement vs `FourModeEffective.LEffAtQ(Q).Evd()`). `IsAnalyticallyDerived = false`; `PendingDerivationNote` records the structural ansätze tried (HD-perm, flip-c, flip-u, flip-both, HD-with-flip), the char-poly evidence (cubic c_3 term), and three promising next directions: (a) approximate factorisation via perturbation in the small probe ↔ SVD coupling → Tier1Candidate with explicit error term; (b) lift to the c=1 PTF framework as a perturbative double-c=1 system; (c) symbolic CAS reduction using the (γ₀, J·g_eff) substitution λ = −4γ₀ + μ·γ₀. **C3 (K-driving pair identification):** `KDrivingPair(Q, bond)` returns `(LamPlus, LamMinus)` — the 2 of 4 eigenvalues whose eigenvectors have largest squared probe overlap |⟨probe \| w_i⟩|². `KDrivingPairIndices(Q, bond)` returns the matching `(IndexPlus, IndexMinus)` ∈ {0, 1, 2, 3} (distinct, ordered so LamPlus has the larger Re). `ProbeOverlapsSquared(Q, bond)` exposes the per-eigenvector overlaps in the same Re-desc order as `Eigenvalues(Q, bond)`. `ProbeProjection` exposes the 4-mode-basis projection of the Dicke probe (cached on construction). **Tier 1 structural sub-fact** (verified by `ProbeProjection_HasZeroSvdTopComponents`): probe ⊥ {\|u_0⟩, \|v_0⟩} at machine precision (< 1e-12) across N=5..8 — the probe lives entirely in the channel-uniform 2D subspace span{\|c_1⟩, \|c_3⟩} per InterChannelSvd's structural finding ("the probe (Dicke state) is orthogonal to \|u_0⟩, \|v_0⟩"). **Tier 2 numerical readout:** the per-(Q, bond) identification of which 2-of-4 eigenvectors are closest to span{\|c_1⟩, \|c_3⟩} depends on the (γ₀, J)-rotation and is computed via `Evd()` + |⟨probe \| w_i⟩|² ranking. The K-driving pair carries dominant probe content (overlap fraction > 0.7 across N=5..8 at Q ∈ {1.0, 1.5, 2.0}). At Q ≈ Q_EP both eigenvalues approach Re(λ) = −4γ₀ per F86a (verified to fall in the [−9γ₀, −γ₀] band across the test grid; the precise EP collapse requires the closed-form g_eff that Stage D will pin). The class-level Tier remains `Tier2Verified` (the C2 char-poly cubic obstruction governs); the C3 structural sub-fact is documented in the `ProbeProjection` XML doc and pinned by its dedicated test |

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
