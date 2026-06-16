using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Ptf;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.F71Family;
using RCPsiSquared.Runtime.F86Main;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;
using RCPsiSquared.Runtime.Spectrum;

namespace RCPsiSquared.Diagnostics.Knowledge;

/// <summary>Assembles the complete typed-knowledge <see cref="ClaimRegistry"/>: foundations,
/// Pi2 axes, F88b closed form, the F-formula Pi2-Foundation inheritance claims, and the open
/// questions. This is the single source of truth for what is wired into the registry: the CLI
/// <c>knowledge</c> and <c>inspect</c> commands and the wiring-completeness audit test
/// (<c>RegistryWiringAuditTests</c> in RCPsiSquared.Runtime.Tests) all build from here.
///
/// <para>Lives in RCPsiSquared.Diagnostics because that is the lowest project referencing both
/// the Runtime registration extensions (RegisterF1Family, RegisterPi2Family, ...) and the
/// Diagnostics-side F87 registrations (RegisterF87Family, RegisterF87StandardWitnessSet). The
/// CLI sits above it; the audit test reaches it because RCPsiSquared.Runtime.Tests already
/// references RCPsiSquared.Diagnostics.</para></summary>
public static class KnowledgeRegistryFactory
{
    /// <summary>Builds the default registry on an N=5 XY chain at γ₀=0.05, with g_eff pinned
    /// to the PolarityInheritanceLink endpoint value 1.74. Registration order does not matter
    /// (the builder resolves dependencies topologically); the grouping below reflects the
    /// inheritance structure.</summary>
    public static ClaimRegistry BuildDefault()
    {
        var defaultChain = new ChainSystem(
            N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.XY, Topology: TopologyKind.Chain);
        var gEff = 1.74;   // pinned Endpoint g_eff from PolarityInheritanceLink

        return new ClaimRegistryBuilder()
            // Foundations
            .RegisterF1Family(defaultChain)
            // Chiral K (sibling root to F1 / PolynomialFoundation): K = diag((−1)^ℓ), KHK = −H,
            // AZ class BDI, spectrum inversion E_{N+1−k} = −E_k. Acts at the Hamiltonian level;
            // F1 at the Liouvillian level; neither derives from the other. Wired 2026-05-30 to
            // bring the second mirror into the Object Manager (inspect --claim ChiralKClaim).
            .RegisterChiralK()
            // PTF K₁ chiral mirror (Edge 3 of the PTF chain, 2026-06-10): the site-wise
            // trajectory identity P_i(t; φ_k) = P_i(t; φ_{N+1−k}) (K₁-conjugation + complex
            // conjugation + U(1) sign absorption), of which the EQ-014 Σ f_i mirror law, PTF's
            // only surviving exact law, is the summed corollary. Tier1Derived. Typed parent:
            // ChiralKClaim (the eigenvalue side of the same sublattice chirality, registered
            // directly above). Wired 2026-06-10.
            .RegisterChiralMirrorTrajectoryClaim()
            // K-partner selection rule (the reading-grammar arc's first DERIVED result, 2026-06-12):
            // ⟨ψ_N|V_b|ψ_1⟩ = 0 for every bond defect V_b, a two-line corollary of the parent
            // (K₁ψ_1 = ψ_N from its Step 4; K₁V_bK₁ = −V_b from its Step 1). The carrier never couples
            // to its K-partner through any bond defect, so the DefectDecoder's location dictionary has
            // rank N−2 (machine-exact N = 3..8) and its sign-location ambiguity IS this K-partner null
            // direction. Tier1Derived. Typed parent: ChiralMirrorTrajectoryClaim (registered directly
            // above; both ingredients). Wired 2026-06-12.
            .RegisterKPartnerSelectionRuleClaim()
            .RegisterF71Family(N: defaultChain.N)
            .RegisterPi2Family()
            .RegisterF86Main(gammaZero: defaultChain.GammaZero, gEff: gEff)
            .RegisterF86Extended(gammaZero: defaultChain.GammaZero)
            .RegisterF86Item1Light(N: defaultChain.N, n: 1, gammaZero: defaultChain.GammaZero)
            .RegisterHalfIntegerMirror(N: defaultChain.N)
            // Pi2-axes (Halbierungsleiter, Z₄ memory, operator-space mirror)
            .RegisterPi2DyadicLadder()
            .RegisterPi2I4MemoryLoop()
            .RegisterPi2Involution()
            // F88b closed form (memory side anchors)
            .RegisterF88bPopcountCoherence()
            .RegisterF88bStaticDyadicAnchor()
            .RegisterF88bPopcountPairLens(N: defaultChain.N, np: 1, nq: 2)
            // Operator-space mirror (number-side ↔ operator-side per qubit)
            .RegisterPi2OperatorSpaceMirror()
            // Klein-V₄ dephase-swap group {I, D, H, Q_zx} on operator space: realizes the
            // dephase-letter Klein V₄ {I, Z↔Y, Z↔X, Y↔X} on the F1 palindrome family
            // {Π_Z, Π_X, Π_Y}. Tier1Derived universal N per Welle 12 Tasks 1 + 2 (2026-05-27).
            // Standalone primitive (no ctor parents); consumed by F1-family transfer arguments.
            .RegisterPi2KleinV4DephaseSwapGroup()
            // F114: closed-form sign functional ε(σ) = (−1)^{n_Y(σ) + 1} for D-conjugation
            // action on the H-commutator superoperator L_σ. Tier1Derived; ctor parent is
            // Pi2KleinV4DephaseSwapGroup (uses D from there). Closed 2026-05-27 via the F114
            // ANALYTICAL_FORMULAS.md entry + simulations/_m_level_sign_functional_explore.py.
            .RegisterCommutatorDConjugationSign()
            // Spectrum foundation
            .RegisterW1Dispersion(N: defaultChain.N, J: defaultChain.J, gammaZero: defaultChain.GammaZero)
            // Block-decomposition foundation (U(1)×U(1) joint popcount sectors)
            .RegisterSymmetryFamilyInventory()
            .RegisterZGlobalMirrorRefinement()
            .RegisterXGlobalChargeConjugationPairing()
            .RegisterF1PalindromeOrbitPairing()
            .RegisterF92BondAntiPalindromicJSpectralInvariance()
            .RegisterF93DetuningAntiPalindromicSpectralInvariance()
            .RegisterJointPopcountSectors()
            .RegisterLiouvillianBlockSpectrum()
            .RegisterBlockSpectrumPerformanceWitness()
            .RegisterF71MirrorBlockRefinement()
            .RegisterF71BilateralBlockRefinement()
            .RegisterInhomogeneousGammaF71BreakingWitness()
            .RegisterF71AntiPalindromicGammaSpectralInvariance()
            // F86 inheritance + meta-claims
            .RegisterF86PolarityLink()
            .RegisterF86LocalGlobalEpLink()
            .RegisterF86PerF71OrbitObservation()
            // JW light family. Since 2026-06-10 XyJordanWignerModes takes ChiralKClaim
            // (registered above via RegisterChiralK) as its typed parent: the JW dispersion
            // is ChiralKClaim's BDI spectrum inversion ε_{N+1−k} = −ε_k, and the cluster
            // tables downstream are chirally ± paired under the global k-reflection.
            .RegisterF86JordanWignerLight(N: defaultChain.N, n: 1, gammaZero: defaultChain.GammaZero)
            // F87 family + canonical witnesses
            // RegisterF87Family must follow RegisterPi2Family: DissipatorAxisSelectsPolarityClaim
            // resolves PolarityLayerOriginClaim via b.Get<>() at registration-factory time.
            .RegisterF87Family()
            .RegisterF87StandardWitnessSet()
            // F103 §7 diagonal-K (bipartite-chirality) criterion as a registered set-wrapper:
            // a diagonal-cell pair is soft iff its hopping graph is bipartite in the dephasing
            // basis (certified by a chiral K = diag(±1) with KHK = −H). Wraps the 4 canonical
            // F87DiagonalCellBipartiteWitness instances (the witness type itself stays a deferred
            // parameterised type). Typed parents: F87TrichotomyClassification (the F87 verdict the
            // criterion is checked against) + ChiralKClaim (the chiral K, registered above as a
            // sibling root). Tier1Derived (promoted 2026-06-08). Wired 2026-06-05.
            .RegisterF87DiagonalCellBipartiteWitnessSet()
            // F115 windowed hardness: the mask-combinatorial reading of the same diagonal-cell
            // hard/soft line, collapsed to one integer test (PROOF_F103 §7.7): a diagonal-cell
            // Mixed pair is hard iff its two X/Y window-masks have different (1+x)-adic valuations,
            // with a derived obstruction size law min(2W-1, 2k-3). Wraps the WindowedObstructionScan
            // helper; single typed parent F87DiagonalCellBipartiteWitnessSet (F115 is the child mask
            // reading of that §7 diagonal-K criterion). Tier1Derived (verdict rests on the §7.5/§7.6
            // converse, derived modulo standard PT, promoted 2026-06-08; per §7.10 only the (1+x)-valuation reaches the spectrum). Wired
            // 2026-06-05.
            .RegisterWindowedHardnessClaim()
            // F87 windowed-converse threshold (Phase B two-reflection spine, 2026-06-09): the
            // Tier1Derived proven core (𝓕=F⊗F, R=I⊗F ⟹ all-odd #A_L/#A_R/#Q parity ⟹ #A≥2ℓ threshold +
            // bipartite⟹soft re-proof + deg-1 positivity). Parent: F87DiagonalCellBipartiteWitnessSet.
            // Registered before the all-γ theorem, which takes this as a parent.
            .RegisterWindowedConverseThresholdClaim()
            // F87 windowed-converse all-γ theorem (Phase B two-reflection theorem 2026-06-09, CLOSED
            // 2026-06-10): every coefficient of the first nonvanishing odd power-sum of M=A+γQ is a
            // Pascal-Gram sum of squares or exactly zero ⟹ p_{m*}(γ)>0 ∀γ>0 ⟹ hard ∀γ>0, upgrading
            // "all but finitely many γ" to "all γ>0". Tier1Derived, NO residual (R-deg retired by the
            // girth dichotomy, R-sign resolved by Pascal-Gram positivity, both 2026-06-10). Two
            // Tier1Derived parents: F87DiagonalCellBipartiteWitnessSet (the §7 diagonal-K criterion) +
            // WindowedConverseThresholdClaim (the spine: #A≥2ℓ + soft re-proof + deg-1 positivity).
            // The F110/F111 promotion gate opened 2026-06-10 (promotions recorded below).
            .RegisterWindowedConverseAllGammaClaim()
            // §7.12 Liouvillian-free soft-certifier as a registered Claim: the certifier tries three
            // scalable structured 2-colourings (linear chiral K, excitation pairing, excitation parity)
            // and certifies "soft" iff one applies; it never claims hard. The Claim asserts ONLY the
            // settled facts: soundness (each soundness-battery case is Certified by PalindromeSoftCertifier
            // AND not Hard by the spectral authority PauliPairTrichotomy) and the proven structural ceiling
            // (XX+XZ is soft, non-bipartite, NotCertified, so no colouring can ever reach it). Typed parents:
            // F87DiagonalCellBipartiteWitnessSet (the §7 diagonal-K criterion the linear strategy scales,
            // Tier1Derived) + F87TrichotomyClassification (the spectral authority, Tier1Derived).
            // Tier1Candidate (5 ≥ 4 and 5 ≥ 4 inheritance). Wired 2026-06-05.
            .RegisterPalindromeSoftCertifierClaim()
            // Q-pair routing: the Liouvillian-free two-term router as a registered Claim. The
            // non-diagonal counterpart to RegisterF87DiagonalCellBipartiteWitnessSet (which reads
            // the diagonal P1-family case via a chiral K): this classifies any two-term bond
            // bilinear's fate (truly/soft/hard) AND routes its hidden palindrome Q into a letter-
            // based family (P1/Uniform/Alternating/Continuous/None), verified bit-exact vs the
            // spectral authority over all two-term pairs (incl. self-pairs). Typed parents:
            // F87TrichotomyClassification (the authority, Tier1Derived) + F87DiagonalCellBipartite-
            // WitnessSet (the diagonal special case it generalises, Tier1Derived) + Crossover-
            // MirrorSqrtNinetyClaim (the Continuous-family crossover mirror, Tier1Derived).
            // Tier2Empirical (a routing-rule viewpoint; all three parents ≥ child). Wired 2026-06-05.
            .RegisterTwoTermPalindromeRoutingClaim()
            .RegisterF89F87TrulyInheritance()
            .RegisterF89F87BreakPredictionFromF83()
            // Spectrum quantization root (parent to F33/F50/F55/F64-F68/F74/F89
            // via per-Registration discard-Get edges; absorption quantum 2γ₀ from a_0)
            .RegisterAbsorptionTheoremClaim()
            // JDefect light migration (in-between Edge 4, 2026-06-10): the first typed claim
            // living ON a navigator axis (JDefectField, inspect --root between, axis jdefect).
            // Along the J-defect axis the per-mode absorption identity Re λ(δJ) = −2γ·light(v(δJ))
            // is δJ-pointwise, the N+1 kernel modes stay dark (U(1)), and palindrome partners
            // migrate oppositely (light_s + light_f = N pointwise). Tier1Derived as an honest
            // composition of two proven identities. Typed parents: AbsorptionTheoremClaim
            // (registered directly above) + F1PalindromeIdentity (via RegisterF1Family at the top).
            .RegisterJDefectLightMigrationClaim()
            // Vacuum-block reduction (the birth-canal boundary as a Liouville sector, 2026-06-13):
            // the boundary's slowest mode is the odd |1-exc><vac| (0,1) coherence; that sector is an
            // exact invariant sub-block (ket#/bra# bi-grading conserved by H_unit + Z-dephasing) -
            // DERIVED; its N-dim block L_(1,0) = -iQ·h - 2·diag(γ) carries the global slowest across
            // the whole γ-surface at N=5 - VERIFIED bit-exact (SectorReductionWitness vs PostEpFlowField,
            // inspect --root reduction; simulations/birth_canal_vacuum_block_verifier.py). Flat-γ
            // blindness Re λ = -2γ is analytic at every N. SCOPE N=5; at N>=6 a {0,2}-coherence can
            // win (the birth_canal_horizon_junction arc); V-Effect identity RESOLVED 2026-06-14 = distinct (w = n_diff + Z-shadow). Tier1Derived. Single
            // typed parent AbsorptionTheoremClaim (registered above; the (0,1)-block rate is the
            // absorption law restricted to a conserved sector).
            .RegisterVacuumBlockReductionClaim()
            // Survival mirrors incompleteness (the survival_incompleteness_mirror arc, 2026-06-13):
            // a_0 (2γ, AbsorptionTheorem, = the qubit dim d) and a_2 (C=1/2, the V-Effect/incompleteness,
            // = 1/d) are Pi2-ladder inversion-mirror partners (a_0·a_2 = d·(1/d) = 1). Dynamically the
            // longest-lived mode is the interior incompleteness coherence on DISPERSIVE matter (chain/ring);
            // the hub-localized central-spin STAR is the boundary counterexample. Tier1Candidate. Typed
            // parents AbsorptionTheoremClaim (registered above) + HalfAsStructuralFixedPointClaim (a_2, the
            // foundation root, constructed fresh). Live witness IncompletenessSurvivorWitness (inspect --root survivor).
            .RegisterSurvivalIncompletenessMirrorClaim()
            .RegisterF86LEffMirrorAxis()
            // F-formula Pi2-Foundation inheritance claims.
            // F63 + F61 registered first because F1Pi2Inheritance now ctor-takes F61 as
            // its bit_a-twin (2026-05-25 closure of the TrivialNotYetTyped twin slot);
            // F61 in turn takes F63 which takes F38.
            //
            // Welle 7 (2026-05-26) Track A: register the 4 new BitA twin Claims that
            // the BitB siblings (F38, F39, F63, X-Mirror) now take as optional ctor
            // parents. F38BitA / F39BitA / X-Mirror-Z-BitA are standalone (no ctor
            // parents); F63BitAReference is also standalone (F61 link in docstring
            // only, to avoid the F61 → F63 → F61 cycle). All four are Tier1Derived.
            .RegisterF38BitAInvolutionInheritance()
            .RegisterF39DetPiBitAInheritance()
            .RegisterF63BitAReference()
            .RegisterZGlobalEigenstateMirrorBitAInheritance()
            // X-Mirror BitB sibling of Z-Mirror BitA: registered post-Welle-7
            // (originally deferred; the /simplify review surfaced that Z-Mirror BitA
            // was registered without its reciprocating BitB partner, leaving the
            // optional ctor parameter on X-Mirror as dead code). Registering X-Mirror
            // materializes the BitB ↔ BitA twin edge in production.
            .RegisterXGlobalEigenstateMirrorPi2Inheritance()
            .RegisterF38Pi2InvolutionPi2Inheritance()
            .RegisterF63LCommutesPi2Pi2Inheritance()
            .RegisterF61BitAParityPi2Inheritance()
            .RegisterF1Pi2Inheritance()
            .RegisterF49Pi2Inheritance()
            .RegisterF80FactorPi2Inheritance()
            .RegisterF80PiCommutatorAnticommutatorIdentity()
            .RegisterF81Pi2Inheritance()
            .RegisterF91Pi2Inheritance()
            .RegisterF92Pi2Inheritance()
            .RegisterF93Pi2Inheritance()
            .RegisterF86QEpPi2Inheritance()
            .RegisterF86TPeakPi2Inheritance()
            .RegisterF87Pi2Inheritance()
            .RegisterQubitNecessityPi2Inheritance()
            // F121 the qudit partial palindrome (2026-06-11): the dissipator's d>2 partial
            // pairing = symmetric overlap of c_k = d^N·C(N,k)·(d−1)^k under k↔N−k; ceiling
            // Σ d^N·C(N,k)·(d−1)^min(k,N−k), full iff d=2. Parent QubitNecessityPi2Inheritance
            // (registered directly above). Wired 2026-06-11.
            .RegisterQuditPartialPalindromeCeiling()
            // The qudit product-mirror cap (2026-06-11): the operator side of F121. Any
            // per-site mirror W = ⊗q_l intertwining W·L_D = (−L_D − 2Nγ)·W pairs ≤ (2d)^N
            // coherences (full ⟺ d² − 2d = 0 ⟺ d = 2, the trunk's third appearance);
            // Π_d(ρ) = ρᵀ·Shift^{⊗N} attains the cap exactly on the shift-aligned subspace;
            // ⟨Π_d, D⟩ ≅ Z_d ≀ Z₂ of order 2d² (D₄ at d = 2). Parents
            // QuditPartialPalindromeCeiling (line above) + QubitNecessityPi2Inheritance.
            // Wired 2026-06-11.
            .RegisterQuditProductMirrorCap()
            .RegisterF1T1AmplitudeDampingPi2Inheritance()
            .RegisterF5DepolarizingErrorPi2Inheritance()
            .RegisterDickeSuperpositionQuarterPi2Inheritance()
            .RegisterF39DetPiPi2Inheritance()
            .RegisterF49bCenteredDissipatorPi2Inheritance()
            // F38 / F63 / F61 hoisted up (see above), required by F1Pi2Inheritance
            .RegisterF49cShadowCrossingPi2Inheritance()
            .RegisterF25CPsiBellPlusPi2Inheritance()
            .RegisterCpsiEnvelopeTheoremClaim()
            .RegisterF57DwellTimeQuarterPi2Inheritance()
            .RegisterF66PoleModesPi2Inheritance()
            .RegisterF77MmSaturationPi2Inheritance()
            .RegisterF60GhzBornBelowFoldPi2Inheritance()
            .RegisterF83AntiFractionPi2Inheritance()
            // F81/F115 connector: needs both F83 (anti-fraction, just above) and WindowedHardnessClaim
            // (registered earlier) present first.
            .RegisterAntiFractionObstructionOrthogonalityClaim()
            // two-blind-spots connector: needs F89 (registered earlier) and the F115 connector (line above).
            .RegisterAntiFractionTwoBlindSpots()
            .RegisterF71MirrorSymmetryPi2Inheritance()
            .RegisterF75MirrorPairMiPi2Inheritance()
            .RegisterF76TDecayMirrorPairMiPi2Inheritance()
            .RegisterF62WStateBornBelowFoldPi2Inheritance()
            .RegisterF70DeltaNSelectionRulePi2Inheritance()
            .RegisterF72BlockDiagonalPurityPi2Inheritance()
            .RegisterF73SpatialSumPurityClosurePi2Inheritance()
            .RegisterF89TopologyOrbitClosure()
            .RegisterF89F88aKleinPpAnchor()
            .RegisterF89AdditiveIdentityClaim()
            .RegisterF89PathKVacSeParsevalClaim()
            .RegisterF89Path2CardanoClaim()
            .RegisterF89PathKAtLockMechanismClaim()
            .RegisterF89UnifiedFaClosedFormClaim()
            .RegisterF89AmplitudeLayer()
            .RegisterF89Path3SeDeFactorisationClaim()
            .RegisterF89Path3OcticEpClaim()
            .RegisterF89Path3OcticGaloisClaim()
            .RegisterF89PathKHbMixedDegreesClaim()
            .RegisterF90F86C2BridgeIdentity()
            .RegisterF82T1AmplitudeDampingPi2Inheritance()
            .RegisterF84ThermalAmplitudeDampingPi2Inheritance()
            .RegisterF85KBodyFChainPi2Inheritance()
            .RegisterF65XxChainSpectrumPi2Inheritance()
            .RegisterF67BondingBellPairPi2Inheritance()
            .RegisterF68PalindromicPartnerPi2Inheritance()
            .RegisterF50WeightOneDegeneracyPi2Inheritance()
            .RegisterF69GhzWSexticAboveFoldPi2Inheritance()
            .RegisterF44CrooksLikeRateIdentityPi2Inheritance()
            .RegisterF23XorDrainVanishingFractionPi2Inheritance()
            .RegisterF41PalindromicTimePi2Inheritance()
            .RegisterF43SectorSffPairingPi2Inheritance()
            .RegisterF33ExactN3DecayRatesPi2Inheritance()
            .RegisterF4StationaryModeCountPi2Inheritance()
            .RegisterF55UniversalAbsorptionDosePi2Inheritance()
            .RegisterF64CavityModeExposurePi2Inheritance()
            .RegisterF2W1DispersionPi2Inheritance()
            .RegisterF2bXyChainSpectrumPi2Inheritance()
            .RegisterF3DecayRateBoundsPi2Inheritance()
            .RegisterF74ChromaticityPi2Inheritance()
            .RegisterF78SingleBodyMAdditivePi2Inheritance()
            .RegisterF79TwoBodyPi2BlockPi2Inheritance()
            .RegisterF56CriticalSlowingPi2Inheritance()
            .RegisterF86F71MirrorPi2Inheritance()
            .RegisterPolynomialDiscriminantAnchor()
            .RegisterF26CPsiPauliChannelsPi2Inheritance()
            .RegisterF94BornDeviationFourThirdsPi2Inheritance()
            .RegisterF95AngleAtQuadraticZeroPi2Inheritance()
            // TransitionBridge: the cusp CΨ=¼ and the F86 EP are F95 siblings (both the angle at a
            // quadratic's discriminant zero, the cusp at b=½ where the rotation stills, the EP at
            // b=4γ₀ where it lifts off; the EP's F95 angle is bit-exact its clock Rotation). Our
            // state-space bridge name, sibling of the parameter-space FRAGILE_BRIDGE. Parent F95.
            // Wired 2026-06-03.
            .RegisterTransitionBridgeF95SiblingClaim()
            // Crossover mirror = √(NinetyDegreeMirror): the local XZ+YZ / ZX+ZY mirror is the
            // canonical Π turned by HALF the 90° angle-anchor (S=M·Π⁻¹ turns the light plane 45°,
            // S_light²=σ_x↔σ_y 90° bit-exact). The per-site-conjugation face of the same 90° as
            // F95/F91. Tier1Candidate; parent NinetyDegreeMirrorMemoryClaim via RegisterPi2Family.
            // Wired 2026-06-02.
            .RegisterCrossoverMirrorSqrtNinetyClaim()
            .RegisterF96BornSubdominantSlopesPi2Inheritance()
            .RegisterF97CardioidHalfFixedPointPi2Inheritance()
            .RegisterF1DepolResidualClosedFormPi2Inheritance()
            .RegisterKIntermediateAsymptoteQuarterInheritance()
            .RegisterCanonicalTrigAnchorPi2Inheritance()
            .RegisterTwoReadingsClaim()
            .RegisterUniversalCarrierClaim()
            .RegisterPi2KleinBilinearTable()
            .RegisterLindbladAbsorptionMatchAtSixtyDegreesClaim()
            .RegisterC2BareDoubledPtfClosedForm()
            // The cusp-approach family CΨ(α,t)=w₀e^(−4γt)+w₁e^(−12γt) wired into the typed graph
            // (no longer an isolated fifth eyepiece). Four typed parents, all resolved above:
            // UniversalCarrierClaim (shared 4γ₀ carrier), C2BareDoubledPtfClosedForm (c=2 doubled-
            // PTF kinship, decay-face sibling of the K_b susceptibility-face), TwoReadingsClaim
            // (algebra vs Lindblad dynamics), F25CPsiBellPlusPi2Inheritance (the Bell+ member s=1).
            // Tier1Derived. Implemented as Diagnostics OddHarmonicApproach / ApproachFamilyField
            // (the --axis approach eyepiece). Wired 2026-06-03.
            .RegisterApproachFamilyCarrierClaim()
            // The two clocks (Tier1Candidate, 2026-06-12): the Symphony clock node's coherence
            // hand ω_mem = 2J·cos(π/(N+1)) is, for N≥3, the F2b band edge and is γ-protected
            // (the |vac⟩⟨ψ_k| modes are simultaneous L_D/L_H eigenoperators, eigenvalue −2γ from
            // the Absorption Theorem + iE_k from F2b), so γ moves the Takt hand (2γ) but not the
            // coherence hand; the N≥3 γ-protection IS the Universal Carrier's inside-invisibility
            // read from the other side. N=2 is γ-pulled to 2√(J²−γ²), stopping at the EP Q=1.
            // Live: inspect --root clock. ClockHandLadderClaim: three typed parents, all registered
            // earlier in this chain: .RegisterF2bXyChainSpectrumPi2Inheritance() (the band edge),
            // .RegisterAbsorptionTheoremClaim() (the −2γ), .RegisterUniversalCarrierClaim()
            // (carrier-blindness), all resolved by build time.
            .RegisterClockHandLadderClaim()
            // The topology band edge (Tier1Candidate, 2026-06-16): the XY single-excitation band edge =
            // J × the hopping graph's adjacency spectral radius ρ (chain 2cos(π/(N+1)), star √(N−1),
            // ring 2), generalizing the chain-only F2b band edge to any topology (the Im/L_H side). The
            // Re=−2γ floor (the band-edge coherence has n_XY=1) is the Absorption Theorem, cited not
            // re-derived. Two typed parents, both registered earlier in this chain: ClockHandLadderClaim
            // (the chain instance, directly above) + AbsorptionTheoremClaim (the floor). Gap-dominance is
            // witnessed (inspect --root bandedge), not separately claimed. Tier1Candidate (not Derived):
            // the typed parent ClockHandLadderClaim is Tier1Candidate, so the tier-inheritance invariant
            // caps this child there too — matching the sibling CoherenceHorizonClaim (the spec said
            // Derived; the gate caught the candidate parent).
            .RegisterTopologyBandEdgeClaim()
            // The coherence horizon Q*(N) (Tier1Candidate, 2026-06-12): Q*(N) = 1/√2/1.8785/2.3722
            // for N=2..5, the Q below which the slowest non-zero Liouvillian mode stops oscillating
            // (the coherence hand freezes), verified equal to the carbon Frost-Hückel coherent↔incoherent
            // threshold under J ↔ |β| (the cross-substrate identity). N=2 (Q*=1) is the EP base the
            // carbon polyene layer cannot reach. The mode that coalesces at Q*(N) is the {0,2}-coherence
            // (population/antisymmetric block) at ALL N=2..5, a genuine √-EP (phase rigidity r→0; corrected
            // 2026-06-13); the band edge 2cos(π/(N+1)) is the co-located γ-protected survivor, not the freezer
            // (no bifurcation at N=4). Closed form OPEN = the {0,2}-block discriminant. Live: inspect --root
            // horizon. Two typed parents,
            // both registered in this chain: ClockHandLadderClaim (the clock's Q-floor made exact,
            // directly above) + F2bXyChainSpectrumPi2Inheritance (the band edge, registered earlier in this chain).
            .RegisterCoherenceHorizonClaim()
            // The handover Q = the F50-floor condition (chain = Q*(N), ring = a distinct (2,2) level
            // crossing). Parents AbsorptionTheoremClaim (the -2g<n_XY> rate, registered above), F50WeightOne-
            // DegeneracyPi2Inheritance (the floor =1, above), CoherenceHorizonClaim (the chain solution, directly
            // above). Live witness: IncompletenessSurvivorWitness handover node (inspect --root survivor).
            .RegisterHandoverFloorClaim()
            // The structural ceiling (Tier1Derived, 2026-06-16, F122): the high-Q gap-rate closed forms
            // g2(K_N)=4/N, g2(star_N)=4/(N−1), g2(K_4)=2−2/√3 — the darkest [H,A]=0 coherence in the largest
            // degenerate single-particle level (the Re-side ceiling companion of TopologyBandEdgeClaim's J·ρ
            // band edge). Single typed parent AbsorptionTheoremClaim (the g2=⟨n_XY⟩ floor, registered above);
            // NOT the Tier1Candidate band-edge claim (the forms are dimensionless, depend only on the
            // Absorption Theorem + commutant linear algebra). Live witness: StructuralCeilingWitness
            // (inspect --root ceiling).
            .RegisterStructuralCeilingClaim()
            // The second clock's regime map (Tier1Candidate, 2026-06-16, "the stitch"): the {0,2}/half-filling
            // coherence (the second clock) is ONE mode whose regime is selected by the single-particle band via
            // two knobs — knob 1 (degeneracy m) → the high-Q ceiling 4/(m+1) (below the −2γ floor iff m≥4), knob 2
            // (dispersion) → the low-Q character (sharp √-EP coherence horizon on a dispersive band, only
            // asymptotic protection on a flat one). It stitches CoherenceHorizonClaim (the EP regime) +
            // StructuralCeilingClaim (the CEILING regime, 4/(m+1) the bridge) + the star-no-horizon (GRADUAL) into
            // one node. Two typed parents = the two regimes, both registered above: CoherenceHorizonClaim (the
            // weaker, Tier1Candidate, caps this child) + StructuralCeilingClaim (Tier1Derived). Gate-verified 15/15
            // (simulations/second_clock_regime_axis.py); live N=4 full-Liouvillian gate (inspect --root secondclock).
            .RegisterSecondClockRegimeClaim()
            .RegisterF86HwhmClosedFormClaim()
            .RegisterIbmBlockCpsiHardwareTable()
            .RegisterPolarityPairQPeakDecompositionClaim()
            // Cubic3-axis (Stage 2b): the 8-cell Z₂³ decomposition (bit_a, bit_b, y_par).
            // First Cubic3-axis Claim; PolarityCubeMap's Cubic3Claims grows from 0 to 1.
            // Structural anchor for F87 Z₂³ refinement work (F103/F105/F106). Must be
            // registered BEFORE every YParity-axis Claim below because Welle 7 wired
            // those Claims to take KleinEightCellClaim as a typed ctor parent (2026-05-26),
            // formalising the Klein2 → Cubic3 → YParity inheritance chain.
            .RegisterKleinEightCellClaim()
            // The mirror group D₄ (2026-06-10): Π_Z = R·D (transpose first, ket reflection
            // by X^⊗N second); ⟨R, D⟩ ≅ D₄ of order 8 closes the mirror inventory, with the
            // windowed-converse spine V₄ {I, 𝓕, R, 𝓕R} as Klein subgroup; the palindrome
            // splits along the generators (D flips L_H, R carries −2Σγ); the polarity cube
            // axes are characters (bit_a/bit_b = Ad_{Z^⊗N}/Ad_{X^⊗N}, y_par = the transpose
            // θ, the antiautomorphism axis). Typed parents: KleinEightCellClaim (registered
            // immediately above), CommutatorDConjugationSign (F114) and
            // Pi2KleinV4DephaseSwapGroup (Welle 12's D), both registered earlier.
            .RegisterMirrorGroupD4Claim()
            // Three dephasing diagonals as one orbit (2026-06-14): {Q_X,Q_Y,Q_Z} is one orbit of the
            // single-qubit Clifford basis-change S3 <h_zx,h_yz>; the three readings are the mirror D4
            // within a diagonal (the structure is S3 |x| D4). Its TWO typed parents are the physics edge
            // welding the mirror-group and absorption clusters (previously joined only at d^2-2d=0):
            // MirrorGroupD4Claim (directly above) + AbsorptionTheoremClaim (registered ~line 174). Anchor
            // simulations/one_diagonal_mirror_group.py (the physics-first gate that corrected the first
            // hypothesis: the Y-transpose, and the permuter is the basis-S3, not <R,D,h>).
            .RegisterThreeDephasingDiagonalsOrbitClaim()
            // YParity-axis seed (F102): standalone Tier1Derived Claim filling the cubic
            // Z₂³ architecture's YParity slot. Must come AFTER every Pi²-Inheritance Claim
            // above but BEFORE RegisterPolarityCubeMap so the cube map's b.Get<...>()
            // dependency resolves.
            .RegisterYParityIndependenceAtK3()
            // YParity-axis F103: Z₂³ refinement of F87 trichotomy at N=4 k=3.
            // First concrete derived class of F87Z2CubedRefinementBase; PolarityCubeMap's
            // YParityClaims grows from 1 to 2.
            .RegisterF87Z2CubedRefinementN4K3()
            // YParity-axis F105: Z₂³ refinement of F87 trichotomy at N=5 k=3.
            // Second concrete derived class of F87Z2CubedRefinementBase; PolarityCubeMap's
            // YParityClaims grows from 2 to 3. Frozen counts captured via SLOW_F105_BATCH
            // tool run; F85-predicted to match F103 bit-exactly (N-stability).
            .RegisterF87Z2CubedRefinementN5K3()
            // YParity-axis F106: Z₂³ refinement of F87 trichotomy at N=4 k=4.
            // Third concrete derived class of F87Z2CubedRefinementBase; PolarityCubeMap's
            // YParityClaims grows from 3 to 4. Frozen counts captured via SLOW_F106_BATCH
            // tool run; tests k-stability where F105 confirmed N-stability.
            .RegisterF87Z2CubedRefinementN4K4()
            // YParity-axis F107: F87 truly classification forces y_par = 0 (closed-form,
            // all dephase letters). First DERIVED-not-EMPIRICAL Claim in the family;
            // direct corollary of F85's k-body truly criterion. PolarityCubeMap's
            // YParityClaims grows from 4 to 5.
            .RegisterTrulyYParityZeroPurity()
            // YParity-axis F109: mother sector Klein (0,0) soft is y_par=1 pure across
            // all dephase letters. Sister to F107; together pin truly + mother-soft
            // y_par signature. Fully unconditional Tier1Derived after F108 Part 1
            // closure (Π_5bilinear, see RegisterF108Part1Pi2EvenAlwaysPalindromic below).
            // PolarityCubeMap's YParityClaims grows from 5 to 6.
            .RegisterMotherSoftYParityOnePurity()
            // YParity-axis F110: F87-hard pairs only in diagonal Klein cells with
            // Y-inversion. Aspect A closed-form via F108 Part 1+2+3 + F87 dissipator-
            // resonance; Aspects B+C derived via the F103 §6 counting rule + §7 mechanism
            // + the closed windowed converse. Tier1Derived (promoted 2026-06-10, gate
            // WindowedConverseAllGammaClaim). PolarityCubeMap's
            // YParityClaims grows from 6 to 7.
            .RegisterHardCellYInversionPattern()
            // YParity-axis F111: Pure-D Template Rule sharpening F110 Aspect B at
            // k = N = 4 in the diagonal Klein cell. Per-pair structural criterion:
            // pair is F87-hard ⟺ at least one term is a pure-D template (uses only
            // dephase letter D and identity I). Implies F110 Aspect B 228:0 Y-
            // inversion as immediate corollary. Tier1Derived (promoted 2026-06-10 via
            // WindowedConverseAllGammaClaim, Pascal-Gram positivity; Mixed+Mixed = soft
            // closed modulo M via PROOF_F103 §7.4). PolarityCubeMap's YParityClaims
            // grows from 7 to 8.
            .RegisterHardCellPureDTemplate()
            // BitA-axis F108 Part 2: Π²_X-even H + X-dephasing admits exact
            // operator-level palindrome (X-deph variant of Π_5bilinear). BitA twin
            // of F108 Part 1; must be registered BEFORE Part 1 (ctor parent).
            // PolarityCubeMap's BitAClaims grows from 1 (F61) to 2.
            .RegisterF108Part2Pi2XEvenAlwaysPalindromic()
            // BitB-axis F108 Part 1: Π²-even H + Z-dephasing admits exact
            // operator-level palindrome (Π_5bilinear). Closes F109's Step 5 for Z-
            // dephasing; together with F108 Part 2 (X-deph) + F108 Part 3 (Y-deph)
            // covers all three dephase letters. PolarityCubeMap's BitBClaims grows
            // by 1.
            .RegisterF108Part1Pi2EvenAlwaysPalindromic()
            // BitB-axis F108 Part 3: Π²_Y-even H + Y-dephasing admits exact
            // operator-level palindrome (Y-deph variant of Π_5bilinear). Y-deph
            // sibling of F108 Part 1 (same bilinear set, different dephase letter);
            // closes F109's Step 5 Y-dephasing branch and promotes F109 to fully
            // unconditional across all three dephase letters. PolarityCubeMap's
            // BitBClaims grows by 1.
            .RegisterF108Part3Pi2YEvenAlwaysPalindromic()
            // BitA-axis F112-X: cross-dephase sibling of F112-Z. Same Π-eigenvalue
            // balance identity with axis_d := bit_a substituted for bit_b: any H +
            // bit_a-homogeneous c_k under X-dephase Π_X polarity gives
            // ‖M_plus_half‖² = ‖M_minus_half‖² bit-exactly. Tier1Derived universal N
            // for both Hermitian and non-Hermitian H via Welle 13 (Route 1 direct
            // axis re-run + Route 2 Hadamard transport from F112-Z via Q_zx). Typed
            // ctor parent: F108 Part 2 (shared bit_a + X-dephase foundation). Welle 15
            // (2026-05-27) wired F112-X as the typed BitA twin of F112-Z, so F112-X
            // MUST be registered BEFORE F112-Z (the F112-Z ctor takes F112-X as the
            // BitATwinClaim parameter). PolarityCubeMap's BitAClaims grows by 1.
            .RegisterLindbladBitAPiBalance()
            // BitB-axis F112: Lindblad Π-eigenvalue balance under bit_b-homogeneous c.
            // Structural identity behind the polarity_coordinates_from_L diagnostic;
            // for any standard Lindblad system (Hermitian H + bit_b-homogeneous c_k)
            // the asymmetry ‖M_plus_half‖² − ‖M_minus_half‖² is exactly 0. Typed ctor
            // parents: F108 Part 1 (shared bit_b axis foundation via F38 / F63 Π²
            // eigenvalue formula on Pauli strings) + LindbladBitAPiBalance (typed BitA
            // twin, Welle 15). Sits on the BitB axis as the "polarity_coordinates
            // diagnostic identity" Claim. Tier1Derived universal N for both Hermitian
            // and non-Hermitian H. PolarityCubeMap's BitBClaims grows by 1; BitA-twin
            // status flipped from BitBSpecific (pre-Welle-15) to Filled.
            .RegisterLindbladBitBPiBalance()
            // BitB-axis F112-Y: Y-dephase sibling of F112-Z on the same bit_b axis
            // (Π_Y² and Π_Z² both grade by bit_b per F38). Same hypothesis on c
            // (bit_b-homogeneous), different polarity axis (Π_Y vs Π_Z). Tier1Derived
            // universal N via Welle 13 Route 1 (axis-direct re-run with d = Y);
            // D-conjugation from F112-Z is NOT available (D lacks Hilbert-space lift,
            // per PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md section (d) Remark). Typed
            // ctor parent: F108 Part 3 (shared bit_b + Y-dephase foundation).
            // PolarityCubeMap's BitBClaims grows by 1.
            .RegisterLindbladBitBPiYBalance()
            // BitB-axis F113: closed-form magnitude for the F112 polarity-asymmetry
            // counterexample. Sister to F112 on the same bit_b axis: F112 says
            // "in-scope → asymmetry = 0", F113 gives the exact magnitude
            // (4^N / 2) · Σ_l ω_l · (γ_T1,l − γ_pump,l) for the out-of-scope
            // Z-drive × amplitude-damping interference regime. Typed ctor parent:
            // F112 (LindbladBitBPiBalance). Tier1Derived at N=2, 3, 4 (bit-exact
            // via constructive parameter sweep); Tier1Candidate general N (universal-
            // N algebraic derivation of the (1/2)·4^N coefficient open). BitBSpecific
            // BitATwin slot (intrinsically about Z-axis single-site drives, no
            // meaningful bit_a analog). PolarityCubeMap's BitBClaims grows by 1.
            .RegisterLindbladBitBPiBreakMagnitude()
            // The antilinear triangle (2026-06-11): θ (transpose), conj (entrywise
            // conjugation), † (adjoint) close with id into one Klein four-group graded by
            // linearity ℓ and multiplicativity m; the transport law
            // μ∘L_H∘μ = ℓ(μ)·m(μ)·L_{μ(H)} is the one engine behind five existing proofs
            // (F114 θ-leg, the girth-ladder reversal kill, F112 Lemmas A+B †-leg, F113/F117
            // Hermitian conjugacy, the K_b mode mirror conj-leg); in the Pauli basis θ = D
            // and † = 𝒦, and ⟨R, D, 𝒦⟩ ≅ D₄ × Z₂ (order 16), the antilinear double of the
            // mirror group. NOT an IZ2AxisClaim (cross-axis structural, like
            // MirrorGroupD4Claim; cube-map counts unchanged). Typed parents:
            // MirrorGroupD4Claim, CommutatorDConjugationSign (F114), and
            // LindbladBitBPiBalance (F112), all registered above; the conj-leg claim
            // (ChiralMirrorTrajectoryClaim, Diagnostics) stays a prose edge.
            .RegisterAntilinearTriangleClaim()
            // The F120 moment-tower pump channel (2026-06-11): amplitude damping is the
            // unique non-unital piece of the standard Lindbladian and pumps along pure
            // local Z (D[σ⁻_l](I) = +Z_l, D[σ⁺_l](I) = −Z_l), so the slope of ⟨H^j⟩ at
            // the maximally mixed state reads the girth-ladder tower t_j(l) = Tr(Z_l·H^j)
            // linearly: slope_j = (1/d)·Σ_l Δγ_l·t_j(l). Dephasing-blind, evolution-blind,
            // closed at detailed balance; rung 1 is F113 (asymmetry = −4^N·slope⟨H⟩
            // exactly); the curvature is exactly affine and fingerprints X/Y-flavored
            // parasites while Z-flavored ones stay invisible (complementary to F113's
            // Z-drive reader); the girth certificate is one-sided. NOT an IZ2AxisClaim
            // (cross-axis structural; cube-map counts unchanged). Typed parents:
            // LindbladBitBPiBreakMagnitude (F113, the bridge) and
            // F84ThermalAmplitudeDampingPi2Inheritance (the Δγ vacuum-rate pump weight),
            // both registered above; the Diagnostics girth-ladder primitive stays a
            // prose edge (compute primitive, not a Claim).
            .RegisterMomentTowerPumpChannelClaim()
            // The three-ladder hinge (Tier1Derived, 2026-06-15): girth(ℓ)/rung(k)/moment(j) are the two
            // factors of one F87-hardness coefficient on M = A + γQ, hinged by Q (spectrum = the rung,
            // action = the girth-moment projector). Parents AbsorptionTheoremClaim (the rung) +
            // MomentTowerPumpChannelClaim (the moments), both registered above. Live: inspect --root ladders.
            .RegisterThreeLadderHingeClaim()
            // Zero-Sector Immunity (Tier1Derived, PROOF_ZERO_IMMUNITY): the palindromic residual M
            // vanishes on the extreme weight blocks (w=0 {I,Z}^⊗N and w=N {X,Y}^⊗N) for every 2-body H
            // + Z-dephasing. Parents F1PalindromeIdentity + F61 + F63 + AbsorptionTheoremClaim, all
            // registered above. Live: inspect --root zeroimmune.
            .RegisterZeroSectorImmunityClaim()
            // Polarity cube map (cubic Z₂³ architecture inventory; aggregates every IZ2AxisClaim)
            // Must come AFTER every Pi²-Inheritance Claim registration above so its b.Get<T>()
            // dependencies resolve. Currently snapshot-only; rebuild registry to refresh.
            .RegisterPolarityCubeMap()
            // Open questions
            .RegisterF1OpenQuestions()
            .RegisterF86OpenQuestions()
            .Build();
    }
}
