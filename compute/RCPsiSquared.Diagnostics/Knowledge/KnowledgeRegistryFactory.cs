using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Diagnostics.F87;
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
            .RegisterF86JordanWignerLight(N: defaultChain.N, n: 1, gammaZero: defaultChain.GammaZero)
            // F87 family + canonical witnesses
            // RegisterF87Family must follow RegisterPi2Family: DissipatorAxisSelectsPolarityClaim
            // resolves PolarityLayerOriginClaim via b.Get<>() at registration-factory time.
            .RegisterF87Family()
            .RegisterF87StandardWitnessSet()
            .RegisterF89F87TrulyInheritance()
            .RegisterF89F87BreakPredictionFromF83()
            // Spectrum quantization root (parent to F33/F50/F55/F64-F68/F74/F89
            // via per-Registration discard-Get edges; absorption quantum 2γ₀ from a_0)
            .RegisterAbsorptionTheoremClaim()
            .RegisterF86LEffMirrorAxis()
            // F-formula Pi2-Foundation inheritance claims
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
            .RegisterF1T1AmplitudeDampingPi2Inheritance()
            .RegisterF5DepolarizingErrorPi2Inheritance()
            .RegisterDickeSuperpositionQuarterPi2Inheritance()
            .RegisterF39DetPiPi2Inheritance()
            .RegisterF49bCenteredDissipatorPi2Inheritance()
            .RegisterF38Pi2InvolutionPi2Inheritance()
            .RegisterF49cShadowCrossingPi2Inheritance()
            .RegisterF25CPsiBellPlusPi2Inheritance()
            .RegisterF57DwellTimeQuarterPi2Inheritance()
            .RegisterF66PoleModesPi2Inheritance()
            .RegisterF63LCommutesPi2Pi2Inheritance()
            .RegisterF61BitAParityPi2Inheritance()
            .RegisterF77MmSaturationPi2Inheritance()
            .RegisterF60GhzBornBelowFoldPi2Inheritance()
            .RegisterF83AntiFractionPi2Inheritance()
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
            .RegisterF86HwhmClosedFormClaim()
            .RegisterIbmBlockCpsiHardwareTable()
            .RegisterPolarityPairQPeakDecompositionClaim()
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
            // y_par signature. Closed-form modulo F108 Part 1 (Π²-even-soft).
            // PolarityCubeMap's YParityClaims grows from 5 to 6.
            .RegisterMotherSoftYParityOnePurity()
            // Cubic3-axis (Stage 2b): the 8-cell Z₂³ decomposition (bit_a, bit_b, y_par).
            // First Cubic3-axis Claim; PolarityCubeMap's Cubic3Claims grows from 0 to 1.
            // Structural anchor for F87 Z₂³ refinement work (F103/F105/F106).
            .RegisterKleinEightCellClaim()
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
