using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Core.SymmetryFamily;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="PolarityCubeMap"/>. Aggregates all
/// IZ2AxisClaim instances visible from the runtime (Schicht-1 registered
/// Pi²-Inheritance Claims and Klein2 Claims) into a single inspectable map.
///
/// <para>The aggregation is performed at registration time; the dependency on each
/// IZ2AxisClaim is declared explicitly via <c>b.Get&lt;T&gt;()</c> so the inheritance
/// graph records PolarityCubeMap as a downstream consumer of every Pi²-Inheritance
/// Claim. The cube map is a snapshot of the typed knowledge graph's Z₂-axis coverage
/// at build time; subsequent typed-knowledge changes require a rebuild to reflect.</para>
///
/// <para>F87Pi2Inheritance lives in RCPsiSquared.Diagnostics and is not depended on
/// here (Runtime does not reference Diagnostics); its absence in the snapshot is
/// intentional and the count reflects the Runtime-visible scope.</para>
///
/// <para>Requires upstream registration of every Pi²-Inheritance Claim referenced
/// below (typically done by <c>KnowledgeRegistryFactory.BuildDefault</c> which
/// composes the registration chain).</para></summary>
public static class PolarityCubeMapRegistration
{
    public static ClaimRegistryBuilder RegisterPolarityCubeMap(
        this ClaimRegistryBuilder builder) =>
        builder.Register<PolarityCubeMap>(b =>
        {
            var z2AxisClaims = new List<IZ2AxisClaim>();
            z2AxisClaims.Add(b.Get<CanonicalTrigAnchorPi2Inheritance>());
            z2AxisClaims.Add(b.Get<DickeSuperpositionQuarterPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F1DepolResidualClosedFormPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F1Pi2Inheritance>());
            z2AxisClaims.Add(b.Get<F1T1AmplitudeDampingPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F23XorDrainVanishingFractionPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F25CPsiBellPlusPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F26CPsiPauliChannelsPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F2W1DispersionPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F2bXyChainSpectrumPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F33ExactN3DecayRatesPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F38Pi2InvolutionPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F39DetPiPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F3DecayRateBoundsPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F41PalindromicTimePi2Inheritance>());
            z2AxisClaims.Add(b.Get<F43SectorSffPairingPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F44CrooksLikeRateIdentityPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F49Pi2Inheritance>());
            z2AxisClaims.Add(b.Get<F49bCenteredDissipatorPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F49cShadowCrossingPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F4StationaryModeCountPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F50WeightOneDegeneracyPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F55UniversalAbsorptionDosePi2Inheritance>());
            z2AxisClaims.Add(b.Get<F56CriticalSlowingPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F57DwellTimeQuarterPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F5DepolarizingErrorPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F60GhzBornBelowFoldPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F61BitAParityPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F62WStateBornBelowFoldPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F63LCommutesPi2Pi2Inheritance>());
            z2AxisClaims.Add(b.Get<F64CavityModeExposurePi2Inheritance>());
            z2AxisClaims.Add(b.Get<F65XxChainSpectrumPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F66PoleModesPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F67BondingBellPairPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F68PalindromicPartnerPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F69GhzWSexticAboveFoldPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F70DeltaNSelectionRulePi2Inheritance>());
            z2AxisClaims.Add(b.Get<F71MirrorSymmetryPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F72BlockDiagonalPurityPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F73SpatialSumPurityClosurePi2Inheritance>());
            z2AxisClaims.Add(b.Get<F74ChromaticityPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F75MirrorPairMiPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F76TDecayMirrorPairMiPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F77MmSaturationPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F78SingleBodyMAdditivePi2Inheritance>());
            z2AxisClaims.Add(b.Get<F79TwoBodyPi2BlockPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F80FactorPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F81Pi2Inheritance>());
            z2AxisClaims.Add(b.Get<F82T1AmplitudeDampingPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F83AntiFractionPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F84ThermalAmplitudeDampingPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F85KBodyFChainPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F86F71MirrorPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F86QEpPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F86TPeakPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F91Pi2Inheritance>());
            z2AxisClaims.Add(b.Get<F92Pi2Inheritance>());
            z2AxisClaims.Add(b.Get<F93Pi2Inheritance>());
            z2AxisClaims.Add(b.Get<F94BornDeviationFourThirdsPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F95AngleAtQuadraticZeroPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F96BornSubdominantSlopesPi2Inheritance>());
            z2AxisClaims.Add(b.Get<F97CardioidHalfFixedPointPi2Inheritance>());
            z2AxisClaims.Add(b.Get<QubitNecessityPi2Inheritance>());
            z2AxisClaims.Add(b.Get<KleinFourCellClaim>());
            z2AxisClaims.Add(b.Get<F89F88aKleinPpAnchor>());
            z2AxisClaims.Add(b.Get<Pi2InvolutionClaim>());
            z2AxisClaims.Add(b.Get<Pi2KleinBilinearTable>());
            // YParity-axis Claims (F102+)
            z2AxisClaims.Add(b.Get<YParityIndependenceAtK3>());
            z2AxisClaims.Add(b.Get<F87Z2CubedRefinementN4K3>());
            z2AxisClaims.Add(b.Get<F87Z2CubedRefinementN5K3>());
            z2AxisClaims.Add(b.Get<F87Z2CubedRefinementN4K4>());
            z2AxisClaims.Add(b.Get<TrulyYParityZeroPurity>());
            z2AxisClaims.Add(b.Get<MotherSoftYParityOnePurity>());
            // Cubic3-axis Claims (Stage 2b+)
            z2AxisClaims.Add(b.Get<KleinEightCellClaim>());
            z2AxisClaims.Add(b.Get<ZGlobalMirrorRefinement>());
            return new PolarityCubeMap(z2AxisClaims);
        });
}
