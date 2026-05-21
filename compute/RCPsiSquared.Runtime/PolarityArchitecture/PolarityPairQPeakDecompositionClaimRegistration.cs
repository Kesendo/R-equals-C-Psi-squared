using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="PolarityPairQPeakDecompositionClaim"/>: the F86
/// schema-level claim that Q_peak decomposes as Q_EP_central ± r_polarity = 2 ± 1/2, so
/// Endpoint orbits sit at the schema value 2.5 and Interior orbits at 1.5. Three parent
/// edges:
///
/// <list type="bullet">
///   <item><see cref="HalfAsStructuralFixedPointClaim"/>: the r_polarity = 1/2 magnitude.</item>
///   <item><see cref="PolarityLayerOriginClaim"/>: the ± pole signs of the polarity pair.</item>
///   <item><see cref="QubitDimensionalAnchorClaim"/>: the d = 2 central anchor.</item>
/// </list>
///
/// <para>Tier consistency: all four Tier1Derived.</para>
///
/// <para>Requires upstream registration: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (registers all three parents).</para></summary>
public static class PolarityPairQPeakDecompositionClaimRegistration
{
    public static ClaimRegistryBuilder RegisterPolarityPairQPeakDecompositionClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<PolarityPairQPeakDecompositionClaim>(b =>
        {
            var half = b.Get<HalfAsStructuralFixedPointClaim>();
            var polarityOrigin = b.Get<PolarityLayerOriginClaim>();
            var dimensionalAnchor = b.Get<QubitDimensionalAnchorClaim>();
            return new PolarityPairQPeakDecompositionClaim(half, polarityOrigin, dimensionalAnchor);
        });
}
