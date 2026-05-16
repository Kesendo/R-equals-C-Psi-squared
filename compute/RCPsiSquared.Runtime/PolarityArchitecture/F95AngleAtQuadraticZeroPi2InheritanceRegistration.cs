using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F95AngleAtQuadraticZeroPi2Inheritance"/>:
/// F95's universal closed form θ(c; b) = arctan(√(c/b² − 1)) brings four typed
/// Pi2-Foundation anchors together via a 4-line polynomial derivation.
///
/// <list type="bullet">
///   <item><see cref="PolynomialFoundationClaim"/>: c = 0 special case is
///         d² − 2d = 0; F95 perturbs c off zero and tracks the angle.</item>
///   <item><see cref="HalfAsStructuralFixedPointClaim"/>: framework b = 1/2.</item>
///   <item><see cref="QuarterAsBilinearMaxvalClaim"/>: framework threshold
///         b² = 1/4.</item>
///   <item><see cref="NinetyDegreeMirrorMemoryClaim"/>: the i in
///         z = b ± i·√(c − b²); the 90°-rotation generator of the angle.</item>
/// </list>
///
/// <para>Tier consistency: F95 Tier 1 derived (bit-exact symbolic). All four
/// parents Tier 1 derived. Tier composition well-formed.</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (which registers PolynomialFoundationClaim, HalfAsStructuralFixedPointClaim,
/// QuarterAsBilinearMaxvalClaim, NinetyDegreeMirrorMemoryClaim).</para></summary>
public static class F95AngleAtQuadraticZeroPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF95AngleAtQuadraticZeroPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F95AngleAtQuadraticZeroPi2Inheritance>(b =>
        {
            var polynomial = b.Get<PolynomialFoundationClaim>();
            var half = b.Get<HalfAsStructuralFixedPointClaim>();
            var quarter = b.Get<QuarterAsBilinearMaxvalClaim>();
            var ninetyDegree = b.Get<NinetyDegreeMirrorMemoryClaim>();
            return new F95AngleAtQuadraticZeroPi2Inheritance(
                polynomial, half, quarter, ninetyDegree);
        });
}
