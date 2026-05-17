using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F97CardioidHalfFixedPointPi2Inheritance"/>:
/// F97 extends F95 from the real-c angle to the full complex-c plane via the
/// Mandelbrot cardioid parametrization c(φ) = b·e^(iφ) − b²·e^(2iφ). The cardioid
/// is the structural curve in complex-c where the period-1 fixed-point magnitude
/// of z² + c equals b = 1/2 (HalfAsStructuralFixedPoint) invariant.
///
/// <list type="bullet">
///   <item><see cref="HalfAsStructuralFixedPointClaim"/>: b = 1/2 is the
///         fixed-point magnitude invariant around the entire cardioid.</item>
///   <item><see cref="QuarterAsBilinearMaxvalClaim"/>: b² = 1/4 is the
///         magnitude of c at the real-axis cusp only (the F95 tangent point).</item>
///   <item><see cref="NinetyDegreeMirrorMemoryClaim"/>: the i in the e^(iφ)
///         parametrization, lifting c from the real axis to the complex plane.</item>
///   <item><see cref="PolynomialFoundationClaim"/>: the c = 0 (origin) special
///         case where the fixed point is z* = 0 (degenerate).</item>
/// </list>
///
/// <para>Tier consistency: F97 Tier 1 derived (bit-exact algebraic identity).
/// All four parents Tier 1 derived. Tier composition well-formed.</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (registers HalfAsStructuralFixedPointClaim, QuarterAsBilinearMaxvalClaim,
/// NinetyDegreeMirrorMemoryClaim, PolynomialFoundationClaim).</para></summary>
public static class F97CardioidHalfFixedPointPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF97CardioidHalfFixedPointPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F97CardioidHalfFixedPointPi2Inheritance>(b =>
        {
            var half = b.Get<HalfAsStructuralFixedPointClaim>();
            var quarter = b.Get<QuarterAsBilinearMaxvalClaim>();
            var ninety = b.Get<NinetyDegreeMirrorMemoryClaim>();
            var polynomial = b.Get<PolynomialFoundationClaim>();
            return new F97CardioidHalfFixedPointPi2Inheritance(
                half, quarter, ninety, polynomial);
        });
}
