using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F56CriticalSlowingPi2Inheritance"/>:
/// F56 critical-slowing iteration count K(ε, tol). Three typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_2 = 1/2</c>
///         (HalfPrefactor), <c>a_{−1} = 4</c> (FourFactor),
///         <c>a_{−3} = 16 = 4²</c> (SixteenFactor), and <c>a_3 = 1/4</c>
///         (CardioidCuspPosition). F56 packs five distinct dyadic-ladder
///         positions in one closed form.</item>
///   <item><see cref="QuarterAsBilinearMaxvalClaim"/>: cardioid cusp at 1/4
///         IS the bilinear-apex maxval; F56's ε measures distance from this
///         maxval position.</item>
///   <item><see cref="HalfAsStructuralFixedPointClaim"/>: the "1/2" prefactor
///         on the leading log; the cardioid argmax is at p = 1/2 (where the
///         bilinear maxval is attained).</item>
/// </list>
///
/// <para>Tier consistency: F56 is Tier 1 with zero fit parameters; verified
/// 0.5–2% accuracy over 5 tol decades and 10 ε decades.</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (registers QuarterAsBilinearMaxvalClaim + HalfAsStructuralFixedPointClaim) +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>.</para></summary>
public static class F56CriticalSlowingPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF56CriticalSlowingPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F56CriticalSlowingPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var quarter = b.Get<QuarterAsBilinearMaxvalClaim>();
            var half = b.Get<HalfAsStructuralFixedPointClaim>();
            return new F56CriticalSlowingPi2Inheritance(ladder, quarter, half);
        });
}
