using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F57DwellTimeQuarterPi2Inheritance"/>:
/// F57's <c>t_dwell(δ) = 2δ / |dCΨ/dt|_{t_cross}</c> closed form's two
/// Pi2-Foundation anchors (CΨ = 1/4 boundary, 2δ-window factor) as inheritance.
/// Two parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_3 = 1/4</c> (the
///         CΨ crossing threshold) and <c>a_0 = 2</c> (the 2δ-window factor =
///         polynomial root d).</item>
///   <item><see cref="QuarterAsBilinearMaxvalClaim"/>: the bilinear-apex maxval
///         anchor identifying 1/4 as max p·(1−p) at p=1/2.</item>
/// </list>
///
/// <para>Tier consistency: F57 is Tier 1 analytical (CRITICAL_SLOWING_AT_THE_CUSP),
/// hardware-verified ibm_kingston 2026-04-16. The Pi2-Foundation anchoring is
/// algebraic-trivial composition. All three claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> in the
/// builder pipeline. <see cref="QuarterAsBilinearMaxvalClaim"/> is registered
/// by <see cref="Pi2FamilyRegistration"/>.</para></summary>
public static class F57DwellTimeQuarterPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF57DwellTimeQuarterPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F57DwellTimeQuarterPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var quarter = b.Get<QuarterAsBilinearMaxvalClaim>();
            return new F57DwellTimeQuarterPi2Inheritance(ladder, quarter);
        });
}
