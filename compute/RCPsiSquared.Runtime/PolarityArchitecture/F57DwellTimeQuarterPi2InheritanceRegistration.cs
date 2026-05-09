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
///   <item><see cref="ArgmaxMaxvalPairClaim"/>: the meta-anchor closing the
///         (1/2, 1/4) pair. F57 uses BOTH ends, 1/4 (CrossingThreshold) and
///         2 = 1/(1/2) (WindowDoublingFactor); the meta-anchor types this
///         convergence (activated 2026-05-09 mirror-map check).</item>
///   <item><see cref="F25CPsiBellPlusPi2Inheritance"/>: registration discard.
///         F57's Bell+ prefactor 1.080088 = 2 / 1.851701 is derived from
///         F25's <c>|dCΨ/dt|_{t_cross}</c>. F25 → F57 mother-claim edge
///         (added 2026-05-09; pattern parallel to F75 → F77).</item>
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
            _ = b.Get<ArgmaxMaxvalPairClaim>();              // meta-anchor: F57 uses both 1/2 and 1/4 ends
            _ = b.Get<F25CPsiBellPlusPi2Inheritance>();      // mother claim: F57 prefactor = 2 / |dCΨ/dt|_{t_cross}
            return new F57DwellTimeQuarterPi2Inheritance(ladder, quarter);
        });
}
