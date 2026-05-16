using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F60GhzBornBelowFoldPi2Inheritance"/>:
/// F60's <c>CΨ(0) = 1/(2^N − 1)</c> closed form's four Pi2-Foundation anchors,
/// all typed ctor parents as of 2026-05-16:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_2 = 1/2</c>
///         (off-diagonal element), <c>a_{1−N}</c> (Hilbert-space dimension),
///         <c>a_3 = 1/4</c> (fold position) via live <c>Term(...)</c> reads.</item>
///   <item><see cref="PolarityLayerOriginClaim"/>: typed parent grounding the
///         "1/2 off-diagonal = ±0.5 polarity pair at d=2" reading. F60 is the
///         first F-formula whose primary anchor sits directly on the 0.5-shift
///         axis (per Tom 2026-05-09 mirror-map check).</item>
///   <item><see cref="QuarterAsBilinearMaxvalClaim"/>: typed parent grounding
///         the "fold = 1/4 = (1/2)² maxval of p·(1−p)" reading. Same anchor as
///         F57 + Dicke superposition.</item>
///   <item><see cref="ArgmaxMaxvalPairClaim"/>: typed meta-parent closing the
///         (1/2, 1/4) pair. F60 uses BOTH: OffDiagonalElement = 1/2 and
///         FoldPosition = 1/4.</item>
/// </list>
///
/// <para>Tier consistency: F60 is Tier 1 geometric corollary; Pi2-Foundation
/// anchoring is algebraic-trivial composition. All four claims Tier1Derived
/// (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (registers PolarityLayerOriginClaim + QuarterAsBilinearMaxvalClaim +
/// ArgmaxMaxvalPairClaim) +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>.</para></summary>
public static class F60GhzBornBelowFoldPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF60GhzBornBelowFoldPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F60GhzBornBelowFoldPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var polarity = b.Get<PolarityLayerOriginClaim>();
            var quarter = b.Get<QuarterAsBilinearMaxvalClaim>();
            var argmaxMaxval = b.Get<ArgmaxMaxvalPairClaim>();
            return new F60GhzBornBelowFoldPi2Inheritance(ladder, polarity, quarter, argmaxMaxval);
        });
}
