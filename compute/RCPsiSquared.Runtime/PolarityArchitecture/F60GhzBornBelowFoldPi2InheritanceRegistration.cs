using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F60GhzBornBelowFoldPi2Inheritance"/>:
/// F60's <c>CΨ(0) = 1/(2^N − 1)</c> closed form's three Pi2-Foundation anchors:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_2 = 1/2</c>
///         (off-diagonal element), <c>a_{1−N}</c> (Hilbert-space dimension),
///         <c>a_3 = 1/4</c> (fold position).</item>
///   <item><see cref="PolarityLayerOriginClaim"/>: registration discard
///         documenting that the "1/2 off-diagonal" IS the ±0.5 polarity pair
///         at d=2. F60 is the first F-formula whose primary anchor sits
///         directly on the 0.5-shift axis (per Tom 2026-05-09 mirror-map
///         check). The 0.5-shift Polarity-Layer-Inheritance is now typed.</item>
///   <item><see cref="QuarterAsBilinearMaxvalClaim"/>: registration discard
///         documenting that the "fold = 1/4" IS the bilinear-apex maxval.
///         Same anchor as F57 + Dicke superposition; F60 inherits the
///         framework's quarter-boundary as its threshold.</item>
/// </list>
///
/// <para>Tier consistency: F60 is Tier 1 geometric corollary; Pi2-Foundation
/// anchoring is algebraic-trivial composition. All four claims Tier1Derived
/// (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (registers PolarityLayerOriginClaim + QuarterAsBilinearMaxvalClaim) +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>.</para></summary>
public static class F60GhzBornBelowFoldPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF60GhzBornBelowFoldPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F60GhzBornBelowFoldPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            _ = b.Get<PolarityLayerOriginClaim>();         // direct: ±0.5 polarity pair
            _ = b.Get<QuarterAsBilinearMaxvalClaim>();     // fold = 1/4 = a_3
            return new F60GhzBornBelowFoldPi2Inheritance(ladder);
        });
}
