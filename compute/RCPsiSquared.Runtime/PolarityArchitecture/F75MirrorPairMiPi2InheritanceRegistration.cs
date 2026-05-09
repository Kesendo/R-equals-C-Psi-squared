using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F75MirrorPairMiPi2Inheritance"/>:
/// F75's MI = 2·h(p) − h(2p) closed form. F75 is the algebraic mother of F77
/// (F77's 1-bit asymptote falls out of F75 by Taylor expansion at p → 0).
/// Three parent edges (one explicit, two as registration discards):
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>
///         (saturation = 2 bits) and <c>a_2 = 1/2</c> (domain upper bound =
///         BilinearApex argmax).</item>
///   <item><see cref="BilinearApexClaim"/>: registration discard. F75's
///         probability domain [0, 1/2] IS the BilinearApex extent; MI saturates
///         at the apex value p = 1/2 (Bell-state mirror-pair). Convex on
///         (0, 1/2) per the F75 structural reading.</item>
///   <item><see cref="F71MirrorSymmetryPi2Inheritance"/>: registration discard.
///         F75 cites F71 as source: "the mirror symmetry that justifies
///         c_{N−1−j} = ±c_j". F71 → F75 is a typed inheritance edge
///         (mother claim).</item>
/// </list>
///
/// <para>Tier consistency: F75 is Tier 1 proven algebraic; verified
/// numerically against C# brecher propagation N=5..13, ~25 data points.
/// All four claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (registers BilinearApexClaim) +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="HalfIntegerMirrorRegistration.RegisterHalfIntegerMirror"/>(N) +
/// <see cref="F71MirrorSymmetryPi2InheritanceRegistration.RegisterF71MirrorSymmetryPi2Inheritance"/>.</para></summary>
public static class F75MirrorPairMiPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF75MirrorPairMiPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F75MirrorPairMiPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            _ = b.Get<BilinearApexClaim>();                          // domain [0, 1/2]
            _ = b.Get<F71MirrorSymmetryPi2Inheritance>();            // mirror symmetry source
            return new F75MirrorPairMiPi2Inheritance(ladder);
        });
}
