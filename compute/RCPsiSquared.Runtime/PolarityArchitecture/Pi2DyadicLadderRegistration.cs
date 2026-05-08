using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="Pi2DyadicLadderClaim"/> — the dyadic halving
/// ladder a_n = 2^(1−n) that organises the Pi2 foundation's three known typed anchors
/// (n=0 d=2 root, n=2 half fixed point, n=3 bilinear maxval) onto a single algebraic
/// continuation.
///
/// <para>Three edges declared, one per anchor on the ladder. The three parents are all
/// Tier1Derived; the ladder itself is Tier1Derived. TierStrength inheritance trivially
/// passes (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> to be called
/// first so all three anchor Claims are registered. <see cref="ClaimRegistryBuilder"/>
/// throws <c>MissingParent</c> for any missing parent.</para></summary>
public static class Pi2DyadicLadderRegistration
{
    public static ClaimRegistryBuilder RegisterPi2DyadicLadder(this ClaimRegistryBuilder builder) =>
        builder.Register<Pi2DyadicLadderClaim>(b =>
        {
            _ = b.Get<QubitDimensionalAnchorClaim>();        // n=0
            _ = b.Get<HalfAsStructuralFixedPointClaim>();    // n=2
            _ = b.Get<QuarterAsBilinearMaxvalClaim>();       // n=3
            return new Pi2DyadicLadderClaim();
        });
}
