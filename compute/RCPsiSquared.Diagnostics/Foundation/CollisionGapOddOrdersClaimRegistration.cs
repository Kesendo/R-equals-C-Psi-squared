using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="CollisionGapOddOrdersClaim"/> (F161; Tier1Derived). ONE typed parent,
/// <see cref="CrackedRingExactCurveClaim"/> (F160, Tier1Derived), whose polynomial is the whole input to the series
/// and whose Theorem G is that series at first order, so the tier rule holds. F129 is the other thing this law stands
/// on and is deliberately NOT a typed edge: its typed home <see cref="CrossTripleOrthogonalityClaim"/> is
/// Tier1Candidate (the code-trust caveat named in that class), and a Tier1Candidate parent would break the tier rule
/// for a Tier1Derived child; F129 is carried by anchor instead, and executably through
/// <c>LevelCollisionCensus.Fires</c> and <c>CollisionFamilyInventory.Count</c>, which the claim calls rather than
/// restates. F2b and F65 are already ancestors through the F160 edge and are not re-added. Resolution is topological,
/// so this registration may sit anywhere after the family that supplies the parent.</summary>
public static class CollisionGapOddOrdersClaimRegistration
{
    public static ClaimRegistryBuilder RegisterCollisionGapOddOrdersClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<CollisionGapOddOrdersClaim>(b =>
            new CollisionGapOddOrdersClaim(b.Get<CrackedRingExactCurveClaim>()));
}
