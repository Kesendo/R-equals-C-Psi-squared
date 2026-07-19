using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="StructuralCeilingClaim"/>: a single typed parent edge to the Tier1Derived
/// <see cref="AbsorptionTheoremClaim"/> (the g2 = ⟨n_XY⟩ floor reading). Deliberately NOT parented on
/// <see cref="TopologyBandEdgeClaim"/> (the same arc's Im side) nor HandoverFloorClaim: the
/// closed forms are dimensionless and depend only on the Absorption Theorem plus self-contained commutant
/// linear algebra, so the claim is genuinely Tier1Derived on that single parent alone.</summary>
public static class StructuralCeilingClaimRegistration
{
    public static ClaimRegistryBuilder RegisterStructuralCeilingClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<StructuralCeilingClaim>(b =>
            new StructuralCeilingClaim(b.Get<AbsorptionTheoremClaim>()));
}
