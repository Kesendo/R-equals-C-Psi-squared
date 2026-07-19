using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class StructuralCeilingClaimTests
{
    [Fact]
    public void Build_WiresAbsorptionParent_AndIsTier1Derived()
    {
        var claim = StructuralCeilingClaim.Build();
        // Tier1Derived: a principal-angle proof of the closed forms + gate-exact verification. The single
        // typed parent is the Tier1Derived AbsorptionTheoremClaim; the claim is NOT parented on
        // TopologyBandEdgeClaim (the forms are dimensionless, no ρ, no gap-dominance dependence).
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.NotNull(claim.Absorption);
        Assert.Equal(Tier.Tier1Derived, claim.Absorption.Tier);
        var kids = ((IInspectable)claim).Children.ToList();
        Assert.Contains(claim.Absorption, kids);
    }

    [Fact]
    public void ClosedForms_MatchTheDerivation()
    {
        Assert.Equal(4.0 / 5.0, StructuralCeilingClaim.CompleteCeiling(5), 12);
        Assert.Equal(4.0 / 7.0, StructuralCeilingClaim.CompleteCeiling(7), 12);
        Assert.Equal(4.0 / 5.0, StructuralCeilingClaim.StarCeiling(6), 12);     // 4/(6−1)
        Assert.Equal(2.0 - 2.0 / System.Math.Sqrt(3.0), StructuralCeilingClaim.K4Ceiling, 12);
    }
}
