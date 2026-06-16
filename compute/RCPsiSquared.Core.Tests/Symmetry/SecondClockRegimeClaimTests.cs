using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class SecondClockRegimeClaimTests
{
    [Fact]
    public void Build_WiresBothRegimeParents_AndIsCappedAtTier1Candidate()
    {
        var claim = SecondClockRegimeClaim.Build();
        // Tier1Candidate: the tier-inheritance invariant caps the child at its weaker parent,
        // CoherenceHorizonClaim (Tier1Candidate); the other parent, StructuralCeilingClaim, is Tier1Derived.
        Assert.Equal(Tier.Tier1Candidate, claim.Tier);
        Assert.NotNull(claim.Horizon);
        Assert.NotNull(claim.Ceiling);
        Assert.Equal(Tier.Tier1Candidate, claim.Horizon.Tier);
        Assert.Equal(Tier.Tier1Derived, claim.Ceiling.Tier);

        var kids = ((IInspectable)claim).Children.ToList();
        Assert.Contains(claim.Horizon, kids);
        Assert.Contains(claim.Ceiling, kids);
    }

    [Fact]
    public void SymmetricManifoldCeiling_Is4OverMPlus1_BelowFloorIffMAtLeast4()
    {
        // the bridge formula 4/(m+1): a structural ceiling (< 1) exactly when m ≥ 4
        Assert.Equal(4.0 / 3.0, SecondClockRegimeClaim.SymmetricManifoldCeiling(2), 12);  // > 1
        Assert.Equal(1.0, SecondClockRegimeClaim.SymmetricManifoldCeiling(3), 12);        // marginal (GRADUAL)
        Assert.Equal(4.0 / 5.0, SecondClockRegimeClaim.SymmetricManifoldCeiling(4), 12);  // < 1 (CEILING)
        Assert.True(SecondClockRegimeClaim.SymmetricManifoldCeiling(4) < 1.0);
        Assert.True(SecondClockRegimeClaim.SymmetricManifoldCeiling(3) >= 1.0 - 1e-12);
    }
}
