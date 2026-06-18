using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class StarFrozenSeamClaimTests
{
    [Fact]
    public void Build_WiresBothParents_AndIsCappedAtTier1Candidate()
    {
        var claim = StarFrozenSeamClaim.Build();
        // Tier1Candidate: the tier-inheritance invariant caps the child at its weaker parent,
        // SecondClockRegimeClaim (Tier1Candidate); the other parent, StructuralCeilingClaim, is Tier1Derived.
        Assert.Equal(Tier.Tier1Candidate, claim.Tier);
        Assert.NotNull(claim.Ceiling);
        Assert.NotNull(claim.Regime);
        Assert.Equal(Tier.Tier1Derived, claim.Ceiling.Tier);
        Assert.Equal(Tier.Tier1Candidate, claim.Regime.Tier);

        var kids = ((IInspectable)claim).Children.ToList();
        Assert.Contains(claim.Ceiling, kids);
        Assert.Contains(claim.Regime, kids);
    }
}
