using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class StarFrozenSeamClaimTests
{
    [Fact]
    public void Build_WiresBothParents_AndStaysTier1CandidateOnOwnStanding()
    {
        var claim = StarFrozenSeamClaim.Build();
        // Tier1Candidate on its OWN standing since 2026-07-19 (no longer parent-capped: the former
        // weaker parent SecondClockRegimeClaim graduated to Tier1Derived): the all-Q survivor
        // statement is gate-verified at N=4..8, not proven for general N.
        Assert.Equal(Tier.Tier1Candidate, claim.Tier);
        Assert.NotNull(claim.Ceiling);
        Assert.NotNull(claim.Regime);
        Assert.Equal(Tier.Tier1Derived, claim.Ceiling.Tier);
        Assert.Equal(Tier.Tier1Derived, claim.Regime.Tier);

        var kids = ((IInspectable)claim).Children.ToList();
        Assert.Contains(claim.Ceiling, kids);
        Assert.Contains(claim.Regime, kids);
    }
}
