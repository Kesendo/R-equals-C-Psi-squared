using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class CpsiEnvelopeTheoremClaimTests
{
    [Fact]
    public void Claim_TierIsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, CpsiEnvelopeTheoremClaim.Shared.Tier);
    }

    [Fact]
    public void Claim_StatementCarriesTheScopeSplitAndTheBoundary()
    {
        var name = CpsiEnvelopeTheoremClaim.Shared.Name;
        Assert.Contains("non-increasing", name);
        Assert.Contains("Proven for 2-qubit", name);
        Assert.Contains("verified N=3-5", name);
        Assert.Contains("reduced", name);
        Assert.Contains("freedom", name);
    }

    [Fact]
    public void Claim_HasBothTypedParents()
    {
        var kids = ((IInspectable)CpsiEnvelopeTheoremClaim.Shared).Children.ToList();
        Assert.Contains(kids, k => k is F25CPsiBellPlusPi2Inheritance);
        Assert.Contains(kids, k => k is QuarterAsBilinearMaxvalClaim);
    }

    [Fact]
    public void Build_SharesOneQuarterInstance_AcrossF25AndTheDirectEdge()
    {
        var c = CpsiEnvelopeTheoremClaim.Shared;
        // Build() threads ONE QuarterAsBilinearMaxvalClaim into both F25 and the ¼-boundary edge.
        Assert.Same(c.Quarter, c.F25.Quarter);
    }
}
