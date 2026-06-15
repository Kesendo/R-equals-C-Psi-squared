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
        Assert.Contains("2-qubit", name);          // proven scope, narrowed from the over-broad "N=3-5"
        Assert.Contains("N=2", name);
        Assert.Contains("OPEN", name);             // the N≥3 full-state envelope is open
        Assert.Contains("RISES", name);            // and genuinely rises at N≥4 strong coupling
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
