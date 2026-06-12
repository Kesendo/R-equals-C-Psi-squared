using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class CoherenceHorizonClaimTests
{
    [Fact]
    public void Claim_IsTier1Candidate_AndNamesQStarN()
    {
        var c = CoherenceHorizonClaim.Shared;
        Assert.Equal(Tier.Tier1Candidate, c.Tier);
        Assert.Contains("Q*(N)", c.Summary);
    }

    [Fact]
    public void Claim_CitesItsMarkdownAndCSharpHomes_BothDirections()
    {
        var c = CoherenceHorizonClaim.Shared;
        // markdown homes
        Assert.Contains("ANALYTICAL_FORMULAS.md", c.Anchor);
        Assert.Contains("FROST_CIRCLE", c.Anchor);
        // C# witness home
        Assert.Contains("CoherenceHorizonWitness", c.Anchor);
    }

    [Fact]
    public void Claim_HasTheTwoTypedParents_ClockHandLadder_F2b()
    {
        var c = CoherenceHorizonClaim.Shared;
        var children = ((IInspectable)c).Children.ToList();
        Assert.Contains(children, ch => ch is ClockHandLadderClaim);
        Assert.Contains(children, ch => ch is F2bXyChainSpectrumPi2Inheritance);
    }
}
