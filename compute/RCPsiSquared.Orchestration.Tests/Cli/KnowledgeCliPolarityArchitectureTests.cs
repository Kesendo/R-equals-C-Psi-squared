using KnowledgeTier = RCPsiSquared.Core.Knowledge.Tier;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Orchestration.Cli;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Orchestration.Tests.Cli;

public class KnowledgeCliPolarityArchitectureTests
{
    private static ClaimRegistry BuildPolarityRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF86PolarityLink()
            .RegisterF88PopcountCoherence()
            .Build();

    [Fact]
    public void Cli_TierT1D_OutputListsAllEightTier1DerivedClaims()
    {
        var cli = new KnowledgeCli(BuildPolarityRegistry());
        var output = cli.Render(new KnowledgeQuery.Tier(KnowledgeTier.Tier1Derived));

        // 7 Pi2 + 1 F88 = 8 Tier1Derived claims (F86 is Tier2Verified, not in this filter)
        Assert.Contains("PolynomialFoundationClaim", output);
        Assert.Contains("QubitDimensionalAnchorClaim", output);
        Assert.Contains("NinetyDegreeMirrorMemoryClaim", output);
        Assert.Contains("PolarityLayerOriginClaim", output);
        Assert.Contains("BilinearApexClaim", output);
        Assert.Contains("HalfAsStructuralFixedPointClaim", output);
        Assert.Contains("KleinFourCellClaim", output);
        Assert.Contains("PopcountCoherenceClaim", output);
        Assert.Contains("8 Claim(s)", output);
    }

    [Fact]
    public void Cli_AncestorsPolarityInheritanceLink_OutputContainsPolynomialAndTrio()
    {
        var cli = new KnowledgeCli(BuildPolarityRegistry());
        var output = cli.Render(new KnowledgeQuery.Ancestors(typeof(PolarityInheritanceLink)));

        Assert.Contains("PolynomialFoundationClaim", output);
        Assert.Contains("PolarityLayerOriginClaim", output);
        Assert.Contains("HalfAsStructuralFixedPointClaim", output);
        Assert.Contains("QubitDimensionalAnchorClaim", output);
    }

    [Fact]
    public void Cli_DescendantsPolynomialFoundation_OutputListsEightDescendants()
    {
        var cli = new KnowledgeCli(BuildPolarityRegistry());
        var output = cli.Render(new KnowledgeQuery.Descendants(typeof(PolynomialFoundationClaim)));

        // 8 descendants: 6 Pi2 (everything except the trunk) + F86 + F88
        Assert.Contains("8 Claim(s)", output);
        Assert.Contains("PolarityInheritanceLink", output);
        Assert.Contains("PopcountCoherenceClaim", output);
    }
}
