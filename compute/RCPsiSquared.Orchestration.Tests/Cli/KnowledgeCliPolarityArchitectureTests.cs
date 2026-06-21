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
            .RegisterF88bPopcountCoherence()
            .Build();

    [Fact]
    public void Cli_TierT1D_OutputListsAllTenTier1DerivedClaims()
    {
        var cli = new KnowledgeCli(BuildPolarityRegistry());
        var output = cli.Render(new KnowledgeQuery.Tier(KnowledgeTier.Tier1Derived));

        // 9 Pi2 + 1 F88b = 10 Tier1Derived claims (F86 PolarityInheritanceLink is Tier2Verified,
        // not in this filter). QuarterAsBilinearMaxval + ArgmaxMaxvalPair joined RegisterPi2Family
        // on 2026-05-08/09 (see Pi2FamilyRegistration), growing the family from 7 to 9.
        Assert.Contains("PolynomialFoundationClaim", output);
        Assert.Contains("QubitDimensionalAnchorClaim", output);
        Assert.Contains("NinetyDegreeMirrorMemoryClaim", output);
        Assert.Contains("PolarityLayerOriginClaim", output);
        Assert.Contains("BilinearApexClaim", output);
        Assert.Contains("HalfAsStructuralFixedPointClaim", output);
        Assert.Contains("QuarterAsBilinearMaxvalClaim", output);
        Assert.Contains("KleinFourCellClaim", output);
        Assert.Contains("ArgmaxMaxvalPairClaim", output);
        Assert.Contains("PopcountCoherenceClaim", output);
        Assert.Contains("10 Claim(s)", output);
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
    public void Cli_DescendantsPolynomialFoundation_OutputListsTenDescendants()
    {
        var cli = new KnowledgeCli(BuildPolarityRegistry());
        var output = cli.Render(new KnowledgeQuery.Descendants(typeof(PolynomialFoundationClaim)));

        // 10 descendants: 8 Pi2 (everything except the trunk) + F86 PolarityInheritanceLink + F88b.
        // Was 8 before QuarterAsBilinearMaxval + ArgmaxMaxvalPair joined RegisterPi2Family (2026-05-08/09).
        Assert.Contains("10 Claim(s)", output);
        Assert.Contains("PolarityInheritanceLink", output);
        Assert.Contains("PopcountCoherenceClaim", output);
    }
}
