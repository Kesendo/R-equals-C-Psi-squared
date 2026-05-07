using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Orchestration.Render;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Orchestration.Tests.Render;

public class KnowledgeRendererTests
{
    [Fact]
    public void Render_EmptyRegistry_ContainsHeaderAndZeroes()
    {
        var registry = new ClaimRegistryBuilder().Build();
        var renderer = new KnowledgeRenderer(registry);

        var md = renderer.Render();

        Assert.Contains("# Knowledge Registry", md);
        Assert.Contains("0 claim(s)", md);
        Assert.Contains("0 edge(s)", md);
    }

    [Fact]
    public void Render_PolarityRegistry_GroupsByTier()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF86PolarityLink()
            .RegisterF88PopcountCoherence()
            .Build();
        var renderer = new KnowledgeRenderer(registry);

        var md = renderer.Render();

        Assert.Contains("## Tier 1 (derived)", md);
        Assert.Contains("## Tier 2 (hardware-verified)", md);  // Polarity link
        Assert.Contains("**PolynomialFoundationClaim**", md);
        Assert.Contains("**KleinFourCellClaim**", md);
        Assert.Contains("**PolarityInheritanceLink**", md);
    }

    [Fact]
    public void Render_PolarityRegistry_IncludesEdgesSection()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF86PolarityLink()
            .RegisterF88PopcountCoherence()
            .Build();
        var renderer = new KnowledgeRenderer(registry);

        var md = renderer.Render();

        Assert.Contains("## Edges", md);
        Assert.Contains("PolynomialFoundationClaim → QubitDimensionalAnchorClaim", md);
        Assert.Contains("KleinFourCellClaim → PopcountCoherenceClaim", md);
    }

    [Fact]
    public void Render_Cached_SecondCallReturnsSameString()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .Build();
        var renderer = new KnowledgeRenderer(registry);

        var first = renderer.Render();
        var second = renderer.Render();

        Assert.Same(first, second);
    }

    [Fact]
    public void Render_EmptyTier_IsOmitted()
    {
        // Empty registry has no Tier sections at all.
        var registry = new ClaimRegistryBuilder().Build();
        var renderer = new KnowledgeRenderer(registry);

        var md = renderer.Render();

        Assert.DoesNotContain("## Tier 1", md);
        Assert.DoesNotContain("## Retracted", md);
    }
}
