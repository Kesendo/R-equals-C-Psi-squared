using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Diagnostics.Ptf;

namespace RCPsiSquared.Diagnostics.Tests.Ptf;

/// <summary>The resolution-limit reading of F124 (the optics/signal facet): κ ~ N², the contrast √κ ~ N, the
/// staggered q=π the diffraction limit. Covers tier (Tier1Derived), the single typed parent
/// (BandEdgeTransitionInvariantClaim), the summary, the anchor, the gate-first battery, and the registry
/// wiring + ancestry.</summary>
public class BandEdgeResolutionLimitClaimTests
{
    private static BandEdgeResolutionLimitClaim Build() => BandEdgeResolutionLimitClaim.Build();

    [Fact]
    public void Claim_TierIsTier1Derived()
    {
        // The single parent F124 is Tier1Derived, so the tier-inheritance invariant leaves the child Tier1Derived.
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void Build_WiresF124Parent_Tier1Derived()
    {
        var claim = Build();
        Assert.NotNull(claim.F124);
        Assert.Equal(Tier.Tier1Derived, claim.F124.Tier);
        var kids = ((IInspectable)claim).Children.Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(BandEdgeTransitionInvariantClaim), kids);
    }

    [Fact]
    public void Summary_NamesResolutionAndKappaAndContrast()
    {
        var summary = Build().Summary;
        Assert.Contains("resolution", summary, System.StringComparison.OrdinalIgnoreCase);
        Assert.Contains("κ", summary);
        Assert.Contains("√κ", summary);
    }

    [Fact]
    public void Anchor_CitesVerifierAndWitness()
    {
        var anchor = Build().Anchor;
        Assert.Contains("f124_inverse_problem_gate.py", anchor);
        Assert.Contains("BandEdgeResolutionLimitWitness.cs", anchor);
    }

    [Fact]
    public void Battery_AllPass()
    {
        var claim = Build();
        Assert.NotEmpty(claim.Cases);
        Assert.Equal(claim.Cases.Count, claim.PassCount);
        foreach (var c in claim.Cases)
            Assert.True(c.Passes, $"battery case '{c.Name}': expected {c.Expected}, got {c.Actual}");
    }

    // --- Registration ---

    [Fact]
    public void BuildDefault_ContainsBandEdgeResolutionLimitClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<BandEdgeResolutionLimitClaim>());
    }

    [Fact]
    public void Claim_Ancestors_ContainF124()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<BandEdgeResolutionLimitClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(BandEdgeTransitionInvariantClaim), ancestors);
    }
}
