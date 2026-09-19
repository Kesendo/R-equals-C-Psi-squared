using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Witnesses for <see cref="TwoReadingsClaim"/>: a discoverable open synthesis
/// whose examples retain their individual scope and do not acquire a common genealogy.</summary>
public class TwoReadingsClaimTests
{
    [Fact]
    public void Build_Tier_IsOpenQuestion()
    {
        var claim = new TwoReadingsClaim();
        Assert.Equal(Tier.OpenQuestion, claim.Tier);
    }

    [Fact]
    public void Build_ConstructorIsParentless()
    {
        var constructor = Assert.Single(typeof(TwoReadingsClaim).GetConstructors());
        Assert.Empty(constructor.GetParameters());
    }

    [Fact]
    public void Build_ExplicitlyDeniesTheUniversalImplication()
    {
        var claim = new TwoReadingsClaim();
        string text = $"{claim.DisplayName} {claim.Summary}".ToLowerInvariant();

        Assert.Contains("does not imply", text);
        Assert.Contains("exactly two", text);
    }

    [Fact]
    public void Build_SummarySeparatesExactExamplesFromInterpretiveComparisons()
    {
        var claim = new TwoReadingsClaim();
        Assert.Contains("bra/ket", claim.Summary);
        Assert.Contains("argmax/maxval", claim.Summary);
        Assert.Contains("interpretive", claim.Summary, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void Build_AnchorCitesTheCatalogueWithoutClaimingDerivation()
    {
        var claim = new TwoReadingsClaim();
        Assert.Contains("ON_BOTH_SIDES_OF_THE_MIRROR", claim.Anchor);
        Assert.Contains("ON_THE_HALF", claim.Anchor);
        Assert.DoesNotContain("parent", claim.Anchor, StringComparison.OrdinalIgnoreCase);
    }
}
