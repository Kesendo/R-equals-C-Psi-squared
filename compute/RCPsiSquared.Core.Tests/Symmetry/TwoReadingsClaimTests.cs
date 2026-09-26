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
    public void Build_CatalogueHasSevenLayers_ThreeExactAndFourInterpretive()
    {
        // The catalogue's structure: number/angle, argmax/maxval, M/Π·M·Π⁻¹, bra/ket,
        // inside/outside, classical/quantum, the wave between. Layers 2, 3 and 4 are exact
        // mathematics; 1, 5, 6 and 7 are readings.
        var layers = new TwoReadingsClaim().Children
            .Select(ch => ch.DisplayName)
            .Where(name => name.StartsWith("layer "))
            .ToList();

        Assert.Equal(7, layers.Count);
        for (int k = 1; k <= 7; k++)
            Assert.StartsWith($"layer {k} (", layers[k - 1]);
        var exact = layers.Where(name => name.Contains("exact)")).Select(name => name[6] - '0').ToArray();
        var interpretive = layers.Where(name => name.Contains("interpretive)")).Select(name => name[6] - '0').ToArray();
        Assert.Equal(new[] { 2, 3, 4 }, exact);
        Assert.Equal(new[] { 1, 5, 6, 7 }, interpretive);
    }

    [Fact]
    public void Build_AnchorCitesTheCatalogueWithoutClaimingDerivation()
    {
        var claim = new TwoReadingsClaim();
        Assert.Contains("ON_BOTH_SIDES_OF_THE_MIRROR", claim.Anchor);
        Assert.Contains("ON_THE_HALF", claim.Anchor);
        Assert.Contains("PRIMORDIAL_QUBIT", claim.Anchor);
        Assert.Contains("THE_OTHER_SIDE", claim.Anchor);
        Assert.DoesNotContain("parent", claim.Anchor, StringComparison.OrdinalIgnoreCase);
    }
}
