using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Tests for the γ-axis polarity-mirror structure (complementary
/// to the α-axis FractionReferenceGraph). Verifies the ±γ pairing per
/// F99 anker and the (1−γ²)/2 folding identity.</summary>
public class PolarityMirrorMapTests
{
    private readonly ITestOutputHelper _out;

    public PolarityMirrorMapTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void CanonicalMap_HasFivePairs()
    {
        var m = new PolarityMirrorMap();
        Assert.Equal(5, m.Pairs.Count);
    }

    [Fact]
    public void FourPairs_AreNonTrivialMirrors()
    {
        // Mirror (γ=±1), depth-3 (γ=±√3/2), silver (γ=±√2/2), KIntermediate (γ=±1/2)
        var m = new PolarityMirrorMap();
        Assert.Equal(4, m.NonTrivialMirrorPairs.Count);
    }

    [Fact]
    public void OnePair_IsSelfMirror_TheGenericGammaZeroCase()
    {
        var m = new PolarityMirrorMap();
        var selfMirrors = m.SelfMirrorPairs;
        Assert.Single(selfMirrors);
        Assert.Equal(0.0, selfMirrors[0].GammaPlus);
        Assert.Equal(1.0 / 2.0, selfMirrors[0].AlphaImage);
    }

    [Fact]
    public void AllPairs_FoldToClaimedAlpha_BitExact()
    {
        var m = new PolarityMirrorMap();
        Assert.True(m.AllPairsFoldToClaimedAlpha(),
            "F86b folding identity: (1−γ²)/2 must equal the claimed AlphaImage for every pair, both ±γ sides");
    }

    [Fact]
    public void PairForAlpha_FindsCorrectGammaPair()
    {
        var m = new PolarityMirrorMap();
        var quarter = m.PairForAlpha(1.0 / 4.0);
        Assert.NotNull(quarter);
        Assert.Equal(System.Math.Sqrt(2) / 2, quarter!.GammaPlus, precision: 12);
        Assert.Equal(-System.Math.Sqrt(2) / 2, quarter.GammaMinus, precision: 12);
    }

    [Fact]
    public void PairForGamma_FindsBothSides()
    {
        var m = new PolarityMirrorMap();
        var fromPlus = m.PairForGamma(System.Math.Sqrt(3) / 2);
        var fromMinus = m.PairForGamma(-System.Math.Sqrt(3) / 2);
        Assert.NotNull(fromPlus);
        Assert.NotNull(fromMinus);
        // Both queries return the same pair (the depth-3 one with α=1/8)
        Assert.Equal(1.0 / 8.0, fromPlus!.AlphaImage);
        Assert.Equal(1.0 / 8.0, fromMinus!.AlphaImage);
    }

    [Fact]
    public void GenericPair_AlphaIsOneHalf_AndIsSelfMirror()
    {
        var m = new PolarityMirrorMap();
        var generic = m.PairForAlpha(1.0 / 2.0);
        Assert.NotNull(generic);
        Assert.True(generic!.IsSelfMirror, "γ=0 Generic is its own polarity mirror");
        Assert.Equal(0.0, generic.GammaPlus);
        Assert.Equal(0.0, generic.GammaMinus);
    }

    [Fact]
    public void MirrorPair_BothSides_GiveZeroAlpha()
    {
        var m = new PolarityMirrorMap();
        var mirror = m.PairForAlpha(0.0);
        Assert.NotNull(mirror);
        Assert.False(mirror!.IsSelfMirror, "γ=±1 Mirror has two distinct realizations");
        Assert.Equal(1.0, mirror.GammaPlus);
        Assert.Equal(-1.0, mirror.GammaMinus);
        Assert.Equal(0.0, mirror.FoldedAlpha(1.0), precision: 14);
        Assert.Equal(0.0, mirror.FoldedAlpha(-1.0), precision: 14);
    }

    [Fact]
    public void FoldingIdentity_HoldsAtOffPairValues()
    {
        // Spot-check the folding formula at a non-anker γ — folding gives a non-anker α
        var pair = new PolarityMirrorMap.PolarityPair(0.3, -0.3, (1 - 0.09) / 2, "test", "test");
        Assert.Equal(0.455, pair.FoldedAlpha(0.3), precision: 12);
        Assert.Equal(0.455, pair.FoldedAlpha(-0.3), precision: 12);
    }

    [Fact]
    public void NullPairs_Throws()
    {
        Assert.Throws<System.ArgumentNullException>(() => new PolarityMirrorMap(null!));
    }

    [Fact]
    public void Render_EmitsTheFullPairTable()
    {
        var m = new PolarityMirrorMap();
        var rendered = m.Render();
        Assert.Contains("Polarity Mirror Map", rendered);
        Assert.Contains("AllPairsFoldToClaimedAlpha = True", rendered);
        Assert.Contains("Non-trivial mirror pairs: 4 of 5", rendered);
        Assert.Contains("Self-mirror pairs:        1 of 5", rendered);
        Assert.Contains("PolarityLayerOriginClaim", rendered);
        _out.WriteLine(rendered);
    }
}
