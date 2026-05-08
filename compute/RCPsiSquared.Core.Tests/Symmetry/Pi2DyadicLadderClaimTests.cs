using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class Pi2DyadicLadderClaimTests
{
    private readonly ITestOutputHelper _out;

    public Pi2DyadicLadderClaimTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var ladder = new Pi2DyadicLadderClaim();
        Assert.Equal(Tier.Tier1Derived, ladder.Tier);
    }

    [Theory]
    [InlineData(0, 2.0)]      // d=2 root
    [InlineData(1, 1.0)]      // identity scale
    [InlineData(2, 0.5)]      // half fixed point
    [InlineData(3, 0.25)]     // bilinear maxval
    [InlineData(4, 0.125)]    // first open prediction
    [InlineData(5, 0.0625)]
    [InlineData(8, 1.0 / 128.0)]    // 2^(1−8) = 2^(−7) = 1/128
    public void Term_ClosedForm(int n, double expected)
    {
        var ladder = new Pi2DyadicLadderClaim();
        Assert.Equal(expected, ladder.Term(n), precision: 12);
    }

    [Fact]
    public void Term_NegativeN_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        Assert.Throws<ArgumentOutOfRangeException>(() => ladder.Term(-1));
    }

    [Fact]
    public void HalfingRule_TermNTimesOneHalfEqualsTermNPlusOne()
    {
        // a_{n+1} = a_n × (1/2). The defining algebraic identity Tom named: every step is
        // the previous halved.
        var ladder = new Pi2DyadicLadderClaim();
        for (int n = 0; n < 8; n++)
            Assert.Equal(ladder.Term(n) * 0.5, ladder.Term(n + 1), precision: 14);
    }

    [Fact]
    public void KnownAnchors_HasThreeEntries()
    {
        var ladder = new Pi2DyadicLadderClaim();
        Assert.Equal(3, ladder.KnownAnchors.Count);
    }

    [Fact]
    public void KnownAnchors_PinsQubitDimAtN0_HalfFixedAtN2_BilinearMaxvalAtN3()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var byN = ladder.KnownAnchors.ToDictionary(a => a.N);

        Assert.Contains(0, byN.Keys);
        Assert.Contains(2, byN.Keys);
        Assert.Contains(3, byN.Keys);

        Assert.Equal(typeof(QubitDimensionalAnchorClaim), byN[0].ClaimType);
        Assert.Equal(2.0, byN[0].Value, precision: 12);

        Assert.Equal(typeof(HalfAsStructuralFixedPointClaim), byN[2].ClaimType);
        Assert.Equal(0.5, byN[2].Value, precision: 12);

        Assert.Equal(typeof(QuarterAsBilinearMaxvalClaim), byN[3].ClaimType);
        Assert.Equal(0.25, byN[3].Value, precision: 12);
    }

    [Fact]
    public void KnownAnchors_ValuesMatchTermAtSameIndex()
    {
        var ladder = new Pi2DyadicLadderClaim();
        foreach (var anchor in ladder.KnownAnchors)
            Assert.Equal(ladder.Term(anchor.N), anchor.Value, precision: 12);
    }

    [Theory]
    [InlineData(0, true)]
    [InlineData(1, false)]   // trivial identity scale; no typed Claim
    [InlineData(2, true)]
    [InlineData(3, true)]
    [InlineData(4, false)]   // open prediction territory
    [InlineData(5, false)]
    public void IsKnownAnchorIndex(int n, bool expected)
    {
        var ladder = new Pi2DyadicLadderClaim();
        Assert.Equal(expected, ladder.IsKnownAnchorIndex(n));
    }

    [Fact]
    public void AnchorAt_KnownIndex_ReturnsAnchor()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var a3 = ladder.AnchorAt(3);
        Assert.NotNull(a3);
        Assert.Equal(typeof(QuarterAsBilinearMaxvalClaim), a3!.ClaimType);
    }

    [Fact]
    public void AnchorAt_UnknownIndex_ReturnsNull()
    {
        var ladder = new Pi2DyadicLadderClaim();
        Assert.Null(ladder.AnchorAt(1));
        Assert.Null(ladder.AnchorAt(4));
    }

    [Fact]
    public void OpenAnchorIndices_ContainsTrivialAndPredictionRange()
    {
        var ladder = new Pi2DyadicLadderClaim();
        Assert.Contains(1, ladder.OpenAnchorIndices);    // identity scale
        Assert.Contains(4, ladder.OpenAnchorIndices);    // first open prediction
        Assert.Contains(5, ladder.OpenAnchorIndices);
        Assert.DoesNotContain(0, ladder.OpenAnchorIndices);
        Assert.DoesNotContain(2, ladder.OpenAnchorIndices);
        Assert.DoesNotContain(3, ladder.OpenAnchorIndices);
    }

    [Fact]
    public void Anchor_References_HalfReflection_AndQuarterProof()
    {
        var ladder = new Pi2DyadicLadderClaim();
        Assert.Contains("ON_THE_HALF.md", ladder.Anchor);
        Assert.Contains("PROOF_BLOCK_CPSI_QUARTER.md", ladder.Anchor);
    }

    [Fact]
    public void Reconnaissance_LadderEmitsFirstNineTerms()
    {
        var ladder = new Pi2DyadicLadderClaim();
        _out.WriteLine("  n |   a_n      | known anchor (Claim)");
        _out.WriteLine("  --|------------|---------------------");
        for (int n = 0; n < 9; n++)
        {
            string anchor = ladder.AnchorAt(n)?.ClaimName ?? "(open)";
            _out.WriteLine($"  {n} | {ladder.Term(n),10:G6} | {anchor}");
        }
    }
}
