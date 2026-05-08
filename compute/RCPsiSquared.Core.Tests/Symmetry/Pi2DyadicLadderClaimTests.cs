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

    [Theory]
    [InlineData(-1, 4.0)]    // upper side: 4 = d² for N=1
    [InlineData(-2, 8.0)]
    [InlineData(-3, 16.0)]   // upper side: 16 = d² for N=2
    [InlineData(-5, 64.0)]   // upper side: 64 = d² for N=3
    public void Term_NegativeN_GivesUpperSideDoubling(int n, double expected)
    {
        // The ladder is bidirectional. n < 0 produces the doubling powers 4, 8, 16, 32, ...
        // which carry the operator-space dimension d² = 4^N at n = −(2N − 1).
        var ladder = new Pi2DyadicLadderClaim();
        Assert.Equal(expected, ladder.Term(n), precision: 12);
    }

    [Fact]
    public void HalfingRule_TermNTimesOneHalfEqualsTermNPlusOne()
    {
        // a_{n+1} = a_n × (1/2). The defining algebraic identity Tom named: every step is
        // the previous halved. Verified across both sides of the ladder.
        var ladder = new Pi2DyadicLadderClaim();
        for (int n = -4; n < 8; n++)
            Assert.Equal(ladder.Term(n) * 0.5, ladder.Term(n + 1), precision: 14);
    }

    [Theory]
    [InlineData(0, 2)]     // a_0 = 2 ↔ a_2 = 1/2
    [InlineData(2, 0)]     // and back
    [InlineData(3, -1)]    // a_3 = 1/4 ↔ a_{-1} = 4
    [InlineData(-1, 3)]    // and back
    [InlineData(4, -2)]    // a_4 = 1/8 ↔ a_{-2} = 8
    [InlineData(1, 1)]     // self-mirror at n=1
    public void MirrorPartnerIndex_ReturnsTwoMinusN(int n, int expected)
    {
        var ladder = new Pi2DyadicLadderClaim();
        Assert.Equal(expected, ladder.MirrorPartnerIndex(n));
    }

    [Theory]
    [InlineData(-3)]
    [InlineData(-1)]
    [InlineData(0)]
    [InlineData(1)]    // self-mirror, a_1 · a_1 = 1
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(7)]
    public void ProductWithMirrorPartner_AlwaysEqualsOne(int n)
    {
        // Tom's structural rule (2026-05-08): "es muss auch eine Regel geben wo es wieder
        // auf 1/1 springt." The rule is the inversion symmetry a_n · a_{2−n} = 1 — every
        // term has a multiplicative mirror; nothing escapes the unit circle.
        var ladder = new Pi2DyadicLadderClaim();
        Assert.Equal(1.0, ladder.ProductWithMirrorPartner(n), precision: 14);
    }

    [Fact]
    public void SelfMirrorIndex_IsOne()
    {
        // a_1 = 1, the unique fixpoint of the inversion. Identity scale, no typed Claim
        // hinged because "1 = 1" is trivial.
        var ladder = new Pi2DyadicLadderClaim();
        Assert.Equal(1, ladder.SelfMirrorIndex);
        Assert.Equal(ladder.SelfMirrorIndex, ladder.MirrorPartnerIndex(ladder.SelfMirrorIndex));
        Assert.Equal(1.0, ladder.Term(ladder.SelfMirrorIndex), precision: 14);
    }

    [Theory]
    [InlineData(1, -1, 4.0)]    // N=1: a_{-1} = 4 = 4^1 = d²
    [InlineData(2, -3, 16.0)]   // N=2: a_{-3} = 16 = 4²
    [InlineData(3, -5, 64.0)]   // N=3: a_{-5} = 64 = 4³
    [InlineData(4, -7, 256.0)]  // N=4: a_{-7} = 256 = 4⁴
    public void OperatorSpaceIndexForN_LandsOnFourToTheN(int N, int expectedIndex, double expectedValue)
    {
        // The upper side carries d² = 4^N for an N-qubit system at index −(2N − 1).
        // Connects the dyadic halving ladder to memory feedback_d2_operator_space.md:
        // d²=4^N is the natural operator-space dimension; it lives on the ladder, not
        // outside it.
        var ladder = new Pi2DyadicLadderClaim();
        Assert.Equal(expectedIndex, ladder.OperatorSpaceIndexForN(N));
        Assert.Equal(expectedValue, ladder.Term(ladder.OperatorSpaceIndexForN(N)), precision: 12);
        Assert.Equal(Math.Pow(4.0, N), ladder.Term(ladder.OperatorSpaceIndexForN(N)), precision: 12);
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
