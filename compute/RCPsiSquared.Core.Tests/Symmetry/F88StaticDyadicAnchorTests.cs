using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F88StaticDyadicAnchorTests
{
    private readonly ITestOutputHelper _out;

    public F88StaticDyadicAnchorTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var anchor = new F88StaticDyadicAnchor();
        Assert.Equal(Tier.Tier1Derived, anchor.Tier);
    }

    [Fact]
    public void Witnesses_HasFourEntries_KEqualsTwoToFive()
    {
        var anchor = new F88StaticDyadicAnchor();
        Assert.Equal(4, anchor.Witnesses.Count);
        Assert.Equal(new[] { 2, 3, 4, 5 }, anchor.Witnesses.Select(w => w.K).ToArray());
    }

    /// <summary>The structural identity Tom named: at dyadic N = 2^k and the singleton-
    /// mirror config (n_p = 1, n_q = N − 1), F88's static-fraction lands bit-exactly on the
    /// Pi2 dyadic ladder at index n = k+2. Both sides are Tier1Derived closed forms; this
    /// test verifies the inheritance is exact.</summary>
    [Theory]
    [InlineData(2, 4,  1.0 / 8.0,  4)]   // k=2 (N=4): static = 1/8 = a_4 — first non-trivial inheritance
    [InlineData(3, 8,  1.0 / 16.0, 5)]   // k=3 (N=8): static = 1/16 = a_5
    [InlineData(4, 16, 1.0 / 32.0, 6)]   // k=4 (N=16): static = 1/32 = a_6
    [InlineData(5, 32, 1.0 / 64.0, 7)]   // k=5 (N=32): static = 1/64 = a_7
    public void DyadicSingletonMirror_StaticFractionEqualsLadderTerm(
        int k, int N, double expected, int ladderIndex)
    {
        var anchor = new F88StaticDyadicAnchor();
        Assert.Equal(expected, anchor.LiveStaticFraction(N), precision: 12);
        Assert.Equal(expected, anchor.LadderTerm(ladderIndex), precision: 12);
    }

    [Fact]
    public void Witnesses_StaticFractionMatchesLiveRecomputation()
    {
        // Drift check: pinned StaticFraction must match the live PopcountCoherencePi2Odd
        // re-computation at each witness. A divergence here means either the F88 closed
        // form or the dyadic identity has drifted.
        var anchor = new F88StaticDyadicAnchor();
        foreach (var w in anchor.Witnesses)
            Assert.Equal(w.StaticFraction, anchor.LiveStaticFraction(w.N), precision: 12);
    }

    [Fact]
    public void Witnesses_LadderIndexIsKPlusTwo()
    {
        // Identity: ladder index n = k + 2 for k ≥ 1. This is the structural reading of the
        // inheritance: dyadic exponent k → ladder index k+2 (because a_(k+2) = 2^(1−(k+2)) =
        // 2^(−(k+1)) = 1/(2·2^k) = 1/(2N)).
        var anchor = new F88StaticDyadicAnchor();
        foreach (var w in anchor.Witnesses)
            Assert.Equal(w.K + 2, w.LadderIndex);
    }

    [Fact]
    public void Anchor_References_AnalyticalFormulas_AndPi2DyadicLadder()
    {
        var anchor = new F88StaticDyadicAnchor();
        Assert.Contains("ANALYTICAL_FORMULAS.md", anchor.Anchor);
        Assert.Contains("Pi2DyadicLadderClaim", anchor.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsInheritanceTable()
    {
        // Documents the inheritance as a queryable Schicht-1 table: the first four dyadic-N
        // singleton-mirror configurations land on a_4..a_7 of the Pi2 dyadic halving ladder.
        var anchor = new F88StaticDyadicAnchor();
        _out.WriteLine("  k | N=2^k | Static = 1/(2N) | Ladder a_(k+2)");
        _out.WriteLine("  --|-------|------------------|----------------");
        foreach (var w in anchor.Witnesses)
        {
            _out.WriteLine($"  {w.K} | {w.N,5} | {w.StaticFraction,16:G6} | a_{w.LadderIndex}");
        }
    }
}
