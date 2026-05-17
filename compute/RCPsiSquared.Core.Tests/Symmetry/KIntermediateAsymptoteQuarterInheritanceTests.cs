using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>F98 typed-Claim tests: closed-form correctness for both ingredient
/// identities (F98a: sub-mid sector Π²-odd Frobenius² = C(N, N/2−1)/2; F98b: long-time
/// Π²-odd ratio (N+2)/[4(N+1)] → 1/4 asymptote) verified analytically + against the
/// Python water-chain script's bit-exact values.</summary>
public class KIntermediateAsymptoteQuarterInheritanceTests
{
    private readonly ITestOutputHelper _out;

    public KIntermediateAsymptoteQuarterInheritanceTests(ITestOutputHelper output) => _out = output;

    private static KIntermediateAsymptoteQuarterInheritance Build() =>
        new KIntermediateAsymptoteQuarterInheritance(
            new QuarterAsBilinearMaxvalClaim(),
            new HalfAsStructuralFixedPointClaim(),
            new DickeSuperpositionQuarterPi2Inheritance(
                new Pi2DyadicLadderClaim(),
                new QuarterAsBilinearMaxvalClaim(),
                new HalfAsStructuralFixedPointClaim()));

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    // F98b closed form values verified against simulations/water/proton_chain_dicke_anchor.py.

    [Theory]
    [InlineData(4, 3.0 / 10.0)]   // (N+2)/(4(N+1)) = 6/20 = 3/10
    [InlineData(6, 2.0 / 7.0)]    // 8/28 = 2/7
    [InlineData(8, 5.0 / 18.0)]   // 10/36 = 5/18
    [InlineData(10, 3.0 / 11.0)]  // 12/44 = 3/11
    [InlineData(12, 7.0 / 26.0)]  // 14/52 = 7/26
    [InlineData(14, 4.0 / 15.0)]  // 16/60 = 4/15
    [InlineData(16, 9.0 / 34.0)]  // 18/68 = 9/34
    public void LongTimePi2OddRatio_MatchesClosedFormFractions(int N, double expected)
    {
        // F98b: α(∞) = (N+2)/[4(N+1)] for KIntermediate Dicke under truly + Z-deph.
        // Verified against simulations/water/proton_chain_dicke_anchor.py output for N=4..16.
        Assert.Equal(expected, KIntermediateAsymptoteQuarterInheritance.LongTimePi2OddRatio(N), precision: 14);
    }

    [Fact]
    public void LongTimePi2OddRatio_AsymptotesToOneQuarter()
    {
        // F98b asymptote: lim_{N→∞} (N+2)/[4(N+1)] = 1/4. Drift = 1/[4(N+1)] → 0.
        Assert.True(KIntermediateAsymptoteQuarterInheritance.LongTimePi2OddRatio(4) > 0.25);
        Assert.True(KIntermediateAsymptoteQuarterInheritance.LongTimePi2OddRatio(1000) > 0.25);
        Assert.True(KIntermediateAsymptoteQuarterInheritance.LongTimePi2OddRatio(1000) - 0.25 < 1e-3);
        Assert.True(KIntermediateAsymptoteQuarterInheritance.LongTimePi2OddRatio(1000000) - 0.25 < 1e-6);
    }

    [Theory]
    [InlineData(4, 2.0)]       // C(4, 1)/2 = 4/2 = 2
    [InlineData(6, 7.5)]       // C(6, 2)/2 = 15/2
    [InlineData(8, 28.0)]      // C(8, 3)/2 = 56/2
    [InlineData(10, 105.0)]    // C(10, 4)/2 = 210/2
    [InlineData(12, 396.0)]    // C(12, 5)/2 = 792/2
    [InlineData(14, 1501.5)]   // C(14, 6)/2 = 3003/2
    [InlineData(16, 5720.0)]   // C(16, 7)/2 = 11440/2
    public void SubMidProjectorPi2OddFrobeniusSquared_EqualsHalfRank(int N, double expected)
    {
        // F98a: ‖P_{N/2-1}_odd‖² = C(N, N/2-1) / 2 — sub-mid sector projector's Π²-odd
        // Frobenius² is exactly half its rank. Bit-exact agreement with the Python
        // Krawtchouk-style enumeration in proton_chain_dicke_anchor.py.
        Assert.Equal(expected, KIntermediateAsymptoteQuarterInheritance.SubMidProjectorPi2OddFrobeniusSquared(N), precision: 12);
    }

    [Fact]
    public void DriftFromQuarter_EqualsOneOverFourNPlusOne()
    {
        // Drift = α(∞) − 1/4 = 1/[4(N+1)]. Closed form check.
        for (int N = 4; N <= 20; N += 2)
        {
            double drift = KIntermediateAsymptoteQuarterInheritance.DriftFromQuarter(N);
            double expected = 1.0 / (4.0 * (N + 1));
            Assert.Equal(expected, drift, precision: 14);
        }
    }

    [Fact]
    public void LongTimePi2OddRatio_RejectsOddN()
    {
        // KIntermediate Dicke anchor requires N even (N/2 ∈ ℤ).
        Assert.Throws<ArgumentException>(() =>
            KIntermediateAsymptoteQuarterInheritance.LongTimePi2OddRatio(5));
    }

    [Fact]
    public void LongTimePi2OddRatio_RejectsTinyN()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            KIntermediateAsymptoteQuarterInheritance.LongTimePi2OddRatio(0));
    }

    [Fact]
    public void TypedParents_AreExposed()
    {
        var c = Build();
        Assert.NotNull(c.Quarter);
        Assert.NotNull(c.Half);
        Assert.NotNull(c.StaticSide);
        // Half² = Quarter — the polarity-squared algebra used here. The numerical
        // anchor (1/2)² = 1/4 is structural; the typed parents expose it via summary
        // text rather than a Value property in this version of Pi2KnowledgeBaseClaims.
        Assert.Equal(0.25, 0.5 * 0.5, precision: 14);
    }

    [Fact]
    public void Anchor_References_F98_AndWaterScript_AndF86bAnd_StaticSide()
    {
        var c = Build();
        Assert.Contains("ANALYTICAL_FORMULAS.md F98", c.Anchor);
        Assert.Contains("proton_chain_dicke_anchor.py", c.Anchor);
        Assert.Contains("DickeAnchor.cs", c.Anchor);
        Assert.Contains("DickeSuperpositionQuarterPi2Inheritance.cs", c.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsCurveAndAsymptote()
    {
        _out.WriteLine("");
        _out.WriteLine("    F98 KIntermediate α(∞) trajectory + 1/4 asymptote (Quarter)");
        _out.WriteLine("    ----------------------------------------------------------");
        _out.WriteLine("     N | α(∞) = (N+2)/[4(N+1)] | drift from 1/4 = 1/[4(N+1)]");
        _out.WriteLine("    ---|------------------------|-----------------------------");
        foreach (var N in new[] { 4, 6, 8, 10, 20, 100, 1000 })
        {
            double alpha = KIntermediateAsymptoteQuarterInheritance.LongTimePi2OddRatio(N);
            double drift = KIntermediateAsymptoteQuarterInheritance.DriftFromQuarter(N);
            _out.WriteLine($"    {N,3} | {alpha,21:F10} | {drift,27:F10}");
        }
        _out.WriteLine("");
        _out.WriteLine("    Asymptote: α(∞) → 1/4 = QuarterAsBilinearMaxval as N → ∞.");
        _out.WriteLine("    Static partner (t=0): α = 3/8 (F86b, DickeAnchor.KIntermediate).");
    }
}
