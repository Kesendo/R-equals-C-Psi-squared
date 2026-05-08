using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F88PopcountPairLensTests
{
    private readonly ITestOutputHelper _out;

    public F88PopcountPairLensTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var lens = new F88PopcountPairLens(N: 6, np: 1, nq: 2);
        Assert.Equal(Tier.Tier1Derived, lens.Tier);
    }

    [Theory]
    [InlineData(2, 0, 2)]    // n_p + n_q = N=2 → mirror
    [InlineData(4, 1, 3)]    // n_p + n_q = N=4 → mirror
    [InlineData(5, 2, 3)]    // n_p + n_q = N=5 → mirror
    [InlineData(6, 3, 3)]    // n_p + n_q = N=6 → mirror (intra-half)
    public void PopcountMirror_AlphaIsZero(int N, int np, int nq)
    {
        var lens = new F88PopcountPairLens(N, np, nq);
        Assert.Equal(ConfigurationKind.PopcountMirror, lens.Kind);
        Assert.Equal(0.0, lens.Alpha, precision: 12);
    }

    [Theory]
    [InlineData(4, 2, 0)]    // even N=4, n_p=N/2=2, n_q=0 ≠ mirror → K-intermediate
    [InlineData(6, 3, 1)]    // even N=6, n_p=N/2=3, n_q=1 ≠ mirror → K-intermediate
    [InlineData(6, 2, 3)]    // even N=6, n_q=N/2=3, n_p=2 ≠ mirror → K-intermediate
    public void KIntermediate_AlphaMatchesClosedForm(int N, int np, int nq)
    {
        var lens = new F88PopcountPairLens(N, np, nq);
        Assert.Equal(ConfigurationKind.KIntermediate, lens.Kind);
        double expected = PopcountCoherencePi2Odd.AlphaKIntermediateClosed(N, np, nq);
        Assert.Equal(expected, lens.Alpha, precision: 12);
    }

    [Theory]
    [InlineData(6, 1, 2)]    // c=2 stratum, generic
    [InlineData(7, 1, 2)]    // c=2 stratum, generic, odd N
    [InlineData(8, 1, 2)]    // c=2 stratum, generic, larger
    [InlineData(5, 1, 2)]    // odd N
    public void Generic_AlphaIsOneHalf(int N, int np, int nq)
    {
        var lens = new F88PopcountPairLens(N, np, nq);
        Assert.Equal(ConfigurationKind.Generic, lens.Kind);
        Assert.Equal(0.5, lens.Alpha, precision: 12);
    }

    [Theory]
    [InlineData(2, 0, 1)]
    [InlineData(3, 1, 2)]
    [InlineData(4, 0, 2)]
    [InlineData(5, 1, 3)]
    [InlineData(6, 1, 2)]
    [InlineData(7, 2, 3)]
    public void AlphaAnchor_BitExactMatchesKrawtchouk(int N, int np, int nq)
    {
        // The closed-form α (anchor) must match the universal Krawtchouk verifier
        // bit-exactly to machine precision; drift between the two indicates a regression
        // in the classification logic.
        var lens = new F88PopcountPairLens(N, np, nq);
        Assert.Equal(lens.AlphaVerifier, lens.Alpha, precision: 14);
    }

    [Fact]
    public void HdEqualsN_Pi2OddIsZero()
    {
        // The Π²-classical anchor: GHZ-like (HD=N) collapses Π²-odd memory to zero.
        var lens = new F88PopcountPairLens(N: 4, np: 1, nq: 3);
        Assert.Equal(0.0, lens.Pi2OddInMemory(hd: 4), precision: 12);
    }

    [Fact]
    public void C2N6_Pi2OddIsHalfForBothHdClasses()
    {
        // The c=2 stratum at N=6 (n_p=1, n_q=2) is generic α=0.5; F88 lens predicts
        // Π²-odd memory ≈ 0.5 for both HD=1 and HD=3 contributions. This is the "facts"
        // surface that downstream Q_peak / orbit work needs. The lens does NOT itself
        // predict bond-position — it pins what F88 says about the c=2 coherence block.
        var lens = new F88PopcountPairLens(N: 6, np: 1, nq: 2);
        double pi2OddHd1 = lens.Pi2OddInMemory(hd: 1);
        double pi2OddHd3 = lens.Pi2OddInMemory(hd: 3);
        Assert.Equal(pi2OddHd1, pi2OddHd3, precision: 12);
        Assert.True(pi2OddHd1 > 0.49 && pi2OddHd1 < 0.51,
            $"c=2 N=6 generic Π²-odd memory expected ≈ 0.5, got {pi2OddHd1}");
    }

    [Theory]
    [InlineData(1)]
    [InlineData(0)]
    public void NLessThanTwo_Throws(int N)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new F88PopcountPairLens(N, np: 0, nq: 0));
    }

    [Theory]
    [InlineData(5, -1, 2)]
    [InlineData(5, 1, 6)]
    public void OutOfRangePopcount_Throws(int N, int np, int nq)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new F88PopcountPairLens(N, np, nq));
    }

    [Fact]
    public void Reconnaissance_C2_Pi2OddAcrossN5To8()
    {
        // Survey the c=2 stratum (n_p=1, n_q=2) at the empirical N range. Documents
        // F88's pin on each (N, hd) entry as a Schicht-1-fact for downstream consumers.
        _out.WriteLine("  N | n_p | n_q | kind     | α       | s       | π²-odd(HD=1) | π²-odd(HD=3)");
        _out.WriteLine("  --|-----|-----|----------|---------|---------|--------------|--------------");
        foreach (int N in new[] { 5, 6, 7, 8 })
        {
            var lens = new F88PopcountPairLens(N, np: 1, nq: 2);
            double hd1 = lens.Pi2OddInMemory(hd: 1);
            double hd3 = lens.Pi2OddInMemory(hd: 3);
            _out.WriteLine($"  {N} |   1 |   2 | {lens.Kind,-8} | {lens.Alpha,7:F4} | {lens.StaticFraction,7:F4} | {hd1,12:F4} | {hd3,12:F4}");
        }
    }
}
