using System;
using RCPsiSquared.Core.F89PathK;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The certificate-grade D-device on the weight-general path
/// (<see cref="FoldResultantCertificate.WeightCertifyDiscMultiplicity"/>): the (1,3)@N=6 scout's
/// deg/v/layers upgraded from two-prime agreement to the lc-divisor-bound certification, the same
/// proof shape as the committed (1,2) device. Controls: the weight path must reproduce the (1,2)
/// certificate values at N=4 (52/24/[20,4], ~90 ms, Category EVEN_DISC13_SCOUT) and N=6
/// (926/536/[124,133], ~11 s per parity, 499 primes, Category SLOW_DISC13_CERT). All pins
/// measured 2026-08-11. Run: dotnet test compute/RCPsiSquared.Core.Tests -c Release --filter
/// "Category=SLOW_DISC13_CERT|Category=EVEN_DISC13_SCOUT" (~3 min; the scout category also
/// carries Disc13LayerScoutTests' fast members).</summary>
public class Disc13CertificateTests
{
    private readonly ITestOutputHelper _out;
    public Disc13CertificateTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(4, 1, 2, false)]
    [InlineData(4, 1, 2, true)]
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void N4_Control_WeightCertificateMatchesCommitted(int n, int wKet, int wBra, bool rOdd)
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var r = FoldResultantCertificate.WeightCertifyDiscMultiplicity(n, wKet, wBra, rOdd);
        _out.WriteLine($"N{n} ({wKet},{wBra}) rOdd={rOdd}: deg={r.TrueDiscriminantDegree} v={r.TrueQValuationD} " +
                       $"layers=[{string.Join(", ", r.DiscLayerDegrees)}] mD={r.InfinityRepeatedD} dBound={r.DiscriminantDegreeBound} " +
                       $"lcBound={r.LcDivisorBoundD} sampled={r.PrimesSampled} certified={r.DiscLayersCertified}  [{sw.Elapsed}]");
        Assert.True(r.DiscLayersCertified);
        Assert.Equal(52, r.TrueDiscriminantDegree);
        Assert.Equal(24, r.TrueQValuationD);
        Assert.Equal(new[] { 20, 4 }, r.DiscLayerDegrees);
    }

    [Theory]
    [InlineData(6, 1, 2, false)]
    [InlineData(6, 1, 2, true)]
    [Trait("Category", "SLOW_DISC13_CERT")]
    public void N6_Control_WeightCertificateMatchesCommitted12(int n, int wKet, int wBra, bool rOdd)
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var r = FoldResultantCertificate.WeightCertifyDiscMultiplicity(n, wKet, wBra, rOdd);
        _out.WriteLine($"N{n} ({wKet},{wBra}) rOdd={rOdd}: deg={r.TrueDiscriminantDegree} v={r.TrueQValuationD} " +
                       $"layers=[{string.Join(", ", r.DiscLayerDegrees)}] mD={r.InfinityRepeatedD} dBound={r.DiscriminantDegreeBound} " +
                       $"lcBound={r.LcDivisorBoundD} sampled={r.PrimesSampled} certified={r.DiscLayersCertified}  [{sw.Elapsed}]");
        Assert.True(r.DiscLayersCertified);
        Assert.Equal(926, r.TrueDiscriminantDegree);
        Assert.Equal(536, r.TrueQValuationD);
        Assert.Equal(new[] { 124, 133 }, r.DiscLayerDegrees);
    }

    /// <summary>The (1,3)@N=6 certificate, both parities (measured 2026-08-11: R-even ~44 s at 924
    /// sampled primes past the lc-divisor bound 923; R-odd ~98 s at 1210 past 1209). The scout's
    /// two-prime deg/v/layers upgraded to certified invariants: TrueDiscriminantDegree and
    /// TrueQValuationD are what the arc's Sturm closing step will lean on, exactly as the (1,2)
    /// gcd certificate leaned on the (1,2) D-device's. The leading-form m_D (51 / 42) sharpens the
    /// interpolation bound to 1620 / 2172 (the scouts use the crude 1722 / 2256); the certified
    /// layer profile reproduces the scout exactly, and MaxDiscMultiplicity = 4 bounds the true
    /// maximum from above (the mod-p lift is one-way upward; "Klein-forced ceiling" is the census
    /// sections' interpretation of the value, not something this method certifies). NOTE the
    /// scope the (1,2) device's consumers must not import: its MaxDiscMultiplicity ≤ 2 IS the
    /// β-exotic exclusion, and at (1,3) that exclusion does NOT fire, a Puiseux-3/2 point being
    /// hideable inside the multiplicity-4 layer.</summary>
    [Theory]
    [InlineData(6, 1, 3, false)]
    [InlineData(6, 1, 3, true)]
    [Trait("Category", "SLOW_DISC13_CERT")]
    public void N6_Certificate_Disc13(int n, int wKet, int wBra, bool rOdd)
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var r = FoldResultantCertificate.WeightCertifyDiscMultiplicity(n, wKet, wBra, rOdd,
            log: s => _out.WriteLine("  " + s));
        _out.WriteLine($"N{n} ({wKet},{wBra}) rOdd={rOdd}: deg={r.TrueDiscriminantDegree} v={r.TrueQValuationD} " +
                       $"layers=[{string.Join(", ", r.DiscLayerDegrees)}] mD={r.InfinityRepeatedD} dBound={r.DiscriminantDegreeBound} " +
                       $"lcBound={r.LcDivisorBoundD} sampled={r.PrimesSampled} certified={r.DiscLayersCertified}  [{sw.Elapsed}]");
        Assert.True(r.DiscLayersCertified);
        Assert.Equal(rOdd ? 2124 : 1572, r.TrueDiscriminantDegree);
        Assert.Equal(rOdd ? 1188 : 876, r.TrueQValuationD);
        Assert.Equal(rOdd ? new[] { 0, 132, 0, 168 } : new[] { 0, 108, 0, 120 }, r.DiscLayerDegrees);
        Assert.Equal(rOdd ? 42 : 51, r.InfinityRepeatedD);
        Assert.Equal(rOdd ? 2172 : 1620, r.DiscriminantDegreeBound);
        Assert.Equal(4, r.MaxDiscMultiplicity);
        Assert.True(r.PrimesSampled > r.LcDivisorBoundD);
    }
}
