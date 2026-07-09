using System.Diagnostics;
using RCPsiSquared.Core.F89PathK;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The D-only β-exotic certificate (<c>CertifyDiscMultiplicity</c>), which drops the corner
/// block, the resultant and the Mignotte lift that make <c>CertifyComplete</c> unaffordable past N = 5.
///
/// <para>The load-bearing test is the AGREEMENT one: at N = 5 the new path must reproduce, exactly, the
/// certified layer degrees, deg_q D, v_q(D) and MaxDiscMultiplicity that the full R1 certificate
/// already reports. A cheaper route to the same statement is worth nothing until it is pinned against
/// the expensive route on the case both can run.</para></summary>
public class DiscMultiplicityCertificateTests
{
    private readonly ITestOutputHelper _out;
    public DiscMultiplicityCertificateTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [Trait("Category", "DISCMULT")]
    [InlineData(false)]
    [InlineData(true)]
    public void N5_DiscOnlyPath_ReproducesTheFullCertificate_ExactlyAndFaster(bool rOdd)
    {
        var sw = Stopwatch.StartNew();
        var full = FoldResultantCertificate.CertifyComplete(5, rOdd);
        long fullMs = sw.ElapsedMilliseconds;

        sw.Restart();
        var disc = FoldResultantCertificate.CertifyDiscMultiplicity(5, rOdd);
        long discMs = sw.ElapsedMilliseconds;

        _out.WriteLine($"N=5 rOdd={rOdd}: full R1 certificate {fullMs} ms, D-only {discMs} ms " +
                       $"({full.PrimesSampled} vs {disc.PrimesSampled} primes sampled, " +
                       $"lc-divisor bound for D = {disc.LcDivisorBoundD})");
        _out.WriteLine($"   layers [{string.Join(", ", disc.DiscLayerDegrees)}], deg_q D = " +
                       $"{disc.TrueDiscriminantDegree}, v_q(D) = {disc.TrueQValuationD}, " +
                       $"max multiplicity {disc.MaxDiscMultiplicity}, certified = {disc.DiscLayersCertified}");

        Assert.True(disc.DiscLayersCertified);
        Assert.Equal(full.ResidualDegree, disc.ResidualDegree);
        Assert.Equal(full.AtDegree, disc.AtDegree);
        Assert.Equal(full.DiscriminantDegreeBound, disc.DiscriminantDegreeBound);
        Assert.Equal(full.TrueDiscriminantDegree, disc.TrueDiscriminantDegree);
        Assert.Equal(full.TrueQValuationD, disc.TrueQValuationD);
        Assert.Equal(full.LcDivisorBoundD, disc.LcDivisorBoundD);
        Assert.Equal(full.DiscLayerDegrees, disc.DiscLayerDegrees);
        Assert.Equal(full.MaxDiscMultiplicity, disc.MaxDiscMultiplicity);

        // the statement itself
        Assert.True(disc.MaxDiscMultiplicity <= 2,
            $"a β-exotic (Puiseux exponent 3/2) would force a disc root of multiplicity ≥ 3; got " +
            $"{disc.MaxDiscMultiplicity}");
    }

    /// <summary>The whole verdict is read off ONE prime, the first the certificate samples. A second,
    /// independent prime must report the same deg_q D, the same v_q(D) and the same layer degrees. A
    /// reduction, degree or single-prime-fluke bug in the layer reading dies here.</summary>
    [Theory]
    [Trait("Category", "DISCMULT")]
    [InlineData(5, false)]
    [InlineData(5, true)]
    public void TheLayerReading_DoesNotDependOnWhichPrime(int n, bool rOdd)
    {
        var first = FoldResultantCertificate.DiscLayersAtNthPrime(n, rOdd, nth: 0);
        var second = FoldResultantCertificate.DiscLayersAtNthPrime(n, rOdd, nth: 1);
        var third = FoldResultantCertificate.DiscLayersAtNthPrime(n, rOdd, nth: 7);

        _out.WriteLine($"n={n} rOdd={rOdd}: prime #0 -> deg {first.DegD}, v_q {first.VD}, " +
                       $"layers [{string.Join(", ", first.LayerDegrees)}]");
        _out.WriteLine($"n={n} rOdd={rOdd}: prime #1 -> deg {second.DegD}, v_q {second.VD}, " +
                       $"layers [{string.Join(", ", second.LayerDegrees)}]");
        _out.WriteLine($"n={n} rOdd={rOdd}: prime #7 -> deg {third.DegD}, v_q {third.VD}, " +
                       $"layers [{string.Join(", ", third.LayerDegrees)}]");

        Assert.Equal(first.DegD, second.DegD);
        Assert.Equal(first.VD, second.VD);
        Assert.Equal(first.LayerDegrees, second.LayerDegrees);
        Assert.Equal(first.DegD, third.DegD);
        Assert.Equal(first.VD, third.VD);
        Assert.Equal(first.LayerDegrees, third.LayerDegrees);
        Assert.True(first.LayerDegrees.Length <= 2);
    }

    /// <summary>The same falsifier at n = 7, where a single prime carries the whole result. Slow only
    /// because the bivariate build is ~15 s; the extra primes are seconds each.</summary>
    [Theory]
    [Trait("Category", "SLOW_DISCMULT_N7")]
    [InlineData(true)]
    [InlineData(false)]
    public void N7_TheLayerReading_DoesNotDependOnWhichPrime(bool rOdd)
    {
        var first = FoldResultantCertificate.DiscLayersAtNthPrime(7, rOdd, nth: 0);
        var second = FoldResultantCertificate.DiscLayersAtNthPrime(7, rOdd, nth: 3);

        _out.WriteLine($"N=7 rOdd={rOdd}: prime #0 -> deg {first.DegD}, v_q {first.VD}, " +
                       $"layers [{string.Join(", ", first.LayerDegrees)}]");
        _out.WriteLine($"N=7 rOdd={rOdd}: prime #3 -> deg {second.DegD}, v_q {second.VD}, " +
                       $"layers [{string.Join(", ", second.LayerDegrees)}]");

        Assert.Equal(first.DegD, second.DegD);
        Assert.Equal(first.VD, second.VD);
        Assert.Equal(first.LayerDegrees, second.LayerDegrees);
        Assert.Equal(2, first.LayerDegrees.Length);
        // the identity that would catch any transcription or aggregation error in the four integers
        Assert.Equal(first.DegD - first.VD, first.LayerDegrees[0] + 2 * first.LayerDegrees[1]);
    }

    /// <summary>Sizing only: read lcDivisorBoundD and dBound at n = 7 without running the prime loop
    /// (maxPrimes = 0), so the cost of a real n = 7 run is known before it is paid.</summary>
    [Theory]
    [Trait("Category", "DISCMULT")]
    [InlineData(false)]
    [InlineData(true)]
    public void N7_Sizing_IsReadBeforeAnyRun(bool rOdd)
    {
        var sw = Stopwatch.StartNew();
        var r = FoldResultantCertificate.CertifyDiscMultiplicity(7, rOdd, maxPrimes: 0);
        _out.WriteLine($"N=7 rOdd={rOdd}: block {r.BlockDimension} = AT {r.AtDegree} ⊎ residual " +
                       $"{r.ResidualDegree}; m_D = {r.InfinityRepeatedD}, dBound = {r.DiscriminantDegreeBound}, " +
                       $"lc-divisor bound for D = {r.LcDivisorBoundD}; bivariate build + sizing in " +
                       $"{sw.Elapsed.TotalSeconds:F1} s");

        Assert.False(r.DiscLayersCertified);        // no primes sampled: must fail closed
        Assert.Equal(0, r.PrimesSampled);
        Assert.True(r.ResidualDegree is 52 or 53);
        Assert.True(r.DiscriminantDegreeBound > 2000);
    }

    /// <summary>The real thing at n = 7: ~2600 interpolation nodes per prime against an lc-divisor bound
    /// of ~1350 primes, no corner block, no resultant, no Mignotte lift. Tagged SLOW because it is minutes,
    /// not seconds; the sizing test above is what a fast gate reads.</summary>
    [Theory]
    [Trait("Category", "SLOW_DISCMULT_N7")]
    [InlineData(true)]
    [InlineData(false)]
    public void N7_DiscOnlyPath_ExcludesTheBetaExotic(bool rOdd)
    {
        var sw = Stopwatch.StartNew();
        var r = FoldResultantCertificate.CertifyDiscMultiplicity(7, rOdd,
            log: m => _out.WriteLine("   " + m));

        _out.WriteLine($"N=7 rOdd={rOdd}: residual degree {r.ResidualDegree}; layers " +
                       $"[{string.Join(", ", r.DiscLayerDegrees)}]; deg_q D = {r.TrueDiscriminantDegree} " +
                       $"(bound {r.DiscriminantDegreeBound}), v_q(D) = {r.TrueQValuationD}; layer prime " +
                       $"{r.LayerPrime}; {r.PrimesSampled} primes vs bound {r.LcDivisorBoundD}; " +
                       $"max multiplicity {r.MaxDiscMultiplicity}; certified = {r.DiscLayersCertified}; " +
                       $"{sw.Elapsed.TotalMinutes:F1} min");

        Assert.True(r.DiscLayersCertified, "the layer prime must attain trueDegD and the minimal v_q(D)");
        Assert.True(r.MaxDiscMultiplicity <= 2,
            $"disc_Λ(F_res) has a root of multiplicity {r.MaxDiscMultiplicity} ≥ 3 off q = 0: the β-exotic " +
            $"is NOT excluded at N = 7, rOdd={rOdd}");
    }
}
