using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The typed registry of (1,2)-block count-change loci used by the shell census
/// (gate RealSeedCensusTests / SLOW_SEEDCENSUS; table in experiments/F89_PATH_K_DIABOLIC.md), with the
/// classification it records recomputed from below at an N=11 locus, where it stands without a
/// certificate.</summary>
public class RealDefectiveSeedsTests
{
    private readonly ITestOutputHelper _out;
    public RealDefectiveSeedsTests(ITestOutputHelper o) => _out = o;

    [Fact]
    public void Counts_MatchTheCountChangeCensus()
    {
        Assert.Equal(4, RealDefectiveSeeds.ForN(5).Count());
        Assert.Equal(6, RealDefectiveSeeds.ForN(7).Count());
        Assert.Equal(7, RealDefectiveSeeds.ForN(9).Count());
        Assert.Equal(9, RealDefectiveSeeds.ForN(11).Count());
    }

    [Fact]
    public void N5Anchors_MatchTheGateValues()
    {
        var n5 = RealDefectiveSeeds.ForN(5).ToList();
        Assert.Contains(n5, s => Math.Abs(s.QStar - 0.620878) < 1e-6 && s.RParity == +1);
        Assert.Contains(n5, s => Math.Abs(s.QStar - 1.077615) < 1e-6 && s.RParity == +1);
        Assert.Contains(n5, s => Math.Abs(s.QStar - 2.804888) < 1e-6 && s.RParity == -1);
        Assert.Contains(n5, s => Math.Abs(s.QStar - 0.643037) < 1e-6 && s.RParity == -1);
    }

    [Fact]
    public void EveryRecordedLocus_HasRealLambdaInTheOneTwoWindow()
    {
        // the (1,2) block's Bendixson window is [-6,-2] for every N (n_diff in {1,3});
        // the window-edge lemma sharpens it to the OPEN interval for certified seeds;
        // every recorded count-change locus also lies inside that same numerical window.
        Assert.All(RealDefectiveSeeds.All, s =>
        {
            Assert.True(s.LambdaA > -6.0 && s.LambdaA < -2.0);
            Assert.True(s.QStar > 0.2 && s.QStar < 3.0);   // the census window
            Assert.True(s.RParity is +1 or -1);
        });
    }

    /// <summary>The N=11 classification, recomputed where it is stored. At each recorded locus the
    /// colliding pair of the locus's R-parity sector (R-odd 300-dimensional, R-even 305: 605 = 5
    /// reflection-fixed + 2·300) is located by bisection on the sign of f(q) = (λ_a − λ_b)² (real at real
    /// q, the sector being self-conjugate: positive where the pair is real, negative where it is a
    /// conjugate pair), then read by the Riesz compression (alg, geo, kind) and by the gap exponent along
    /// real q. The record must hold what the computation returns, the stored exponent ½ included.
    ///
    /// <para>The exponent window is the one tolerance here and it has an error law: f is analytic with a
    /// simple zero at q*, so gap(δ) = √|f′δ|·(1 + (f″/4f′)·δ + O(δ²)), and the slope fitted over
    /// δ ∈ [1e-6, 1e-3] departs from ½ by the O(δ) term (at most ~10⁻³ here) plus the √ε eigenvalue floor
    /// over the smallest gap (~10⁻⁴). ±0.01 sits well above both and well below the distance to the
    /// alternatives a count drop could otherwise show, ⅓ (a 3×3 Jordan block) and 1 (a semisimple
    /// crossing, linear splitting).</para>
    ///
    /// <para>The fast fact runs the R-odd locus q* = 2.700712 (~10 s); the theory runs all nine
    /// (~90 s, Category SLOW_SEEDCHARACTER).</para></summary>
    [Fact]
    [Trait("Category", "SEEDCHARACTER")]
    public void N11_ROddLocus_2_700712_IsClassifiedDefectiveEp2() => AssertClassifiedDefectiveEp2(2.700712, -1);

    [Theory]
    [Trait("Category", "SLOW_SEEDCHARACTER")]
    [InlineData(0.502123, +1)]
    [InlineData(0.598253, +1)]
    [InlineData(0.812457, +1)]
    [InlineData(1.902993, +1)]
    [InlineData(0.587252, -1)]
    [InlineData(0.631587, -1)]
    [InlineData(0.727968, -1)]
    [InlineData(1.010462, -1)]
    [InlineData(2.700712, -1)]
    public void N11_EveryLocus_IsClassifiedDefectiveEp2(double recordedQ, int rParity) =>
        AssertClassifiedDefectiveEp2(recordedQ, rParity);

    [Fact]
    public void N11_TheoryCoversEveryRecordedLocus()
    {
        // the theory's rows are the registry's N=11 rows, no more and no fewer
        var rows = typeof(RealDefectiveSeedsTests).GetMethod(nameof(N11_EveryLocus_IsClassifiedDefectiveEp2))!
            .GetCustomAttributes(typeof(InlineDataAttribute), false).Cast<InlineDataAttribute>()
            .Select(a => a.GetData(null!).Single()).Select(d => ((double)d[0], (int)d[1])).ToHashSet();
        Assert.Equal(RealDefectiveSeeds.ForN(11).Select(s => (s.QStar, s.RParity)).ToHashSet(), rows);
    }

    private void AssertClassifiedDefectiveEp2(double recordedQ, int rParity)
    {
        var record = RealDefectiveSeeds.ForN(11).Single(s => s.QStar == recordedQ && s.RParity == rParity);
        Assert.NotNull(record.Classified);
        var c = record.Classified!.Value;
        var anchor = new Complex(record.LambdaA, 0.0);
        bool odd = rParity < 0;

        var probe = Sector(record.QStar, odd);
        Assert.Equal(odd ? 300 : 305, probe.RowCount);

        // bisection on sign f, bracketed by the recorded 6-decimal q* ± 2e-6
        double lo = record.QStar - 2e-6, hi = record.QStar + 2e-6;
        double fLo = F(lo, anchor, odd), fHi = F(hi, anchor, odd);
        Assert.True(fLo * fHi < 0, $"f must change sign across the recorded locus: f({lo})={fLo:E3}, f({hi})={fHi:E3}");
        for (int it = 0; it < 80 && hi - lo > 0; it++)
        {
            double mid = 0.5 * (lo + hi);
            if (mid <= lo || mid >= hi) break;
            double fm = F(mid, anchor, odd);
            if (Math.Sign(fm) == Math.Sign(fLo)) { lo = mid; fLo = fm; } else hi = mid;
        }
        double qStar = 0.5 * (lo + hi);
        var m = Sector(qStar, odd);
        var (a, b) = Pair(m, anchor);
        Complex lambdaStar = 0.5 * (a + b);
        _out.WriteLine($"q* = {qStar:R}, lambda* = {lambdaStar.Real:R}{lambdaStar.Imaginary:+0.0E0;-0.0E0}i, pair gap {(a - b).Magnitude:E2}");

        // the recorded values are this locus to their last printed digit (half a unit in that place)
        Assert.True(Math.Abs(qStar - record.QStar) <= 5e-7);
        Assert.True(Math.Abs(lambdaStar.Real - record.LambdaA) <= 5e-5);

        // character: one Jordan pair, alg 2 and geo 1, Riesz circle inside the gap to the third eigenvalue
        double third = m.Evd().EigenValues.Select(z => (z - lambdaStar).Magnitude).OrderBy(x => x).ElementAt(2);
        var reading = EpCharacter.Characterize(m, lambdaStar, 0.4 * third);
        _out.WriteLine($"third eigenvalue at {third:F3}; kind {reading.Kind}, alg {reading.Algebraic}, geo {reading.Geometric}, dep {reading.Departure:F3}");
        Assert.Equal(EpCharacter.EpKind.Defective, reading.Kind);
        Assert.Equal(2, reading.Algebraic);
        Assert.Equal(1, reading.Geometric);

        // Puiseux exponent on both sides of q*
        foreach (int side in new[] { -1, +1 })
        {
            double g3 = Gap(qStar + side * 1e-3, anchor, odd), g6 = Gap(qStar + side * 1e-6, anchor, odd);
            double p = Math.Log10(g3 / g6) / 3.0;
            _out.WriteLine($"side {side:+0;-0}: gap(1e-3) = {g3:E3}, gap(1e-6) = {g6:E3}, exponent {p:F5}");
            Assert.InRange(p, c.PuiseuxExponent - 0.01, c.PuiseuxExponent + 0.01);
        }

        Assert.Equal(reading.Kind, c.Kind);
        Assert.Equal(reading.Algebraic, c.Algebraic);
        Assert.Equal(reading.Geometric, c.Geometric);
    }

    private static Matrix<Complex> Sector(double q, bool odd)
    {
        var (a, d) = WeightCoherenceBlock.BuildReflectionSectorColumnMajor(11, 1, 2, new Complex(q, 0.0), odd);
        return Matrix<Complex>.Build.Dense(d, d, a);   // column-major storage
    }

    private static (Complex A, Complex B) Pair(Matrix<Complex> m, Complex anchor)
    {
        var two = m.Evd().EigenValues.OrderBy(z => (z - anchor).Magnitude).Take(2).ToArray();
        return (two[0], two[1]);
    }

    private static double F(double q, Complex anchor, bool odd)
    {
        var (a, b) = Pair(Sector(q, odd), anchor);
        return ((a - b) * (a - b)).Real;
    }

    private static double Gap(double q, Complex anchor, bool odd)
    {
        var (a, b) = Pair(Sector(q, odd), anchor);
        return (a - b).Magnitude;
    }
}
