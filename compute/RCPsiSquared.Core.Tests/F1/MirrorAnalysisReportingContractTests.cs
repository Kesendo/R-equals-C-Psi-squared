using RCPsiSquared.Compute;
using RCPsiSquared.Core.F1;
using Xunit;

namespace RCPsiSquared.Core.Tests.F1;

/// <summary><see cref="MirrorAnalysis"/> is a thin reporting wrapper: the pairing distance itself is
/// computed by <see cref="F1SpectrumStatistics.MaxF1RatePairingDistance"/> and is pinned beside the
/// canon in <see cref="F1SpectrumStatisticsTests"/>, including the λ = −d rate embedding. What the
/// wrapper adds on top of that call, and what had no C# test until 2026-08-06, is exactly three
/// things: the eps ratio, the rate radius it is divided by, and the empty-input policy. Those are
/// what is pinned here, from below and without a tolerance anywhere.
///
/// <para>The Python producer for the same surface is
/// <c>simulations/f1_rate_embedding_check.py</c>; it covers the embedding equality, which is the
/// part that lives in Core.</para></summary>
public class MirrorAnalysisReportingContractTests
{
    /// <summary>Machine epsilon reached WITHOUT reading the constant under test: the gap above 1.0
    /// is 2^−52 by definition of the double format. If <c>MirrorAnalysis.Eps</c> is ever edited to
    /// something else, the ratio pin below fails rather than silently rescaling every reported
    /// number. Not asserted to equal 2.220446049250313e-16 anywhere: that identity is a property
    /// of the double format and holds on every conforming runtime, so it could not fail.</summary>
    static readonly double MachineEps = Math.BitIncrement(1.0) - 1.0;

    // The centre is σ, the SUM of the per-site rates, and the reflection is d ↦ 2σ − d. Around
    // σ = 1 the multiset {0.5, 1.5} is its own reflection, so the distance is 0, and 0 exactly:
    // the reflection of a dyadic rate about a dyadic centre is dyadic, so there is nothing here
    // for the arithmetic to round and a tolerance would only hide a real defect.
    [Fact]
    public void APalindromicRateMultiset_ScoresExactlyZero()
    {
        var result = MirrorAnalysis.CheckSymmetry(new List<double> { 0.5, 1.5 }, center: 1.0);
        Assert.True(result.MaxPairingDistance == 0.0,
            $"{{0.5, 1.5}} is closed under d -> 2-d, so the distance must be exactly 0.0; "
            + $"got {result.MaxPairingDistance:E17}");
        Assert.True(result.DistanceInEps == 0.0, $"eps ratio {result.DistanceInEps:E17}");
    }

    // The break-input for the row above, and the fixture the ratio pins read: {0.5, 1.0} about
    // σ = 1 reflects to {1.5, 1.0}, whose nearest-neighbour matching tops out at 0.5.
    // Hand-computed, not read off a run.
    //
    // NOT a removal discriminator, and it would be easy to write it up as one: without removal
    // the second rate pairs onto 1.0 at distance 0, but the method returns only the MAX and that
    // is 0.5 either way. Removal is pinned where it belongs, on the matcher itself, by
    // F1SpectrumStatisticsTests.MaxF1PairingDistance_DetectsAMultiplicityDefect.
    [Fact]
    public void ANonPalindromicRateMultiset_ScoresTheForcedDistance()
    {
        var result = MirrorAnalysis.CheckSymmetry(new List<double> { 0.5, 1.0 }, center: 1.0);
        Assert.True(result.MaxPairingDistance == 0.5,
            $"{{0.5, 1.0}} reflects to {{1.5, 1.0}}, whose worst pairing is 0.5; "
            + $"got {result.MaxPairingDistance:E17}");
        Assert.Equal(2, result.Count);
        Assert.Equal(1.0, result.Center);
    }

    // The ratio is distance / (eps · radius), and on the fixture above every factor is a power of
    // two: 0.5 / (2^-52 · 1.0) = 2^51, which is an integer a double holds exactly. Both the
    // independently derived eps and the literal are asserted, so the pin fails whether
    // MirrorAnalysis.Eps drifts off machine epsilon or the radius stops being max|rate|; neither
    // route recomputes the expression under test.
    [Fact]
    public void TheEpsRatio_IsTheDistanceInUnitsOfEpsTimesTheRadius()
    {
        var result = MirrorAnalysis.CheckSymmetry(new List<double> { 0.5, 1.0 }, center: 1.0);
        Assert.True(result.RateRadius == 1.0, $"radius {result.RateRadius:E17}");
        Assert.True(result.DistanceInEps == 0.5 / MachineEps,
            $"the unit must be machine epsilon; got {result.DistanceInEps:E17}");
        Assert.True(result.DistanceInEps == 2251799813685248.0,
            $"0.5 / (2^-52 * 1.0) = 2^51 = 2251799813685248; got {result.DistanceInEps:E17}");
    }

    // The radius is the largest MAGNITUDE, not the largest value. Rates come out of
    // Liouvillian.ExtractRates as d = −Re(λ) and are non-negative there, so a sign error in this
    // line would never show up in the suite's own runs; the pin is what makes it visible.
    [Fact]
    public void TheRadius_IsTheLargestMagnitude_NotTheLargestValue()
    {
        var result = MirrorAnalysis.CheckSymmetry(new List<double> { -3.0, 1.0 }, center: 0.0);
        Assert.True(result.RateRadius == 3.0, $"radius {result.RateRadius:E17}");
    }

    // A radius of zero has no eps·ρ to divide by. The wrapper must report NaN rather than an
    // infinity, which would print as a number and read as an enormous violation.
    [Fact]
    public void AZeroRadius_LeavesTheRatioUndefined_RatherThanInfinite()
    {
        var result = MirrorAnalysis.CheckSymmetry(new List<double> { 0.0, 0.0 }, center: 0.0);
        Assert.True(result.RateRadius == 0.0, $"radius {result.RateRadius:E17}");
        Assert.True(result.MaxPairingDistance == 0.0, $"distance {result.MaxPairingDistance:E17}");
        Assert.True(double.IsNaN(result.DistanceInEps),
            $"a zero radius has no eps unit, so the ratio must be NaN and not infinite; "
            + $"got {result.DistanceInEps:E17}");
    }

    // The empty-input contract, and the reason it is worth a test of its own: the ONE answer no
    // route may give is 0.0, the best possible reading for a multiset that was never scored. Core
    // refuses by throwing; this wrapper reports NaN because its caller prints a table row and
    // should see a hole rather than a crash. Both halves are asserted together so the difference
    // is on the record as a deliberate split rather than an inconsistency waiting to be
    // "harmonised". (The third value in the family, fw.max_f1_pairing_distance returning 0.0, is
    // the open gap; it is Python and cannot be reached from here.)
    [Fact]
    public void AnEmptyRateMultiset_ReportsAHole_WhileCoreRefuses()
    {
        var result = MirrorAnalysis.CheckSymmetry(new List<double>(), center: 0.25);
        Assert.True(double.IsNaN(result.MaxPairingDistance),
            $"an unscored multiset must not report a distance; got {result.MaxPairingDistance:E17}");
        Assert.True(double.IsNaN(result.DistanceInEps), $"eps ratio {result.DistanceInEps:E17}");
        Assert.True(result.RateRadius == 0.0, $"radius {result.RateRadius:E17}");
        Assert.Equal(0, result.Count);
        Assert.Equal(0.25, result.Center);

        Assert.Throws<ArgumentException>(() =>
            F1SpectrumStatistics.MaxF1RatePairingDistance(new List<double>(), sigma: 0.25));
    }

    [Fact]
    public void ANullRateMultiset_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => MirrorAnalysis.CheckSymmetry(null!, center: 0.0));
    }
}
