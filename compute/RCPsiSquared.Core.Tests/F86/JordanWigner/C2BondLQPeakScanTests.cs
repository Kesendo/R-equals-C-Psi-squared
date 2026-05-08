using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.JordanWigner;

public class C2BondLQPeakScanTests
{
    private readonly ITestOutputHelper _out;

    public C2BondLQPeakScanTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void F71MirrorInvariance_QPeak_BitExact(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var scan = C2BondLQPeakScan.Build(block);

        Assert.True(scan.MaxF71MirrorDeviationQPeak < 1e-8,
            $"N={N}: F71 Q_peak deviation {scan.MaxF71MirrorDeviationQPeak:G3}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void F71MirrorInvariance_KMax_BitExact(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var scan = C2BondLQPeakScan.Build(block);

        Assert.True(scan.MaxF71MirrorDeviationKMax < 1e-8,
            $"N={N}: F71 KMax deviation {scan.MaxF71MirrorDeviationKMax:G3}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void EndpointKMax_LessThan_InnermostKMax(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var scan = C2BondLQPeakScan.Build(block);

        var endpoint = scan.Bonds.First(b => b.BondClass == BondClass.Endpoint);
        var innermost = ClosestToCenter(scan.Bonds, N);
        if (innermost is null || innermost.Bond == endpoint.Bond) return;

        Assert.True(endpoint.XbNormAtPeak < innermost.XbNormAtPeak,
            $"N={N}: Endpoint ‖xB‖_max={endpoint.XbNormAtPeak:F4} should be < Innermost ‖xB‖_max={innermost.XbNormAtPeak:F4}");
    }

    private static BondLQPeakWitness? ClosestToCenter(IReadOnlyList<BondLQPeakWitness> bonds, int N)
    {
        double center = (N - 2) / 2.0;
        BondLQPeakWitness? best = null;
        double bestDist = double.PositiveInfinity;
        foreach (var b in bonds)
        {
            if (b.BondClass != BondClass.Interior) continue;
            double d = Math.Abs(b.Bond - center);
            if (d < bestDist) { bestDist = d; best = b; }
        }
        return best;
    }

    [Theory]
    [InlineData(5)]
    [InlineData(8)]
    public void QPeak_WithinDefaultGridRange(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var scan = C2BondLQPeakScan.Build(block);

        foreach (var b in scan.Bonds)
        {
            Assert.InRange(b.QPeak, scan.QGrid[0], scan.QGrid[^1]);
        }
    }

    [Fact]
    public void Build_RejectsNonC2Block()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);
        Assert.Throws<ArgumentException>(() => C2BondLQPeakScan.Build(block));
    }

    [Fact]
    public void Build_RejectsTooSmallQGrid()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        Assert.Throws<ArgumentException>(() => C2BondLQPeakScan.Build(block, new[] { 1.0, 2.0 }));
    }

    [Fact]
    public void Tier_IsTier2Verified()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var scan = C2BondLQPeakScan.Build(block);
        Assert.Equal(Tier.Tier2Verified, scan.Tier);
    }

    [Fact]
    public void Anchor_References_PROOF_F86_QPEAK_AndDirectionB()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var scan = C2BondLQPeakScan.Build(block);
        Assert.Contains("PROOF_F86_QPEAK", scan.Anchor);
        Assert.Contains("Direction (b'')", scan.Anchor);
    }

    // The big reconnaissance: per bond, Q_peak / KMax / HWHM_left / HWHM_left/Q_peak
    // across N=5..8. Then comparison anchor: C2HwhmRatio's empirical Q_peak / HWHM
    // are around Q_peak ≈ 2.2 with HwhmLeft/Qpeak ≈ 0.75 (Endpoint) / 0.75 (Interior).
    [Fact]
    public void Reconnaissance_EmitsPerBondQPeakHwhm_AcrossN5To8()
    {
        _out.WriteLine("  N | b | class    | Q_peak  | ‖xB‖_max | HWHM_left | HWHM_left/Q_peak");
        _out.WriteLine("  --|---|----------|---------|-----------|-----------|-------------------");
        foreach (int N in new[] { 5, 6, 7, 8 })
        {
            var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
            var scan = C2BondLQPeakScan.Build(block);
            foreach (var b in scan.Bonds)
            {
                string hwhmLeftStr = b.HwhmLeft?.ToString("F4") ?? "n/a";
                string hwhmRatioStr = b.HwhmLeftOverQPeak?.ToString("F4") ?? "n/a";
                _out.WriteLine(
                    $"  {N} | {b.Bond} | {b.BondClass,-8} | {b.QPeak,7:F4} | " +
                    $"{b.XbNormAtPeak,9:F4} | {hwhmLeftStr,9} | {hwhmRatioStr,17}");
            }
            _out.WriteLine($"  -- N={N}: F71 mirror Q_peak dev = {scan.MaxF71MirrorDeviationQPeak:G3}, " +
                           $"KMax dev = {scan.MaxF71MirrorDeviationKMax:G3}");
            _out.WriteLine("");
        }
    }
}
