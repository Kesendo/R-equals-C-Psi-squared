using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.JordanWigner;

public class PerF71OrbitLQPeakTableTests
{
    private readonly ITestOutputHelper _out;

    public PerF71OrbitLQPeakTableTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void F71MirrorInvariance_QPeak_BelowTolerance(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitLQPeakTable.Build(block);

        Assert.True(table.MaxQPeakWithinOrbitDeviation < PerF71OrbitLQPeakTable.F71MirrorQPeakTolerance,
            $"N={N}: Q_peak within-orbit deviation {table.MaxQPeakWithinOrbitDeviation:G3} " +
            $"exceeds tolerance {PerF71OrbitLQPeakTable.F71MirrorQPeakTolerance:G3}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void F71MirrorInvariance_KMax_BelowTolerance(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitLQPeakTable.Build(block);

        Assert.True(table.MaxKMaxWithinOrbitDeviation < PerF71OrbitLQPeakTable.F71MirrorKMaxTolerance,
            $"N={N}: KMax within-orbit deviation {table.MaxKMaxWithinOrbitDeviation:G3} " +
            $"exceeds tolerance {PerF71OrbitLQPeakTable.F71MirrorKMaxTolerance:G3}");
    }

    [Theory]
    [InlineData(5, 2)]   // 4 bonds → 2 orbit pairs (no self-paired)
    [InlineData(6, 3)]   // 5 bonds → 2 orbit pairs + 1 self-paired
    [InlineData(7, 3)]   // 6 bonds → 3 orbit pairs (no self-paired)
    [InlineData(8, 4)]   // 7 bonds → 3 orbit pairs + 1 self-paired
    public void OrbitWitnesses_Count_MatchesF71Decomposition(int N, int expectedOrbitCount)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitLQPeakTable.Build(block);

        Assert.Equal(expectedOrbitCount, table.OrbitWitnesses.Count);
    }

    [Theory]
    [InlineData(5, false)]  // 4 bonds (even) → no self-paired
    [InlineData(6, true)]   // 5 bonds (odd) → self-paired at b=2
    [InlineData(7, false)]  // 6 bonds (even) → no self-paired
    [InlineData(8, true)]   // 7 bonds (odd) → self-paired at b=3
    public void SelfPairedCentralOrbit_PresentForOddNumBonds(int N, bool expectedSelfPaired)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitLQPeakTable.Build(block);
        bool hasSelfPaired = table.OrbitWitnesses.Any(w => w.Orbit.IsSelfPaired);

        Assert.Equal(expectedSelfPaired, hasSelfPaired);
    }

    [Fact]
    public void Build_RejectsNonC2Block()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);  // c=3
        Assert.Throws<ArgumentException>(() => PerF71OrbitLQPeakTable.Build(block));
    }

    [Fact]
    public void Tier_IsTier2Verified()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitLQPeakTable.Build(block);
        Assert.Equal(Tier.Tier2Verified, table.Tier);
    }

    [Fact]
    public void Anchor_References_PROOF_F86_QPEAK_AndC1Mirror()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitLQPeakTable.Build(block);
        Assert.Contains("PROOF_F86_QPEAK", table.Anchor);
        Assert.Contains("PROOF_C1_MIRROR_SYMMETRY", table.Anchor);
    }

    // The orbit-level Q_peak / KMax / HWHM structure across N=5..8 is the EP-resonance
    // counterpart of PerF71OrbitKModeTable's K-mode reconnaissance. Emit per-N orbit
    // Q_peak / KMax / HWHM_left so the Endpoint vs Interior split (load-bearing for the
    // BareDoubledPtfXPeak ≈ 2.197 floor + g_eff_probe ratio) is visible across the
    // transition.
    [Fact]
    public void Reconnaissance_EmitsPerNOrbitLQPeak_AcrossN5To8()
    {
        _out.WriteLine("  N | orbit                       | Q_peak  | KMax    | HWHM_L  | HWHM_L/Q_peak");
        _out.WriteLine("  --|-----------------------------|---------|---------|---------|---------------");
        foreach (int N in new[] { 5, 6, 7, 8 })
        {
            var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
            var table = PerF71OrbitLQPeakTable.Build(block);
            Assert.Equal((N - 1 + 1) / 2, table.OrbitWitnesses.Count);
            foreach (var w in table.OrbitWitnesses)
            {
                string hwhmL = w.HwhmLeft?.ToString("F4") ?? "n/a";
                string hwhmRatio = w.HwhmLeftOverQPeak?.ToString("F4") ?? "n/a";
                _out.WriteLine($"  {N} | {w.Orbit.DisplayName,-27} | {w.QPeak,7:F4} | {w.KMax,7:F4} | {hwhmL,7} | {hwhmRatio,13}");
            }
            _out.WriteLine($"  -- N={N}: max-Δ Q_peak within orbit = {table.MaxQPeakWithinOrbitDeviation:G3}, " +
                           $"max-Δ KMax within orbit = {table.MaxKMaxWithinOrbitDeviation:G3}");
            _out.WriteLine("");
        }
    }
}
