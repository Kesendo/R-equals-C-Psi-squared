using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.JordanWigner;

public class PerF71OrbitKModeTableTests
{
    private readonly ITestOutputHelper _out;

    public PerF71OrbitKModeTableTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    [InlineData(10)]
    public void F71MirrorInvariance_K90_BitExact(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKModeTable.Build(block);

        Assert.Equal(0, table.MaxK90WithinOrbitDeviation);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    [InlineData(10)]
    public void F71MirrorInvariance_K99_BitExact(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKModeTable.Build(block);

        Assert.Equal(0, table.MaxK99WithinOrbitDeviation);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    [InlineData(10)]
    public void F71MirrorInvariance_RowL1_BelowFloatPrecisionTolerance(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKModeTable.Build(block);

        Assert.True(table.MaxRowL1WithinOrbitDeviation < PerF71OrbitKModeTable.F71MirrorRowL1Tolerance,
            $"N={N}: row-L1 within-orbit deviation {table.MaxRowL1WithinOrbitDeviation:G3} " +
            $"exceeds tolerance {PerF71OrbitKModeTable.F71MirrorRowL1Tolerance:G3}");
    }

    [Theory]
    [InlineData(5, 2)]   // 4 bonds → 2 orbit pairs (no self-paired)
    [InlineData(6, 3)]   // 5 bonds → 2 orbit pairs + 1 self-paired
    [InlineData(7, 3)]   // 6 bonds → 3 orbit pairs (no self-paired)
    [InlineData(8, 4)]   // 7 bonds → 3 orbit pairs + 1 self-paired
    [InlineData(9, 4)]   // 8 bonds → 4 orbit pairs (no self-paired)
    [InlineData(10, 5)]  // 9 bonds → 4 orbit pairs + 1 self-paired
    public void OrbitWitnesses_Count_MatchesF71Decomposition(int N, int expectedOrbitCount)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKModeTable.Build(block);

        Assert.Equal(expectedOrbitCount, table.OrbitWitnesses.Count);
    }

    [Theory]
    [InlineData(5, false)]  // 4 bonds (even) → no self-paired
    [InlineData(6, true)]   // 5 bonds (odd) → self-paired at b=2
    [InlineData(7, false)]  // 6 bonds (even) → no self-paired
    [InlineData(8, true)]   // 7 bonds (odd) → self-paired at b=3
    [InlineData(9, false)]  // 8 bonds (even) → no self-paired
    [InlineData(10, true)]  // 9 bonds (odd) → self-paired at b=4
    public void SelfPairedCentralOrbit_PresentForOddNumBonds(int N, bool expectedSelfPaired)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKModeTable.Build(block);
        bool hasSelfPaired = table.OrbitWitnesses.Any(w => w.Orbit.IsSelfPaired);

        Assert.Equal(expectedSelfPaired, hasSelfPaired);
    }

    [Fact]
    public void Build_RejectsNonC2Block()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);  // c=3
        Assert.Throws<ArgumentException>(() => PerF71OrbitKModeTable.Build(block));
    }

    [Fact]
    public void Tier_IsTier2Verified()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKModeTable.Build(block);
        Assert.Equal(Tier.Tier2Verified, table.Tier);
    }

    [Fact]
    public void Anchor_References_PROOF_F86_QPEAK_AndC1Mirror()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var table = PerF71OrbitKModeTable.Build(block);
        Assert.Contains("PROOF_F86_QPEAK", table.Anchor);
        Assert.Contains("PROOF_C1_MIRROR_SYMMETRY", table.Anchor);
    }

    // The orbit-level K_90 reorganisation N=5..8 is the precursor of the F71-orbit
    // Q_peak escape that becomes loud at N≥9 in the slow K-resonance lens. Emit per-N
    // orbit K_90 / K_99 / TopThreeK so the structural reorganisation is visible across
    // the transition.
    [Fact]
    public void Reconnaissance_EmitsPerNOrbitKMode_AcrossN5To10()
    {
        _out.WriteLine("  N | orbit                       | K_90 | K_99 | top-3 k indices");
        _out.WriteLine("  --|-----------------------------|------|------|----------------");
        foreach (int N in new[] { 5, 6, 7, 8, 9, 10 })
        {
            var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
            var table = PerF71OrbitKModeTable.Build(block);
            Assert.Equal((N - 1 + 1) / 2, table.OrbitWitnesses.Count);
            foreach (var w in table.OrbitWitnesses)
            {
                string topKStr = string.Join(",", w.TopThreeKIndices);
                _out.WriteLine($"  {N} | {w.Orbit.DisplayName,-27} | {w.K90,4} | {w.K99,4} | {topKStr,15}");
            }
            _out.WriteLine($"  -- N={N}: max-Δ K_90 within orbit = {table.MaxK90WithinOrbitDeviation}, " +
                           $"max-Δ row-L1 within orbit = {table.MaxRowL1WithinOrbitDeviation:G3}");
            _out.WriteLine("");
        }
    }
}
