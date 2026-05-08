using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.JordanWigner;

public class JwBondClusterPairAffinityTests
{
    private readonly ITestOutputHelper _out;

    public JwBondClusterPairAffinityTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void NumBondsMatchesBlockNumBonds(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var aff = JwBondClusterPairAffinity.Build(block);
        Assert.Equal(N - 1, aff.Bonds.Count);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void EveryBond_HasNonEmptyTopPair(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var aff = JwBondClusterPairAffinity.Build(block);
        foreach (var b in aff.Bonds)
        {
            Assert.NotNull(b.TopPair);
            Assert.True(b.TopPair!.FrobeniusSquared > 0);
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void RankedPairs_AreSortedByFrobeniusSquaredDescending(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var aff = JwBondClusterPairAffinity.Build(block);
        foreach (var b in aff.Bonds)
            for (int i = 1; i < b.RankedPairs.Count; i++)
                Assert.True(b.RankedPairs[i].FrobeniusSquared <= b.RankedPairs[i - 1].FrobeniusSquared + 1e-12,
                    $"bond {b.Bond}: rankings not descending at index {i}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void F71MirrorBonds_HaveSameTopPairFrobenius(int N)
    {
        // Bonds b and N-2-b are F71-mirror partners; their cluster-pair affinity rankings
        // must match (cluster pair → mirror cluster pair, same Frobenius²).
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var aff = JwBondClusterPairAffinity.Build(block);
        int numBonds = aff.Bonds.Count;
        for (int b = 0; b < numBonds; b++)
        {
            int mirror = numBonds - 1 - b;
            if (mirror <= b) continue;
            Assert.Equal(aff.Bonds[b].TopPair!.FrobeniusSquared,
                         aff.Bonds[mirror].TopPair!.FrobeniusSquared,
                         precision: 8);
        }
    }

    // The bond-class structural finding: Endpoint-bonds prefer cluster-pairs with smaller
    // |Δδ| than Innermost-Interior bonds. This explains the empirical Q_peak ordering
    // (Endpoint Q_peak > Innermost Q_peak) via Q_EP ∝ 1/|Δδ|.
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void EndpointBonds_PreferSmaller_AbsoluteDeltaδ_ThanInnermost(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var aff = JwBondClusterPairAffinity.Build(block);
        var endpoint = aff.Bonds.First(b => b.BondClass == BondClass.Endpoint);
        var innermost = ClosestToCenter(aff.Bonds, N);
        if (innermost is null || innermost.Bond == endpoint.Bond) return;

        Assert.True(endpoint.TopPair!.AbsoluteDeltaδ <= innermost.TopPair!.AbsoluteDeltaδ + 1e-10,
            $"N={N}: Endpoint top |Δδ| ({endpoint.TopPair.AbsoluteDeltaδ:F4}) should be ≤ " +
            $"Innermost top |Δδ| ({innermost.TopPair.AbsoluteDeltaδ:F4})");
    }

    private static BondClusterPairAffinity? ClosestToCenter(IReadOnlyList<BondClusterPairAffinity> bonds, int N)
    {
        double center = (N - 2) / 2.0;
        BondClusterPairAffinity? best = null;
        double bestDist = double.PositiveInfinity;
        foreach (var b in bonds)
        {
            if (b.BondClass != BondClass.Interior) continue;
            double d = Math.Abs(b.Bond - center);
            if (d < bestDist) { bestDist = d; best = b; }
        }
        return best;
    }

    [Fact]
    public void Build_RejectsNonC2Block()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);
        Assert.Throws<ArgumentException>(() => JwBondClusterPairAffinity.Build(block));
    }

    [Fact]
    public void Tier_IsTier2Verified()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var aff = JwBondClusterPairAffinity.Build(block);
        Assert.Equal(Tier.Tier2Verified, aff.Tier);
    }

    [Fact]
    public void Anchor_References_PROOF_F86_QPEAK_AndDirectionB()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var aff = JwBondClusterPairAffinity.Build(block);
        Assert.Contains("PROOF_F86_QPEAK", aff.Anchor);
        Assert.Contains("Direction (b'')", aff.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsTopPairsPerBond_AcrossN5To7()
    {
        _out.WriteLine("  N | bond | class    | top-pair (c1, c2) | δ_c1   δ_c2   |Δδ|   | frob²");
        _out.WriteLine("  --|------|----------|-------------------|----------------------|--------");
        foreach (int N in new[] { 5, 6, 7 })
        {
            var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
            var aff = JwBondClusterPairAffinity.Build(block);
            foreach (var b in aff.Bonds)
            {
                var t = b.TopPair!;
                _out.WriteLine($"  {N} | {b.Bond,4} | {b.BondClass,-8} | ({t.Cluster1Index,2}, {t.Cluster2Index,2}) sz=({t.Cluster1Size},{t.Cluster2Size}) | " +
                               $"{t.Cluster1Delta,6:F3} {t.Cluster2Delta,6:F3} {t.AbsoluteDeltaδ,6:F3} | {t.FrobeniusSquared:F4}");
            }
            _out.WriteLine("");
        }
    }
}
