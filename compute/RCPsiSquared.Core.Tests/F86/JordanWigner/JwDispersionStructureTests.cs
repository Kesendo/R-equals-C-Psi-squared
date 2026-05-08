using System;
using System.Linq;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.JordanWigner;

public class JwDispersionStructureTests
{
    private readonly ITestOutputHelper _out;

    public JwDispersionStructureTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(4, 24)]
    [InlineData(5, 50)]
    [InlineData(6, 90)]
    [InlineData(7, 147)]
    [InlineData(8, 224)]
    public void TotalTriples_EqualsNTimesCN2(int N, int expected)
    {
        var s = JwDispersionStructure.Build(N);
        Assert.Equal(expected, s.TotalTriples);
        Assert.Equal(N * N * (N - 1) / 2, s.TotalTriples);
    }

    [Theory]
    [InlineData(4, 4)]
    [InlineData(5, 6)]
    [InlineData(6, 7)]
    [InlineData(7, 9)]
    [InlineData(8, 10)]
    public void LargestClusterSize_MatchesEnumerationCount(int N, int expected)
    {
        var s = JwDispersionStructure.Build(N);
        Assert.Equal(expected, s.LargestClusterSize);
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Clusters_AreDeltaSorted(int N)
    {
        var s = JwDispersionStructure.Build(N);
        for (int i = 1; i < s.Clusters.Count; i++)
            Assert.True(s.Clusters[i].Delta >= s.Clusters[i - 1].Delta - JwDispersionStructure.Tolerance,
                $"N={N}: clusters not δ-sorted at index {i}: {s.Clusters[i - 1].Delta} > {s.Clusters[i].Delta}");
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void DeltaSymmetry_AroundZero_FromCosineIdentity(int N)
    {
        // Cosine identity ε_{N+1−k} = −ε_k implies that for every cluster at δ there is a
        // mirror cluster at −δ with the same size. Tested via cluster-size multiset.
        var s = JwDispersionStructure.Build(N);
        var positive = s.Clusters.Where(c => c.Delta > JwDispersionStructure.Tolerance).Select(c => c.Triples.Count).OrderBy(x => x).ToArray();
        var negative = s.Clusters.Where(c => c.Delta < -JwDispersionStructure.Tolerance).Select(c => c.Triples.Count).OrderBy(x => x).ToArray();
        Assert.Equal(positive, negative);
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void EveryTriple_AppearsExactlyOnce_AcrossClusters(int N)
    {
        var s = JwDispersionStructure.Build(N);
        var allTriples = s.Clusters.SelectMany(c => c.Triples).ToList();
        Assert.Equal(s.TotalTriples, allTriples.Count);
        Assert.Equal(s.TotalTriples, allTriples.Distinct().Count());
    }

    [Theory]
    [InlineData(2)]   // smallest valid N (2 sites, 1 popcount-2 state, 2 popcount-1 states)
    public void SmallestValidN_BuildsCleanly(int N)
    {
        var s = JwDispersionStructure.Build(N);
        Assert.Equal(N * N * (N - 1) / 2, s.TotalTriples);
        Assert.True(s.Clusters.Count >= 1);
    }

    [Fact]
    public void Build_RejectsTooSmallN()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => JwDispersionStructure.Build(1));
        Assert.Throws<ArgumentOutOfRangeException>(() => JwDispersionStructure.Build(0));
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var s = JwDispersionStructure.Build(5);
        Assert.Equal(Tier.Tier1Derived, s.Tier);
    }

    [Fact]
    public void Anchor_References_PROOF_F86_QPEAK_AndJwTrack()
    {
        var s = JwDispersionStructure.Build(5);
        Assert.Contains("PROOF_F86_QPEAK", s.Anchor);
        Assert.Contains("JW track", s.Anchor);
    }

    // Verify the three algebraic patterns A/B/C generate the expected cluster members at
    // δ = -ε_1 for N=4. Tests the structural understanding of the cluster origin.
    [Fact]
    public void PatternA_TrivialInversion_AtN4DeltaMinusEpsilon1()
    {
        // At N=4 the cluster at δ = -ε_1 = -2cos(π/5) ≈ -1.6180 contains:
        //   Pattern A: (2,1,2), (3,1,3), (4,1,4) — (k, m=1, k) for k > 1
        //   Pattern B: (4,2,3) — k=N+1-m=4, k1+k2=N+1=5
        var s = JwDispersionStructure.Build(4);
        double targetDelta = -2.0 * Math.Cos(Math.PI / 5.0);
        var cluster = s.Clusters.FirstOrDefault(c => Math.Abs(c.Delta - targetDelta) < 1e-6);
        Assert.NotNull(cluster);
        Assert.Equal(4, cluster!.Triples.Count);
        var memberSet = new HashSet<(int, int, int)>(cluster.Triples.Select(t => (t.K, t.K1, t.K2)));
        Assert.Contains((2, 1, 2), memberSet);
        Assert.Contains((3, 1, 3), memberSet);
        Assert.Contains((4, 1, 4), memberSet);
        Assert.Contains((4, 2, 3), memberSet);
    }

    [Fact]
    public void Reconnaissance_EmitsClusterStats_AcrossN4To8()
    {
        _out.WriteLine("  N | Triples | Clusters | Multi (≥2) | Max | % in clusters");
        _out.WriteLine("  --|---------|----------|------------|-----|---------------");
        foreach (int N in new[] { 4, 5, 6, 7, 8 })
        {
            var s = JwDispersionStructure.Build(N);
            double pct = 100.0 * s.TriplesInDegenerateClusters / s.TotalTriples;
            _out.WriteLine($"  {N} | {s.TotalTriples,7} | {s.Clusters.Count,8} | {s.DegenerateClusterCount,10} | {s.LargestClusterSize,3} | {pct,12:F1}%");
        }
        _out.WriteLine("");
    }
}
