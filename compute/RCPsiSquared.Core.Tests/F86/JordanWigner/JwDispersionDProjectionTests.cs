using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.JordanWigner;

public class JwDispersionDProjectionTests
{
    private readonly ITestOutputHelper _out;

    public JwDispersionDProjectionTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void TotalFrobeniusSquared_MatchesDFrobeniusOnComputationalBasis(int N)
    {
        // Unitary basis change preserves Frobenius² (the JW transform U is unitary per T9).
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var proj = JwDispersionDProjection.Build(block);
        double dCompFrobSq = block.Decomposition.D.FrobeniusNorm();
        dCompFrobSq *= dCompFrobSq;
        Assert.Equal(dCompFrobSq, proj.TotalFrobeniusSquared, precision: 8);
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void BetweenClusterFrobenius_DominatesWithinCluster(int N)
    {
        // EP-relevance: D mixes JW dispersion-clusters, so the between-cluster off-diagonal
        // Frobenius² is larger than the within-cluster off-diagonal.
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var proj = JwDispersionDProjection.Build(block);
        Assert.True(proj.BetweenClusterFrobeniusSquared > proj.WithinClusterFrobeniusSquared,
            $"N={N}: between-cluster ({proj.BetweenClusterFrobeniusSquared:F4}) should exceed " +
            $"within-cluster ({proj.WithinClusterFrobeniusSquared:F4})");
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void DiagonalDominates_AbovePercent60(int N)
    {
        // The bulk of D's structure is captured by its JW-basis diagonal (the HD-class
        // eigenvalue projected onto each JW vector). This justifies treating D's
        // off-diagonal mixing as a perturbation in the closed-form derivation.
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var proj = JwDispersionDProjection.Build(block);
        double diagFrac = proj.DiagonalFrobeniusSquared / proj.TotalFrobeniusSquared;
        Assert.True(diagFrac > 0.60,
            $"N={N}: diagonal Frobenius² fraction {100 * diagFrac:F1}% should exceed 60%");
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void ClusterCardinalitySymmetry_BitExact(int N)
    {
        // Clusters with the same size have bit-identical intra-cluster D-Frobenius².
        // F71-mirror invariance + δ ↔ -δ symmetry rolled into one empirical fact.
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var proj = JwDispersionDProjection.Build(block);
        Assert.True(proj.MaxIntraClusterFrobeniusDeviationForSameSize < 1e-10,
            $"N={N}: max intra-cluster Frobenius² deviation for same-size clusters " +
            $"= {proj.MaxIntraClusterFrobeniusDeviationForSameSize:G3}");
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void ClusterProjections_CountMatchesDispersionStructure(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var proj = JwDispersionDProjection.Build(block);
        Assert.Equal(proj.DispersionStructure.Clusters.Count, proj.ClusterProjections.Count);
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void ClusterFrobeniusSums_MatchAggregates(int N)
    {
        // Σ_c diagonal_c = total diagonal; Σ_c within-cluster_off_c = total within-cluster.
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var proj = JwDispersionDProjection.Build(block);
        double sumDiag = proj.ClusterProjections.Sum(c => c.DiagonalFrobeniusSquared);
        double sumOff = proj.ClusterProjections.Sum(c => c.OffDiagonalFrobeniusSquared);
        Assert.Equal(proj.DiagonalFrobeniusSquared, sumDiag, precision: 8);
        Assert.Equal(proj.WithinClusterFrobeniusSquared, sumOff, precision: 8);
    }

    [Fact]
    public void Build_RejectsNonC2Block()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);  // c=3
        Assert.Throws<ArgumentException>(() => JwDispersionDProjection.Build(block));
    }

    [Fact]
    public void Tier_IsTier2Verified()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var proj = JwDispersionDProjection.Build(block);
        Assert.Equal(Tier.Tier2Verified, proj.Tier);
    }

    [Fact]
    public void Anchor_References_PROOF_F86_QPEAK_AndDirectionB()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var proj = JwDispersionDProjection.Build(block);
        Assert.Contains("PROOF_F86_QPEAK", proj.Anchor);
        Assert.Contains("Direction (b'')", proj.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsDecomposition_AcrossN4To6()
    {
        _out.WriteLine("  N | total | diag (%) | within (%) | between (%) | same-size dev");
        _out.WriteLine("  --|-------|----------|------------|-------------|--------------");
        foreach (int N in new[] { 4, 5, 6 })
        {
            var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
            var proj = JwDispersionDProjection.Build(block);
            double total = proj.TotalFrobeniusSquared;
            _out.WriteLine($"  {N} | {total,5:F2} | {100 * proj.DiagonalFrobeniusSquared / total,7:F1} | " +
                           $"{100 * proj.WithinClusterFrobeniusSquared / total,9:F1} | " +
                           $"{100 * proj.BetweenClusterFrobeniusSquared / total,10:F1} | " +
                           $"{proj.MaxIntraClusterFrobeniusDeviationForSameSize,12:G3}");
        }
    }
}
