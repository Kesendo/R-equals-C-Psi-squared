using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.JordanWigner;

public class JwClusterDEigenstructureTests
{
    private readonly ITestOutputHelper _out;

    public JwClusterDEigenstructureTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void ClusterCardinalitySymmetry_AtEigenvalueLevel_BitExact(int N)
    {
        // Theorem: clusters of the same size have IDENTICAL eigenvalue spectra (W_c
        // matrices are unitarily equivalent across same-size clusters). Stronger than
        // T11's Frobenius²-level symmetry.
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var es = JwClusterDEigenstructure.Build(block);
        Assert.True(es.MaxIntraClusterEigenvalueDeviationForSameSize < JwClusterDEigenstructure.SymmetryTolerance,
            $"N={N}: max same-size eigenvalue deviation = {es.MaxIntraClusterEigenvalueDeviationForSameSize:G3} " +
            $"≥ tolerance {JwClusterDEigenstructure.SymmetryTolerance:G3}");
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void EigenvalueSpectraSorted_PerCluster(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var es = JwClusterDEigenstructure.Build(block);
        foreach (var w in es.Clusters)
            for (int i = 1; i < w.Eigenvalues.Count; i++)
                Assert.True(w.Eigenvalues[i] >= w.Eigenvalues[i - 1] - 1e-10,
                    $"cluster δ={w.Cluster.Delta:F4}: eigenvalues not sorted at index {i}");
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void EigenvalueCount_MatchesClusterSize(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var es = JwClusterDEigenstructure.Build(block);
        foreach (var w in es.Clusters)
            Assert.Equal(w.Cluster.Triples.Count, w.Eigenvalues.Count);
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void EigenvalueRange_BoundedBy_HD1AndHD3DiagonalValues(int N)
    {
        // The W_c eigenvalues of any cluster lie within [-6γ, -2γ] (the HD=3 and HD=1
        // diagonal D values), since W_c is the projection of D (which is diagonal in
        // computational basis with these two values) onto the cluster sub-space.
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        double γ = block.GammaZero;
        var es = JwClusterDEigenstructure.Build(block);
        foreach (var w in es.Clusters)
            foreach (var λ in w.Eigenvalues)
                Assert.InRange(λ, -6 * γ - 1e-10, -2 * γ + 1e-10);
    }

    [Fact]
    public void N4_AllSize4Clusters_HaveExpectedSpectrum()
    {
        // Pinned: at N=4 γ=0.05 all 4 size-4 clusters share the spectrum
        // {-6γ, -4.4γ, -3.6γ, -2γ} = {-0.30, -0.22, -0.18, -0.10}.
        var block = new CoherenceBlock(N: 4, n: 1, gammaZero: 0.05);
        var es = JwClusterDEigenstructure.Build(block);
        var size4 = es.Clusters.Where(c => c.Cluster.Triples.Count == 4).ToList();
        Assert.Equal(4, size4.Count);

        double[] expected = { -0.30, -0.22, -0.18, -0.10 };
        foreach (var w in size4)
            for (int i = 0; i < 4; i++)
                Assert.Equal(expected[i], w.Eigenvalues[i], precision: 6);
    }

    [Fact]
    public void Build_RejectsNonC2Block()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);
        Assert.Throws<ArgumentException>(() => JwClusterDEigenstructure.Build(block));
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var es = JwClusterDEigenstructure.Build(block);
        Assert.Equal(Tier.Tier1Derived, es.Tier);
    }

    [Fact]
    public void Anchor_References_PROOF_F86_QPEAK_AndDirectionB()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var es = JwClusterDEigenstructure.Build(block);
        Assert.Contains("PROOF_F86_QPEAK", es.Anchor);
        Assert.Contains("Direction (b'')", es.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsClusterSpectra_AcrossN4To6()
    {
        _out.WriteLine("  N | cluster δ | size | eigenvalues");
        _out.WriteLine("  --|-----------|------|------------");
        foreach (int N in new[] { 4, 5, 6 })
        {
            var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
            var es = JwClusterDEigenstructure.Build(block);
            // Emit one representative cluster per size
            var bySize = es.Clusters.GroupBy(c => c.Cluster.Triples.Count).OrderBy(g => g.Key);
            foreach (var group in bySize)
            {
                var rep = group.First();
                string evals = string.Join(", ", rep.Eigenvalues.Select(e => e.ToString("F4")));
                _out.WriteLine($"  {N} | {rep.Cluster.Delta,9:F4} | {rep.Cluster.Triples.Count,4} | [{evals}]");
            }
            _out.WriteLine($"  -- N={N}: max same-size eigenvalue deviation = {es.MaxIntraClusterEigenvalueDeviationForSameSize:G3}");
            _out.WriteLine("");
        }
    }
}
