using System;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.JordanWigner;

public class JwBlockBasisTests
{
    private readonly ITestOutputHelper _out;

    public JwBlockBasisTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void OrthonormalityResidual_BelowTolerance(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var jw = JwBlockBasis.Build(block);
        Assert.True(jw.OrthonormalityResidual < JwBlockBasis.Tolerance,
            $"N={N}: orthonormality residual {jw.OrthonormalityResidual:G3} ≥ tolerance");
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void MhTotalDiagonalityResidual_BelowTolerance(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var jw = JwBlockBasis.Build(block);
        Assert.True(jw.MhTotalDiagonalityResidual < JwBlockBasis.Tolerance,
            $"N={N}: MhTotal diagonality residual {jw.MhTotalDiagonalityResidual:G3} ≥ tolerance");
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void MhTotalEigenvalueMatchResidual_BelowTolerance(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var jw = JwBlockBasis.Build(block);
        Assert.True(jw.MhTotalEigenvalueMatchResidual < JwBlockBasis.Tolerance,
            $"N={N}: MhTotal eigenvalue-match residual {jw.MhTotalEigenvalueMatchResidual:G3} ≥ tolerance");
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void Triples_LengthEquals_Mtot(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var jw = JwBlockBasis.Build(block);
        int expected = N * N * (N - 1) / 2;  // N · C(N, 2)
        Assert.Equal(expected, jw.Triples.Count);
        Assert.Equal(block.Basis.MTotal, jw.Triples.Count);
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void Triples_AreOrdered_K1LessThanK2(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var jw = JwBlockBasis.Build(block);
        foreach (var triple in jw.Triples)
        {
            Assert.True(triple.K >= 1 && triple.K <= N,
                $"k={triple.K} outside [1, {N}]");
            Assert.True(triple.K1 >= 1 && triple.K1 <= N,
                $"k₁={triple.K1} outside [1, {N}]");
            Assert.True(triple.K2 >= 1 && triple.K2 <= N,
                $"k₂={triple.K2} outside [1, {N}]");
            Assert.True(triple.K1 < triple.K2,
                $"triple ({triple.K}, {triple.K1}, {triple.K2}): k₁ < k₂ violated");
        }
    }

    [Fact]
    public void Build_RejectsNonC2Block()
    {
        // For N=5, n=2: c = min(2, 5-1-2) + 1 = min(2, 2) + 1 = 3.
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);
        Assert.Equal(3, block.C);
        Assert.Throws<ArgumentException>(() => JwBlockBasis.Build(block));
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var jw = JwBlockBasis.Build(block);
        Assert.Equal(Tier.Tier1Derived, jw.Tier);
    }

    [Fact]
    public void Anchor_References_PROOF_F86_QPEAK_AndJWBasis()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var jw = JwBlockBasis.Build(block);
        Assert.Contains("PROOF_F86_QPEAK", jw.Anchor);
        Assert.Contains("Jordan-Wigner", jw.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsResiduals_AcrossN4To7()
    {
        _out.WriteLine("  N | Mtot | ortho-resid | diag-resid  | eig-match-resid");
        _out.WriteLine("  --|------|-------------|-------------|-----------------");
        foreach (int N in new[] { 4, 5, 6, 7 })
        {
            var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
            var jw = JwBlockBasis.Build(block);
            _out.WriteLine(
                $"  {N} | {jw.Triples.Count,4} | {jw.OrthonormalityResidual,11:G3} | " +
                $"{jw.MhTotalDiagonalityResidual,11:G3} | {jw.MhTotalEigenvalueMatchResidual,15:G3}");
        }
    }
}
