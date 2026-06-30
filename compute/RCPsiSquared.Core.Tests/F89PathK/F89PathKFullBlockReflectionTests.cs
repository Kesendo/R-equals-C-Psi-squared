using System;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The FULL (SE,DE) block (R-even AND R-odd, not just the symmetric sector the scout uses) and the
/// site-reflection permutation R: i ↦ nBlock−1−i. This is the substrate for the odd-N real-q diabolic
/// mechanism: R splits the block into an R-even and an R-odd sector, the reflection-fixed central site
/// (only at odd nBlock) makes dim(even) − dim(odd) = (nBlock−1)/2 ≠ 0, and that dimension mismatch is what
/// forces the realness antiunitarity to act within each sector at odd N (self-conjugate, real eigenvalues
/// that can cross) instead of across the sectors at even N. Structural pins here; the EVD-based
/// self/cross-conjugacy census and the diabolic reproduction live in the Diagnostics witness/tests.</summary>
public class F89PathKFullBlockReflectionTests
{
    [Theory]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    public void BuildFullBlock_HasFullCoherenceDimension(int nBlock)
    {
        int pairs = nBlock * (nBlock - 1) / 2;
        int full = nBlock * pairs;                                  // SE index × DE pair, no symmetrization
        var block = F89PathKSeDeBlock.BuildFullBlock(nBlock, new Complex(1.0, 0.0));
        Assert.Equal(full, block.GetLength(0));
        Assert.Equal(full, block.GetLength(1));
    }

    [Theory]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    public void ReflectionPermutation_IsInvolution(int nBlock)
    {
        var perm = F89PathKSeDeBlock.ReflectionPermutation(nBlock);
        for (int t = 0; t < perm.Length; t++)
            Assert.Equal(t, perm[perm[t]]);
    }

    [Theory]
    [InlineData(7)]
    [InlineData(8)]
    public void FullBlock_CommutesWithReflection(int nBlock)
    {
        var block = F89PathKSeDeBlock.BuildFullBlock(nBlock, new Complex(0.7, 0.3));
        var perm = F89PathKSeDeBlock.ReflectionPermutation(nBlock);
        int n = block.GetLength(0);
        double maxDiff = 0;
        for (int r = 0; r < n; r++)                                 // [R, block] = 0 ⟺ block[perm r, perm c] = block[r, c]
            for (int c = 0; c < n; c++)
                maxDiff = Math.Max(maxDiff, (block[perm[r], perm[c]] - block[r, c]).Magnitude);
        Assert.True(maxDiff < 1e-12, $"[R, block] != 0: maxDiff={maxDiff}");
    }

    [Theory]
    [InlineData(7, 3)]      // odd nBlock: (nBlock−1)/2 reflection-fixed singletons (center SE × self-mirror DE)
    [InlineData(9, 4)]
    [InlineData(6, 0)]      // even nBlock: no fixed central site, no singletons
    [InlineData(8, 0)]
    public void ReflectionPermutation_FixedSingletonCount_IsCenterSiteFingerprint(int nBlock, int expectedFixed)
    {
        var perm = F89PathKSeDeBlock.ReflectionPermutation(nBlock);
        int fixedCount = 0;
        for (int t = 0; t < perm.Length; t++) if (perm[t] == t) fixedCount++;
        Assert.Equal(expectedFixed, fixedCount);
    }
}
