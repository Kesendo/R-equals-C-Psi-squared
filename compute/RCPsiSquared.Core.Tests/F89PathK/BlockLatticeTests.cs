using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The block-lattice enumeration primitive behind the step-3 shell census: the D₄ fold-lattice
/// orbit images with their spectral cocycle (PROOF_CODIM1_BY_ADDITIVITY §7, gate BlockLatticeFoldGroupTests),
/// the {p ≤ w, p + w ≤ N} quotient region, and the Bendixson rate windows of the window-shell lemma
/// (gate WindowShellLemmaTests). The windows are pinned FROM BELOW here: against the enumerated n_diff
/// values of the actual block bases, not against the closed-form combinatorics.</summary>
public class BlockLatticeTests
{
    [Theory]
    [InlineData(5)]
    [InlineData(9)]
    public void EveryBlock_TransfersIntoFundamentalDomain(int n)
    {
        foreach (var (p, w) in BlockLattice.AllBlocks(n))
            Assert.Contains(BlockLattice.OrbitImages(n, p, w),
                img => BlockLattice.InFundamentalDomain(n, img.P, img.W));
    }

    [Fact]
    public void OrbitImages_AreTheEightD4Legs_WithFoldParity()
    {
        int n = 5, p = 1, w = 2;
        var imgs = BlockLattice.OrbitImages(n, p, w);
        Assert.Equal(8, imgs.Count);
        // ker chi (same-lambda legs): identity, transpose, Klein, transpose*Klein
        Assert.Contains((1, 2, 0), imgs);
        Assert.Contains((2, 1, 0), imgs);
        Assert.Contains((4, 3, 0), imgs);
        Assert.Contains((3, 4, 0), imgs);
        // fold coset (lambda -> -lambda-2N legs)
        Assert.Contains((1, 3, 1), imgs);
        Assert.Contains((4, 2, 1), imgs);
        Assert.Contains((3, 1, 1), imgs);
        Assert.Contains((2, 4, 1), imgs);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void RateWindow_MatchesEnumeratedNdiff_FromBelow(int n)
    {
        foreach (var (p, w) in BlockLattice.AllBlocks(n))
        {
            var kets = WeightCoherenceBlock.Configs(n, p);
            var bras = WeightCoherenceBlock.Configs(n, w);
            var diffs = kets.SelectMany(k => bras.Select(b =>
                BitOperations.PopCount((uint)(k ^ b)))).ToList();
            var (nMin, nMax) = BlockLattice.NdiffRange(n, p, w);
            Assert.Equal(diffs.Min(), nMin);
            Assert.Equal(diffs.Max(), nMax);
            var (lo, hi) = BlockLattice.RateWindow(n, p, w);
            Assert.Equal(-2.0 * diffs.Max(), lo, 12);
            Assert.Equal(-2.0 * diffs.Min(), hi, 12);
        }
    }

    [Fact]
    public void WindowDistance_ZeroInside_PositiveOutside()
    {
        // (1,2) at N=5: window [-6,-2]
        Assert.Equal(0.0, BlockLattice.WindowDistance(5, 1, 2, -4.0), 12);
        Assert.Equal(1.0, BlockLattice.WindowDistance(5, 1, 2, -7.0), 12);
        Assert.Equal(0.5, BlockLattice.WindowDistance(5, 1, 2, -1.5), 12);
        // (1,6) at N=9: window [-14,-10]; the scout's control margin 5.619 at lambda_A=-4.3807
        Assert.Equal(5.6193, BlockLattice.WindowDistance(9, 1, 6, -4.3807), 4);
    }

    [Fact]
    public void Dim_IsBinomialProduct()
    {
        Assert.Equal(15876, BlockLattice.Dim(9, 4, 5));  // C(9,4)*C(9,5) = 126*126
        Assert.Equal(7056, BlockLattice.Dim(9, 3, 3));   // 84*84
        Assert.Equal(100, BlockLattice.Dim(5, 2, 3));    // 10*10
        Assert.Equal(213444, BlockLattice.Dim(11, 5, 5)); // 462*462, the deferred N=11 core
    }
}
