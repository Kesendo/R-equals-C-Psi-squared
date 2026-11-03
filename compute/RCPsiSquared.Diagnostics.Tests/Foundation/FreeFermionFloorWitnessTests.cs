using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The free-fermion floor of the number-changing (n, n+1) blocks
/// (docs/proofs/PROOF_CHAIN_GAP_DOMINANCE.md §2-§3, read block by block): every block of the
/// uniform XY chain holds exactly N modes at Re λ = −2γ₀, the Jordan-Wigner family
/// P_n·c_k·P_{n+1}, and the uniform Dicke probe sees only the mirror-even half of them,
/// which in the (1,2) block is F89's F_a set.</summary>
public class FreeFermionFloorWitnessTests
{
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void EveryBlock_HoldsExactlyN_FloorModes_CertifiedExactly(int n)
    {
        var w = new FreeFermionFloorWitness(n);
        Assert.Equal(n, w.Blocks.Count);
        foreach (var b in w.Blocks)
        {
            // the lower bound, in integers: the N Jordan-Wigner vectors sit on HD = 1 cells and the
            // hopping maps their span onto itself as the chain's adjacency matrix (one global sign)
            Assert.True(b.FamilyOnHd1Cells, $"N={n} n={b.LowerPopcount}: family leaves the HD=1 cells");
            Assert.True(b.FamilyClosedUnderHopping, $"N={n} n={b.LowerPopcount}: family not hopping-invariant");
            Assert.True(b.HoppingIsChainAdjacency, $"N={n} n={b.LowerPopcount}: family's hopping is not ±adjacency");
            // the upper bound, mod p: the largest hopping-invariant subspace inside the HD = 1 cells
            Assert.Equal(n, b.FloorUpperBoundModP);
            Assert.True(b.FloorCertified);
        }
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void Eigensolver_MeetsTheExactCount_AndTheModesAboveSitOffTheFloor(int n)
    {
        var w = new FreeFermionFloorWitness(n);
        foreach (var b in w.Blocks)
        {
            Assert.Equal(n, b.EigensolverFloorCount);
            // the floor frequencies are the single-particle comb, read against the rounding law
            Assert.True(b.MaxFloorFrequencyErrorOverRounding <= 1.0,
                $"N={n} n={b.LowerPopcount}: frequency error {b.MaxFloorFrequencyErrorOverRounding} rounding units");
            if (b.BlockDimension > n)
            {
                // two-sided: everything that is not the family sits strictly above the floor
                Assert.False(double.IsNaN(b.NextRateAboveFloor));
                Assert.True(b.NextRateAboveFloor > 2.0 * (1.0 + 1e-6),
                    $"N={n} n={b.LowerPopcount}: next rate {b.NextRateAboveFloor} γ₀ not above the floor");
            }
            else
            {
                // the pinned blocks (F153): the family is the whole block
                Assert.True(b.LowerPopcount == 0 || b.LowerPopcount == n - 1);
                Assert.True(double.IsNaN(b.NextRateAboveFloor));
            }
        }
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void DickeProbe_SeesExactlyTheMirrorEvenHalf(int n)
    {
        var w = new FreeFermionFloorWitness(n);
        foreach (var b in w.Blocks)
        {
            Assert.True(b.ReflectionRelationExact, $"N={n} n={b.LowerPopcount}: integer reflection relation fails");
            Assert.Equal(n, b.DickeVisibleK.Count + b.DickeDarkK.Count);
            // the exact parity route and the eigensolver's biorthogonal overlap agree on the dark set
            Assert.Equal(b.PredictedDarkK.OrderBy(k => k), b.DickeDarkK.OrderBy(k => k));
            // two-sided: the dark overlaps are rounding, the visible ones are not
            Assert.True(b.MaxDarkOverlapSquared < 1e-20, $"dark overlap {b.MaxDarkOverlapSquared}");
            Assert.True(b.MinVisibleOverlapSquared > 1e-6, $"visible overlap {b.MinVisibleOverlapSquared}");
        }
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void TheC2Block_VisibleHalf_IsF89sFaSet(int n)
    {
        var b = new FreeFermionFloorWitness(n).Blocks[1];
        Assert.Equal(1, b.LowerPopcount);
        Assert.Equal(F89PathKAtLockMechanismClaim.FaCount(n), b.DickeVisibleK.Count);
        Assert.Equal(F89PathKAtLockMechanismClaim.SeAntiBlochOrbit(n), b.DickeVisibleK.OrderBy(k => k).ToList());
        // the other half exists on the floor all the same, dark only by mirror parity
        Assert.Equal(n - F89PathKAtLockMechanismClaim.FaCount(n), b.DickeDarkK.Count);
    }

    [Fact]
    public void TheFloorIsQFlat_TheChromaticityDoesNotMatter()
    {
        // sterile at every chromaticity: the family's count does not move with Q, at c = 1, 2, 3
        foreach (double q in new[] { 0.3, 2.0, 40.0 })
        {
            var w = new FreeFermionFloorWitness(5, q);
            Assert.Equal(new[] { 1, 2, 3, 2, 1 }, w.Blocks.Select(b => b.Chromaticity).ToArray());
            Assert.All(w.Blocks, b => Assert.Equal(5, b.EigensolverFloorCount));
        }
    }

    [Fact]
    public void Constructor_RejectsOutOfRangeArguments()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new FreeFermionFloorWitness(1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new FreeFermionFloorWitness(8));
        Assert.Throws<ArgumentOutOfRangeException>(() => new FreeFermionFloorWitness(5, 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new FreeFermionFloorWitness(5, 2.0, -0.05));
    }

    [Fact]
    public void Inspect_ShowsOneNodePerBlock()
    {
        IInspectable w = new FreeFermionFloorWitness(4);
        var kids = w.Children.ToList();
        Assert.Equal(1 + 4, kids.Count); // the gate node, then one node per (n, n+1) block
        Assert.Contains("certified", kids[0].Summary);
    }
}
