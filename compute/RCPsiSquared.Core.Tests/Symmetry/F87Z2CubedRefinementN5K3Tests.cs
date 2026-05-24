using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F87Z2CubedRefinementN5K3Tests
{
    [Fact]
    public void Z2Axis_IsYParity()
    {
        var claim = new F87Z2CubedRefinementN5K3();
        Assert.Equal(Z2Axis.YParity, claim.Z2Axis);
    }

    [Fact]
    public void BitATwin_IsNull()
    {
        var claim = new F87Z2CubedRefinementN5K3();
        Assert.Null(claim.BitATwin);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var claim = new F87Z2CubedRefinementN5K3();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void Regime_NkTotalPairs_AreN5K3With294Pairs()
    {
        // F105's anchor is N=5, k=3, same 294 Z2-cubed-homogeneous pairs as F103.
        // The pair count is N-independent at k=3 (letter enumeration); only the
        // chain on which we classify changes (N=4 to N=5). F85 predicts the
        // Pi2-class trichotomy is N-stable, so per-cell counts should match F103
        // bit-exactly. This was CONFIRMED by the Task 7 SLOW_F105_BATCH run.
        var claim = new F87Z2CubedRefinementN5K3();
        Assert.Equal(5, claim.N);
        Assert.Equal(3, claim.K);
        Assert.Equal(294, claim.TotalPairs);
    }

    [Fact]
    public void TrulyIsYParityZeroPure()
    {
        // F87 truly classifications across all 12 (Klein x dephase) cells at N=5 k=3
        // have y_par=1 count exactly 0; total truly across the grid is 300.
        // Bit-exact match to F103 at N=4: F85 N-stability prediction confirmed.
        var claim = new F87Z2CubedRefinementN5K3();
        Assert.Equal(0, claim.TrulyPurity.YParityOneCount);
        Assert.Equal(300, claim.TrulyPurity.TotalTrulyClassifications);
    }

    [Fact]
    public void HardDiagonalSplit42To8WithYInversion()
    {
        // Hard appears only in the diagonal Klein cells (Klein matches dephase letter).
        // Z and X dephase split 42:8 (y_par=0 dominant); Y dephase inverts to 8:42
        // because Y carries y_par=1. Each sum equals 50. Bit-exact match to F103.
        var claim = new F87Z2CubedRefinementN5K3();
        Assert.Equal((42, 8), claim.HardDiagonal.ZDephKlein01);
        Assert.Equal((42, 8), claim.HardDiagonal.XDephKlein10);
        Assert.Equal((8, 42), claim.HardDiagonal.YDephKlein11);

        // Y-inversion: Y-deph (8, 42) is the swap of Z/X-deph (42, 8).
        Assert.Equal(50, claim.HardDiagonal.ZDephKlein01.YPar0 + claim.HardDiagonal.ZDephKlein01.YPar1);
        Assert.Equal(50, claim.HardDiagonal.XDephKlein10.YPar0 + claim.HardDiagonal.XDephKlein10.YPar1);
        Assert.Equal(50, claim.HardDiagonal.YDephKlein11.YPar0 + claim.HardDiagonal.YDephKlein11.YPar1);

        Assert.Equal(claim.HardDiagonal.ZDephKlein01.YPar0, claim.HardDiagonal.YDephKlein11.YPar1);
        Assert.Equal(claim.HardDiagonal.ZDephKlein01.YPar1, claim.HardDiagonal.YDephKlein11.YPar0);
    }

    [Fact]
    public void DiagonalSoftSplit13To13Universal()
    {
        // The same 3 diagonal cells that host hard 42:8 also host soft 13:13 (sum 26).
        // Unlike the hard 42:8 asymmetry, soft is y_par-symmetric in these cells.
        // Bit-exact match to F103.
        var claim = new F87Z2CubedRefinementN5K3();
        Assert.Equal((13, 13), claim.DiagonalSoft.ZDephKlein01);
        Assert.Equal((13, 13), claim.DiagonalSoft.XDephKlein10);
        Assert.Equal((13, 13), claim.DiagonalSoft.YDephKlein11);

        Assert.Equal(26, claim.DiagonalSoft.ZDephKlein01.YPar0 + claim.DiagonalSoft.ZDephKlein01.YPar1);
        Assert.Equal(claim.DiagonalSoft.ZDephKlein01.YPar0, claim.DiagonalSoft.ZDephKlein01.YPar1);
        Assert.Equal(claim.DiagonalSoft.XDephKlein10.YPar0, claim.DiagonalSoft.XDephKlein10.YPar1);
        Assert.Equal(claim.DiagonalSoft.YDephKlein11.YPar0, claim.DiagonalSoft.YDephKlein11.YPar1);
    }

    [Fact]
    public void MotherSoftIsYParityOnePure()
    {
        // Klein (0,0) soft under any dephase letter is (0, 21): zero y_par=0 pairs,
        // 21 y_par=1 pairs. Cell sum 21. Bit-exact match to F103.
        var claim = new F87Z2CubedRefinementN5K3();
        Assert.Equal(0, claim.MotherSoft.ZDephCounts.YPar0);
        Assert.Equal(21, claim.MotherSoft.ZDephCounts.YPar1);
        Assert.Equal(0, claim.MotherSoft.XDephCounts.YPar0);
        Assert.Equal(21, claim.MotherSoft.XDephCounts.YPar1);
        Assert.Equal(0, claim.MotherSoft.YDephCounts.YPar0);
        Assert.Equal(21, claim.MotherSoft.YDephCounts.YPar1);
    }

    [Fact]
    public void OffDiagonalSoftPatternsHaveSixCells()
    {
        var claim = new F87Z2CubedRefinementN5K3();
        Assert.Equal(6, claim.OffDiagonalSoft.Cells.Count);

        // Pattern B (3 cells): proportional to (Klein, y_par) enum breakdown.
        // (0,1) and (1,0) Klein cells have enum split (55, 21); (1,1) has inverted (21, 55).
        // Bit-exact match to F103.
        var expectedPatternB = new Dictionary<(int KleinA, int KleinB, char Dephase), (int YPar0, int YPar1)>
        {
            { (0, 1, 'Y'), (55, 21) },
            { (1, 1, 'Z'), (21, 55) },
            { (1, 1, 'X'), (21, 55) },
        };
        foreach (var (key, expected) in expectedPatternB)
        {
            var actual = claim.OffDiagonalSoft.Cells[key];
            Assert.Equal(expected, actual);
            Assert.Equal(76, actual.YPar0 + actual.YPar1);
        }

        // Pattern C (3 cells): sum 21 each, YPar0 == 0. Bit-exact match to F103.
        var patternC = new (int KleinA, int KleinB, char Dephase)[]
        {
            (0, 1, 'X'),
            (1, 0, 'Z'),
            (1, 0, 'Y'),
        };
        foreach (var key in patternC)
        {
            var counts = claim.OffDiagonalSoft.Cells[key];
            Assert.Equal(0, counts.YPar0);
            Assert.Equal(21, counts.YPar1);
        }
    }
}
