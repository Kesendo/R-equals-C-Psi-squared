using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F87Z2CubedRefinementN4K4Tests
{
    [Fact]
    public void Z2Axis_IsYParity()
    {
        var claim = new F87Z2CubedRefinementN4K4(new KleinEightCellClaim(new KleinFourCellClaim()));
        Assert.Equal(Z2Axis.YParity, claim.Z2Axis);
    }

    [Fact]
    public void BitATwin_IsNull()
    {
        var claim = new F87Z2CubedRefinementN4K4(new KleinEightCellClaim(new KleinFourCellClaim()));
        Assert.Null(claim.BitATwin);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var claim = new F87Z2CubedRefinementN4K4(new KleinEightCellClaim(new KleinFourCellClaim()));
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void Regime_NkTotalPairs_AreN4K4With4248Pairs()
    {
        // F106's anchor is N=4, k=4, 4248 Z2-cubed-homogeneous pairs.
        var claim = new F87Z2CubedRefinementN4K4(new KleinEightCellClaim(new KleinFourCellClaim()));
        Assert.Equal(4, claim.N);
        Assert.Equal(4, claim.K);
        Assert.Equal(4248, claim.TotalPairs);
    }

    [Fact]
    public void TrulyCounts_MatchFrozen()
    {
        // At k=4, total truly classifications across the 12 (Klein × dephase) grid is
        // 3924; y_par=1 truly count is 0. Y_par=0-purity HELD at k=4 (not predicted by F85;
        // empirically confirmed). All truly entries in the JSON have y_par=0.
        var claim = new F87Z2CubedRefinementN4K4(new KleinEightCellClaim(new KleinFourCellClaim()));
        Assert.Equal(3924, claim.TrulyPurity.TotalTrulyClassifications);
        Assert.Equal(0, claim.TrulyPurity.YParityOneCount);
    }

    [Fact]
    public void HardDiagonalSplit_MatchesFrozen()
    {
        // Diagonal cells (Klein matches dephase). At k=4 the 42:8 ratio from k=3
        // (mixed y_par) BROKE: now polarized to (228, 0) for Z/X and (0, 228) for Y.
        // Y-inversion structure SURVIVED qualitatively (the y_par split flips for
        // Y-dephase); the ratio became 100% pure per cell, not 42:8 mixed.
        var claim = new F87Z2CubedRefinementN4K4(new KleinEightCellClaim(new KleinFourCellClaim()));
        Assert.Equal((228, 0), claim.HardDiagonal.ZDephKlein01);
        Assert.Equal((228, 0), claim.HardDiagonal.XDephKlein10);
        Assert.Equal((0, 228), claim.HardDiagonal.YDephKlein11);
    }

    [Fact]
    public void DiagonalSoftSplit_MatchesFrozen()
    {
        // Same 3 diagonal cells, soft classifications. At k=4 the y_par-symmetry from
        // k=3 (13:13) BROKE: now asymmetric (300, 528) for Z/X and (528, 300) for Y
        // (Y-inverted). Total soft per diagonal cell = 828; the 13+13 = 26 at k=3
        // proportionally maps to 828 = 300+528 (the larger of the two halves is the
        // 528 count, matching the mother-soft sibling structure).
        var claim = new F87Z2CubedRefinementN4K4(new KleinEightCellClaim(new KleinFourCellClaim()));
        Assert.Equal((300, 528), claim.DiagonalSoft.ZDephKlein01);
        Assert.Equal((300, 528), claim.DiagonalSoft.XDephKlein10);
        Assert.Equal((528, 300), claim.DiagonalSoft.YDephKlein11);
    }

    [Fact]
    public void MotherSoft_MatchesFrozen()
    {
        // Klein (0,0) under each dephase letter. At k=4 the (0,0) enum balance is 780/300
        // (vs 45/21 at k=3). DESPITE the structural worry that the y_par=1-purity might
        // break, it SURVIVED: still (0, 300) i.e. y_par=1-pure across all 3 dephase
        // letters. The truly is now 780 (y_par=0) and the soft is now 300 (y_par=1).
        var claim = new F87Z2CubedRefinementN4K4(new KleinEightCellClaim(new KleinFourCellClaim()));
        Assert.Equal((0, 300), claim.MotherSoft.ZDephCounts);
        Assert.Equal((0, 300), claim.MotherSoft.XDephCounts);
        Assert.Equal((0, 300), claim.MotherSoft.YDephCounts);
    }

    [Fact]
    public void OffDiagonalSoftPatterns_HasExpectedCells()
    {
        // 6 off-diagonal soft cells at k=4 (Klein non-mother, Klein != dephase Klein):
        // 3 are y_par=1-pure (0, 528) (Pattern C analog: F105 had (0, 21) at the same
        // positions); 3 are fully y_par-symmetric (528, 528) (new pattern at k=4;
        // F105's same positions had (55, 21)/(21, 55) asymmetric Pattern B).
        var claim = new F87Z2CubedRefinementN4K4(new KleinEightCellClaim(new KleinFourCellClaim()));
        Assert.Equal(6, claim.OffDiagonalSoft.Cells.Count);
        Assert.Equal((0, 528), claim.OffDiagonalSoft.Cells[(0, 1, 'X')]);
        Assert.Equal((528, 528), claim.OffDiagonalSoft.Cells[(0, 1, 'Y')]);
        Assert.Equal((0, 528), claim.OffDiagonalSoft.Cells[(1, 0, 'Y')]);
        Assert.Equal((0, 528), claim.OffDiagonalSoft.Cells[(1, 0, 'Z')]);
        Assert.Equal((528, 528), claim.OffDiagonalSoft.Cells[(1, 1, 'X')]);
        Assert.Equal((528, 528), claim.OffDiagonalSoft.Cells[(1, 1, 'Z')]);
    }
}
