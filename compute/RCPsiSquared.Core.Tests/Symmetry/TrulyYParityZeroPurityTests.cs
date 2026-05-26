using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class TrulyYParityZeroPurityTests
{
    [Fact]
    public void Z2Axis_IsYParity()
    {
        var claim = new TrulyYParityZeroPurity(new KleinEightCellClaim(new KleinFourCellClaim()));
        Assert.Equal(Z2Axis.YParity, claim.Z2Axis);
    }

    [Fact]
    public void BitATwin_IsNull()
    {
        var claim = new TrulyYParityZeroPurity(new KleinEightCellClaim(new KleinFourCellClaim()));
        Assert.Null(claim.BitATwin);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var claim = new TrulyYParityZeroPurity(new KleinEightCellClaim(new KleinFourCellClaim()));
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void TrulyCriterionHolds_ZDephase_RequiresNYEvenAndNZEven()
    {
        // F85 Z-dephase truly: #Y even AND #Z even.
        // Examples: XXI (k_body=2, #Y=0, #Z=0) ⟹ truly under Z.
        var xxi = new PauliTerm(new[] { PauliLetter.X, PauliLetter.X, PauliLetter.I }, Complex.One);
        Assert.True(TrulyYParityZeroPurity.TrulyCriterionHolds(xxi, PauliLetter.Z));

        // XYI (#Y=1 odd) ⟹ NOT truly under Z.
        var xyi = new PauliTerm(new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.I }, Complex.One);
        Assert.False(TrulyYParityZeroPurity.TrulyCriterionHolds(xyi, PauliLetter.Z));

        // XZI (#Z=1 odd) ⟹ NOT truly under Z.
        var xzi = new PauliTerm(new[] { PauliLetter.X, PauliLetter.Z, PauliLetter.I }, Complex.One);
        Assert.False(TrulyYParityZeroPurity.TrulyCriterionHolds(xzi, PauliLetter.Z));
    }

    [Fact]
    public void TrulyCriterionHolds_XDephase_RequiresNXEvenAndNYEven()
    {
        // F85 SU(2)-X truly: #X even AND #Y even.
        // ZZI (#X=0, #Y=0) ⟹ truly under X.
        var zzi = new PauliTerm(new[] { PauliLetter.Z, PauliLetter.Z, PauliLetter.I }, Complex.One);
        Assert.True(TrulyYParityZeroPurity.TrulyCriterionHolds(zzi, PauliLetter.X));

        // XZI (#X=1 odd) ⟹ NOT truly under X.
        var xzi = new PauliTerm(new[] { PauliLetter.X, PauliLetter.Z, PauliLetter.I }, Complex.One);
        Assert.False(TrulyYParityZeroPurity.TrulyCriterionHolds(xzi, PauliLetter.X));
    }

    [Fact]
    public void TrulyCriterionHolds_YDephase_RequiresNYEvenAndNZEven()
    {
        // F85 SU(2)-Y truly: #Y even AND #Z even (same combined form as Z-dephase
        // because Π_Y has same per-letter swap as Π_Z, only the phase differs).
        // XII (#Y=0, #Z=0) ⟹ truly under Y.
        var xii = new PauliTerm(new[] { PauliLetter.X, PauliLetter.I, PauliLetter.I }, Complex.One);
        Assert.True(TrulyYParityZeroPurity.TrulyCriterionHolds(xii, PauliLetter.Y));

        // XYI (#Y=1 odd) ⟹ NOT truly under Y.
        var xyi = new PauliTerm(new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.I }, Complex.One);
        Assert.False(TrulyYParityZeroPurity.TrulyCriterionHolds(xyi, PauliLetter.Y));
    }

    [Fact]
    public void VerifyOnTerm_TrulyTermsAllHaveYParityZero()
    {
        // Enumerate all 64 k=3 letter sequences. For each (sequence, dephase letter),
        // if the truly criterion holds, the term must have y_par = 0. This is F107's
        // central claim.
        var letters = new[] { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z };
        var dephases = new[] { PauliLetter.Z, PauliLetter.X, PauliLetter.Y };
        int trulyCount = 0;
        int verifiedCount = 0;
        foreach (var a in letters)
            foreach (var b in letters)
                foreach (var c in letters)
                {
                    var term = new PauliTerm(new[] { a, b, c }, Complex.One);
                    foreach (var d in dephases)
                    {
                        Assert.True(TrulyYParityZeroPurity.VerifyOnTerm(term, d),
                            $"F107 violated: term [{a},{b},{c}] under {d}-dephase has truly criterion " +
                            $"but y_par = {term.YParity} (expected 0)");
                        if (TrulyYParityZeroPurity.TrulyCriterionHolds(term, d))
                        {
                            trulyCount++;
                            if (term.YParity == 0) verifiedCount++;
                        }
                    }
                }
        // All truly cases have y_par = 0 ⟺ trulyCount == verifiedCount.
        Assert.Equal(trulyCount, verifiedCount);
        // Sanity: some truly cases exist (not vacuous).
        Assert.True(trulyCount > 0, "Expected at least some truly classifications across 64×3 cases");
    }

    [Fact]
    public void VerifyOnTerm_AtK4_StillForcesYParityZero()
    {
        // Same enumeration at k=4 (256 letter sequences) to confirm F107 extends to k=4
        // (where F106 anchored empirically). Total truly count across 256×3 = 768 should
        // match the 9 cells × non-zero values in F106's TrulyCounts (total = 3924 from F106).
        // Here we just spot-check: each truly term has y_par = 0.
        var letters = new[] { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z };
        var dephases = new[] { PauliLetter.Z, PauliLetter.X, PauliLetter.Y };
        int trulyCount = 0;
        foreach (var a in letters)
            foreach (var b in letters)
                foreach (var c in letters)
                    foreach (var e in letters)
                    {
                        var term = new PauliTerm(new[] { a, b, c, e }, Complex.One);
                        foreach (var d in dephases)
                        {
                            Assert.True(TrulyYParityZeroPurity.VerifyOnTerm(term, d),
                                $"F107 violated at k=4: term [{a},{b},{c},{e}] under {d}-dephase " +
                                $"has truly criterion but y_par = {term.YParity}");
                            if (TrulyYParityZeroPurity.TrulyCriterionHolds(term, d))
                                trulyCount++;
                        }
                    }
        // Sanity: many truly cases at k=4 (more than k=3 since 256 sequences vs 64).
        Assert.True(trulyCount > 100, $"Expected substantial truly count at k=4; got {trulyCount}");
    }
}
