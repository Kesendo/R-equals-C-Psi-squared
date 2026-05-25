using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public sealed class HardCellPureDTemplateTests
{
    private readonly HardCellPureDTemplate _claim = new();

    [Fact]
    public void Tier_IsTier1Candidate() =>
        Assert.Equal(Tier.Tier1Candidate, _claim.Tier);

    [Fact]
    public void Z2Axis_IsYParity() =>
        Assert.Equal(Z2Axis.YParity, _claim.Z2Axis);

    [Fact]
    public void BitATwin_IsNull_ForYParityAxis() =>
        Assert.Null(_claim.BitATwin);

    [Fact]
    public void BitATwinStatus_IsNotApplicableForThisAxis() =>
        Assert.Equal(BitATwinClassification.NotApplicableForThisAxis,
            ((IZ2AxisClaim)_claim).BitATwinStatus);

    [Fact]
    public void Theorem_MentionsPureDTemplateAndDiagonalCell()
    {
        Assert.Contains("pure-D template", _claim.Theorem);
        Assert.Contains("diagonal Klein cell", _claim.Theorem);
        Assert.Contains("F87-hard", _claim.Theorem);
    }

    [Fact]
    public void AnchorFile_References_PROOF_F111() =>
        Assert.Contains("PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md", _claim.Anchor);

    [Fact]
    public void AnchorFile_References_F110_Parent_Observation() =>
        Assert.Contains("PROOF_F110_HARD_CELL_Y_INVERSION.md", _claim.Anchor);

    // ============================================================
    // IsPureDTemplate (static helper)
    // ============================================================

    [Fact]
    public void IsPureDTemplate_DetectsOnlyDAndI()
    {
        // Pure-Z template: only Z and I letters.
        var pureZ = new PauliTerm(
            new[] { PauliLetter.Z, PauliLetter.I, PauliLetter.Z, PauliLetter.Z },
            Complex.One);
        Assert.True(HardCellPureDTemplate.IsPureDTemplate(pureZ, PauliLetter.Z));

        // Mixed (contains X) is NOT a pure-Z template.
        var mixed = new PauliTerm(
            new[] { PauliLetter.X, PauliLetter.Z, PauliLetter.Z, PauliLetter.I },
            Complex.One);
        Assert.False(HardCellPureDTemplate.IsPureDTemplate(mixed, PauliLetter.Z));

        // Pure-Y template under Y-deph.
        var pureY = new PauliTerm(
            new[] { PauliLetter.Y, PauliLetter.I, PauliLetter.Y, PauliLetter.Y },
            Complex.One);
        Assert.True(HardCellPureDTemplate.IsPureDTemplate(pureY, PauliLetter.Y));

        // Pure-Y template under Z-deph: contains Y (not Z), so NOT pure-Z.
        Assert.False(HardCellPureDTemplate.IsPureDTemplate(pureY, PauliLetter.Z));

        // All-I term is pure-D for any D in {X, Y, Z} (no non-I letter to violate).
        var allI = new PauliTerm(
            new[] { PauliLetter.I, PauliLetter.I, PauliLetter.I, PauliLetter.I },
            Complex.One);
        Assert.True(HardCellPureDTemplate.IsPureDTemplate(allI, PauliLetter.Z));
        Assert.True(HardCellPureDTemplate.IsPureDTemplate(allI, PauliLetter.X));
        Assert.True(HardCellPureDTemplate.IsPureDTemplate(allI, PauliLetter.Y));
    }

    [Fact]
    public void IsPureDTemplate_RejectsDephaseI() =>
        Assert.Throws<ArgumentException>(() =>
            HardCellPureDTemplate.IsPureDTemplate(
                new PauliTerm(
                    new[] { PauliLetter.I, PauliLetter.I, PauliLetter.I, PauliLetter.Z },
                    Complex.One),
                PauliLetter.I));

    // ============================================================
    // IsInDiagonalCellAtK4N4 (static helper)
    // ============================================================

    [Fact]
    public void IsInDiagonalCellAtK4N4_RequiresLength4AndMatchingKlein()
    {
        // Z-deph diagonal cell is Klein (0, 1). A pure-Z 4-letter term lives there.
        var pureZ4 = new PauliTerm(
            new[] { PauliLetter.Z, PauliLetter.I, PauliLetter.I, PauliLetter.I },
            Complex.One);
        Assert.Equal((0, 1), pureZ4.KleinIndex);
        Assert.True(HardCellPureDTemplate.IsInDiagonalCellAtK4N4(pureZ4, pureZ4, PauliLetter.Z));

        // Same term at length 3 fails the k=N=4 scope check.
        var pureZ3 = new PauliTerm(
            new[] { PauliLetter.Z, PauliLetter.I, PauliLetter.I },
            Complex.One);
        Assert.False(HardCellPureDTemplate.IsInDiagonalCellAtK4N4(pureZ3, pureZ3, PauliLetter.Z));

        // Pure-Z 4-letter under X-deph: X-deph diagonal is Klein (1, 0), so the
        // Z-Klein term (0, 1) is NOT in the X-deph diagonal cell.
        Assert.False(HardCellPureDTemplate.IsInDiagonalCellAtK4N4(pureZ4, pureZ4, PauliLetter.X));
    }

    [Fact]
    public void IsInDiagonalCellAtK4N4_RejectsDephaseI()
    {
        var t = new PauliTerm(
            new[] { PauliLetter.Z, PauliLetter.I, PauliLetter.I, PauliLetter.I },
            Complex.One);
        Assert.Throws<ArgumentException>(() =>
            HardCellPureDTemplate.IsInDiagonalCellAtK4N4(t, t, PauliLetter.I));
    }

    // ============================================================
    // IsPredictedHardAtK4N4 (static helper)
    // ============================================================

    [Fact]
    public void IsPredictedHardAtK4N4_TrueWhenAtLeastOneIsPureDTemplate()
    {
        // Z-deph diagonal cell is Klein (0, 1) = (bit_a=0, bit_b=1).
        // pureZ: one Z, rest I. bit_a = 0, bit_b = 1 ⇒ in Z-deph diagonal cell.
        var pureZ = new PauliTerm(
            new[] { PauliLetter.Z, PauliLetter.I, PauliLetter.I, PauliLetter.I },
            Complex.One);
        Assert.Equal((0, 1), pureZ.KleinIndex);

        // mixedInDiag: two X and one Z, rest I. bit_a = 2 mod 2 = 0, bit_b = 1 ⇒ in diagonal.
        // contains X, so NOT a pure-Z template.
        var mixedInDiag = new PauliTerm(
            new[] { PauliLetter.X, PauliLetter.X, PauliLetter.I, PauliLetter.Z },
            Complex.One);
        Assert.Equal((0, 1), mixedInDiag.KleinIndex);
        Assert.False(HardCellPureDTemplate.IsPureDTemplate(mixedInDiag, PauliLetter.Z));

        // (pure-Z, mixed-in-diag) is the Pure-Mixed cell case, predicted HARD.
        Assert.True(HardCellPureDTemplate.IsPredictedHardAtK4N4(pureZ, mixedInDiag, PauliLetter.Z));
        // Symmetric.
        Assert.True(HardCellPureDTemplate.IsPredictedHardAtK4N4(mixedInDiag, pureZ, PauliLetter.Z));
        // (pure-Z, pure-Z) is the Pure-Pure cell case, predicted HARD.
        Assert.True(HardCellPureDTemplate.IsPredictedHardAtK4N4(pureZ, pureZ, PauliLetter.Z));
    }

    [Fact]
    public void IsPredictedHardAtK4N4_FalseForTwoMixedTerms()
    {
        var mixed1 = new PauliTerm(
            new[] { PauliLetter.X, PauliLetter.X, PauliLetter.I, PauliLetter.Z },
            Complex.One);
        var mixed2 = new PauliTerm(
            new[] { PauliLetter.X, PauliLetter.X, PauliLetter.Z, PauliLetter.I },
            Complex.One);
        // Both Klein (0, 1) ⇒ in Z-deph diagonal cell.
        Assert.Equal((0, 1), mixed1.KleinIndex);
        Assert.Equal((0, 1), mixed2.KleinIndex);
        // Both mixed ⇒ Pure-D Template Rule predicts SOFT (open subclaim d).
        Assert.False(HardCellPureDTemplate.IsPredictedHardAtK4N4(mixed1, mixed2, PauliLetter.Z));
    }

    [Fact]
    public void IsPredictedHardAtK4N4_FalseWhenOutsideDiagonalCell()
    {
        // Pure-Z term under X-deph: NOT in the X-deph diagonal cell, so rule's
        // scope excludes it ⇒ predicted-hard = false (scope check).
        var pureZ = new PauliTerm(
            new[] { PauliLetter.Z, PauliLetter.I, PauliLetter.I, PauliLetter.I },
            Complex.One);
        Assert.False(HardCellPureDTemplate.IsPredictedHardAtK4N4(pureZ, pureZ, PauliLetter.X));
    }

    // ============================================================
    // VerifyYInversionCorollaryAtK4N4 (static helper)
    // ============================================================

    [Fact]
    public void YInversionCorollary_PureDHasYParEqualsYParOfD()
    {
        // Pure-Z template, length 4, in Z-deph diagonal cell.
        var pureZ = new PauliTerm(
            new[] { PauliLetter.Z, PauliLetter.I, PauliLetter.I, PauliLetter.I },
            Complex.One);
        // pure-Z has #Y=0 ⇒ y_par(term) = 0; y_par(Z) = Z.BitA() & Z.BitB() = 0 & 1 = 0.
        Assert.Equal(0, pureZ.YParity);
        Assert.True(HardCellPureDTemplate.VerifyYInversionCorollaryAtK4N4(pureZ, pureZ, PauliLetter.Z));

        // Pure-Y template, length 4, in Y-deph diagonal cell.
        var pureY = new PauliTerm(
            new[] { PauliLetter.Y, PauliLetter.I, PauliLetter.I, PauliLetter.I },
            Complex.One);
        // pure-Y has #Y=1 ⇒ y_par(term) = 1; y_par(Y) = Y.BitA() & Y.BitB() = 1 & 1 = 1.
        Assert.Equal(1, pureY.YParity);
        Assert.True(HardCellPureDTemplate.VerifyYInversionCorollaryAtK4N4(pureY, pureY, PauliLetter.Y));
    }

    [Fact]
    public void YInversionCorollary_VacuousWhenRuleDoesNotPredictHard()
    {
        // Two mixed terms in Z-deph diagonal cell: rule predicts NOT hard,
        // so corollary is vacuously satisfied.
        var mixed1 = new PauliTerm(
            new[] { PauliLetter.X, PauliLetter.X, PauliLetter.I, PauliLetter.Z },
            Complex.One);
        var mixed2 = new PauliTerm(
            new[] { PauliLetter.X, PauliLetter.X, PauliLetter.Z, PauliLetter.I },
            Complex.One);
        Assert.True(HardCellPureDTemplate.VerifyYInversionCorollaryAtK4N4(mixed1, mixed2, PauliLetter.Z));
    }

    // ============================================================
    // Subclaim properties
    // ============================================================

    [Fact]
    public void SubclaimA_MentionsHardAndPureD()
    {
        Assert.Contains("HARD", _claim.SubclaimA_PureDSingleTermHard);
        Assert.Contains("Pure-D", _claim.SubclaimA_PureDSingleTermHard);
    }

    [Fact]
    public void SubclaimB_MentionsSoftAndMixed()
    {
        Assert.Contains("SOFT", _claim.SubclaimB_MixedSingleTermSoft);
        Assert.Contains("Mixed", _claim.SubclaimB_MixedSingleTermSoft);
    }

    [Fact]
    public void SubclaimC_MentionsHardAndPureMixedPair()
    {
        Assert.Contains("HARD", _claim.SubclaimC_PureMixedPairHard);
        Assert.Contains("Pure-D", _claim.SubclaimC_PureMixedPairHard);
        Assert.Contains("Mixed", _claim.SubclaimC_PureMixedPairHard);
    }

    [Fact]
    public void SubclaimD_MentionsBlocking()
    {
        Assert.Contains("SOFT", _claim.SubclaimD_MixedMixedPairSoft_OPEN);
        Assert.True(
            _claim.SubclaimD_MixedMixedPairSoft_OPEN.Contains("open")
            || _claim.SubclaimD_MixedMixedPairSoft_OPEN.Contains("BLOCKED")
            || _claim.SubclaimD_MixedMixedPairSoft_OPEN.Contains("Tier1Derived"));
    }

    [Fact]
    public void YInversionCorollary_MentionsYParAndDephaseLetter()
    {
        Assert.Contains("y_par", _claim.YInversionCorollary);
        Assert.Contains("F87-hard pair", _claim.YInversionCorollary);
    }

    [Fact]
    public void DecompositionPerCell_Lists36_192_300_228()
    {
        Assert.Contains("36", _claim.DecompositionPerCell);
        Assert.Contains("192", _claim.DecompositionPerCell);
        Assert.Contains("300", _claim.DecompositionPerCell);
        Assert.Contains("228", _claim.DecompositionPerCell);
    }
}
