using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>F112 LindbladBitBPiBalance Claim metadata and helper-predicate tests.
/// Covers the typed-knowledge surface (Z2Axis, BitATwin classification, Tier,
/// Theorem / Summary strings) plus the static helpers
/// (<see cref="LindbladBitBPiBalance.BitBParity"/>,
/// <see cref="LindbladBitBPiBalance.IsBitBHomogeneous"/>,
/// <see cref="LindbladBitBPiBalance.IsHermitianHamiltonian"/>).
///
/// <para>The numerical bit-exact verification of the theorem itself
/// (‖M_plus_half‖² = ‖M_minus_half‖² across 14 probes) lives in the Python
/// framework / probe scripts referenced by the Claim's anchor; this test file
/// covers the C# typed-knowledge surface and the helper predicates.</para></summary>
public class LindbladBitBPiBalanceTests
{
    private static LindbladBitBPiBalance Make()
    {
        var part2 = new F108Part2Pi2XEvenAlwaysPalindromic();
        return new LindbladBitBPiBalance(
            new F108Part1Pi2EvenAlwaysPalindromic(part2),
            new LindbladBitAPiBalance(part2));
    }

    // ============================================================
    // Claim metadata
    // ============================================================

    [Fact]
    public void Z2Axis_IsBitB() =>
        Assert.Equal(Z2Axis.BitB, Make().Z2Axis);

    [Fact]
    public void BitATwin_IsLindbladBitAPiBalance()
    {
        // Welle 15 (2026-05-27) wired F112-X (LindbladBitAPiBalance) as the typed
        // BitA twin per PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md section (f) item 4.
        // Pre-Welle-15, BitATwin was null and BitATwinStatus = BitBSpecific.
        var claim = Make();
        Assert.NotNull(claim.BitATwin);
        Assert.IsType<LindbladBitAPiBalance>(claim.BitATwin);
    }

    [Fact]
    public void BitATwinStatus_IsFilled() =>
        // Flipped from BitBSpecific to Filled on Welle 15 wiring of F112-X as the
        // typed BitA twin.
        Assert.Equal(BitATwinClassification.Filled, Make().BitATwinStatus);

    [Fact]
    public void Tier_IsTier1Derived()
    {
        // The non-Hermitian H extension is also Tier1Derived universal N (Welle 11,
        // 2026-05-27, structural proof in
        // `PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md`); both Hermitian and non-Hermitian
        // H now sit under this single Tier1Derived Claim.
        Assert.Equal(Tier.Tier1Derived, Make().Tier);
    }

    [Fact]
    public void Theorem_Mentions_bitB_homogeneous_and_Hermitian()
    {
        var theorem = Make().Theorem;
        Assert.Contains("bit_b-homogeneous", theorem);
        Assert.Contains("Hermitian H", theorem);
        Assert.Contains("‖M_plus_half‖² = ‖M_minus_half‖²", theorem);
    }

    [Fact]
    public void Summary_NotEmpty()
    {
        var summary = Make().Summary;
        Assert.False(string.IsNullOrWhiteSpace(summary));
    }

    // ============================================================
    // BitBParity static helper
    // ============================================================

    [Fact]
    public void BitBParity_OfSingleY_Is1()
    {
        var term = PauliTerm.SingleSite(3, 0, PauliLetter.Y, Complex.One);
        Assert.Equal(1, LindbladBitBPiBalance.BitBParity(term));
    }

    [Fact]
    public void BitBParity_OfSingleZ_Is1()
    {
        var term = PauliTerm.SingleSite(3, 0, PauliLetter.Z, Complex.One);
        Assert.Equal(1, LindbladBitBPiBalance.BitBParity(term));
    }

    [Fact]
    public void BitBParity_OfSingleX_Is0()
    {
        var term = PauliTerm.SingleSite(3, 0, PauliLetter.X, Complex.One);
        Assert.Equal(0, LindbladBitBPiBalance.BitBParity(term));
    }

    [Fact]
    public void BitBParity_OfSingleI_Is0()
    {
        // All-I term has bit_b sum = 0 at any N.
        var term = new PauliTerm(new[] { PauliLetter.I, PauliLetter.I, PauliLetter.I },
                                 Complex.One);
        Assert.Equal(0, LindbladBitBPiBalance.BitBParity(term));
    }

    [Fact]
    public void BitBParity_OfYZ_Is0()
    {
        // bit_b(Y) + bit_b(Z) = 1 + 1 = 2 ≡ 0 mod 2.
        var term = PauliTerm.TwoSite(3, 0, PauliLetter.Y, 1, PauliLetter.Z, Complex.One);
        Assert.Equal(0, LindbladBitBPiBalance.BitBParity(term));
    }

    [Fact]
    public void BitBParity_OfXY_Is1()
    {
        // bit_b(X) + bit_b(Y) = 0 + 1 = 1 mod 2.
        var term = PauliTerm.TwoSite(3, 0, PauliLetter.X, 1, PauliLetter.Y, Complex.One);
        Assert.Equal(1, LindbladBitBPiBalance.BitBParity(term));
    }

    // ============================================================
    // IsBitBHomogeneous static helper
    // ============================================================

    [Fact]
    public void IsBitBHomogeneous_SingleTerm_IsTrue()
    {
        // Single-term list is trivially homogeneous: a bit_b = 0 term alone.
        var zeroOnly = new[] {
            PauliTerm.SingleSite(3, 0, PauliLetter.X, Complex.One),
        };
        Assert.True(LindbladBitBPiBalance.IsBitBHomogeneous(zeroOnly));

        // A bit_b = 1 term alone.
        var oneOnly = new[] {
            PauliTerm.SingleSite(3, 0, PauliLetter.Z, Complex.One),
        };
        Assert.True(LindbladBitBPiBalance.IsBitBHomogeneous(oneOnly));
    }

    [Fact]
    public void IsBitBHomogeneous_AllSameBitB_IsTrue()
    {
        // Multi-term, all bit_b = 0: {X⊗I⊗I, I⊗X⊗I, X⊗X⊗I}. bit_b values 0, 0, 0.
        var allZero = new[] {
            PauliTerm.SingleSite(3, 0, PauliLetter.X, Complex.One),
            PauliTerm.SingleSite(3, 1, PauliLetter.X, Complex.One),
            PauliTerm.TwoSite(3, 0, PauliLetter.X, 1, PauliLetter.X, Complex.One),
        };
        Assert.True(LindbladBitBPiBalance.IsBitBHomogeneous(allZero));

        // Multi-term, all bit_b = 1: {Z⊗I⊗I, I⊗Y⊗I, Y⊗Z⊗Z}. bit_b values 1, 1, 1+1+1 = 1.
        var allOne = new[] {
            PauliTerm.SingleSite(3, 0, PauliLetter.Z, Complex.One),
            PauliTerm.SingleSite(3, 1, PauliLetter.Y, Complex.One),
            new PauliTerm(new[] { PauliLetter.Y, PauliLetter.Z, PauliLetter.Z }, Complex.One),
        };
        Assert.True(LindbladBitBPiBalance.IsBitBHomogeneous(allOne));
    }

    [Fact]
    public void IsBitBHomogeneous_MixedBitB_IsFalse()
    {
        // One bit_b = 0 (X⊗I⊗I) + one bit_b = 1 (Z⊗I⊗I) breaks homogeneity.
        var mixed = new[] {
            PauliTerm.SingleSite(3, 0, PauliLetter.X, Complex.One),
            PauliTerm.SingleSite(3, 0, PauliLetter.Z, Complex.One),
        };
        Assert.False(LindbladBitBPiBalance.IsBitBHomogeneous(mixed));
    }

    [Fact]
    public void IsBitBHomogeneous_EmptyList_IsTrueVacuously()
    {
        // Empty list has no terms to disagree; vacuous true per the docstring.
        var empty = Array.Empty<PauliTerm>();
        Assert.True(LindbladBitBPiBalance.IsBitBHomogeneous(empty));
    }

    // ============================================================
    // IsHermitianHamiltonian static helper (scope assertion)
    // ============================================================

    [Fact]
    public void IsHermitianHamiltonian_AnyRealCoefficients_IsTrue()
    {
        // Pauli strings are Hermitian matrices; any real-coefficient linear
        // combination is Hermitian. The method always returns true; its purpose is
        // to make the F112 Tier1Derived scope assertion explicit in code.
        Assert.True(LindbladBitBPiBalance.IsHermitianHamiltonian(new[] { 1.0, 2.0, -0.5 }));
        Assert.True(LindbladBitBPiBalance.IsHermitianHamiltonian(new[] { 0.0 }));
        Assert.True(LindbladBitBPiBalance.IsHermitianHamiltonian(Array.Empty<double>()));
    }
}
