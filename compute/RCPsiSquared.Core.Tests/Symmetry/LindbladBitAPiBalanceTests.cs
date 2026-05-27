using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>F112-X Claim structure tests: BitA axis, Tier1Derived universal N, and the
/// static helpers BitAParity / IsBitAHomogeneous that encode the bit_a-homogeneity
/// hypothesis on the dissipator operators c_k. Mirrors
/// <see cref="LindbladBitBPiBalanceTests"/> (if present) on the BitB axis.</summary>
public class LindbladBitAPiBalanceTests
{
    private static LindbladBitAPiBalance Make() =>
        new LindbladBitAPiBalance(new F108Part2Pi2XEvenAlwaysPalindromic());

    // ============================================================
    // Claim metadata
    // ============================================================

    [Fact]
    public void Z2Axis_Is_BitA() =>
        Assert.Equal(Z2Axis.BitA, Make().Z2Axis);

    [Fact]
    public void Tier_Is_Tier1Derived() =>
        Assert.Equal(Tier.Tier1Derived, Make().Tier);

    [Fact]
    public void BitATwin_Is_Null_On_BitA_Axis() =>
        Assert.Null(Make().BitATwin);

    [Fact]
    public void BitATwinStatus_Is_NotApplicableForThisAxis() =>
        Assert.Equal(BitATwinClassification.NotApplicableForThisAxis, Make().BitATwinStatus);

    [Fact]
    public void Theorem_Mentions_BitA_Homogeneity_And_PiX_Polarity()
    {
        string theorem = Make().Theorem;
        Assert.Contains("bit_a-homogeneous", theorem);
        Assert.Contains("Π_X", theorem);
        Assert.Contains("M_plus_half", theorem);
        Assert.Contains("M_minus_half", theorem);
    }

    [Fact]
    public void Part2_Parent_Is_Wired()
    {
        var claim = Make();
        Assert.NotNull(claim.Part2);
        Assert.IsType<F108Part2Pi2XEvenAlwaysPalindromic>(claim.Part2);
    }

    // ============================================================
    // BitAParity helper
    // ============================================================

    [Fact]
    public void BitAParity_SingleSite_X_Is_One()
    {
        // X has bit_a = 1, single-site at N=2: bit_a sum = 1.
        var x = PauliTerm.SingleSite(2, 0, PauliLetter.X, Complex.One);
        Assert.Equal(1, LindbladBitAPiBalance.BitAParity(x));
    }

    [Fact]
    public void BitAParity_SingleSite_Z_Is_Zero()
    {
        // Z has bit_a = 0, single-site at N=2: bit_a sum = 0.
        var z = PauliTerm.SingleSite(2, 0, PauliLetter.Z, Complex.One);
        Assert.Equal(0, LindbladBitAPiBalance.BitAParity(z));
    }

    [Fact]
    public void BitAParity_TwoSite_XX_Is_Zero()
    {
        // XX: bit_a sum = 1 + 1 = 2 → 0 mod 2.
        var xx = PauliTerm.TwoSite(2, 0, PauliLetter.X, 1, PauliLetter.X, Complex.One);
        Assert.Equal(0, LindbladBitAPiBalance.BitAParity(xx));
    }

    [Fact]
    public void BitAParity_TwoSite_XZ_Is_One()
    {
        // XZ: bit_a sum = 1 + 0 = 1.
        var xz = PauliTerm.TwoSite(2, 0, PauliLetter.X, 1, PauliLetter.Z, Complex.One);
        Assert.Equal(1, LindbladBitAPiBalance.BitAParity(xz));
    }

    [Fact]
    public void BitAParity_Null_Throws() =>
        Assert.Throws<ArgumentNullException>(() => LindbladBitAPiBalance.BitAParity(null!));

    // ============================================================
    // IsBitAHomogeneous helper
    // ============================================================

    [Fact]
    public void IsBitAHomogeneous_EmptyList_IsVacuouslyTrue()
    {
        Assert.True(LindbladBitAPiBalance.IsBitAHomogeneous(Array.Empty<PauliTerm>()));
    }

    [Fact]
    public void IsBitAHomogeneous_SingleTerm_IsTrue()
    {
        var single = new[] { PauliTerm.TwoSite(2, 0, PauliLetter.X, 1, PauliLetter.Y, Complex.One) };
        Assert.True(LindbladBitAPiBalance.IsBitAHomogeneous(single));
    }

    [Fact]
    public void IsBitAHomogeneous_TwoTerms_BothBitAEven_IsTrue()
    {
        // XX (bit_a=0) and ZZ (bit_a=0): both bit_a-even.
        var terms = new[]
        {
            PauliTerm.TwoSite(2, 0, PauliLetter.X, 1, PauliLetter.X, Complex.One),
            PauliTerm.TwoSite(2, 0, PauliLetter.Z, 1, PauliLetter.Z, Complex.One),
        };
        Assert.True(LindbladBitAPiBalance.IsBitAHomogeneous(terms));
    }

    [Fact]
    public void IsBitAHomogeneous_TwoTerms_BothBitAOdd_IsTrue()
    {
        // XZ (bit_a=1) and ZX (bit_a=1): both bit_a-odd.
        var terms = new[]
        {
            PauliTerm.TwoSite(2, 0, PauliLetter.X, 1, PauliLetter.Z, Complex.One),
            PauliTerm.TwoSite(2, 0, PauliLetter.Z, 1, PauliLetter.X, Complex.One),
        };
        Assert.True(LindbladBitAPiBalance.IsBitAHomogeneous(terms));
    }

    [Fact]
    public void IsBitAHomogeneous_MixedParity_IsFalse()
    {
        // XX (bit_a=0) and XZ (bit_a=1): mixed parity.
        var terms = new[]
        {
            PauliTerm.TwoSite(2, 0, PauliLetter.X, 1, PauliLetter.X, Complex.One),
            PauliTerm.TwoSite(2, 0, PauliLetter.X, 1, PauliLetter.Z, Complex.One),
        };
        Assert.False(LindbladBitAPiBalance.IsBitAHomogeneous(terms));
    }

    [Fact]
    public void IsBitAHomogeneous_Null_Throws() =>
        Assert.Throws<ArgumentNullException>(() => LindbladBitAPiBalance.IsBitAHomogeneous(null!));
}
