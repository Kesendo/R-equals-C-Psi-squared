using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class YParityIndependenceAtK3Tests
{
    [Fact]
    public void Z2Axis_IsYParity()
    {
        var claim = new YParityIndependenceAtK3();
        Assert.Equal(Z2Axis.YParity, claim.Z2Axis);
    }

    [Fact]
    public void BitATwin_IsNull()
    {
        var claim = new YParityIndependenceAtK3();
        Assert.Null(claim.BitATwin);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var claim = new YParityIndependenceAtK3();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void K2CollapseIdentity_HoldsForAllKBody2Bilinears()
    {
        // At k_body=2, all 9 Pauli bilinears (XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ)
        // satisfy y_par = bit_a XOR bit_b. The identity also holds trivially at
        // k_body=0 (II). At k_body=1 it strictly fails (the 6 k=2-position k_body=1
        // pairs IX, IY, IZ, XI, YI, ZI all violate it). At odd k_body in general,
        // y_par and (bit_a XOR bit_b) differ by 1 mod 2.
        var nonIdentityLetters = new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z };
        foreach (var l1 in nonIdentityLetters)
        {
            foreach (var l2 in nonIdentityLetters)
            {
                var term = new PauliTerm(new[] { l1, l2 }, Complex.One);
                Assert.True(YParityIndependenceAtK3.K2CollapseIdentityHolds(term),
                    $"Identity y_par = bit_a XOR bit_b failed for k_body=2 bilinear ({l1}, {l2}): " +
                    $"bit_a={term.TotalBitA & 1}, bit_b={term.Pi2Parity}, y_par={term.YParity}");
            }
        }

        // k_body=0 (II) trivially satisfies the identity (0 = 0).
        var ii = new PauliTerm(new[] { PauliLetter.I, PauliLetter.I }, Complex.One);
        Assert.True(YParityIndependenceAtK3.K2CollapseIdentityHolds(ii));
    }

    [Fact]
    public void K2CollapseIdentity_FailsAtKBody1()
    {
        // At k_body=1 (single non-I letter embedded in a multi-position tensor),
        // y_par = n_Y mod 2 while bit_a XOR bit_b = (n_X + n_Z) mod 2 = 1 - n_Y mod 2.
        // They differ by 1, so the identity always fails. Parity is position-independent,
        // so both right-position (IX, IY, IZ) and left-position (XI, YI, ZI) k_body=1
        // pairs must violate the identity.
        var nonIdentityLetters = new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z };
        foreach (var letter in nonIdentityLetters)
        {
            var rightPos = new PauliTerm(new[] { PauliLetter.I, letter }, Complex.One);
            Assert.False(YParityIndependenceAtK3.K2CollapseIdentityHolds(rightPos),
                $"k_body=1 right-position (I, {letter}) should violate the k_body=2 collapse identity");

            var leftPos = new PauliTerm(new[] { letter, PauliLetter.I }, Complex.One);
            Assert.False(YParityIndependenceAtK3.K2CollapseIdentityHolds(leftPos),
                $"k_body=1 left-position ({letter}, I) should violate the k_body=2 collapse identity");
        }
    }

    [Fact]
    public void K3_XYZ_vs_III_Counterexample_DemonstratesIndependence()
    {
        // XYZ: (bit_a, bit_b, y_par) = (1+1, 1+1, 1) mod 2 = (0, 0, 1)
        // III: (0, 0, 0)
        // Both Klein (0, 0); y_par differs (1 vs 0).
        var xyz = new PauliTerm(
            new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z }, Complex.One);
        var iii = new PauliTerm(
            new[] { PauliLetter.I, PauliLetter.I, PauliLetter.I }, Complex.One);

        Assert.Equal((0, 0), xyz.KleinIndex);
        Assert.Equal((0, 0), iii.KleinIndex);
        Assert.NotEqual(xyz.YParity, iii.YParity);
        Assert.Equal(1, xyz.YParity);
        Assert.Equal(0, iii.YParity);

        // K2 collapse identity must fail for XYZ (since k_body=3 is odd).
        Assert.False(YParityIndependenceAtK3.K2CollapseIdentityHolds(xyz),
            "XYZ at k_body=3 should violate the k_body=2 collapse identity (y_par=1, bit_a XOR bit_b=0)");
    }
}
