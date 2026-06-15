using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live Π-level witness (inspect --root mirrorgroup): the palindromizer factorization
/// Π_Z = R·D, the dihedral inversion D·Π_Z·D = Π_Y (reflections/D_PI_Z_EQUALS_PI_Y), the D₄ closure,
/// the §3 palindrome split WITH the −2Σγ shift, and the palindrome product Π·L·Π⁻¹ = −L − 2Σγ.
/// The first two are the GATE: they fire if Π_coh = M·PiOperator·M⁻¹ and R·D are built in mismatched
/// vec conventions (the exact (−1)^{n_Y} twist of the origin reflection).</summary>
public class MirrorGroupWitnessTests
{
    // ---- the gate: the convention bridge (vec_F/vec_R). Pass = the representation is consistent. ----
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void factorization_Pi_Z_equals_R_times_D(int n)
    {
        var w = new MirrorGroupWitness(n);
        Assert.True(w.FactorizationDev() < 1e-9, $"Π_Z = R·D failed at N={n}: dev {w.FactorizationDev()}");
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void dihedral_inversion_D_Pi_Z_D_equals_Pi_Y(int n)
    {
        // reflections/D_PI_Z_EQUALS_PI_Y, live: D·Π_Z·D = Π_Y = Π_Z⁻¹ (the dihedral inversion s·r·s = r⁻¹)
        var w = new MirrorGroupWitness(n);
        Assert.True(w.DihedralInversionDev() < 1e-9, $"D·Π_Z·D = Π_Y failed at N={n}: dev {w.DihedralInversionDev()}");
    }

    // ---- the group ----
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void mirror_group_closes_to_order_8(int n)
    {
        var w = new MirrorGroupWitness(n);
        Assert.Equal(8, w.GroupOrder());   // ⟨R, D⟩ ≅ D₄
    }

    // ---- the §3 palindrome split, generator-by-generator, WITH the −2Σγ shift on the dissipator rows ----
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void palindrome_splits_along_generators_with_shift(int n)
    {
        var w = new MirrorGroupWitness(n);
        var (dLh, dLdiss, rLh, rLdiss) = w.PalindromeSplitDevs();
        Assert.True(dLh < 1e-9, $"D·L_H·D = −L_H failed: dev {dLh}");
        Assert.True(dLdiss < 1e-9, $"D·L_diss·D = +L_diss failed: dev {dLdiss}");
        Assert.True(rLh < 1e-9, $"R·L_H·R = +L_H failed: dev {rLh}");
        Assert.True(rLdiss < 1e-9, $"R·L_diss·R = −L_diss − 2Σγ·I failed: dev {rLdiss}");
    }

    // ---- the palindrome as the closure of the factorization (the spectral form lives in MirrorSystem) ----
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void palindrome_product_is_minus_L_minus_2_sigma_gamma(int n)
    {
        var w = new MirrorGroupWitness(n);
        Assert.True(w.PalindromeProductDev() < 1e-9,
            $"Π_Z·L·Π_Z⁻¹ = −L − 2Σγ·I failed at N={n}: dev {w.PalindromeProductDev()}");
    }

    [Fact]
    public void guards_against_too_large_N()
    {
        // 4^6 = 4096 > MaxDim (1024)
        Assert.Throws<ArgumentOutOfRangeException>(() => new MirrorGroupWitness(6));
    }
}
