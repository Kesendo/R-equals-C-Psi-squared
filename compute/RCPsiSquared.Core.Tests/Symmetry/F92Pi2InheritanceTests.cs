using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F92Pi2InheritanceTests
{
    private static F92Pi2Inheritance Build() =>
        new F92Pi2Inheritance(new Pi2I4MemoryLoopClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void Z4ClosureOrder_IsExactlyFour()
    {
        // F92's 90°-rotation R_{90} on J closes at order 4 (i⁴ = 1 on the Z₄ memory loop).
        Assert.Equal(4, Build().Z4ClosureOrder);
    }

    [Fact]
    public void MemoryLoopClosure_IsExactlyOne()
    {
        // i⁴ = 1 exactly; drift indicator on the Z₄ generator parent.
        var f = Build();
        Assert.Equal(1.0, f.MemoryLoopClosure.Real, precision: 14);
        Assert.Equal(0.0, f.MemoryLoopClosure.Imaginary, precision: 14);
    }

    [Fact]
    public void ParameterAxis_NamesJBondCoupling()
    {
        // F92 lives on the per-bond XY-coupling axis J_b (distinct from F91 γ_l and F93 h_l).
        var f = Build();
        Assert.Contains("J_b", f.ParameterAxis);
        Assert.Contains("bond", f.ParameterAxis);
    }

    [Fact]
    public void DisplayNameAndSummary_AreNonEmpty()
    {
        var f = Build();
        Assert.False(string.IsNullOrWhiteSpace(f.DisplayName));
        Assert.False(string.IsNullOrWhiteSpace(f.Summary));
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F92Pi2Inheritance(null!));
    }

    [Fact]
    public void Anchor_References_F92Proof_AndPi2MemoryLoop()
    {
        var f = Build();
        Assert.Contains("PROOF_F92_BOND_ANTI_PALINDROMIC_J.md", f.Anchor);
        Assert.Contains("F92BondAntiPalindromicJSpectralInvariance.cs", f.Anchor);
        Assert.Contains("Pi2I4MemoryLoopClaim.cs", f.Anchor);
        Assert.Contains("NinetyDegreeMirrorMemoryClaim", f.Anchor);
    }
}
