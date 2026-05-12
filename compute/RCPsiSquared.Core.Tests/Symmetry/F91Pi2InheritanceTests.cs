using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F91Pi2InheritanceTests
{
    private static F91Pi2Inheritance Build() =>
        new F91Pi2Inheritance(new Pi2I4MemoryLoopClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void Z4ClosureOrder_IsExactlyFour()
    {
        // F91's 90°-rotation R_{90} on γ closes at order 4 (i⁴ = 1 on the Z₄ memory loop).
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
    public void ParameterAxis_NamesGammaZDephasing()
    {
        // F91 lives on the per-site Z-dephasing axis γ_l (distinct from F92 J_b and F93 h_l).
        var f = Build();
        Assert.Contains("γ_l", f.ParameterAxis);
        Assert.Contains("Z-dephasing", f.ParameterAxis);
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
            new F91Pi2Inheritance(null!));
    }

    [Fact]
    public void Anchor_References_F91Proof_AndPi2MemoryLoop()
    {
        var f = Build();
        Assert.Contains("PROOF_F91_GAMMA_NINETY_DEGREES.md", f.Anchor);
        Assert.Contains("F71AntiPalindromicGammaSpectralInvariance.cs", f.Anchor);
        Assert.Contains("Pi2I4MemoryLoopClaim.cs", f.Anchor);
        Assert.Contains("NinetyDegreeMirrorMemoryClaim", f.Anchor);
    }
}
