using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F93Pi2InheritanceTests
{
    private static F93Pi2Inheritance Build() =>
        new F93Pi2Inheritance(new Pi2I4MemoryLoopClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void Z4ClosureOrder_IsExactlyFour()
    {
        // F93's 90°-rotation R_{90} on h closes at order 4 (i⁴ = 1 on the Z₄ memory loop).
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
    public void ParameterAxis_NamesHLongitudinalDetuning()
    {
        // F93 lives on the per-site longitudinal Z-detuning axis h_l (distinct from F91 γ_l and F92 J_b).
        var f = Build();
        Assert.Contains("h_l", f.ParameterAxis);
        Assert.Contains("longitudinal", f.ParameterAxis);
    }

    [Fact]
    public void Summary_NotesLongitudinalScope()
    {
        // Scope note: only longitudinal h_l Z_l is in scope; transverse breaks joint-popcount.
        var f = Build();
        Assert.False(string.IsNullOrWhiteSpace(f.DisplayName));
        Assert.False(string.IsNullOrWhiteSpace(f.Summary));
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F93Pi2Inheritance(null!));
    }

    [Fact]
    public void Anchor_References_F93Proof_AndPi2MemoryLoop()
    {
        var f = Build();
        Assert.Contains("PROOF_F93_DETUNING_ANTI_PALINDROMIC.md", f.Anchor);
        Assert.Contains("F93DetuningAntiPalindromicSpectralInvariance.cs", f.Anchor);
        Assert.Contains("Pi2I4MemoryLoopClaim.cs", f.Anchor);
        Assert.Contains("NinetyDegreeMirrorMemoryClaim", f.Anchor);
    }
}
