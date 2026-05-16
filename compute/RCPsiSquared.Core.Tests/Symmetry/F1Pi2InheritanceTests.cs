using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F1Pi2InheritanceTests
{
    private static F1Pi2Inheritance Build() =>
        new F1Pi2Inheritance(
            new RCPsiSquared.Core.F1.F1PalindromeIdentity(),
            new Pi2DyadicLadderClaim(),
            new Pi2I4MemoryLoopClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void TwoFactor_IsExactlyTwo_FromLadderTermZero()
    {
        // F1's "2" coefficient in -L - 2Σγ·I is exactly a_0 = 2 on the Pi2 dyadic
        // ladder (= d, qubit dimension). Not a free parameter.
        var f = Build();
        Assert.Equal(2.0, f.TwoFactor, precision: 14);
    }

    [Fact]
    public void SignFlipFromZ4_IsExactlyMinusOne()
    {
        // F1's "−L" sign flip = i² = −1 on the Z₄ memory loop (Pi2I4MemoryLoop Layer 1).
        var f = Build();
        Assert.Equal(-1.0, f.SignFlipFromZ4.Real, precision: 14);
        Assert.Equal(0.0, f.SignFlipFromZ4.Imaginary, precision: 14);
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var f1 = new RCPsiSquared.Core.F1.F1PalindromeIdentity();
        var ladder = new Pi2DyadicLadderClaim();
        var loop = new Pi2I4MemoryLoopClaim();
        Assert.Throws<ArgumentNullException>(() =>
            new F1Pi2Inheritance(null!, ladder, loop));
        Assert.Throws<ArgumentNullException>(() =>
            new F1Pi2Inheritance(f1, null!, loop));
        Assert.Throws<ArgumentNullException>(() =>
            new F1Pi2Inheritance(f1, ladder, null!));
    }

    [Fact]
    public void TypedParents_F1IsExposed()
    {
        Assert.NotNull(Build().F1);
    }

    [Fact]
    public void Anchor_References_F1Proof_AndPi2Ladder_AndMemoryLoop()
    {
        var f = Build();
        Assert.Contains("MIRROR_SYMMETRY_PROOF.md", f.Anchor);
        Assert.Contains("F1PalindromeIdentity.cs", f.Anchor);
        Assert.Contains("Pi2DyadicLadderClaim.cs", f.Anchor);
        Assert.Contains("Pi2I4MemoryLoopClaim.cs", f.Anchor);
    }
}
