using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F1Pi2InheritanceTests
{
    private static F1Pi2Inheritance Build() =>
        new F1Pi2Inheritance(new Pi2DyadicLadderClaim());

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
    public void Constructor_RejectsNullLadder()
    {
        Assert.Throws<ArgumentNullException>(() => new F1Pi2Inheritance(null!));
    }

    [Fact]
    public void Anchor_References_F1Proof_AndPi2Ladder()
    {
        var f = Build();
        Assert.Contains("MIRROR_SYMMETRY_PROOF.md", f.Anchor);
        Assert.Contains("F1PalindromeIdentity.cs", f.Anchor);
        Assert.Contains("Pi2DyadicLadderClaim.cs", f.Anchor);
    }
}
