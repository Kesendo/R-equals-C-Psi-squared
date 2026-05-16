using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F81Pi2InheritanceTests
{
    private static F81Pi2Inheritance Build() =>
        new F81Pi2Inheritance(
            new RCPsiSquared.Core.F1.F1PalindromeIdentity(),
            new Pi2DyadicLadderClaim(),
            new Pi2OperatorSpaceMirrorClaim(),
            new Pi2I4MemoryLoopClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void TwoFactor_IsExactlyTwo_FromLadderTermZero()
    {
        // F81's "2" in M − 2·L_{H_odd} is exactly a_0 = 2 = d.
        var f = Build();
        Assert.Equal(2.0, f.TwoFactor, precision: 14);
    }

    [Fact]
    public void HalfFactor_IsExactlyOneHalf_FromLadderTermTwo()
    {
        // F81's "1/2" in the 50/50 split is exactly a_2 = 1/2 = 1/d
        // (= HalfAsStructuralFixedPoint).
        var f = Build();
        Assert.Equal(0.5, f.HalfFactor, precision: 14);
    }

    [Fact]
    public void TwoTimesHalf_IsExactlyOne_InversionSymmetry()
    {
        // The two F81 coefficients are mirror partners on the Pi2 dyadic ladder via
        // a_n · a_{2−n} = 1: the "2" times the "1/2" gives unity exactly.
        var f = Build();
        Assert.Equal(1.0, f.TwoTimesHalf, precision: 14);
    }

    [Theory]
    [InlineData(1, 4.0)]    // d² for 1 qubit
    [InlineData(2, 16.0)]
    [InlineData(3, 64.0)]
    [InlineData(6, 4096.0)]
    public void OperatorSpaceDimension_AgreesWithPi2OperatorSpaceMirrorPinnedTable(int N, double expected)
    {
        // F81's M lives in operator-space d² = 4^N. The dimension is pulled from
        // Pi2OperatorSpaceMirror's pinned table (Mirror Space connection Tom asked about).
        var f = Build();
        Assert.Equal(expected, f.OperatorSpaceDimension(N), precision: 12);
    }

    [Fact]
    public void OperatorSpaceDimension_OutsidePinnedRange_Throws()
    {
        var f = Build();
        Assert.Throws<ArgumentOutOfRangeException>(() => f.OperatorSpaceDimension(7));
    }

    [Fact]
    public void Z4OrderOfPi_IsFour()
    {
        Assert.Equal(4, Build().Z4OrderOfPi);
    }

    [Fact]
    public void MemoryLoopClosure_IsExactlyOne()
    {
        // i⁴ = 1 exactly; drift indicator on the Z₄ generator anchor.
        var f = Build();
        Assert.Equal(1.0, f.MemoryLoopClosure.Real, precision: 14);
        Assert.Equal(0.0, f.MemoryLoopClosure.Imaginary, precision: 14);
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var f1 = new RCPsiSquared.Core.F1.F1PalindromeIdentity();
        var ladder = new Pi2DyadicLadderClaim();
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        Assert.Throws<ArgumentNullException>(() => new F81Pi2Inheritance(null!, ladder, mirror, memoryLoop));
        Assert.Throws<ArgumentNullException>(() => new F81Pi2Inheritance(f1, null!, mirror, memoryLoop));
        Assert.Throws<ArgumentNullException>(() => new F81Pi2Inheritance(f1, ladder, null!, memoryLoop));
        Assert.Throws<ArgumentNullException>(() => new F81Pi2Inheritance(f1, ladder, mirror, null!));
    }

    [Fact]
    public void Anchor_References_F81Proof_AndAllFourPi2Foundation()
    {
        var f = Build();
        Assert.Contains("PROOF_F81_PI_CONJUGATION_OF_M.md", f.Anchor);
        Assert.Contains("F1PalindromeIdentity.cs", f.Anchor);
        Assert.Contains("Pi2DyadicLadderClaim.cs", f.Anchor);
        Assert.Contains("Pi2KnowledgeBaseClaims.cs", f.Anchor);
        Assert.Contains("Pi2OperatorSpaceMirrorClaim.cs", f.Anchor);
        Assert.Contains("Pi2I4MemoryLoopClaim.cs", f.Anchor);
    }
}
