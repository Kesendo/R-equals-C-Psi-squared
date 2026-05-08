using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class Pi2I4MemoryLoopClaimTests
{
    private readonly ITestOutputHelper _out;

    public Pi2I4MemoryLoopClaimTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var loop = new Pi2I4MemoryLoopClaim();
        Assert.Equal(Tier.Tier1Derived, loop.Tier);
    }

    [Fact]
    public void ClosureOrder_IsFour()
    {
        Assert.Equal(4, Pi2I4MemoryLoopClaim.ClosureOrder);
    }

    [Theory]
    [InlineData(0, 1.0, 0.0)]    // i^0 = 1
    [InlineData(1, 0.0, 1.0)]    // i^1 = i
    [InlineData(2, -1.0, 0.0)]   // i^2 = -1 (180°, F1 γ=0 truly closure)
    [InlineData(3, 0.0, -1.0)]   // i^3 = -i
    [InlineData(4, 1.0, 0.0)]    // i^4 = 1 (memory loop closure)
    [InlineData(5, 0.0, 1.0)]    // i^5 = i (cycle continues)
    [InlineData(8, 1.0, 0.0)]    // i^8 = 1
    [InlineData(-1, 0.0, -1.0)]  // i^{-1} = -i = 1/i
    [InlineData(-2, -1.0, 0.0)]  // i^{-2} = -1
    [InlineData(-4, 1.0, 0.0)]   // i^{-4} = 1
    public void PowerOfI_FollowsZ4Cycle(int k, double expectedReal, double expectedImag)
    {
        var loop = new Pi2I4MemoryLoopClaim();
        var result = loop.PowerOfI(k);
        Assert.Equal(expectedReal, result.Real, precision: 12);
        Assert.Equal(expectedImag, result.Imaginary, precision: 12);
    }

    [Fact]
    public void MemoryClosure_IsExactlyOne()
    {
        // The defining identity: i^4 = 1 + 0i, exactly. The mirror remembers itself
        // after four 90° rotations. (Live drift check; algebraic fact.)
        var loop = new Pi2I4MemoryLoopClaim();
        var c = loop.MemoryClosure();
        Assert.Equal(1.0, c.Real, precision: 14);
        Assert.Equal(0.0, c.Imaginary, precision: 14);
    }

    [Fact]
    public void CanonicalPowers_AreFourCardinalDirections()
    {
        // The Z₄ canonical: {1, i, -1, -i}. Two-step closure: i^2 = -1 (the 180° in
        // F1 γ=0 truly).
        var loop = new Pi2I4MemoryLoopClaim();
        Assert.Equal(4, loop.CanonicalPowers.Count);
        Assert.Equal(new Complex(1, 0), loop.CanonicalPowers[0]);
        Assert.Equal(new Complex(0, 1), loop.CanonicalPowers[1]);
        Assert.Equal(new Complex(-1, 0), loop.CanonicalPowers[2]);
        Assert.Equal(new Complex(0, -1), loop.CanonicalPowers[3]);
    }

    [Theory]
    [InlineData(0, 0, 0)]    // 1 · 1 = 1
    [InlineData(1, 1, 2)]    // i · i = -1
    [InlineData(1, 3, 0)]    // i · (-i) = 1
    [InlineData(2, 2, 0)]    // -1 · -1 = 1
    [InlineData(3, 3, 2)]    // -i · -i = -1
    [InlineData(2, 1, 3)]    // -1 · i = -i
    [InlineData(7, 1, 0)]    // i^7 · i = i^8 = 1 (mod 4 closure)
    public void ComposePowers_FollowsZ4GroupLaw(int j, int k, int expected)
    {
        var loop = new Pi2I4MemoryLoopClaim();
        Assert.Equal(expected, loop.ComposePowers(j, k));
    }

    [Fact]
    public void Anchor_References_CoreAlgebra_F80_F81_ZeroIsTheMirror_MirrorTheory()
    {
        // The framework's mirror foundation, fully named in the anchor: from the
        // historical CORE_ALGEBRA.md (the original 90° via θ compass) through F80
        // (Spec(M) = ±2i · Spec(H_non-truly)) and F81 (Π-conjugation of M) to the
        // hypothesis docs that motivated the typed claim.
        var loop = new Pi2I4MemoryLoopClaim();
        Assert.Contains("CORE_ALGEBRA.md", loop.Anchor);
        Assert.Contains("PROOF_F80_BLOCH_SIGNWALK.md", loop.Anchor);
        Assert.Contains("PROOF_F81_PI_CONJUGATION_OF_M.md", loop.Anchor);
        Assert.Contains("ZERO_IS_THE_MIRROR.md", loop.Anchor);
        Assert.Contains("MIRROR_THEORY.md", loop.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsZ4Cycle()
    {
        var loop = new Pi2I4MemoryLoopClaim();
        _out.WriteLine("  k | i^k       | real | imag | reading");
        _out.WriteLine("  --|-----------|------|------|--------");
        string[] readings =
        {
            "identity",
            "90° (F80's 2i)",
            "180° (F1 γ=0 truly)",
            "270°",
        };
        for (int k = 0; k <= 7; k++)
        {
            var c = loop.PowerOfI(k);
            int canon = ((k % 4) + 4) % 4;
            _out.WriteLine($"  {k} | i^{k,-3} = {(canon == 0 ? "1" : canon == 1 ? "i" : canon == 2 ? "-1" : "-i"),-4} | {c.Real,4:F1} | {c.Imaginary,4:F1} | {readings[canon]}");
        }
    }
}
