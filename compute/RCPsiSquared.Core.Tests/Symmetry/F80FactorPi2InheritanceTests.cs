using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F80FactorPi2InheritanceTests
{
    private readonly ITestOutputHelper _out;

    public F80FactorPi2InheritanceTests(ITestOutputHelper output) => _out = output;

    private static F80FactorPi2Inheritance Build() =>
        new F80FactorPi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), new Pi2DyadicLadderClaim(), new Pi2I4MemoryLoopClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void TwoFactor_IsExactlyTwo_FromLadderTermZero()
    {
        // The "2" in F80's "±2i" is exactly a_0 = 2 on the Pi2 dyadic ladder
        // (= d, the qubit dimension; the number-anchor side of d=2).
        var f = Build();
        Assert.Equal(2.0, f.TwoFactor, precision: 14);
    }

    [Fact]
    public void IFactor_IsExactlyImaginaryUnit_FromI4LoopPowerOne()
    {
        // The "i" in F80's "±2i" is exactly i^1 on the Z₄ memory loop
        // (= 90° rotation; the angle-anchor side of d=2).
        var f = Build();
        Assert.Equal(new Complex(0, 1), f.IFactor);
    }

    [Fact]
    public void MinusIFactor_IsExactlyNegativeImaginaryUnit_FromI4LoopPowerThree()
    {
        var f = Build();
        Assert.Equal(new Complex(0, -1), f.MinusIFactor);
    }

    [Fact]
    public void PlusTwoIFactor_IsExactlyTwoI()
    {
        // Composition: 2 · i = +2i = (0, 2). Live composition of TwoFactor (Pi2 ladder)
        // and IFactor (Z₄ loop); drift in either parent surfaces here.
        var f = Build();
        var plusTwoI = f.PlusTwoIFactor;
        Assert.Equal(0.0, plusTwoI.Real, precision: 14);
        Assert.Equal(2.0, plusTwoI.Imaginary, precision: 14);
    }

    [Fact]
    public void MinusTwoIFactor_IsExactlyMinusTwoI()
    {
        var f = Build();
        var minusTwoI = f.MinusTwoIFactor;
        Assert.Equal(0.0, minusTwoI.Real, precision: 14);
        Assert.Equal(-2.0, minusTwoI.Imaginary, precision: 14);
    }

    [Fact]
    public void PlusAndMinusTwoI_AreMirrorPartnersInZ4()
    {
        // The ± in F80's "±2i": +2i and −2i are reflections across the real axis,
        // which on the Z₄ loop is the i ↔ i^3 = -i partnership.
        var f = Build();
        Assert.Equal(f.PlusTwoIFactor.Real, f.MinusTwoIFactor.Real, precision: 14);
        Assert.Equal(f.PlusTwoIFactor.Imaginary, -f.MinusTwoIFactor.Imaginary, precision: 14);
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var f1 = new RCPsiSquared.Core.F1.F1PalindromeIdentity();
        var ladder = new Pi2DyadicLadderClaim();
        var loop = new Pi2I4MemoryLoopClaim();
        Assert.Throws<ArgumentNullException>(() => new F80FactorPi2Inheritance(null!, ladder, loop));
        Assert.Throws<ArgumentNullException>(() => new F80FactorPi2Inheritance(f1, null!, loop));
        Assert.Throws<ArgumentNullException>(() => new F80FactorPi2Inheritance(f1, ladder, null!));
    }

    [Fact]
    public void Anchor_References_F80Proof_AnalyticalFormulas_AndBothPi2Parents()
    {
        var f = Build();
        Assert.Contains("PROOF_F80_BLOCH_SIGNWALK.md", f.Anchor);
        Assert.Contains("ANALYTICAL_FORMULAS.md", f.Anchor);
        Assert.Contains("Pi2DyadicLadderClaim", f.Anchor);
        Assert.Contains("Pi2I4MemoryLoopClaim", f.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsFactorDecomposition()
    {
        var f = Build();
        _out.WriteLine("");
        _out.WriteLine("    F80: Spec(M) = ±2i · Spec(H_non-truly)");
        _out.WriteLine("    decomposition via Pi2-Foundation:");
        _out.WriteLine($"      2 = Pi2DyadicLadder.Term(0) = {f.TwoFactor}");
        _out.WriteLine($"      i = Pi2I4MemoryLoop.PowerOfI(1) = ({f.IFactor.Real}, {f.IFactor.Imaginary})");
        _out.WriteLine($"      +2i = ({f.PlusTwoIFactor.Real}, {f.PlusTwoIFactor.Imaginary})");
        _out.WriteLine($"      -2i = ({f.MinusTwoIFactor.Real}, {f.MinusTwoIFactor.Imaginary})");
        _out.WriteLine("");
        _out.WriteLine("    one algebra, two Pi2 axes participating");
    }
}
