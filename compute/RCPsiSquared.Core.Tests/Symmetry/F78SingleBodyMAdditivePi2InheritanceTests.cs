using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F78SingleBodyMAdditivePi2InheritanceTests
{
    private static F78SingleBodyMAdditivePi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var loop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(ladder, loop);
        return new F78SingleBodyMAdditivePi2Inheritance(ladder, loop, f1);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void EigenvalueCoefficient_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().EigenvalueCoefficient, precision: 14);
    }

    [Fact]
    public void PerSiteDimension_IsExactlyFour()
    {
        Assert.Equal(4.0, BuildClaim().PerSiteDimension, precision: 14);
    }

    [Fact]
    public void ImaginaryUnit_IsI()
    {
        Assert.Equal(0.0, BuildClaim().ImaginaryUnit.Real, precision: 14);
        Assert.Equal(1.0, BuildClaim().ImaginaryUnit.Imaginary, precision: 14);
    }

    [Fact]
    public void TwoIFactor_IsExactly2I()
    {
        Assert.Equal(0.0, BuildClaim().TwoIFactor.Real, precision: 14);
        Assert.Equal(2.0, BuildClaim().TwoIFactor.Imaginary, precision: 14);
    }

    [Theory]
    [InlineData('X')]
    [InlineData('x')]
    public void PerSiteEigenvalues_ForX_AllZero(char p)
    {
        var eig = BuildClaim().PerSiteEigenvalues(p, cl: 1.0, gammaZero: 0.05);
        Assert.Equal(4, eig.Count);
        Assert.All(eig, e => Assert.Equal(Complex.Zero, e));
    }

    [Theory]
    [InlineData('Y')]
    [InlineData('Z')]
    public void PerSiteEigenvalues_ForYorZ_PlusMinus2ClGammaI(char p)
    {
        var eig = BuildClaim().PerSiteEigenvalues(p, cl: 1.0, gammaZero: 0.05);
        Assert.Equal(4, eig.Count);
        // Two +0.1·i, two −0.1·i
        Assert.Equal(2, eig.Count(e => Math.Abs(e.Real) < 1e-12 && Math.Abs(e.Imaginary - 0.1) < 1e-12));
        Assert.Equal(2, eig.Count(e => Math.Abs(e.Real) < 1e-12 && Math.Abs(e.Imaginary + 0.1) < 1e-12));
    }

    [Fact]
    public void IsTruly_TrueOnlyForX()
    {
        var f = BuildClaim();
        Assert.True(f.IsTruly('X'));
        Assert.False(f.IsTruly('Y'));
        Assert.False(f.IsTruly('Z'));
    }

    [Fact]
    public void YEqualsZSvdSpectrum_HoldsAcrossWeights()
    {
        var f = BuildClaim();
        Assert.True(f.YEqualsZSvdSpectrum(cl: 1.0, gammaZero: 0.05));
        Assert.True(f.YEqualsZSvdSpectrum(cl: 2.5, gammaZero: 0.1));
    }

    [Fact]
    public void PerSiteEigenvalues_InvalidPauliLetter_Throws()
    {
        Assert.Throws<ArgumentException>(() => BuildClaim().PerSiteEigenvalues('I', 1.0, 0.05));
        Assert.Throws<ArgumentException>(() => BuildClaim().PerSiteEigenvalues('a', 1.0, 0.05));
    }

    [Fact]
    public void PerSiteEigenvalues_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().PerSiteEigenvalues('Y', 1.0, -0.05));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var loop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(ladder, loop);
        Assert.Throws<ArgumentNullException>(() =>
            new F78SingleBodyMAdditivePi2Inheritance(null!, loop, f1));
    }

    [Fact]
    public void Constructor_NullLoop_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var loop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(ladder, loop);
        Assert.Throws<ArgumentNullException>(() =>
            new F78SingleBodyMAdditivePi2Inheritance(ladder, null!, f1));
    }

    [Fact]
    public void Constructor_NullF1_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var loop = new Pi2I4MemoryLoopClaim();
        Assert.Throws<ArgumentNullException>(() =>
            new F78SingleBodyMAdditivePi2Inheritance(ladder, loop, null!));
    }
}
