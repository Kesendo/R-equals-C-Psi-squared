using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F78SingleBodyMAdditivePi2InheritanceTests
{
    private static F78SingleBodyMAdditivePi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var loop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, loop);
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
        var eig = BuildClaim().PerSiteEigenvalues(p, cl: 1.0);
        Assert.Equal(4, eig.Count);
        Assert.All(eig, e => Assert.Equal(Complex.Zero, e));
    }

    [Theory]
    [InlineData('Y')]
    [InlineData('Z')]
    public void PerSiteEigenvalues_ForYorZ_PlusMinus2ClI(char p)
    {
        var eig = BuildClaim().PerSiteEigenvalues(p, cl: 1.5);
        Assert.Equal(4, eig.Count);
        // Two +3·i, two −3·i: 2·c_l with c_l = 1.5, exact in floating point.
        Assert.Equal(2, eig.Count(e => e.Real == 0.0 && e.Imaginary == 3.0));
        Assert.Equal(2, eig.Count(e => e.Real == 0.0 && e.Imaginary == -3.0));
    }

    [Theory]
    [InlineData('X', 0.05)]
    [InlineData('Y', 0.05)]
    [InlineData('Z', 0.05)]
    [InlineData('X', 0.37)]
    [InlineData('Y', 0.37)]
    [InlineData('Z', 0.37)]
    public void PerSiteEigenvalues_MatchTheBuiltResidual_AtTwoRates(char p, double gamma)
    {
        // One site, H = c·P under Z-dephasing at rate γ: build M = Π·L·Π⁻¹ + L + 2σ·I and compare its
        // spectrum with the closed form, without an eigensolver: for a 4×4 matrix the power sums
        // tr(M^k), k = 1..4, fix the eigenvalue multiset (Newton's identities). M does not depend on γ
        // (Master Lemma), so the closed form, which takes no rate, must match it at both rates. The
        // comparison is exact: every product in this one-site construction multiplies a γ or c term by
        // 0, ±1, ±i or a power of two, and no sum adds a γ term to a c term within one real or
        // imaginary part, so the γ entries cancel to 0.0; with the dyadic c = 1.5 every power and sum
        // of M is an integer below 2^53. This rests on the conventional four-multiplication complex
        // product (a 3M kernel would leave rounding on the Z rows).
        const double c = 1.5;
        var letter = p switch { 'X' => PauliLetter.X, 'Y' => PauliLetter.Y, _ => PauliLetter.Z };
        var H = PauliString.SiteOp(1, 0, letter).Multiply(new Complex(c, 0));
        var L = PauliDephasingDissipator.BuildZ(H, new[] { gamma });
        var M = PalindromeResidual.Build(L, 1, gamma, PauliLetter.Z);
        if (p == 'X')
            Assert.True(M.Enumerate().All(z => z == Complex.Zero),
                $"X at γ = {gamma}: M_l = 0 (the truly case) expected entry by entry");
        var predicted = BuildClaim().PerSiteEigenvalues(p, c);
        var power = M;
        for (int k = 1; k <= 4; k++)
        {
            Complex measured = power.Trace();
            Complex expected = Complex.Zero;
            foreach (var lambda in predicted)
            {
                Complex term = Complex.One;
                for (int j = 0; j < k; j++) term *= lambda;
                expected += term;
            }
            Assert.True(measured == expected,
                $"{p} at γ = {gamma}: tr(M^{k}) = {measured}, closed form Σλ^{k} = {expected}");
            power = power * M;
        }
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
        Assert.True(f.YEqualsZSvdSpectrum(cl: 1.0));
        Assert.True(f.YEqualsZSvdSpectrum(cl: 2.5));
    }

    [Fact]
    public void PerSiteEigenvalues_InvalidPauliLetter_Throws()
    {
        Assert.Throws<ArgumentException>(() => BuildClaim().PerSiteEigenvalues('I', 1.0));
        Assert.Throws<ArgumentException>(() => BuildClaim().PerSiteEigenvalues('a', 1.0));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var loop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, loop);
        Assert.Throws<ArgumentNullException>(() =>
            new F78SingleBodyMAdditivePi2Inheritance(null!, loop, f1));
    }

    [Fact]
    public void Constructor_NullLoop_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var loop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, loop);
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
