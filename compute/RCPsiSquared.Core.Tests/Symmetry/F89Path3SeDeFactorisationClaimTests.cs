using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F89Path3SeDeFactorisationClaimTests
{
    private static F89Path3SeDeFactorisationClaim BuildClaim()
    {
        var f89 = new F89TopologyOrbitClosure(new Pi2DyadicLadderClaim());
        var atLock = new F89PathKAtLockMechanismClaim(f89);
        return new F89Path3SeDeFactorisationClaim(f89, atLock);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void S2SymSubBlockDimension_Is12()
    {
        Assert.Equal(12, F89Path3SeDeFactorisationClaim.S2SymSubBlockDimension);
    }

    [Fact]
    public void FactorDegrees_Are2_2_8()
    {
        var degs = F89Path3SeDeFactorisationClaim.FactorDegrees;
        Assert.Equal(new[] { 2, 2, 8 }, degs);
        Assert.Equal(12, degs.Sum());
    }

    [Theory]
    [InlineData(0.05, 0.075)]   // γ=0.05, J=0.075 (q=1.5)
    [InlineData(1.0, 1.5)]      // γ=1, J=1.5
    [InlineData(0.1, 0.05)]     // γ=0.1, J=0.05
    public void FaRoots_Match_NegativeTwoGammaPlusIJTimesMinusOnePlusMinusSqrt5(double gamma, double j)
    {
        var roots = F89Path3SeDeFactorisationClaim.FaRoots(gamma, j);
        Assert.Equal(2, roots.Length);

        double sqrt5 = Math.Sqrt(5);
        // Sort by imaginary part to make ordering deterministic
        var sorted = roots.OrderBy(r => r.Imaginary).ToArray();
        Assert.Equal(-2 * gamma, sorted[0].Real, precision: 12);
        Assert.Equal(j * (-1 - sqrt5), sorted[0].Imaginary, precision: 12);
        Assert.Equal(-2 * gamma, sorted[1].Real, precision: 12);
        Assert.Equal(j * (-1 + sqrt5), sorted[1].Imaginary, precision: 12);
    }

    [Theory]
    [InlineData(0.05, 0.075)]
    [InlineData(1.0, 1.5)]
    public void FbRoots_Match_NegativeSixGammaPlusIJTimesMinusOnePlusMinusSqrt5(double gamma, double j)
    {
        var roots = F89Path3SeDeFactorisationClaim.FbRoots(gamma, j);
        Assert.Equal(2, roots.Length);

        double sqrt5 = Math.Sqrt(5);
        var sorted = roots.OrderBy(r => r.Imaginary).ToArray();
        Assert.Equal(-6 * gamma, sorted[0].Real, precision: 12);
        Assert.Equal(j * (-1 - sqrt5), sorted[0].Imaginary, precision: 12);
        Assert.Equal(-6 * gamma, sorted[1].Real, precision: 12);
        Assert.Equal(j * (-1 + sqrt5), sorted[1].Imaginary, precision: 12);
    }

    [Fact]
    public void OcticIsIrreducibleOverQiSqrt5_AsDocumentedFact()
    {
        // The octic F_8 is irreducible over Q[i, √5] per
        // simulations/_f89_path3_octic_factor_test.py (sympy verification).
        // Documented external fact (Tier 1 derived).
        Assert.True(F89Path3SeDeFactorisationClaim.OcticIsIrreducibleOverQiSqrt5);
    }

    [Fact]
    public void FaRoots_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89Path3SeDeFactorisationClaim.FaRoots(-0.05, 0.075));
    }

    [Fact]
    public void FaRoots_NegativeJ_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89Path3SeDeFactorisationClaim.FaRoots(0.05, -0.075));
    }

    [Fact]
    public void FbRoots_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89Path3SeDeFactorisationClaim.FbRoots(-0.05, 0.075));
    }

    [Fact]
    public void FbRoots_NegativeJ_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89Path3SeDeFactorisationClaim.FbRoots(0.05, -0.075));
    }
}
