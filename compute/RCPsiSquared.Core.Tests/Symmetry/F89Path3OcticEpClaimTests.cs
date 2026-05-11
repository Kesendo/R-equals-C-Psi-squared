using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F89Path3OcticEpClaimTests
{
    private static F89Path3OcticEpClaim BuildClaim()
    {
        var f89 = new F89TopologyOrbitClosure(new Pi2DyadicLadderClaim());
        var atLock = new F89PathKAtLockMechanismClaim(f89);
        return new F89Path3OcticEpClaim(f89, atLock);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void QEpSquared_IsMinusOnePlusSqrt13Over6()
    {
        // q² = (-1 + √13) / 6 from the (3q⁴+q²-1)² discriminant factor zero
        double expected = (-1.0 + Math.Sqrt(13)) / 6.0;
        Assert.Equal(expected, F89Path3OcticEpClaim.QEpSquared, precision: 14);
    }

    [Fact]
    public void QEp_IsApproximately_0_658983()
    {
        Assert.Equal(0.658983, F89Path3OcticEpClaim.QEp, precision: 5);
    }

    [Fact]
    public void DiscriminantFactor_AtQEp_IsZero()
    {
        // 3q⁴ + q² - 1 = 0 at q_EP
        double q = F89Path3OcticEpClaim.QEp;
        double residual = 3.0 * q * q * q * q + q * q - 1.0;
        Assert.True(Math.Abs(residual) < 1e-14, $"3q⁴+q²-1 at q_EP: {residual} should be ~0");
    }

    [Theory]
    [InlineData(0.05, 0.05 * 0.658983)]
    [InlineData(1.0, 0.658983)]
    [InlineData(1.0, 0.0)]   // J=0 degenerate: λ_EP = -4γ + 0i (no oscillation)
    public void MergedEigenvalueLambdaEp_IsMinusFourGammaPlusTwoIJ(double gamma, double j)
    {
        var lam = F89Path3OcticEpClaim.MergedEigenvalue(gamma, j);
        Assert.Equal(-4 * gamma, lam.Real, precision: 12);
        Assert.Equal(2 * j, lam.Imaginary, precision: 12);
    }

    [Fact]
    public void MergedEigenvalueRate_IsFourGamma_AtAtSpectralMidpoint()
    {
        // Re(λ_EP)/γ = -4 sits at the AT-spectral midpoint of rates 2γ (overlap) and 6γ (no-overlap)
        // (-2 + -6) / 2 = -4 ✓
        var lam = F89Path3OcticEpClaim.MergedEigenvalue(1.0, 0.658983);
        Assert.Equal(-4.0, lam.Real, precision: 12);
        double midpoint = (-2.0 + -6.0) / 2.0;
        Assert.Equal(midpoint, lam.Real, precision: 12);
    }

    [Fact]
    public void MergedEigenvalueAtQEp_OmegaIsTwoTimesQEp()
    {
        // At J = QEp (the EP location in q = J/γ units), Im(λ_EP) should equal 2·QEp.
        // This tests consistency between the QEp constant and the MergedEigenvalue formula
        // (catches any drift between the analytical EP location and the merged-eigenvalue claim).
        var lam = F89Path3OcticEpClaim.MergedEigenvalue(1.0, F89Path3OcticEpClaim.QEp);
        Assert.Equal(2.0 * F89Path3OcticEpClaim.QEp, lam.Imaginary, precision: 12);
    }

    [Fact]
    public void MergedEigenvalue_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89Path3OcticEpClaim.MergedEigenvalue(-0.05, 0.075));
    }

    [Fact]
    public void MergedEigenvalue_NegativeJ_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89Path3OcticEpClaim.MergedEigenvalue(0.05, -0.075));
    }
}
