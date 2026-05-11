using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F89PathKVacSeParsevalClaimTests
{
    private static F89PathKVacSeParsevalClaim BuildClaim() =>
        new F89PathKVacSeParsevalClaim(new F89TopologyOrbitClosure(new Pi2DyadicLadderClaim()));

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Theory]
    [InlineData(1, 7, 2.0 * 25.0 / (49.0 * 6.0))]   // (k+1)·(N-k-1)²/(N²(N-1)) = 2·25/(49·6) = 50/294
    [InlineData(2, 7, 3.0 * 16.0 / (49.0 * 6.0))]   // 3·16/(49·6) = 48/294
    [InlineData(3, 7, 4.0 * 9.0 / (49.0 * 6.0))]    // 4·9/(49·6) = 36/294
    [InlineData(1, 5, 2.0 * 9.0 / (25.0 * 4.0))]    // k=1, N=5: 2·9/(25·4) = 18/100
    [InlineData(2, 5, 3.0 * 4.0 / (25.0 * 4.0))]    // k=2, N=5: 3·4/(25·4) = 12/100
    public void Coefficient_MatchesClosedForm(int k, int n, double expected)
    {
        Assert.Equal(expected, F89PathKVacSeParsevalClaim.Coefficient(k, n), precision: 14);
    }

    [Theory]
    [InlineData(1, 7, 0.05, 0.0)]
    [InlineData(2, 7, 0.05, 5.0)]
    [InlineData(3, 7, 0.05, 10.0)]
    public void VacSeBlockClosedForm_Decays4Gamma(int k, int n, double gamma, double t)
    {
        double val_t = F89PathKVacSeParsevalClaim.VacSeBlockClosedForm(k, n, gamma, t);
        double val_0 = F89PathKVacSeParsevalClaim.VacSeBlockClosedForm(k, n, gamma, 0);
        // S(t) / S(0) = exp(-4γt)
        Assert.Equal(Math.Exp(-4.0 * gamma * t), val_t / val_0, precision: 12);
    }

    [Fact]
    public void Coefficient_AtN_EqualsKPlus1_GivesZero()
    {
        // N = k+1 means N_E = 0 (no bare sites; block IS the whole chain).
        // Coefficient = (k+1)·0²/((k+1)²·k) = 0
        for (int k = 1; k <= 5; k++)
            Assert.Equal(0.0, F89PathKVacSeParsevalClaim.Coefficient(k, k + 1), precision: 14);
    }

    [Fact]
    public void Coefficient_KLessThan1_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89PathKVacSeParsevalClaim.Coefficient(0, 5));
    }

    [Fact]
    public void Coefficient_NLessThanKPlus1_Throws()
    {
        // N=4, k=4 → N < k+1 = 5
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89PathKVacSeParsevalClaim.Coefficient(4, 4));
    }

    [Fact]
    public void VacSeBlockClosedForm_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89PathKVacSeParsevalClaim.VacSeBlockClosedForm(1, 5, -0.05, 0));
    }
}
