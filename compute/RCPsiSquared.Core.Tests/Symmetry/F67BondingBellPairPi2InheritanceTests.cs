using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F67BondingBellPairPi2InheritanceTests
{
    private static F67BondingBellPairPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var qubitAnchor = new QubitDimensionalAnchorClaim();
        var f66 = new F66PoleModesPi2Inheritance(ladder, qubitAnchor);
        var f65 = new F65XxChainSpectrumPi2Inheritance(ladder, f66);
        return new F67BondingBellPairPi2Inheritance(ladder, f65);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void NumeratorCoefficient_IsExactlyFour()
    {
        // 4 = a_{-1} on the dyadic ladder, transitively via F65.
        Assert.Equal(4.0, BuildClaim().NumeratorCoefficient, precision: 14);
    }

    [Theory]
    // F67's α_1 = F65's single-excitation rate at k=1.
    // N=3: α_1/γ₀ = (4/4)·sin²(π/4) = 1·(1/2) = 1/2 ✓
    // N=5: α_1/γ₀ = (4/6)·sin²(π/6) = (2/3)·(1/4) = 1/6 ✓
    [InlineData(3, 1.0, 0.5)]
    [InlineData(5, 1.0, 1.0 / 6.0)]
    [InlineData(3, 0.5, 0.25)]      // half γ
    [InlineData(5, 0.1, 1.0 / 60.0)]
    public void BondingModeDecayRate_MatchesClosedForm(int N, double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().BondingModeDecayRate(N, gammaZero), precision: 12);
    }

    [Fact]
    public void T2_AtN3GammaOne_IsTwo()
    {
        // T_2 = 1/α_1 at N=3, γ=1: α_1 = 1/2 → T_2 = 2.
        Assert.Equal(2.0, BuildClaim().T2(N: 3, gammaZero: 1.0), precision: 12);
    }

    [Fact]
    public void T2_AtN5GammaOne_IsSix()
    {
        // T_2 = 1/α_1 at N=5, γ=1: α_1 = 1/6 → T_2 = 6.
        Assert.Equal(6.0, BuildClaim().T2(N: 5, gammaZero: 1.0), precision: 12);
    }

    [Fact]
    public void AsymptoticT2_LargeN_ApproachesActualT2()
    {
        // T_2 → (N+1)³/(4π²γ) for large N. At N=100, asymptotic and exact differ by ~1%.
        var f = BuildClaim();
        double exact = f.T2(N: 100, gammaZero: 1.0);
        double asympt = f.AsymptoticT2(N: 100, gammaZero: 1.0);
        // Exact = (N+1)/(4γ·sin²(π/(N+1))). At large N, sin(π/(N+1)) ≈ π/(N+1) so
        // exact ≈ (N+1)/(4γ·π²/(N+1)²) = (N+1)³/(4π²γ) = asymptotic. Convergence ~1/N².
        Assert.True(Math.Abs(exact - asympt) / asympt < 0.01,
            $"asymptotic-exact mismatch: exact={exact}, asympt={asympt}");
    }

    [Theory]
    // F65's bonding-mode amplitude |ψ_1(j)|² = (2/(N+1))·sin²(π(j+1)/(N+1)) at k=1.
    // N=3 j=0: (2/4)·sin²(π/4) = (1/2)·(1/2) = 1/4
    // N=3 j=1: (2/4)·sin²(π/2) = (1/2)·1 = 1/2
    // N=3 j=2: (2/4)·sin²(3π/4) = (1/2)·(1/2) = 1/4 (palindromic equal to j=0)
    [InlineData(3, 0, 0.25)]
    [InlineData(3, 1, 0.5)]
    [InlineData(3, 2, 0.25)]
    public void BondingModeAmplitude_MatchesClosedForm(int N, int site, double expected)
    {
        Assert.Equal(expected, BuildClaim().BondingModeAmplitude(N, site), precision: 12);
    }

    [Fact]
    public void PalindromicAEqualsCAtBondingMode_HoldsForN3()
    {
        Assert.True(BuildClaim().PalindromicAEqualsCAtBondingMode(N: 3));
    }

    [Fact]
    public void PalindromicAEqualsCAtBondingMode_HoldsForN5()
    {
        Assert.True(BuildClaim().PalindromicAEqualsCAtBondingMode(N: 5));
    }

    [Fact]
    public void PalindromicAEqualsCAtBondingMode_HoldsForN7()
    {
        Assert.True(BuildClaim().PalindromicAEqualsCAtBondingMode(N: 7));
    }

    [Fact]
    public void BondingModeDecayRate_NLessThanTwo_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            BuildClaim().BondingModeDecayRate(N: 1, gammaZero: 1.0));
    }

    [Fact]
    public void T2_NonPositiveGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            BuildClaim().T2(N: 3, gammaZero: 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            BuildClaim().T2(N: 3, gammaZero: -0.1));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f66 = new F66PoleModesPi2Inheritance(ladder, new QubitDimensionalAnchorClaim());
        var f65 = new F65XxChainSpectrumPi2Inheritance(ladder, f66);
        Assert.Throws<ArgumentNullException>(() =>
            new F67BondingBellPairPi2Inheritance(null!, f65));
    }

    [Fact]
    public void Constructor_NullF65_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F67BondingBellPairPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }
}
