using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F44AlgebraicPairRateLogIdentityPi2InheritanceTests
{
    private static F44AlgebraicPairRateLogIdentityPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, memoryLoop);
        return new F44AlgebraicPairRateLogIdentityPi2Inheritance(ladder, f1);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void SumCoefficient_IsExactlyTwo()
    {
        // 2 = a_0 on dyadic ladder (palindromic-pair sum coefficient).
        Assert.Equal(2.0, BuildClaim().SumCoefficient, precision: 14);
    }

    [Fact]
    public void ArTanhCoefficient_IsExactlyTwo()
    {
        // 2 in 2·artanh(...) — from ln((1+x)/(1-x)) = 2·artanh(x) identity.
        Assert.Equal(2.0, BuildClaim().ArTanhCoefficient, precision: 14);
    }

    [Theory]
    // F44 closed form: ln(d_fast/d_slow) = 2·artanh(Δd/(2·Σγ))
    // Test: F68 N=3 γ=0.05 example: d_fast = 0.075, d_slow = 0.025, Δd = 0.05, Σγ = 0.05
    // → arg = 0.05/(2·0.05) = 0.5; 2·artanh(0.5) = 2·0.5493 = 1.0986; ln(0.075/0.025) = ln(3) ≈ 1.0986
    [InlineData(0.05, 0.05, 1.0986122886681098)]   // ln(3) at F68 N=3 example
    [InlineData(0.0, 0.05, 0.0)]                   // Δd=0 → ratio = 1 → log = 0
    [InlineData(0.05, 0.10, 0.5108256237659907)]   // arg = 0.25; 2·artanh(0.25) ≈ 0.5108
    public void LogRatio_MatchesClosedForm(double deltaD, double totalGamma, double expected)
    {
        Assert.Equal(expected, BuildClaim().LogRatio(deltaD, totalGamma), precision: 12);
    }

    [Fact]
    public void LogRatio_ConsistentWithDirectComputation()
    {
        // Direct verification: pick d_fast, d_slow, compute LogRatio, compare to ln(d_fast/d_slow).
        // For d_fast = 0.075, d_slow = 0.025: Σγ = (0.075 + 0.025)/2 = 0.05; Δd = 0.05.
        var f = BuildClaim();
        double dFast = 0.075;
        double dSlow = 0.025;
        double totalGamma = (dFast + dSlow) / 2.0;   // palindromic pair sum
        double deltaD = dFast - dSlow;
        double directLog = Math.Log(dFast / dSlow);
        double formulaLog = f.LogRatio(deltaD, totalGamma);
        Assert.Equal(directLog, formulaLog, precision: 12);
    }

    [Theory]
    // Linear approximation: ln(d_fast/d_slow) ≈ Δd/Σγ for small Δd/Σγ.
    [InlineData(0.001, 0.05, 0.02)]   // small ratio: linear close to exact
    [InlineData(0.01, 0.05, 0.2)]
    public void LogRatioLinearApproximation_GivesDeltaOverGamma(double deltaD, double totalGamma, double expected)
    {
        Assert.Equal(expected, BuildClaim().LogRatioLinearApproximation(deltaD, totalGamma), precision: 12);
    }

    [Fact]
    public void LinearApproximation_ApproachesExactAtSmallDelta()
    {
        // For Δd/Σγ ≪ 1, the linear approximation should match LogRatio to leading order.
        // Δd/Σγ = 0.001/0.05 = 0.02. Linear gives 0.02; exact gives 2·artanh(0.01) ≈ 0.020001.
        var f = BuildClaim();
        double linear = f.LogRatioLinearApproximation(0.001, 0.05);
        double exact = f.LogRatio(0.001, 0.05);
        Assert.True(Math.Abs(exact - linear) / Math.Abs(linear) < 0.001,
            $"linear ({linear}) and exact ({exact}) should agree to <0.1% at small Δd/Σγ");
    }

    [Theory]
    [InlineData(0.05, 20.0)]
    [InlineData(0.1, 10.0)]
    [InlineData(1.0, 1.0)]
    public void InverseRateLinearCoefficient_IsOneOverSumGamma(double totalGamma, double expected)
    {
        Assert.Equal(expected, BuildClaim().InverseRateLinearCoefficient(totalGamma), precision: 12);
    }

    [Fact]
    public void PalindromicSumHolds_ForF68N3Example()
    {
        // F68 N=3, γ=0.05: d_fast = 0.075, d_slow = 0.025, sum = 0.1 = 2·γ ✓
        Assert.True(BuildClaim().PalindromicSumHolds(dFast: 0.075, dSlow: 0.025, totalGamma: 0.05));
    }

    [Fact]
    public void PalindromicSumHolds_FailsForBrokenPair()
    {
        // Non-palindromic pair: d_fast = 0.10, d_slow = 0.025, sum = 0.125 ≠ 0.1 = 2·0.05
        Assert.False(BuildClaim().PalindromicSumHolds(dFast: 0.10, dSlow: 0.025, totalGamma: 0.05));
    }

    [Fact]
    public void PublicTypedSurface_ContainsNoThermodynamicApi()
    {
        string[] names = typeof(F44AlgebraicPairRateLogIdentityPi2Inheritance)
            .GetMembers().Select(member => member.Name).ToArray();
        Assert.DoesNotContain(names, name => name.Contains("Temperature", StringComparison.OrdinalIgnoreCase));
        Assert.DoesNotContain(names, name => name.Contains("Crooks", StringComparison.OrdinalIgnoreCase));
        Assert.DoesNotContain(names, name => name.Contains("Jarzynski", StringComparison.OrdinalIgnoreCase));
    }

    [Fact]
    public void LogRatio_NonPositiveTotalGamma_Throws()
    {
        foreach (double totalGamma in new[] { -0.05, 0.0, double.NaN, double.PositiveInfinity, double.NegativeInfinity })
            Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().LogRatio(0.05, totalGamma));
    }

    [Fact]
    public void LogRatio_DeltaOutsideFiniteOrderedPairDomain_Throws()
    {
        // Δd = d_fast-d_slow is finite and non-negative; Δd = 2·Σγ is the d_slow=0 singular endpoint.
        foreach (double deltaD in new[] { -0.01, 0.10, 0.20, double.NaN, double.PositiveInfinity, double.NegativeInfinity })
            Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().LogRatio(deltaD, 0.05));
    }

    [Fact]
    public void LinearApproximation_UsesSameFiniteOrderedPairDomainAsExactIdentity()
    {
        var claim = BuildClaim();
        foreach (double deltaD in new[] { -0.01, 0.10, 0.20, double.NaN, double.PositiveInfinity, double.NegativeInfinity })
            Assert.Throws<ArgumentOutOfRangeException>(() => claim.LogRatioLinearApproximation(deltaD, 0.05));
        foreach (double totalGamma in new[] { -0.05, 0.0, double.NaN, double.PositiveInfinity, double.NegativeInfinity })
            Assert.Throws<ArgumentOutOfRangeException>(() => claim.LogRatioLinearApproximation(0.01, totalGamma));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, memoryLoop);
        Assert.Throws<ArgumentNullException>(() =>
            new F44AlgebraicPairRateLogIdentityPi2Inheritance(null!, f1));
    }

    [Fact]
    public void Constructor_NullF1_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F44AlgebraicPairRateLogIdentityPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }
}
