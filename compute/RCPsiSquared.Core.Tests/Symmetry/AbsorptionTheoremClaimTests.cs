using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class AbsorptionTheoremClaimTests
{
    private static AbsorptionTheoremClaim BuildClaim() =>
        new AbsorptionTheoremClaim(new Pi2DyadicLadderClaim());

    [Fact]
    public void PublicSurface_SeparatesBasisPairCostFromEigenmodeExpectation()
    {
        var methods = typeof(AbsorptionTheoremClaim).GetMethods()
            .Select(m => m.Name).ToHashSet();
        Assert.Contains("BasisPairDissipatorCost", methods);
        Assert.Contains("EigenmodeDecayRate", methods);
        Assert.Contains("AverageNXyFromEigenmodeDecayRate", methods);
        Assert.DoesNotContain("Rate", methods);
        Assert.DoesNotContain("NXyFromRate", methods);

        var surface = $"{BuildClaim().Name} {BuildClaim().DisplayName} {BuildClaim().Summary}";
        Assert.DoesNotContain("spectrum quantized", surface, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("rate-quantization", surface, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("basis-pair", surface, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("eigenmode", surface, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void DissipatorCoefficient_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().DissipatorCoefficient, precision: 14);
    }

    [Fact]
    public void DissipatorCoefficientMatchesLiteral_HoldsExactly()
    {
        Assert.True(BuildClaim().DissipatorCoefficientMatchesLiteral());
    }

    [Theory]
    [InlineData(0.0, 0.0)]
    [InlineData(0.05, 0.1)]
    [InlineData(0.5, 1.0)]
    [InlineData(1.0, 2.0)]
    [InlineData(2.5, 5.0)]
    public void SingleDisagreementCellCost_IsTwoTimesGammaZero(double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().SingleDisagreementCellCost(gammaZero), precision: 12);
    }

    [Theory]
    [InlineData(0, 1.0, 0.0)]
    [InlineData(1, 1.0, 2.0)]
    [InlineData(3, 1.0, 6.0)]
    [InlineData(5, 1.0, 10.0)]
    [InlineData(3, 0.05, 0.30)]
    public void EigenmodeDecayRate_IsTwoGammaTimesAverageNXY(int nXY, double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().EigenmodeDecayRate(nXY, gammaZero), precision: 12);
    }

    [Fact]
    public void EigenmodeDecayRate_AcceptsNonIntegerLightExpectation()
    {
        Assert.Equal(8.0 / 3.0,
            BuildClaim().EigenmodeDecayRate(averageNXy: 4.0 / 3.0, gammaZero: 1.0),
            precision: 12);
    }

    [Theory]
    [InlineData(0, 1.0, 0.0)]
    [InlineData(1, 1.0, 2.0)]
    [InlineData(3, 1.0, 6.0)]
    public void BasisPairDissipatorCost_EqualsEigenmodeFormulaAtSameNumericLight(int nDiff, double gammaZero, double expected)
    {
        var c = BuildClaim();
        Assert.Equal(expected, c.BasisPairDissipatorCost(nDiff, gammaZero), precision: 12);
        Assert.Equal(c.EigenmodeDecayRate(nDiff, gammaZero), c.BasisPairDissipatorCost(nDiff, gammaZero), precision: 14);
    }

    [Theory]
    [InlineData(1, 1.0, 2.0)]
    [InlineData(3, 1.0, 6.0)]
    [InlineData(5, 1.0, 10.0)]
    [InlineData(7, 0.05, 0.7)]
    public void EigenmodeDecayRateCeiling_IsTwoGammaTimesN(int n, double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().EigenmodeDecayRateCeiling(n, gammaZero), precision: 12);
    }

    [Theory]
    [InlineData(0.0, 1.0, 0.0)]
    [InlineData(2.0, 1.0, 1.0)]
    [InlineData(6.0, 1.0, 3.0)]
    [InlineData(0.10, 0.05, 1.0)]
    [InlineData(0.30, 0.05, 3.0)]
    public void AverageNXyFromEigenmodeDecayRate_InvertsRate(double rate, double gammaZero, double expectedNXy)
    {
        Assert.Equal(expectedNXy, BuildClaim().AverageNXyFromEigenmodeDecayRate(rate, gammaZero), precision: 12);
    }

    [Theory]
    [InlineData(2, 1.0, 4.0)]
    [InlineData(3, 1.0, 6.0)]   // F89c path-2 anchor: pair-sum = 6γ for 3-qubit block
    [InlineData(4, 1.0, 8.0)]   // F89c path-3 prediction: pair-sum = 8γ for 4-qubit block
    [InlineData(3, 0.05, 0.30)]
    public void HammingComplementCellCostSum_IsTwoGammaTimesBlockSize(int blockSize, double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().HammingComplementCellCostSum(blockSize, gammaZero), precision: 12);
    }

    [Fact]
    public void HammingComplementCellCostSum_AtN3_MatchesF89cEmpiricalAnchor()
    {
        // F89c bit-exact verification: (SE,SE) + (SE,DE) eigenvalue pairs sum to 6γ at path-2.
        Assert.Equal(6.0, BuildClaim().HammingComplementCellCostSum(3, 1.0), precision: 12);
    }

    [Fact]
    public void EigenmodeDecayRate_NegativeAverageNXY_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().EigenmodeDecayRate(-1, 1.0));
    }

    [Fact]
    public void EigenmodeDecayRate_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().EigenmodeDecayRate(1, -0.1));
    }

    [Fact]
    public void AverageNXyFromEigenmodeDecayRate_NegativeRate_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().AverageNXyFromEigenmodeDecayRate(-1.0, 1.0));
    }

    [Fact]
    public void AverageNXyFromEigenmodeDecayRate_ZeroGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().AverageNXyFromEigenmodeDecayRate(1.0, 0.0));
    }

    [Fact]
    public void EigenmodeDecayRateCeiling_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().EigenmodeDecayRateCeiling(0, 1.0));
    }

    [Fact]
    public void HammingComplementCellCostSum_BlockSizeZero_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().HammingComplementCellCostSum(0, 1.0));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => new AbsorptionTheoremClaim(null!));
    }
}
