using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class AbsorptionTheoremClaimTests
{
    private static AbsorptionTheoremClaim BuildClaim() =>
        new AbsorptionTheoremClaim(new Pi2DyadicLadderClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void AbsorptionQuantumCoefficient_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().AbsorptionQuantumCoefficient, precision: 14);
    }

    [Fact]
    public void AbsorptionQuantumMatchesLiteral_HoldsExactly()
    {
        Assert.True(BuildClaim().AbsorptionQuantumMatchesLiteral());
    }

    [Theory]
    [InlineData(0.0, 0.0)]
    [InlineData(0.05, 0.1)]
    [InlineData(0.5, 1.0)]
    [InlineData(1.0, 2.0)]
    [InlineData(2.5, 5.0)]
    public void AbsorptionQuantum_IsTwoTimesGammaZero(double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().AbsorptionQuantum(gammaZero), precision: 12);
    }

    [Theory]
    [InlineData(0, 1.0, 0.0)]
    [InlineData(1, 1.0, 2.0)]
    [InlineData(3, 1.0, 6.0)]
    [InlineData(5, 1.0, 10.0)]
    [InlineData(3, 0.05, 0.30)]
    public void Rate_IsTwoGammaTimesNXY(int nXY, double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().Rate(nXY, gammaZero), precision: 12);
    }

    [Theory]
    [InlineData(0, 1.0, 0.0)]
    [InlineData(1, 1.0, 2.0)]
    [InlineData(3, 1.0, 6.0)]
    public void PerCoherenceRateComputationalBasis_EqualsRateAtNXyEqualsNDiff(int nDiff, double gammaZero, double expected)
    {
        var c = BuildClaim();
        Assert.Equal(expected, c.PerCoherenceRateComputationalBasis(nDiff, gammaZero), precision: 12);
        Assert.Equal(c.Rate(nDiff, gammaZero), c.PerCoherenceRateComputationalBasis(nDiff, gammaZero), precision: 14);
    }

    [Theory]
    [InlineData(1, 1.0, 2.0)]
    [InlineData(3, 1.0, 6.0)]
    [InlineData(5, 1.0, 10.0)]
    [InlineData(7, 0.05, 0.7)]
    public void MaxRate_IsTwoGammaTimesN(int n, double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().MaxRate(n, gammaZero), precision: 12);
    }

    [Theory]
    [InlineData(0.0, 1.0, 0.0)]
    [InlineData(2.0, 1.0, 1.0)]
    [InlineData(6.0, 1.0, 3.0)]
    [InlineData(0.10, 0.05, 1.0)]
    [InlineData(0.30, 0.05, 3.0)]
    public void NXyFromRate_InvertsRate(double rate, double gammaZero, double expectedNXy)
    {
        Assert.Equal(expectedNXy, BuildClaim().NXyFromRate(rate, gammaZero), precision: 12);
    }

    [Theory]
    [InlineData(2, 1.0, 4.0)]
    [InlineData(3, 1.0, 6.0)]   // F89c path-2 anchor: pair-sum = 6γ for 3-qubit block
    [InlineData(4, 1.0, 8.0)]   // F89c path-3 prediction: pair-sum = 8γ for 4-qubit block
    [InlineData(3, 0.05, 0.30)]
    public void HammingComplementPairSum_IsTwoGammaTimesBlockSize(int blockSize, double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().HammingComplementPairSum(blockSize, gammaZero), precision: 12);
    }

    [Fact]
    public void HammingComplementPairSum_AtN3_MatchesF89cEmpiricalAnchor()
    {
        // F89c bit-exact verification: (SE,SE) + (SE,DE) eigenvalue pairs sum to 6γ at path-2.
        Assert.Equal(6.0, BuildClaim().HammingComplementPairSum(3, 1.0), precision: 12);
    }

    [Fact]
    public void Rate_NegativeNXY_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().Rate(-1, 1.0));
    }

    [Fact]
    public void Rate_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().Rate(1, -0.1));
    }

    [Fact]
    public void NXyFromRate_NegativeRate_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().NXyFromRate(-1.0, 1.0));
    }

    [Fact]
    public void NXyFromRate_ZeroGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().NXyFromRate(1.0, 0.0));
    }

    [Fact]
    public void MaxRate_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().MaxRate(0, 1.0));
    }

    [Fact]
    public void HammingComplementPairSum_BlockSizeZero_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().HammingComplementPairSum(0, 1.0));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => new AbsorptionTheoremClaim(null!));
    }
}
