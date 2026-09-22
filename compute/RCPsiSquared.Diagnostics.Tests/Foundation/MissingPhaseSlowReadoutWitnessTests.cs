using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class MissingPhaseSlowReadoutWitnessTests
{
    [Fact]
    public void CompletePreparedResidueHasThePhysicalEndpointReadouts()
    {
        var reading = new MissingPhaseSlowReadoutWitness().Reading;
        Assert.Equal(new BigRational(-1, 2), reading.EndZzResidue);
        Assert.Equal(new BigRational(-1, 4), reading.EndXxResidue);
        Assert.Equal(new BigRational(1, 4), reading.ResidueNormSquared);
        Assert.All(reading.Checks, check => Assert.True(check.Passes, check.Detail));
    }

    [Fact]
    public void RemovingAClusterMemberAndChangingPreparationAreDetected()
    {
        var reading = new MissingPhaseSlowReadoutWitness().Reading;
        Assert.Equal(new BigRational(-1, 4), reading.RankOneEndZzResidue);
        Assert.NotEqual(reading.EndZzResidue, reading.RankOneEndZzResidue);
        Assert.True(reading.StationaryResidueNormSquared.IsZero);
        Assert.Equal(new BigRational(1, 16), reading.SingleCopySinePairNormSquared);
        Assert.True(reading.PairedCopySinePairNormSquared.IsZero);
    }

    [Fact]
    public void PhysicalPairMapAndDecoderKeepDifferentQuadratures()
    {
        var reading = new MissingPhaseSlowReadoutWitness().Reading;
        Assert.Equal(96 * 21, reading.ComplexBasisPairCases);
        Assert.Equal(0, reading.PhysicalPairMismatchCount);
        Assert.True(reading.AllPairSineNormSquared.IsZero);
        Assert.Equal(new BigRational(1, 4), reading.DecodedHalfTraceNormSquared);
        Assert.Equal(3, reading.OddPairZeros);
    }

    [Theory]
    [InlineData(1, 100)]
    [InlineData(1, 7)]
    [InlineData(3, 10)]
    [InlineData(2, 1)]
    [InlineData(100, 1)]
    public void LeakageGermIsRecomputedFromTheBorderedDerivative(int numerator, int denominator)
    {
        var gamma = new BigRational(numerator, denominator);
        var reading = new MissingPhaseSlowReadoutWitness(gamma).Reading;
        Assert.Equal(new BigRational(-1, 16), reading.LeakageRealCoefficient);
        Assert.Equal(gamma / 16, reading.LeakageImaginaryRootTwoCoefficient);
        Assert.All(reading.Checks, check => Assert.True(check.Passes, check.Detail));
    }

    [Fact]
    public void WitnessRejectsNonphysicalGamma()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new MissingPhaseSlowReadoutWitness(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new MissingPhaseSlowReadoutWitness(-1));
    }
}
