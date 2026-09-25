using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class MissingPhaseSlowReadoutWitnessTests
{
    private static readonly Lazy<MissingPhaseSlowReadoutWitness.Snapshot> Canonical =
        new(() => new MissingPhaseSlowReadoutWitness().Reading);

    [Fact]
    public void CompletePreparedResidueHasThePhysicalEndpointReadouts()
    {
        var reading = Canonical.Value;
        Assert.Equal(new BigRational(-1, 2), reading.EndZzResidue);
        Assert.Equal(new BigRational(-1, 4), reading.EndXxResidue);
        Assert.Equal(new BigRational(1, 4), reading.ResidueNormSquared);
        Assert.All(reading.Checks, check => Assert.True(check.Passes, check.Detail));
    }

    [Fact]
    public void TheTwoDyadsExhaustTheClusterWithoutAJordanPartner()
    {
        var reading = Canonical.Value;
        Assert.Equal(47, reading.ClusterRank);
        Assert.Equal(47, reading.ClusterSquareRank);
    }

    [Fact]
    public void RemovingAClusterMemberAndChangingPreparationAreDetected()
    {
        var reading = Canonical.Value;
        Assert.Equal(new BigRational(-1, 4), reading.RankOneEndZzResidue);
        Assert.NotEqual(reading.EndZzResidue, reading.RankOneEndZzResidue);
        Assert.True(reading.StationaryResidueNormSquared.IsZero);
        Assert.Equal(new BigRational(1, 16), reading.SingleCopySinePairNormSquared);
        Assert.True(reading.PairedCopySinePairNormSquared.IsZero);
    }

    [Fact]
    public void PhysicalPairMapAndDecoderKeepDifferentQuadratures()
    {
        var reading = Canonical.Value;
        Assert.Equal(96 * 21, reading.ComplexBasisPairCases);
        Assert.Equal(0, reading.PhysicalPairMismatchCount);
        Assert.Equal(168, reading.RealPartMutationMismatchCount);
        Assert.True(reading.AllPairSineNormSquared.IsZero);
        Assert.Equal(new BigRational(1, 4), reading.DecodedHalfTraceNormSquared);
    }

    [Theory]
    [InlineData(1, 100)]
    [InlineData(3, 10)]
    [InlineData(100, 1)]
    public void SevenPairsMissTheUniformResidueAndFourOfThemSwitchOnAtFirstOrder(int numerator, int denominator)
    {
        var gamma = new BigRational(numerator, denominator);
        var reading = new MissingPhaseSlowReadoutWitness(gamma).Reading;
        Assert.Equal(MissingPhaseSlowReadoutWitness.UniformNullPairs, reading.NullPairsAtUniformPoint);
        Assert.Equal(MissingPhaseSlowReadoutWitness.OddNullPairs, reading.NullPairsAtFirstOrder);
        Assert.Equal(new BigRational(11, 128), reading.FirstOrderPairNormSquared[(0, 2)]);
        Assert.Equal(new BigRational(3, 128), reading.FirstOrderPairNormSquared[(4, 6)]);
        var evenCross = (3 + 6 * gamma * gamma) / 128;
        Assert.Equal(evenCross, reading.FirstOrderPairNormSquared[(0, 4)]);
        Assert.Equal(evenCross, reading.FirstOrderPairNormSquared[(2, 6)]);
    }

    [Fact]
    public void ThreeSitesReadTheSineQuadratureThatEveryPairLoses()
    {
        var reading = Canonical.Value;
        Assert.Equal(new BigRational(1, 8), reading.ThreeSpinResidueImaginaryRootTwoCoefficient);
        Assert.Equal(new BigRational(1, 4), reading.ThreeSpinSineRootTwoCoefficient);
        Assert.Equal(24, reading.SeparatingThreeSitePages);
    }

    [Theory]
    [InlineData(1, 1)]
    [InlineData(1, 7)]
    public void ACentreOnsiteEnergyBreaksThePairingGatesButNotTheOddTier(int numerator, int denominator)
    {
        // Breaks C h C = -h and nothing the uniform dyads see. The two gates that rest on the chiral
        // pairing must go red; the odd-pair null, which rests on the blind ray's support, must not.
        var energy = new BigRational(numerator, denominator);
        var mutant = new MissingPhaseSlowReadoutWitness(new BigRational(3, 10), energy).Reading;
        MissingPhaseSlowReadoutWitness.Check Named(string name) => mutant.Checks.Single(check => check.Name == name);
        Assert.True(Named("first dyad exact A eigenoperator").Passes);
        Assert.True(Named("second dyad exact A eigenoperator").Passes);
        Assert.True(Named("physical endpoint ZZ residue").Passes);
        Assert.False(Named("chiral pairing at first order: e-' = C conj(e+') from two independent solves").Passes);
        Assert.False(Named("first-order mixed even/odd cells are antisymmetric").Passes);
        Assert.True(Named("of the seven uniform null pairs, those still null at first order").Passes);
    }

    private static GaussianRational[,] DefectHopping() =>
        MissingPhaseOnsetWitness.Hopping(MissingPhaseSlowReadoutWitness.SiteCount, new BigRational(1, 3));

    private static readonly BigRational DefectRatio = new(4, 3);

    [Fact]
    public void TheBlindRayGatesPassOnThePhysicalChain()
    {
        Assert.All(MissingPhaseSlowReadoutWitness.BlindRayChecks(DefectHopping(), DefectRatio),
            check => Assert.True(check.Passes, check.Detail));
    }

    [Fact]
    public void ABondBetweenTwoEvenSitesRemovesTheBlindRay()
    {
        var hopping = DefectHopping();
        hopping[0, 2] = new BigRational(1);
        hopping[2, 0] = new BigRational(1);
        var checks = MissingPhaseSlowReadoutWitness.BlindRayChecks(hopping, DefectRatio);
        Assert.False(checks[0].Passes);
        Assert.Equal("0", checks[0].Actual);
    }

    [Fact]
    public void TwoTunedEvenBondsKeepOneRayWithOddEntriesAndTheOddSiteGateRejectsIt()
    {
        // Bonds (0,2) = 1 and (4,6) = -3/4 at r = 4/3 leave a one-dimensional kernel,
        // (-3/4, -3/8, 1, 3/4, -1, -3/8, 1): the dimension gate stays green and the odd-site gate fires.
        var hopping = DefectHopping();
        hopping[0, 2] = new BigRational(1);
        hopping[2, 0] = new BigRational(1);
        hopping[4, 6] = new BigRational(-3, 4);
        hopping[6, 4] = new BigRational(-3, 4);
        var checks = MissingPhaseSlowReadoutWitness.BlindRayChecks(hopping, DefectRatio);
        Assert.True(checks[0].Passes, checks[0].Detail);
        Assert.False(checks[1].Passes);
        Assert.Equal("3", checks[1].Actual);
    }

    // Five distinct gamma make every first-order statement a statement at every gamma: both bordered
    // determinants are gamma-free, so each derivative, each first-order pair entry and the leakage germ is
    // affine in gamma, and two agreeing values already fix an affine function.
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
