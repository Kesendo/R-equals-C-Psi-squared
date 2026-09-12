using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class MissingPhaseRelaxationScaleWitnessTests
{
    [Theory]
    [InlineData(0.0)]
    [InlineData(0.11)]
    [InlineData(-0.11)]
    [InlineData(1e-20)]
    [InlineData(-1e-20)]
    [InlineData(1e-15)]
    [InlineData(-1e-15)]
    [InlineData(double.NaN)]
    [InlineData(double.PositiveInfinity)]
    [InlineData(double.NegativeInfinity)]
    public void Constructor_RejectsNonlocalOrNonfiniteEpsilon(double epsilon)
    {
        Assert.Equal("epsilon", Assert.Throws<ArgumentOutOfRangeException>(
            () => new MissingPhaseRelaxationScaleWitness(epsilon, 0.3)).ParamName);
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(-0.3)]
    [InlineData(double.NaN)]
    [InlineData(double.PositiveInfinity)]
    [InlineData(double.NegativeInfinity)]
    public void Constructor_RejectsNonpositiveOrNonfiniteGamma(double gamma)
    {
        Assert.Equal("gamma", Assert.Throws<ArgumentOutOfRangeException>(
            () => new MissingPhaseRelaxationScaleWitness(0.01, gamma)).ParamName);
    }

    [Fact]
    public void Constructor_RejectsPositiveGammaWhenThePredictedRateIsBelowTheSpectralErrorModel()
    {
        Assert.Equal("gamma", Assert.Throws<ArgumentOutOfRangeException>(
            () => new MissingPhaseRelaxationScaleWitness(0.01, 1e-20)).ParamName);
    }

    [Fact]
    public void DefaultWitness_IsTheFixedN7LocalA11Companion()
    {
        var witness = new MissingPhaseRelaxationScaleWitness();

        Assert.Equal(7, MissingPhaseRelaxationScaleWitness.SiteCount);
        Assert.Equal(3, MissingPhaseRelaxationScaleWitness.WatchedSeat);
        Assert.Equal(0.01, witness.Epsilon);
        Assert.Equal(0.3, witness.Gamma);
        Assert.Contains("N = 7", witness.Scope);
        Assert.Contains("A/(1,1)", witness.Scope);
        Assert.Contains("fixed gamma > 0", witness.Scope);
        Assert.Contains("punctured", witness.Scope);
        Assert.Contains("not a punctured-neighbourhood radius proved uniformly in gamma", witness.Scope);
        Assert.Contains("not the short-time onset", witness.Scope);
        Assert.Contains("not an all-N theorem", witness.Scope);
        Assert.Contains("not the full 4^7 Liouvillian gap", witness.Scope);
        Assert.Contains("not an observable lifetime", witness.Scope);
        Assert.Contains("independently", witness.ConstructionProvenance);
        Assert.Contains("Evd(A)", witness.ConstructionProvenance);
        Assert.Contains("Evd(A†)", witness.ConstructionProvenance);
        Assert.Contains("row-major", witness.ConstructionProvenance);
        Assert.Contains("I_E", witness.FullSpaceOnlyFence);
        Assert.Contains("3 gamma/2", witness.FullSpaceOnlyFence);
    }

    [Fact]
    public void SpectralRoute_ResolvesBothConjugateRankTwoClustersOnBothSides()
    {
        var witness = new MissingPhaseRelaxationScaleWitness();

        Assert.Equal(new[] { -1, +1 }, witness.Clusters.Select(cluster => cluster.FrequencySign).ToArray());
        foreach (var cluster in witness.Clusters)
        {
            Assert.Equal(2, cluster.Right.Rank); // a rank-one mutation turns this red
            Assert.Equal(2, cluster.Adjoint.Rank);
            Assert.Equal(2.0, cluster.Right.ProjectorTrace, 11);
            Assert.Equal(2.0, cluster.Adjoint.ProjectorTrace, 11);
            Assert.True(cluster.Right.SeedOverlap > 1.9);
            Assert.True(cluster.Adjoint.SeedOverlap > 1.9);
            Assert.True(cluster.IndependentProjectorDistance > witness.ProjectorErrorModel,
                "copying the right projector into the A† route must be observable");
        }
    }

    [Fact]
    public void SpectralRoute_ProjectorsAndAbsorptionCloseUnderTheirNumericalErrorLaw()
    {
        var witness = new MissingPhaseRelaxationScaleWitness();

        Assert.Equal(MissingPhaseRelaxationScaleWitness.MachineEpsilon
                     * MissingPhaseRelaxationScaleWitness.GeneratorDimension
                     * MissingPhaseRelaxationScaleWitness.GeneratorDimension,
            witness.ProjectorErrorModel);
        Assert.Equal(MissingPhaseRelaxationScaleWitness.MachineEpsilon
                     * Math.Max(1.0, witness.GeneratorFrobeniusNorm)
                     * MissingPhaseRelaxationScaleWitness.GeneratorDimension
                     * MissingPhaseRelaxationScaleWitness.GeneratorDimension,
            witness.SpectralErrorModel);

        foreach (var side in witness.Clusters.SelectMany(cluster => new[] { cluster.Right, cluster.Adjoint }))
        {
            Assert.Equal(2, side.Rank);
            Assert.True(side.HermitianResidual <= witness.ProjectorErrorModel,
                $"{side.Source}: Hermitian residual/model = {side.HermitianResidual / witness.ProjectorErrorModel:E3}");
            Assert.True(side.IdempotenceResidual <= witness.ProjectorErrorModel,
                $"{side.Source}: idempotence residual/model = {side.IdempotenceResidual / witness.ProjectorErrorModel:E3}");
            Assert.True(side.RitzResidual <= witness.SpectralErrorModel,
                $"{side.Source}: Ritz residual/model = {side.RitzResidual / witness.SpectralErrorModel:E3}");
            Assert.True(side.AbsorptionResidual <= witness.SpectralErrorModel,
                $"{side.Source}: absorption residual/model = {side.AbsorptionResidual / witness.SpectralErrorModel:E3}");
            Assert.True(side.MeanDecay > 0.0);
            Assert.True(side.CentreLight > 0.0);
        }
    }

    [Fact]
    public void SpectralRoute_LargeGammaNegativeEpsilonKeepsEachCompleteRankTwoClusterTogether()
    {
        // Regression: selecting the two best individual vectors mixed two nearby doubled
        // eigenvalues here, while the conjugate branch happened to remain intact.
        var witness = new MissingPhaseRelaxationScaleWitness(-0.1, 100.0);

        foreach (var cluster in witness.Clusters)
        {
            Assert.Equal(2, cluster.Right.Rank);
            Assert.Equal(2, cluster.Adjoint.Rank);
            Assert.True(Math.Abs(cluster.Right.MeanDecay - cluster.Adjoint.MeanDecay)
                        <= witness.SpectralErrorModel,
                $"sign {cluster.FrequencySign}: right={cluster.Right.MeanDecay:R}, " +
                $"adjoint={cluster.Adjoint.MeanDecay:R}, model={witness.SpectralErrorModel:R}");
        }
        Assert.True(Math.Abs(witness.Clusters[0].Right.MeanDecay - witness.Clusters[1].Right.MeanDecay)
                    <= witness.SpectralErrorModel,
            "the conjugate Evd(A) clusters must have the same decay; a mixed pair breaks this");
        Assert.True(Math.Abs(witness.Clusters[0].Adjoint.MeanDecay - witness.Clusters[1].Adjoint.MeanDecay)
                    <= witness.SpectralErrorModel,
            "the conjugate Evd(A†) clusters must have the same decay; a mixed pair breaks this");
        const double independentlyContinuedDecay = 0.005523959421952;
        foreach (var side in witness.Clusters.SelectMany(cluster => new[] { cluster.Right, cluster.Adjoint }))
            Assert.True(Math.Abs(side.MeanDecay - independentlyContinuedDecay)
                        <= witness.SpectralErrorModel,
                $"{side.Source}: decay={side.MeanDecay:R}, continued anchor={independentlyContinuedDecay:R}");
    }

    [Fact]
    public void SpectralRoute_FailsClosedWhenStrongDephasingMakesTheSelectedRateUnresolved()
    {
        var witness = new MissingPhaseRelaxationScaleWitness(-0.1, 1e8);

        var error = Assert.Throws<InvalidOperationException>(() => witness.Clusters);
        Assert.Contains("unresolved", error.Message);
        Assert.Contains("spectral error model", error.Message);
    }

    [Fact]
    public void DirectGenerator_WatchesSeatThreeAndTheMovedSeatControlChangesTheSameCell()
    {
        const double epsilon = 0.01;
        const double gamma = 0.3;
        int chargedCell = MissingPhaseRelaxationScaleWitness.RowMajorIndex(3, 0);

        Matrix<Complex> centre = MissingPhaseRelaxationScaleWitness.BuildAGenerator(epsilon, gamma, 3);
        Matrix<Complex> moved = MissingPhaseRelaxationScaleWitness.BuildAGenerator(epsilon, gamma, 2);

        Assert.Equal(-2.0 * gamma, centre[chargedCell, chargedCell].Real, 13);
        Assert.Equal(0.0, moved[chargedCell, chargedCell].Real, 13);
        Assert.True((centre - moved).FrobeniusNorm() > 10.0 * double.Epsilon);

        int x00 = MissingPhaseRelaxationScaleWitness.RowMajorIndex(0, 0);
        int x10 = MissingPhaseRelaxationScaleWitness.RowMajorIndex(1, 0);
        int x01 = MissingPhaseRelaxationScaleWitness.RowMajorIndex(0, 1);
        double endHopping = 2.0 * (1.0 + epsilon);
        Assert.Equal(-Complex.ImaginaryOne * endHopping, centre[x00, x10]);
        Assert.Equal(+Complex.ImaginaryOne * endHopping, centre[x00, x01]);
    }

    [Fact]
    public void HilbertRoute_ContinuesThreeBlindModesAndReadsTheCubicLightGerm()
    {
        var witness = new MissingPhaseRelaxationScaleWitness();
        var modes = witness.HilbertModes.ToDictionary(mode => mode.Label);

        Assert.Equal(new[] { "d-", "d0", "d+" }, witness.HilbertModes.Select(mode => mode.Label).ToArray());
        Assert.True(modes["d0"].CentreWeight <= witness.HilbertErrorModel);
        Assert.True(modes["d-"].SeedOverlap > 0.99);
        Assert.True(modes["d0"].SeedOverlap > 0.99);
        Assert.True(modes["d+"].SeedOverlap > 0.99);
        Assert.Equal(modes["d-"].CentreWeight, modes["d+"].CentreWeight, 12);

        double predicted = Math.Pow(witness.Epsilon, 2) / 4.0
                         + Math.Pow(witness.Epsilon, 3) / 4.0;
        double fourthOrderScale = Math.Pow(Math.Abs(witness.Epsilon), 4);
        Assert.True(Math.Abs(modes["d-"].CentreWeight - predicted) <= fourthOrderScale);
        Assert.True(Math.Abs(modes["d+"].CentreWeight - predicted) <= fourthOrderScale);
    }

    [Fact]
    public void HilbertDyads_ReproduceTheFiveNonIESecondOrderCoefficients()
    {
        var witness = new MissingPhaseRelaxationScaleWitness();

        Assert.Equal(5, witness.DyadCoefficients.Count);
        Assert.Equal(new[]
        {
            "-2 sqrt(2)i: |d0><d-|",
            "+2 sqrt(2)i: |d-><d0|",
            "-4 sqrt(2)i: |d+><d-|",
            "+4 sqrt(2)i: |d-><d+|",
            "0: |d+><d+| - |d-><d-|",
        }, witness.DyadCoefficients.Select(reading => reading.Label).ToArray());
        Assert.Equal(new[] { 0.5, 0.5, 1.0, 1.0, 1.0 },
            witness.DyadCoefficients.Select(reading => reading.ExpectedCoefficientOverGamma).ToArray());
        foreach (var reading in witness.DyadCoefficients)
        {
            Assert.False(reading.IncludesOuterComplementIdentity);
            Assert.True(Math.Abs(reading.MeasuredCoefficientOverGamma - reading.CubicPredictionOverGamma)
                        <= reading.FourthOrderScaleOverGamma,
                $"{reading.Label}: measured/gamma {reading.MeasuredCoefficientOverGamma:R}, " +
                $"cubic/gamma {reading.CubicPredictionOverGamma:R}, " +
                $"scale/gamma {reading.FourthOrderScaleOverGamma:R}");
            Assert.Equal(witness.Gamma * reading.MeasuredCoefficientOverGamma,
                reading.MeasuredRateCoefficient, 14);
            Assert.Equal(witness.Gamma * reading.ExpectedCoefficientOverGamma,
                reading.ExpectedRateCoefficient, 14);
        }
        Assert.DoesNotContain(witness.DyadCoefficients,
            reading => reading.ExpectedCoefficientOverGamma == 1.5);

        string summary = witness.Children.Single(child => child.DisplayName.Contains("five non-I_E")).Summary;
        Assert.Contains("gamma-normalized", summary);
        Assert.Contains("rate coefficient", summary);
    }

    [Fact]
    public void InspectableTree_CarriesBothIndependentSpectralSidesAndTheHilbertControl()
    {
        var witness = new MissingPhaseRelaxationScaleWitness();
        var labels = witness.Children.Select(child => child.DisplayName).ToList();

        Assert.Contains(labels, label => label.Contains("Evd(A)"));
        Assert.Contains(labels, label => label.Contains("Evd(A†)"));
        Assert.Contains(labels, label => label.Contains("Hilbert"));
        Assert.Contains(labels, label => label.Contains("five non-I_E"));
        Assert.Contains(labels, label => label.Contains("scope and non-claims"));
    }
}
