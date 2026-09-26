using RCPsiSquared.Core.Calibration;

namespace RCPsiSquared.Core.Tests.Calibration;

/// <summary>Tests the calibration-level normalized-purity proxy and its neutral
/// r = T2/(2*T1) classification. The labels describe only a side of the named
/// R* threshold in the stated free-single-transmon |+> model.</summary>
public class QubitRegimeTests
{
    private const long RStarBits = 0x3fcb3b8c72b39bd9;
    private static Lazy<IReadOnlyList<QubitData>> Marrakesh20260425 =>
        CalibrationFixtures.Marrakesh20260425;

    [Fact]
    public void RStar_PinsBinary64ValueAndDefiningTouch()
    {
        const double yStar = 0.613819895538536;

        Assert.Equal(RStarBits, BitConverter.DoubleToInt64Bits(QubitRegime.R_STAR));
        Assert.InRange(QubitRegime.R_STAR, 0.2127547798220052, 0.2127547798220055);
        Assert.InRange(
            Math.Abs(QubitRegime.NormalizedPurityProxyStationarity(yStar, QubitRegime.R_STAR)),
            0.0,
            1e-13);
        Assert.Equal(
            0.25,
            QubitRegime.NormalizedPurityProxy(yStar, QubitRegime.R_STAR),
            precision: 13);
    }

    [Fact]
    public void RStar_RoundedSixDecimalReplacementFailsTheExactPin()
    {
        Assert.NotEqual(
            BitConverter.DoubleToInt64Bits(0.212755),
            BitConverter.DoubleToInt64Bits(QubitRegime.R_STAR));
    }

    [Fact]
    public void NormalizedPurityProxy_AtPointFourDoesNotTrackActiveCpsiQuarterCrossing()
    {
        const double r = 0.4;
        const double yAtProxyMinimum = 0.5247887155814954;
        const double tauAtActiveCrossing = 0.8503020378519266;

        double proxyMinimum = QubitRegime.NormalizedPurityProxy(yAtProxyMinimum, r);
        double activeCpsi = ActiveCpsiForFreePlus(tauAtActiveCrossing, r);

        Assert.InRange(
            Math.Abs(QubitRegime.NormalizedPurityProxyStationarity(yAtProxyMinimum, r)),
            0.0,
            1e-13);
        Assert.Equal(0.4253341805025406, proxyMinimum, precision: 14);
        Assert.True(proxyMinimum > 0.25);
        Assert.Equal(0.25, activeCpsi, precision: 13);
    }

    [Fact]
    public void RParam_ComputesT2OverTwoT1()
    {
        Assert.Equal(0.5, QubitRegime.RParam(t1Us: 100, t2Us: 100), precision: 6);
        Assert.Equal(1.0, QubitRegime.RParam(t1Us: 100, t2Us: 200), precision: 6);
        Assert.Equal(0.1, QubitRegime.RParam(t1Us: 100, t2Us: 20), precision: 6);
    }

    [Fact]
    public void RParam_NonOperationalT1_IsNotMisreadAsBelowThreshold()
    {
        Assert.Equal(double.PositiveInfinity, QubitRegime.RParam(t1Us: 0, t2Us: 100));
        Assert.Equal(double.PositiveInfinity, QubitRegime.RParam(t1Us: -1, t2Us: 100));
        Assert.Equal(Regime.AtOrAboveRStar, QubitRegime.Classify(t1Us: 0, t2Us: 100));
        Assert.False(QubitRegime.IsBelowRStar(t1Us: 0, t2Us: 100));
    }

    [Fact]
    public void Classify_BinaryDefault_AssignsEqualityToAtOrAbove()
    {
        Assert.Equal(Regime.BelowRStar, QubitRegime.Classify(t1Us: 100, t2Us: 40));
        Assert.Equal(Regime.AtOrAboveRStar, QubitRegime.Classify(t1Us: 100, t2Us: 50));
        Assert.Equal(
            Regime.AtOrAboveRStar,
            QubitRegime.Classify(t1Us: 0.5, t2Us: QubitRegime.R_STAR));
    }

    [Fact]
    public void Classify_FourUlpBandIsClosedAndFifthUlpIsOutside()
    {
        double lowerFour = BitConverter.Int64BitsToDouble(RStarBits - 4);
        double upperFour = BitConverter.Int64BitsToDouble(RStarBits + 4);
        double lowerFive = BitConverter.Int64BitsToDouble(RStarBits - 5);
        double upperFive = BitConverter.Int64BitsToDouble(RStarBits + 5);
        double epsilon = QubitRegime.R_STAR - lowerFour;

        Assert.Equal(Regime.NearRStar, QubitRegime.Classify(0.5, lowerFour, epsilon));
        Assert.Equal(Regime.NearRStar, QubitRegime.Classify(0.5, upperFour, epsilon));
        Assert.Equal(Regime.BelowRStar, QubitRegime.Classify(0.5, lowerFive, epsilon));
        Assert.Equal(Regime.AtOrAboveRStar, QubitRegime.Classify(0.5, upperFive, epsilon));
    }

    [Fact]
    public void MarrakeshRows_Q0IsBelowAndQ1IsAtOrAboveRStar()
    {
        var qubits = Marrakesh20260425.Value;
        var q0 = qubits.Single(q => q.Qubit == 0);
        var q1 = qubits.Single(q => q.Qubit == 1);

        Assert.True(q0.IsBelowRStar);
        Assert.Equal(Regime.BelowRStar, q0.RStarBand);
        Assert.False(q1.IsBelowRStar);
        Assert.Equal(Regime.AtOrAboveRStar, q1.RStarBand);
    }

    [Fact]
    public void PathComposition_UsesThreeNeutralBandCounts()
    {
        var qubits = Marrakesh20260425.Value;

        var mixed = QubitRegime.PathComposition(qubits, new[] { 0, 1, 2 });
        Assert.Equal((1, 0, 2), mixed);

        var oneBand = QubitRegime.PathComposition(qubits, new[] { 48, 49, 50 });
        Assert.Equal((0, 0, 3), oneBand);

        var apr25Best5 = QubitRegime.PathComposition(qubits, new[] { 1, 2, 3, 4, 5 });
        Assert.Equal((0, 0, 5), apr25Best5);
    }

    [Fact]
    public void Classify_NegativeEpsilon_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => QubitRegime.Classify(100, 30, epsilon: -1e-15));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            QubitRegime.PathComposition(Marrakesh20260425.Value, new[] { 0 }, epsilon: -1.0));
    }

    [Fact]
    public void PathComposition_RejectsUnknownQubit()
    {
        Assert.Throws<ArgumentException>(() =>
            QubitRegime.PathComposition(Marrakesh20260425.Value, new[] { 0, 9999 }));
    }

    [Fact]
    public void QubitDataDerivedValues_AgreeWithStaticClassifier()
    {
        foreach (var q in Marrakesh20260425.Value)
        {
            Assert.Equal(QubitRegime.Classify(q.T1Us, q.T2Us), q.RStarBand);
            Assert.Equal(QubitRegime.RParam(q.T1Us, q.T2Us), q.RParam);
            Assert.Equal(QubitRegime.IsBelowRStar(q.T1Us, q.T2Us), q.IsBelowRStar);
        }
    }

    private static double ActiveCpsiForFreePlus(double tau, double r)
    {
        double y = Math.Exp(-tau);
        double purity = 1.0 - y + 0.5 * y * y + 0.5 * Math.Exp(-tau / r);
        return purity * Math.Exp(-tau / (2.0 * r));
    }
}
