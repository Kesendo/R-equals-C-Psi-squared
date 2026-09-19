using RCPsiSquared.Core.Calibration;

namespace RCPsiSquared.Core.Tests.Calibration;

/// <summary>Tests the addressability, empirical score, and neutral R* band
/// composition as separate readings of a calibration snapshot.</summary>
public class RegimeSummaryTests
{
    private static Lazy<IReadOnlyList<QubitData>> Marrakesh20260425 =>
        CalibrationFixtures.Marrakesh20260425;

    [Fact]
    public void FrameworkSnapshotsPath_IsBandMixedAndAddressable()
    {
        var s = RegimeSummary.For(Marrakesh20260425.Value, new[] { 0, 1, 2 });

        Assert.Equal(RegimeVerdict.RStarBandMixed, s.RStarBandVerdict);
        Assert.False(s.HasSingleNonBoundaryRStarBand);
        Assert.True(s.IsAddressable);
        Assert.False(s.IsAddressableAndSingleRStarBand);
        Assert.Equal(1, s.BelowRStarCount);
        Assert.Equal(0, s.NearRStarCount);
        Assert.Equal(2, s.AtOrAboveRStarCount);
        Assert.Equal(597.27, s.Score, precision: 1);
    }

    [Fact]
    public void SoftBreakPath_IsAllAtOrAbove()
    {
        var s = RegimeSummary.For(Marrakesh20260425.Value, new[] { 48, 49, 50 });

        Assert.Equal(RegimeVerdict.AllAtOrAboveRStar, s.RStarBandVerdict);
        Assert.True(s.HasSingleNonBoundaryRStarBand);
        Assert.True(s.IsAddressableAndSingleRStarBand);
        Assert.Equal(3, s.AtOrAboveRStarCount);
        Assert.Equal(682.50, s.Score, precision: 1);
    }

    [Fact]
    public void PerQubitRow_CarriesRawCalibrationAndDerivedBand()
    {
        var s = RegimeSummary.For(Marrakesh20260425.Value, new[] { 0, 1, 2 });
        var q0 = s.Qubits[0];
        var q1 = s.Qubits[1];

        Assert.Equal(0, q0.Qubit);
        Assert.Equal(Regime.BelowRStar, q0.RStarBand);
        Assert.True(q0.RParam < QubitRegime.R_STAR);
        Assert.Equal(1, q1.Qubit);
        Assert.Equal(Regime.AtOrAboveRStar, q1.RStarBand);
        Assert.True(q1.RParam > QubitRegime.R_STAR);
    }

    [Fact]
    public void ExplicitNearBand_PreventsSingleNonBoundaryBandVerdict()
    {
        var one = new QubitData(
            Qubit: 7,
            T1Us: 0.5,
            T2Us: QubitRegime.R_STAR,
            ReadoutError: 0,
            SxError: 0,
            PauliXError: 0,
            Operational: true,
            CzNeighbours: new Dictionary<int, double>(),
            RzzNeighbours: new Dictionary<int, double>());

        var s = RegimeSummary.For(new[] { one }, new[] { 7 }, epsilon: 1e-15);

        Assert.Equal(1, s.NearRStarCount);
        Assert.False(s.HasSingleNonBoundaryRStarBand);
        Assert.Equal(RegimeVerdict.RStarBandMixed, s.RStarBandVerdict);
    }

    [Fact]
    public void NonCzCoupledPath_IsNotAddressableAndVerdictTakesPrecedence()
    {
        var s = RegimeSummary.For(Marrakesh20260425.Value, new[] { 0, 50 });

        Assert.False(s.AllCzCoupled);
        Assert.True(s.AllOperational);
        Assert.False(s.IsAddressable);
        Assert.False(s.IsAddressableAndSingleRStarBand);
        Assert.Equal(RegimeVerdict.NotAddressable, s.RStarBandVerdict);
    }

    [Fact]
    public void InvalidPaths_AreRejected()
    {
        Assert.Throws<ArgumentException>(() =>
            RegimeSummary.For(Marrakesh20260425.Value, new[] { 0, 9999 }));
        Assert.Throws<ArgumentException>(() =>
            RegimeSummary.For(Marrakesh20260425.Value, Array.Empty<int>()));
    }

    [Fact]
    public void Headline_UsesNeutralLabelsAndExplicitCounts()
    {
        string h = RegimeSummary.For(Marrakesh20260425.Value, new[] { 0, 1, 2 }).ToHeadline();

        Assert.Contains("[0, 1, 2]", h);
        Assert.Contains("r-star-band-mixed", h);
        Assert.Contains("1 below, 0 near, 2 at-or-above", h);
        Assert.Contains("addressable", h);
        Assert.Contains("empirical score 597", h);
    }

    [Fact]
    public void BestChain_ComposesThroughSameNeutralSurface()
    {
        var qubits = Marrakesh20260425.Value;
        var calChain = IbmCalibration.SelectBestChain(qubits, length: 5);
        var s = RegimeSummary.For(qubits, calChain.QubitIds);

        Assert.True(s.IsAddressable);
        Assert.Equal(RegimeVerdict.AllAtOrAboveRStar, s.RStarBandVerdict);
        Assert.Equal(calChain.Score, s.Score, precision: 1);
    }

    [Fact]
    public void RelatedConfirmations_RemainsAReadOnlyLookup()
    {
        var s = RegimeSummary.For(Marrakesh20260425.Value, new[] { 0, 1, 2 });
        var hits = s.RelatedConfirmations().Select(c => c.Name).ToList();

        Assert.Contains("pi_protected_xiz_yzzy", hits);
        Assert.Contains("lebensader_skeleton_trace_decoupling", hits);
        Assert.All(
            s.RelatedConfirmations(machine: "ibm_marrakesh"),
            c => Assert.Equal("ibm_marrakesh", c.Machine));
    }

    [Fact]
    public void Labels_AreRoundTripStableAtTheLoggingBoundary()
    {
        Assert.Equal("all-below-r-star", RegimeVerdict.AllBelowRStar.Label());
        Assert.Equal("all-at-or-above-r-star", RegimeVerdict.AllAtOrAboveRStar.Label());
        Assert.Equal("r-star-band-mixed", RegimeVerdict.RStarBandMixed.Label());
        Assert.Equal("not-addressable", RegimeVerdict.NotAddressable.Label());
    }
}
