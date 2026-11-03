using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.F86.Item1Derivation;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

/// <summary>The live f86b2 alpha extraction on the envelope anchors: per BondSubClass the line
/// lift = alpha · g_eff + beta, and what the anchors can say about it.</summary>
public class F86HwhmAlphaExtractionTests
{
    [Fact]
    public void OrbitAnchors_AreThePerF71OrbitObservationValues()
    {
        // One source of numbers: every anchor's Q_peak is the frozen orbit table's literal.
        var obs = new PerF71OrbitObservation();
        foreach (var (n, b, q, _) in F86HwhmAlphaExtraction.OrbitAnchors)
        {
            var w = obs.Witnesses.Single(x => x.Chromaticity == 2 && x.N == n);
            Assert.Equal(w.QPeakPerOrbit[b], q);
        }
    }

    [Fact]
    public void Anchors_CoverEveryBondOfN5To8()
    {
        Assert.Equal(4 + 5 + 6 + 7, F86HwhmAlphaExtraction.Anchors.Count);
    }

    [Fact]
    public void Verdicts_FollowTheCountOfDistinctGEff()
    {
        var e = F86HwhmAlphaExtraction.Extract();
        Assert.Equal(AlphaVerdict.Overdetermined, e[BondSubClass.Endpoint].Verdict);
        Assert.Equal(4, e[BondSubClass.Endpoint].DistinctG);
        Assert.Equal(AlphaVerdict.Overdetermined, e[BondSubClass.Flanking].Verdict);
        Assert.Equal(4, e[BondSubClass.Flanking].DistinctG);
        Assert.Equal(AlphaVerdict.TwoPoint, e[BondSubClass.Mid].Verdict);
        Assert.Equal(AlphaVerdict.TwoPoint, e[BondSubClass.CentralSelfPaired].Verdict);
        // On the canonical reading no peak lies past Q = 4: the far-maximum classes are empty.
        Assert.False(e.ContainsKey(BondSubClass.Orbit2Escape));
        Assert.False(e.ContainsKey(BondSubClass.CentralEscapeOrbit3));
        Assert.Throws<InvalidOperationException>(() => F86HwhmAlphaExtraction.FitParams(BondSubClass.Orbit2Escape));
    }

    [Fact]
    public void EndpointAndFlanking_MissALine_FarAboveTheAnchorResolution()
    {
        // The four anchors of each class (N=5..8) are resolved to ~1e-10; a least-squares line
        // misses them by 1.28e-3 (Endpoint) and 4.15e-4 (Flanking), so the linear form is not
        // exact in either. The expected values are Python's polyfit on the second route's anchors,
        // which differ from these literals by up to 1e-10 (the scan resolution); the pins' 1e-9
        // covers that, five orders below the residuals themselves.
        var e = F86HwhmAlphaExtraction.Extract();
        Assert.Equal(0.0012811015073008747, e[BondSubClass.Endpoint].MaxResidual, 1e-9);
        Assert.Equal(0.00041510993185450706, e[BondSubClass.Flanking].MaxResidual, 1e-9);
    }

    [Fact]
    public void TwoPointClasses_ReproduceTheirAnchors()
    {
        // A line through two points reproduces them up to the round-off of the fit: alpha·g and
        // beta are of size 0.2..0.3, so a few ulp of those (~5e-17 each).
        var e = F86HwhmAlphaExtraction.Extract();
        foreach (var sub in new[] { BondSubClass.Mid, BondSubClass.CentralSelfPaired })
            Assert.True(e[sub].MaxResidual < 1e-15, $"{sub}: {e[sub].MaxResidual}");
    }

    [Fact]
    public void BareFloor_IsTheEnvelopeFloor()
    {
        Assert.Equal(C2BareDoubledPtfClosedForm.EnvelopeHwhmLeftOverXPeak, F86HwhmAlphaExtraction.BareFloor);
    }
}
