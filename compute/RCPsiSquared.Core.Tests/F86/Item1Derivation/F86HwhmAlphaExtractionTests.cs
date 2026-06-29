using RCPsiSquared.Core.F86.Item1Derivation;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

/// <summary>The honest, live f86b2 alpha-extraction (the f86b2_robust_extraction arc landed
/// 2026-06-29, "honest and phased" stance). The per-sub-class slope alpha is recomputed live
/// from the 22 F90-bridge anchors with its grid-noise uncertainty, and the Phase-1 verdict
/// (which alpha is a real slope vs extraction noise) is pinned here.</summary>
public class F86HwhmAlphaExtractionTests
{
    [Fact]
    public void Extract_RecomputesEndpointSlope_ReproducesTheHistoricFit()
    {
        // The live polyfit reproduces the historically-frozen Endpoint slope (-0.129110),
        // proving the recompute is the same extraction, not a new convention.
        var e = F86HwhmAlphaExtraction.Extract()[BondSubClass.Endpoint];
        Assert.Equal(-0.129110, e.FittedAlpha, precision: 4);
    }

    [Fact]
    public void Verdicts_MatchThePhase1HonestJudgement()
    {
        var e = F86HwhmAlphaExtraction.Extract();
        Assert.Equal(AlphaVerdict.ResolvedSlope, e[BondSubClass.Endpoint].Verdict);
        Assert.Equal(AlphaVerdict.ResolvedSlope, e[BondSubClass.Flanking].Verdict);
        Assert.Equal(AlphaVerdict.ConstantNoise, e[BondSubClass.Mid].Verdict);
        Assert.Equal(AlphaVerdict.SinglePoint, e[BondSubClass.CentralSelfPaired].Verdict);
        Assert.Equal(AlphaVerdict.EscapeArtefact, e[BondSubClass.Orbit2Escape].Verdict);
        Assert.Equal(AlphaVerdict.SinglePoint, e[BondSubClass.CentralEscapeOrbit3].Verdict);
    }

    [Fact]
    public void ReshapedParams_ZeroSlope_ForNoiseAndSinglePoint()
    {
        // The indefensible classes collapse to a per-class constant lift (alpha = 0).
        Assert.Equal(0.0, F86HwhmAlphaExtraction.ReshapedParams(BondSubClass.Mid).Alpha);
        Assert.Equal(0.0, F86HwhmAlphaExtraction.ReshapedParams(BondSubClass.CentralSelfPaired).Alpha);
        Assert.Equal(0.0, F86HwhmAlphaExtraction.ReshapedParams(BondSubClass.CentralEscapeOrbit3).Alpha);
    }

    [Fact]
    public void ReshapedParams_KeepSlope_ForResolvedAndFlaggedEscape()
    {
        // Endpoint/Flanking keep their resolved slope; Orbit2Escape keeps its (flagged)
        // slope only because its two grid-edge lift levels differ by more than the 0.005 floor.
        Assert.NotEqual(0.0, F86HwhmAlphaExtraction.ReshapedParams(BondSubClass.Endpoint).Alpha);
        Assert.NotEqual(0.0, F86HwhmAlphaExtraction.ReshapedParams(BondSubClass.Flanking).Alpha);
        Assert.NotEqual(0.0, F86HwhmAlphaExtraction.ReshapedParams(BondSubClass.Orbit2Escape).Alpha);
    }

    [Fact]
    public void GridNoiseSigma_FlagsMidAsTheLeastDetermined()
    {
        // Mid's microscopic g_eff lever (~0.011) gives by far the largest slope uncertainty.
        var e = F86HwhmAlphaExtraction.Extract();
        Assert.True(e[BondSubClass.Mid].SigmaAlpha > e[BondSubClass.Endpoint].SigmaAlpha);
        Assert.True(e[BondSubClass.Mid].SigmaAlpha > e[BondSubClass.Flanking].SigmaAlpha);
    }

    [Fact]
    public void ReshapedConstantLift_PredictsItsAnchorsWithinGridFloor()
    {
        // A constant lift (alpha=0) for Mid still reproduces its anchors within 0.005.
        var (alpha, beta) = F86HwhmAlphaExtraction.ReshapedParams(BondSubClass.Mid);
        double predicted = F86HwhmAlphaExtraction.BareFloor + alpha * F86HwhmAlphaExtraction.GEff(1.54) + beta;
        Assert.True(System.Math.Abs(predicted - 0.7469) <= 0.005);
    }
}
