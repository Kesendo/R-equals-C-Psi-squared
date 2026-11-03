using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class F86HwhmClosedFormClaimTests
{
    [Fact]
    public void Tier_IsTier2Empirical()
    {
        // The floor is derived (C2BareDoubledPtfClosedForm's envelope ratio); the (alpha, beta)
        // per sub-class are a line through the anchors, which Endpoint and Flanking miss, so the
        // lift is observed, not held by an exact law.
        var claim = new F86HwhmClosedFormClaim();
        Assert.Equal(Tier.Tier2Empirical, claim.Tier);
    }

    [Fact]
    public void PredictHwhmRatio_ReproducesEveryAnchor_WithinItsClassResidual()
    {
        // Each anchor is reproduced to its class's own MaxResidual: 1.28e-3 for Endpoint and
        // 4.2e-4 for Flanking, the round-off of the fit for the two-point classes.
        var claim = new F86HwhmClosedFormClaim();
        var e = F86HwhmAlphaExtraction.Extract();
        foreach (var (n, b, q, hwhm) in F86HwhmAlphaExtraction.Anchors)
        {
            var sub = BondSubClassExtensions.Classify(n, b, q);
            double residual = Math.Abs(claim.PredictHwhmRatio(n, b, q) - hwhm);
            Assert.True(residual <= e[sub].MaxResidual + 1e-15,
                $"N={n} b={b} ({sub}): residual {residual:G6} above the class's {e[sub].MaxResidual:G6}");
        }
    }

    [Fact]
    public void BareFloor_IsTheEnvelopeFloor()
    {
        var claim = new F86HwhmClosedFormClaim();
        Assert.Equal(C2BareDoubledPtfClosedForm.EnvelopeHwhmLeftOverXPeak, claim.BareFloor);
    }
}
