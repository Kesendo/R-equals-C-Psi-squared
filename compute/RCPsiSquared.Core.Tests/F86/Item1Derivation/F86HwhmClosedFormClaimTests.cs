using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class F86HwhmClosedFormClaimTests
{
    [Fact]
    public void Tier_IsTier1Candidate()
    {
        // Tier-reviewed 2026-05-16: was Tier1Derived, downgraded to Tier1Candidate because
        // the (alpha, beta) per sub-class are fitted via polyfit on N=5..8 anchors, not
        // derived from F89/F90 structure. Bare floor 0.671535 IS derived (C2BareDoubledPtfClosedForm),
        // and the linear-in-g_eff form is the candidate structure, but the 12 fit values
        // remain phenomenological. Tier1Derived requires analytical derivation of (alpha, beta).
        var claim = new F86HwhmClosedFormClaim();
        Assert.Equal(Tier.Tier1Candidate, claim.Tier);
    }

    [Theory]
    [InlineData(5, 0, 2.40, 0.7700)]   // Endpoint
    [InlineData(5, 1, 1.50, 0.7454)]   // Flanking (N=5 b=1: classifier returns Flanking, not Mid)
    [InlineData(6, 0, 2.52, 0.7737)]   // Endpoint
    [InlineData(6, 2, 1.44, 0.7449)]   // CentralSelfPaired
    [InlineData(7, 0, 2.53, 0.7738)]   // Endpoint
    [InlineData(7, 2, 1.54, 0.7469)]   // Mid (N=7 b=2: q-peak in Mid range, no escape)
    [InlineData(7, 1, 7.27, 0.9162)]   // Orbit2Escape
    [InlineData(8, 3, 16.79, 0.5778)]  // CentralEscapeOrbit3
    public void PredictHwhmRatio_MatchesEmpiricalWithin0p005(int n, int b, double qPeak, double expected)
    {
        var claim = new F86HwhmClosedFormClaim();
        double predicted = claim.PredictHwhmRatio(n, b, qPeak);
        Assert.True(Math.Abs(predicted - expected) <= 0.005,
            $"Predicted {predicted}, empirical {expected}, residual {Math.Abs(predicted - expected):G6}");
    }

    [Fact]
    public void BareFloor_Equals0p671535()
    {
        var claim = new F86HwhmClosedFormClaim();
        Assert.Equal(0.671535, claim.BareFloor, precision: 6);
    }
}
