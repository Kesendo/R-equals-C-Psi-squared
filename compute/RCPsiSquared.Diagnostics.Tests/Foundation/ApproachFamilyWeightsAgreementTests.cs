using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The typed claim (<see cref="ApproachFamilyCarrierClaim"/>) and the live primitive
/// (<see cref="OddHarmonicApproach"/>) both carry the approach family's weights. They are the same
/// doubles: w₀ = s(1 − s²/2)/3 is computed as s·(1 − s·s/2)/3 in one and s·(1 − 0.5·s·s)/3 in the
/// other, and halving is exact in binary64 while s²/2 stays normal (s ≳ 2.1e-154), so the two
/// agree bit for bit; w₁ = s·s·s/6 is the same
/// expression in both.</summary>
public class ApproachFamilyWeightsAgreementTests
{
    [Theory]
    [InlineData(0.0)]
    [InlineData(0.1)]
    [InlineData(0.3)]
    [InlineData(0.5)]
    [InlineData(0.75)]
    [InlineData(0.9)]
    [InlineData(0.987654321)]
    [InlineData(1.0)]
    public void CoreClaimAndDiagnosticsPrimitive_CarryTheSameWeights(double s)
    {
        var (w0, w1) = ApproachFamilyCarrierClaim.Weights(s);
        var (carrier, harmonic) = OddHarmonicApproach.Weights(s);
        Assert.True(w0 == carrier, $"w₀: {w0:R} vs {carrier:R}");
        Assert.True(w1 == harmonic, $"w₁: {w1:R} vs {harmonic:R}");
    }
}
