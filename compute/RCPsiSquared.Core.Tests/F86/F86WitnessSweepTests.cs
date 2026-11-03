using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Resonance;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>End-to-end sweep: for (c, N) witnesses of the F86 Statement 2 band, run a real
/// <see cref="ResonanceScan"/> on the default Q grid (the grid the band was read on) and check that
/// each class ratio lies in the band, widened only by the scan's reproducibility
/// (<see cref="UniversalShapeBand.Resolution"/>). The band is the span of exactly these readings, so
/// a failure here means the scan no longer reproduces them.</summary>
public class F86WitnessSweepTests(ITestOutputHelper output)
{
    [Theory]
    [InlineData(2, 5, 1)]
    [InlineData(2, 6, 1)]
    [InlineData(2, 7, 1)]
    [InlineData(2, 8, 1)]
    [InlineData(3, 5, 2)]
    [InlineData(3, 6, 2)]
    public void Sweep_F86_EveryWitnessLiesInTheBand_OnTheDefaultGrid(int chromaticity, int N, int n)
    {
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        Assert.Equal(chromaticity, block.C);

        var curve = new ResonanceScan(block).ComputeKCurve();
        var kb = new F86KnowledgeBase(block);
        var matches = kb.CompareTo(curve);

        var interior = matches[0];
        var endpoint = matches[1];

        output.WriteLine(
            $"c={chromaticity} N={N} n={n} | " +
            $"Interior {interior.Actual:R} (band {UniversalShapeBand.InteriorLow:F4}..{UniversalShapeBand.InteriorHigh:F4}) | " +
            $"Endpoint {endpoint.Actual:R} (band {UniversalShapeBand.EndpointLow:F4}..{UniversalShapeBand.EndpointHigh:F4})");

        Assert.True(interior.Within, interior.Description);
        Assert.True(endpoint.Within, endpoint.Description);
    }

    [Fact]
    public void Sweep_QEpConvergence_CrossesOneOverSqrt2_BetweenN7AndN8_AtC2()
    {
        // For c=2 the inter-channel singular value g_eff = σ_0 has 1/√2 = 2/(2√2) as the
        // relevant scale. Numerical observation 2026-05-02: Q_EP crosses 1/√2 between
        // N=7 (≈ 0.7071, on it bit-exactly) and N=8 (≈ 0.7044, below).
        var qEpValues = new (int N, double Value)[4];
        for (int i = 0; i < 4; i++)
        {
            int N = 5 + i;
            var block = new CoherenceBlock(N, 1, 0.05);
            var kb = new F86KnowledgeBase(block);
            Assert.NotNull(kb.QEp);
            qEpValues[i] = (N, kb.QEp!.Value);
            output.WriteLine($"N={N} c=2 → Q_EP = {kb.QEp.Value:F6}");
        }

        // All four are within 0.02 of 1/√2.
        const double oneOverSqrt2 = 0.70710678;
        foreach (var (N, value) in qEpValues)
            Assert.True(Math.Abs(value - oneOverSqrt2) < 0.02,
                $"N={N} c=2 Q_EP={value:F6} should be within 0.02 of 1/√2");

        // Monotone decreasing all the way through (passes through 1/√2 between N=7 and N=8).
        for (int i = 0; i < 3; i++)
            Assert.True(qEpValues[i].Value > qEpValues[i + 1].Value,
                $"Q_EP should decrease with N: N={qEpValues[i].N} → {qEpValues[i].Value} vs N={qEpValues[i + 1].N} → {qEpValues[i + 1].Value}");

        // N=7 lands ON 1/√2 to within 0.001 (the structural sweet spot).
        Assert.InRange(qEpValues[2].Value, oneOverSqrt2 - 0.001, oneOverSqrt2 + 0.001);
        // N=8 has crossed below.
        Assert.True(qEpValues[3].Value < oneOverSqrt2,
            $"N=8 Q_EP={qEpValues[3].Value:F6} should be below 1/√2 = {oneOverSqrt2:F6}");
    }

}
