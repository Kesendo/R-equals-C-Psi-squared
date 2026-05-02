using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Diagnostics.F84;
using F82Predict = RCPsiSquared.Diagnostics.F82.T1DissipatorPrediction;

namespace RCPsiSquared.Diagnostics.Tests.F84;

public class AmplitudeDampingPredictionTests
{
    [Fact]
    public void DetailedBalance_GivesZeroViolation()
    {
        // γ_↓ = γ_↑ → Δγ = 0 → violation = 0 (thermal equilibrium).
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        double[] gT1 = { 0.05, 0.05, 0.05 };
        double[] gPump = { 0.05, 0.05, 0.05 };
        Assert.Equal(0.0, AmplitudeDampingPrediction.PredictViolation(chain, gT1, gPump), 12);
    }

    [Fact]
    public void PureCooling_ReducesToF82()
    {
        // γ_↑ = 0 → F84 = F82.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        double[] gT1 = { 0.05, 0.05, 0.05 };
        double f84 = AmplitudeDampingPrediction.PredictViolation(chain, gT1);  // pump = null
        double f82 = F82Predict.PredictViolation(chain, gT1);
        Assert.Equal(f82, f84, 12);
    }

    [Fact]
    public void NetCooling_ScalesWithDifference()
    {
        // Δγ uniform = 0.02 → ‖D_odd‖ = 0.02 · √N · 2^(N−1).
        var chain = new ChainSystem(N: 4, J: 1.0, GammaZero: 0.05);
        double[] gT1 = { 0.10, 0.10, 0.10, 0.10 };
        double[] gPump = { 0.08, 0.08, 0.08, 0.08 };
        double v = AmplitudeDampingPrediction.PredictViolation(chain, gT1, gPump);
        double expected = 0.02 * Math.Sqrt(4) * Math.Pow(2, 3); // 0.02 · 2 · 8 = 0.32
        Assert.Equal(expected, v, 10);
    }
}
