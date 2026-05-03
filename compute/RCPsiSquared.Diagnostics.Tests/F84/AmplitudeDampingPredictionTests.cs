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

    [Fact]
    public void NetHeating_GivesSameMagnitudeAsCoolingByPCC()
    {
        // (γ_↓ − γ_↑)² is symmetric under sign flip; PCC Lemma in PROOF_F84 says cooling-only and
        // heating-only at equal |Δγ| produce identical violation magnitudes.
        var chain = new ChainSystem(N: 4, J: 1.0, GammaZero: 0.05);
        double[] cooling_gT1 = { 0.10, 0.10, 0.10, 0.10 };
        double[] cooling_gPump = { 0.08, 0.08, 0.08, 0.08 };
        double[] heating_gT1 = { 0.08, 0.08, 0.08, 0.08 };
        double[] heating_gPump = { 0.10, 0.10, 0.10, 0.10 };

        double vCool = AmplitudeDampingPrediction.PredictViolation(chain, cooling_gT1, cooling_gPump);
        double vHeat = AmplitudeDampingPrediction.PredictViolation(chain, heating_gT1, heating_gPump);
        Assert.Equal(vCool, vHeat, 12);
    }

    [Fact]
    public void NonUniformDelta_AggregatesAsRootSumOfSquares()
    {
        // Per-site (γ_↓ − γ_↑) = (0.01, 0.02, 0.03) → ‖D_odd‖ = √(Σ(Δγ)²) · 2^(N−1).
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        double[] gT1 = { 0.10, 0.10, 0.10 };
        double[] gPump = { 0.09, 0.08, 0.07 };  // Δγ = (0.01, 0.02, 0.03)
        double v = AmplitudeDampingPrediction.PredictViolation(chain, gT1, gPump);
        double sumSq = 0.01 * 0.01 + 0.02 * 0.02 + 0.03 * 0.03; // 0.0014
        double expected = Math.Sqrt(sumSq) * Math.Pow(2, chain.N - 1);
        Assert.Equal(expected, v, 12);
    }

    [Fact]
    public void Estimate_InvertsPrediction()
    {
        // Forward then inverse: violation → RMS Δγ, with |Δγ|_RMS = √(Σ(Δγ)²/N).
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        double[] gT1 = { 0.10, 0.10, 0.10 };
        double[] gPump = { 0.09, 0.08, 0.07 };
        double v = AmplitudeDampingPrediction.PredictViolation(chain, gT1, gPump);

        double rms = AmplitudeDampingPrediction.EstimateRmsNetCoolingFromViolation(chain, v);
        double sumSq = 0.01 * 0.01 + 0.02 * 0.02 + 0.03 * 0.03;
        double expectedRms = Math.Sqrt(sumSq / chain.N);
        Assert.Equal(expectedRms, rms, 12);
    }
}
