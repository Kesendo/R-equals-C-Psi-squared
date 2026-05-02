using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F49;
using RCPsiSquared.Diagnostics.F81;
using RCPsiSquared.Diagnostics.F82;

namespace RCPsiSquared.Diagnostics.Tests.F82;

public class T1DissipatorPredictionTests
{
    [Theory]
    [InlineData(2, 0.05, 0.141421356)]   // 0.05·√2·2 = 0.05·2.828 = 0.14142
    [InlineData(3, 0.05, 0.346410162)]   // 0.05·√3·4 = 0.05·1.7321·4 = 0.34641
    [InlineData(3, 0.10, 0.692820323)]   // 0.10·√3·4 = 0.6928
    public void PredictViolation_UniformT1_MatchesClosedForm(int N, double gammaT1, double expected)
    {
        // ‖D_{T1, odd}‖_F = γ_T1 · √N · 2^(N−1)
        var chain = new ChainSystem(N, J: 1.0, GammaZero: 0.05);
        double v = T1DissipatorPrediction.PredictViolation(chain, Enumerable.Repeat(gammaT1, N).ToArray());
        Assert.Equal(expected, v, 6);
    }

    [Fact]
    public void Predict_MatchesActualF81Violation()
    {
        // Cross-check: F82 closed form should match the F81-violation actually computed by
        // PiDecomposition for a truly Hamiltonian under T1.
        var chain = new ChainSystem(N: 2, J: 1.0, GammaZero: 0.05);
        double[] gT1 = { 0.03, 0.05 };
        var truly = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
        };
        var d = PiDecomposition.Decompose(chain, truly, gT1);
        double predicted = T1DissipatorPrediction.PredictViolation(chain, gT1);
        // The F81 violation for "truly H + T1" should match F82 within ~3% (truly means
        // M_anti = L_{H_odd} = 0, so violation = ‖M_anti‖ = ‖D_{T1, odd}‖).
        Assert.InRange(d.F81Violation, 0.97 * predicted, 1.03 * predicted);
    }

    [Fact]
    public void Estimate_InvertsPrediction()
    {
        var chain = new ChainSystem(N: 4, J: 1.0, GammaZero: 0.05);
        double gammaT1 = 0.07;
        double v = T1DissipatorPrediction.PredictViolation(chain, Enumerable.Repeat(gammaT1, 4).ToArray());
        double recovered = T1DissipatorPrediction.EstimateRmsT1FromViolation(chain, v);
        Assert.Equal(gammaT1, recovered, 10);
    }
}
