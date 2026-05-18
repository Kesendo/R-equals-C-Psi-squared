using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Tests.F1;

public class F1T1ResidualPi2DecompositionTests
{
    [Fact]
    public void T1ResidualPi2Decomposition_IsTier1Derived_AnchorsProofAndParent()
    {
        var claim = new F1T1ResidualPi2Decomposition();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        // Anchor should point to both the proof file (Step 7 derivation) and the parent total claim.
        Assert.Contains("PROOF_F1_T1_RESIDUAL_CLOSED_FORM", claim.Anchor);
        Assert.Contains("F1T1ResidualClosedForm", claim.Anchor);
    }

    [Fact]
    public void T1ResidualPi2Decomposition_Coefficients_AntiIsLocalOnly_SymCarriesCooperative()
    {
        // Anti side: entire content is the F82/F84 (Z, I)-channel per-site contribution → (1, 0).
        Assert.Equal(1.0, F1T1ResidualPi2Decomposition.AntisymmetricLocalCoefficient);
        Assert.Equal(0.0, F1T1ResidualPi2Decomposition.AntisymmetricCrossCoefficient);

        // Sym side: cooperative cross-site piece + remainder of per-site kernel → (2, 4).
        Assert.Equal(2.0, F1T1ResidualPi2Decomposition.SymmetricLocalCoefficient);
        Assert.Equal(4.0, F1T1ResidualPi2Decomposition.SymmetricCrossCoefficient);

        // Coefficient sums must reproduce the parent F1T1ResidualClosedForm total (3, 4).
        Assert.Equal(
            F1T1ResidualClosedForm.LocalCoefficient,
            F1T1ResidualPi2Decomposition.AntisymmetricLocalCoefficient
                + F1T1ResidualPi2Decomposition.SymmetricLocalCoefficient);
        Assert.Equal(
            F1T1ResidualClosedForm.CrossSiteCoefficient,
            F1T1ResidualPi2Decomposition.AntisymmetricCrossCoefficient
                + F1T1ResidualPi2Decomposition.SymmetricCrossCoefficient);
    }

    [Theory]
    // Pythagorean closure ‖M‖² = ‖M_anti‖² + ‖M_sym‖² must hold bit-exact (Π²-orthogonal split).
    // Uniform γ_T1 = 0.1 at N = 3, 4, 5: parent totals from F1T1ResidualClosedFormTests numerics.
    [InlineData(3, 0.1)]
    [InlineData(4, 0.1)]
    [InlineData(5, 0.1)]
    public void PythagoreanSum_Uniform_RecoversParentTotal(int N, double gammaT1)
    {
        var gammas = Enumerable.Repeat(gammaT1, N).ToArray();
        double anti = F1T1ResidualPi2Decomposition.PredictAntisymmetric(N, gammas);
        double sym = F1T1ResidualPi2Decomposition.PredictSymmetric(N, gammas);
        double total = F1T1ResidualClosedForm.Predict(N, gammas);

        Assert.True(total > 0, $"sanity: parent total at N={N} should be positive");
        Assert.Equal(total, anti + sym, 12);
        // Tight relative tolerance: the split is algebraic, no rounding beyond Math.Pow.
        Assert.True(Math.Abs(anti + sym - total) / total < 1e-12,
            $"Pythagorean closure violated at N={N}, uniform γ={gammaT1}: anti+sym={anti + sym}, total={total}");
    }

    [Fact]
    public void PythagoreanSum_NonUniformGamma_RecoversParentTotal()
    {
        // Non-uniform γ_T1 = [0.1, 0.2, 0.3] at N = 3: the (Σγ)² cooperative term differs from Σγ².
        var gammas = new[] { 0.1, 0.2, 0.3 };
        double anti = F1T1ResidualPi2Decomposition.PredictAntisymmetric(3, gammas);
        double sym = F1T1ResidualPi2Decomposition.PredictSymmetric(3, gammas);
        double total = F1T1ResidualClosedForm.Predict(3, gammas);

        // Sanity: anti = 4² · (0.01 + 0.04 + 0.09) = 16 · 0.14 = 2.24
        Assert.Equal(2.24, anti, 12);
        // Sanity: sym = 4² · [2·0.14 + 4·0.36] = 16 · [0.28 + 1.44] = 16 · 1.72 = 27.52
        Assert.Equal(27.52, sym, 12);
        // Pythagoras: 2.24 + 27.52 = 29.76; parent = 4² · [3·0.14 + 4·0.36] = 16 · 1.86 = 29.76
        Assert.Equal(29.76, total, 12);
        Assert.Equal(total, anti + sym, 12);
    }

    [Fact]
    public void PredictAntisymmetric_UniformAndExplicitListAgree()
    {
        const int N = 4;
        const double gamma = 0.1;
        var explicitList = Enumerable.Repeat(gamma, N).ToArray();
        Assert.Equal(
            F1T1ResidualPi2Decomposition.PredictAntisymmetricUniform(N, gamma),
            F1T1ResidualPi2Decomposition.PredictAntisymmetric(N, explicitList), 12);
    }

    [Fact]
    public void PredictSymmetric_UniformAndExplicitListAgree()
    {
        const int N = 4;
        const double gamma = 0.1;
        var explicitList = Enumerable.Repeat(gamma, N).ToArray();
        Assert.Equal(
            F1T1ResidualPi2Decomposition.PredictSymmetricUniform(N, gamma),
            F1T1ResidualPi2Decomposition.PredictSymmetric(N, explicitList), 12);
    }

    [Fact]
    public void T1ResidualPi2Decomposition_AppearsInF1KnowledgeBaseAtTier1()
    {
        var kb = new F1KnowledgeBase(N: 5);
        Assert.NotNull(kb.T1ResidualPi2Decomposition);
        Assert.Equal(Tier.Tier1Derived, kb.T1ResidualPi2Decomposition.Tier);
        // Discoverable via the tier-1 query so the inspector picks it up.
        var tier1 = kb.ClaimsAtTier(Tier.Tier1Derived).ToList();
        Assert.Contains(tier1, c => c is F1T1ResidualPi2Decomposition);
        // Parent total must also still be there; we are an addition, not a replacement.
        Assert.Contains(tier1, c => c is F1T1ResidualClosedForm);
    }

    [Fact]
    public void Predict_RejectsTooSmallN()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F1T1ResidualPi2Decomposition.PredictAntisymmetricUniform(1, 0.1));
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F1T1ResidualPi2Decomposition.PredictSymmetricUniform(1, 0.1));
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F1T1ResidualPi2Decomposition.PredictAntisymmetric(1, new[] { 0.1 }));
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F1T1ResidualPi2Decomposition.PredictSymmetric(1, new[] { 0.1 }));
    }

    [Fact]
    public void Predict_RejectsGammaLengthMismatch()
    {
        Assert.Throws<ArgumentException>(
            () => F1T1ResidualPi2Decomposition.PredictAntisymmetric(N: 3, new[] { 0.1, 0.2 }));
        Assert.Throws<ArgumentException>(
            () => F1T1ResidualPi2Decomposition.PredictSymmetric(N: 3, new[] { 0.1, 0.2 }));
    }
}
