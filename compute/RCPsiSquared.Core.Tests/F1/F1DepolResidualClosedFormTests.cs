using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Tests.F1;

public class F1DepolResidualClosedFormTests
{
    [Fact]
    public void DepolResidualClosedForm_IsTier1Derived_WithProofAnchor()
    {
        var claim = new F1DepolResidualClosedForm();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.Contains("PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM", claim.Anchor);
    }

    [Fact]
    public void DepolResidualClosedForm_Constants_MatchProofDerivation()
    {
        // Per-site kernel: ‖M_l‖²_F = 160/9, |tr(M_l)|² = 64 (Step 3 of the proof).
        Assert.Equal(160.0 / 9.0, F1DepolResidualClosedForm.PerSiteFrobeniusSquared);
        Assert.Equal(64.0, F1DepolResidualClosedForm.PerSiteTraceSquared);

        // Closed-form coefficients: c1 = ‖M_l‖² − |tr(M_l)|²/4 = 160/9 − 16 = 16/9, c2 = |tr(M_l)|²/4 = 16.
        Assert.Equal(16.0 / 9.0, F1DepolResidualClosedForm.LocalCoefficient);
        Assert.Equal(16.0, F1DepolResidualClosedForm.CrossSiteCoefficient);

        // Compare as rationals × 9 to dodge IEEE 754 representation of 16/9.
        Assert.Equal(
            9.0 * (F1DepolResidualClosedForm.PerSiteFrobeniusSquared - F1DepolResidualClosedForm.PerSiteTraceSquared / 4.0),
            9.0 * F1DepolResidualClosedForm.LocalCoefficient,
            12);
        Assert.Equal(
            F1DepolResidualClosedForm.PerSiteTraceSquared / 4.0,
            F1DepolResidualClosedForm.CrossSiteCoefficient);
    }

    [Theory]
    // Uniform γ = 0.1, mirrors simulations/_f1_depol_residual_verify.py section 2 numerics
    // (bit-exact match against the Python framework's palindrome_residual at N=2..5).
    // 4^(N−1) · 0.01 · ((16/9)·N + 16·N²).
    [InlineData(2, 0.1, 2.7022222222)]    // 4^1 · 0.01 · ((16/9)·2 + 16·4) = 4 · 0.01 · (32/9 + 64) = 4 · 0.01 · 67.5555…
    [InlineData(3, 0.1, 23.8933333333)]   // 4^2 · 0.01 · ((16/9)·3 + 16·9) = 16 · 0.01 · (16/3 + 144)
    [InlineData(4, 0.1, 168.3911111111)]  // 4^3 · 0.01 · ((16/9)·4 + 16·16) = 64 · 0.01 · (64/9 + 256)
    [InlineData(5, 0.1, 1046.7555555556)] // 4^4 · 0.01 · ((16/9)·5 + 16·25) = 256 · 0.01 · (80/9 + 400)
    public void PredictUniform_MatchesVerificationNumerics(int N, double gamma, double expected)
    {
        double predicted = F1DepolResidualClosedForm.PredictUniform(N, gamma);
        // 16/9 is irrational in IEEE 754 ⟹ absolute tolerance, not exact equality.
        Assert.Equal(expected, predicted, 9);
    }

    [Theory]
    // Non-uniform γ = [0.05·(k+1)], mirrors verification script numerics.
    [InlineData(2, 1.5288888889)]    // γ=(0.05,0.10): Σγ²=0.0125, (Σγ)²=0.0225 → 4·((16/9)·0.0125 + 16·0.0225)
    [InlineData(3, 24.0355555556)]   // γ=(0.05,0.10,0.15): Σγ²=0.035, (Σγ)²=0.09 → 16·((16/9)·0.035 + 16·0.09)
    public void Predict_NonUniformGamma_MatchesVerificationNumerics(int N, double expected)
    {
        var gammas = new double[N];
        for (int l = 0; l < N; l++) gammas[l] = 0.05 * (l + 1);
        double predicted = F1DepolResidualClosedForm.Predict(N, gammas);
        Assert.Equal(expected, predicted, 9);
    }

    [Fact]
    public void PredictUniform_MatchesPredict_ExplicitList()
    {
        const int N = 4;
        const double gamma = 0.1;
        var explicitList = Enumerable.Repeat(gamma, N).ToArray();
        Assert.Equal(F1DepolResidualClosedForm.PredictUniform(N, gamma),
                     F1DepolResidualClosedForm.Predict(N, explicitList), 12);
    }

    [Fact]
    public void Predict_RejectsTooSmallN()
    {
        var exUniform = Assert.Throws<ArgumentOutOfRangeException>(
            () => F1DepolResidualClosedForm.PredictUniform(1, 0.1));
        Assert.Contains("N must be ≥ 2", exUniform.Message);

        var exList = Assert.Throws<ArgumentOutOfRangeException>(
            () => F1DepolResidualClosedForm.Predict(1, new[] { 0.1 }));
        Assert.Contains("N must be ≥ 2", exList.Message);
    }

    [Fact]
    public void Predict_RejectsGammaLengthMismatch()
    {
        var ex = Assert.Throws<ArgumentException>(
            () => F1DepolResidualClosedForm.Predict(N: 3, new[] { 0.1, 0.2 }));
        Assert.Contains("must equal N", ex.Message);
    }

    [Fact]
    public void DepolResidualClosedForm_ExtraChildren_ExposeKernelConstants()
    {
        var claim = new F1DepolResidualClosedForm();
        var names = claim.Children.Select(c => c.DisplayName).ToList();
        Assert.Contains("tier", names);
        Assert.Contains("anchor", names);
        Assert.Contains("statement", names);
        Assert.Contains("per-site ‖M_l‖²_F (γ=1)", names);
        Assert.Contains("per-site |tr(M_l)|² (γ=1)", names);
        Assert.Contains("local coefficient (Σγ²)", names);
        Assert.Contains("cross-site coefficient ((Σγ)²)", names);
        Assert.Contains("derivation", names);
        Assert.Contains("orthogonality", names);
        // The two structural surprises must be surfaced as inspectable children.
        Assert.Contains("Π²-decomposition (trivial)", names);
        Assert.Contains("F1 σ-shift = 0", names);
        Assert.Contains("verification", names);
    }
}
