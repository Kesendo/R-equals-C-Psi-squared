using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Tests.F1;

public class F1T1ResidualClosedFormTests
{
    [Fact]
    public void T1ResidualClosedForm_IsTier1Derived_WithProofAnchor()
    {
        var claim = new F1T1ResidualClosedForm();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.Contains("PROOF_F1_T1_RESIDUAL_CLOSED_FORM", claim.Anchor);
    }

    [Fact]
    public void T1ResidualClosedForm_Constants_MatchProofDerivation()
    {
        // Per-site kernel: ‖M_l‖²_F = 7, |tr(M_l)|² = 16 (Step 3 of the proof).
        Assert.Equal(7.0, F1T1ResidualClosedForm.PerSiteFrobeniusSquared);
        Assert.Equal(16.0, F1T1ResidualClosedForm.PerSiteTraceSquared);

        // Closed-form coefficients: c1 = ‖M_l‖² − |tr(M_l)|²/4 = 3, c2 = |tr(M_l)|²/4 = 4.
        Assert.Equal(3.0, F1T1ResidualClosedForm.LocalCoefficient);
        Assert.Equal(4.0, F1T1ResidualClosedForm.CrossSiteCoefficient);
        Assert.Equal(
            F1T1ResidualClosedForm.PerSiteFrobeniusSquared - F1T1ResidualClosedForm.PerSiteTraceSquared / 4.0,
            F1T1ResidualClosedForm.LocalCoefficient);
        Assert.Equal(
            F1T1ResidualClosedForm.PerSiteTraceSquared / 4.0,
            F1T1ResidualClosedForm.CrossSiteCoefficient);
    }

    [Theory]
    // Uniform γ_T1 = 0.1, mirrors simulations/f1_t1_residual_verify.py section 2 numerics
    // (bit-exact match against the Python framework's palindrome_residual at N=2..5).
    [InlineData(2, 0.1, 0.880000)]    // 4^1 · 0.01 · (3·2 + 4·4) = 4 · 0.01 · 22 = 0.88
    [InlineData(3, 0.1, 7.200000)]    // 4^2 · 0.01 · (3·3 + 4·9) = 16 · 0.01 · 45 = 7.2
    [InlineData(4, 0.1, 48.640000)]   // 4^3 · 0.01 · (3·4 + 4·16) = 64 · 0.01 · 76 = 48.64
    [InlineData(5, 0.1, 294.400000)]  // 4^4 · 0.01 · (3·5 + 4·25) = 256 · 0.01 · 115 = 294.4
    public void PredictUniform_MatchesVerificationNumerics(int N, double gammaT1, double expected)
    {
        double predicted = F1T1ResidualClosedForm.PredictUniform(N, gammaT1);
        Assert.Equal(expected, predicted, 9);
    }

    [Theory]
    // Non-uniform γ_T1 = [0.05, 0.10, ..., 0.05·N], mirrors verification script numerics.
    [InlineData(2, 0.510000)]    // γ=(0.05,0.10): Σγ²=0.0125, (Σγ)²=0.0225 → 4·(3·0.0125+4·0.0225) = 0.51
    [InlineData(3, 7.440000)]    // γ=(0.05,0.10,0.15): Σγ²=0.035, (Σγ)²=0.09 → 16·(0.105+0.36) = 7.44
    [InlineData(4, 78.400000)]   // γ=(0.05,0.10,0.15,0.20): Σγ²=0.075, (Σγ)²=0.25 → 64·(0.225+1.0) = 78.4
    [InlineData(5, 681.600000)]  // γ=(0.05,…,0.25): Σγ²=0.1375, (Σγ)²=0.5625 → 256·(0.4125+2.25) = 681.6
    public void Predict_NonUniformGamma_MatchesVerificationNumerics(int N, double expected)
    {
        var gammas = new double[N];
        for (int l = 0; l < N; l++) gammas[l] = 0.05 + 0.05 * l;
        double predicted = F1T1ResidualClosedForm.Predict(N, gammas);
        Assert.Equal(expected, predicted, 9);
    }

    [Fact]
    public void Predict_UniformAndExplicitListAgree()
    {
        const int N = 4;
        const double gamma = 0.1;
        var explicitList = Enumerable.Repeat(gamma, N).ToArray();
        Assert.Equal(F1T1ResidualClosedForm.PredictUniform(N, gamma),
                     F1T1ResidualClosedForm.Predict(N, explicitList), 12);
    }

    [Fact]
    public void Predict_RejectsTooSmallN()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F1T1ResidualClosedForm.PredictUniform(1, 0.1));
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F1T1ResidualClosedForm.Predict(1, new[] { 0.1 }));
    }

    [Fact]
    public void Predict_RejectsGammaLengthMismatch()
    {
        Assert.Throws<ArgumentException>(
            () => F1T1ResidualClosedForm.Predict(N: 3, new[] { 0.1, 0.2 }));
    }

    [Fact]
    public void T1ResidualClosedForm_AppearsInF1KnowledgeBaseAtTier1()
    {
        var kb = new F1KnowledgeBase(N: 5);
        Assert.NotNull(kb.T1ResidualClosedForm);
        Assert.Equal(Tier.Tier1Derived, kb.T1ResidualClosedForm.Tier);
        // Discoverable via the tier-1 query so the inspector picks it up.
        var tier1 = kb.ClaimsAtTier(Tier.Tier1Derived).ToList();
        Assert.Contains(tier1, c => c is F1T1ResidualClosedForm);
    }

    [Fact]
    public void T1ResidualClosedForm_ExtraChildren_ExposeKernelConstants()
    {
        var claim = new F1T1ResidualClosedForm();
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
        Assert.Contains("verification", names);
    }
}
