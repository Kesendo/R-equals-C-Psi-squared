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

        // Centering with the F1 shift sigma = Sum gamma removes the trace and
        // therefore every cross-site term.
        Assert.Equal(16.0 / 9.0, F1DepolResidualClosedForm.LocalCoefficient);
        Assert.Equal(0.0, F1DepolResidualClosedForm.CrossSiteCoefficient);
    }

    [Theory]
    // Uniform γ = 0.1, mirrors simulations/f1_depol_residual_verify.py section 2 numerics
    // (machine-precision match against the Python framework's palindrome_residual at N=2..5).
    // 4^(N−1) · 0.01 · (16/9)·N.
    [InlineData(2, 0.1, 0.1422222222)]
    [InlineData(3, 0.1, 0.8533333333)]
    [InlineData(4, 0.1, 4.5511111111)]
    [InlineData(5, 0.1, 22.7555555556)]
    public void PredictUniform_MatchesVerificationNumerics(int N, double gamma, double expected)
    {
        double predicted = F1DepolResidualClosedForm.PredictUniform(N, gamma);
        // 16/9 is not exactly representable in binary floating point ⟹ tolerance, not exact equality.
        Assert.Equal(expected, predicted, 9);
    }

    [Theory]
    // Non-uniform γ = [0.05·(k+1)], mirrors verification script numerics.
    [InlineData(2, 0.0888888889)]
    [InlineData(3, 0.9955555556)]
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
        Assert.Contains("centered cross-site coefficient ((Σγ)²)", names);
        Assert.Contains("derivation", names);
        Assert.Contains("orthogonality", names);
        // The two structural surprises must be surfaced as inspectable children.
        Assert.Contains("Π²-decomposition (trivial)", names);
        Assert.Contains("F1-centering shift σ = Σγ", names);
        Assert.Contains("verification", names);
    }
}
