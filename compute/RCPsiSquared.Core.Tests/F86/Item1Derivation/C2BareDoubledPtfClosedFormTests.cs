using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using SC = System.Numerics.Complex;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

/// <summary>Bit-exact verification of <see cref="C2BareDoubledPtfClosedForm.EvaluateKb"/>
/// against a brute-force Duhamel evaluator that diagonalises L_2(x) numerically and
/// computes the K_b observable via eigendecomposition + Duhamel integral. The closed
/// form must reproduce the brute evaluator at machine precision across both regimes
/// (pre-EP, post-EP, and the EP limit itself).</summary>
public class C2BareDoubledPtfClosedFormTests
{
    private readonly ITestOutputHelper _out;
    public C2BareDoubledPtfClosedFormTests(ITestOutputHelper output) => _out = output;

    private const double BitExactTolerance = 1e-12;

    private static double BruteKb(double x, double gamma0 = 1.0, double tPeak = 0.25)
    {
        // L_2(Q) with J = Q·γ₀, dimensionless x = Q·g_eff/2 → so Q·g_eff = 2x.
        // Choose g_eff = 2 for convenience; then Q = x, J = x.
        const double g = 2.0;
        double J = x;
        var L = Matrix<SC>.Build.DenseOfArray(new SC[,]
        {
            { new SC(-2 * gamma0, 0), new SC(0,  J * g) },
            { new SC(0,  J * g),       new SC(-6 * gamma0, 0) },
        });
        var Vb = Matrix<SC>.Build.DenseOfArray(new SC[,]
        {
            { SC.Zero, new SC(0, g) },
            { new SC(0, g), SC.Zero },
        });

        var evd = L.Evd();
        var evals = evd.EigenValues.ToArray();
        var R = evd.EigenVectors;
        var Rinv = R.Inverse();
        var rho0 = MathNet.Numerics.LinearAlgebra.Vector<SC>.Build.Dense(new SC[] { SC.One, SC.Zero });

        int n = 2;
        var expLam = new SC[n];
        for (int i = 0; i < n; i++) expLam[i] = SC.Exp(evals[i] * tPeak);
        var I_mat = Matrix<SC>.Build.Dense(n, n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                var diff = evals[j] - evals[i];
                I_mat[i, j] = diff.Magnitude < 1e-10
                    ? tPeak * expLam[i]
                    : (expLam[j] - expLam[i]) / diff;
            }
        var V_eig = Rinv * Vb * R;
        var c0 = Rinv * rho0;
        var drho_eig = MathNet.Numerics.LinearAlgebra.Vector<SC>.Build.Dense(n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                drho_eig[i] += V_eig[i, j] * I_mat[i, j] * c0[j];
        var drho = R * drho_eig;
        var rho_t = R * MathNet.Numerics.LinearAlgebra.Vector<SC>.Build.Dense(new SC[] { expLam[0] * c0[0], expLam[1] * c0[1] });

        SC inner = SC.Zero;
        for (int i = 0; i < n; i++) inner += SC.Conjugate(rho_t[i]) * drho[i];
        // For g = 2 the brute V_b = ∂L/∂J and the closed-form V_x = ∂L/∂x are the
        // SAME matrix (since x = J·g/2 = J when g = 2), so K_b_brute = K_b_closed
        // without further normalisation.
        return 2.0 * inner.Real;
    }

    [Theory]
    [InlineData(0.5)]
    [InlineData(0.721607)]   // x_half (pre-EP)
    [InlineData(0.95)]
    [InlineData(1.5)]
    [InlineData(2.196910)]   // x_peak (post-EP)
    [InlineData(3.0)]
    [InlineData(5.0)]
    public void EvaluateKb_AtSampleX_MatchesBruteDuhamel(double x)
    {
        double closed = C2BareDoubledPtfClosedForm.EvaluateKb(x);
        double brute = BruteKb(x);
        double residual = Math.Abs(closed - brute);
        _out.WriteLine($"x = {x}: closed = {closed:G14}, brute = {brute:G14}, residual = {residual:E2}");
        Assert.True(residual < BitExactTolerance,
            $"x = {x}: closed = {closed:G14}, brute = {brute:G14}, residual = {residual:E2} exceeds tolerance");
    }

    [Fact]
    public void EvaluateKb_AtEpLimit_MatchesMinusFiveOverTwelveTimesEMinus2()
    {
        double closed = C2BareDoubledPtfClosedForm.EvaluateKb(1.0);
        double expected = -5.0 * Math.Exp(-2.0) / 12.0;
        Assert.Equal(expected, closed, precision: 12);
        Assert.Equal(C2BareDoubledPtfClosedForm.KbAtEp, closed, precision: 12);
    }

    [Theory]
    // The Taylor expansion near EP must match the canonical closed form within tolerance
    [InlineData(1.001)]    // post-EP, very close to EP
    [InlineData(0.999)]    // pre-EP, very close to EP
    [InlineData(1.0001)]
    [InlineData(0.9999)]
    public void EvaluateKb_NearEp_StableAcrossSingularity(double x)
    {
        double k = C2BareDoubledPtfClosedForm.EvaluateKb(x);
        double expectedAtEp = -5.0 * Math.Exp(-2.0) / 12.0;
        double residual = Math.Abs(k - expectedAtEp);
        _out.WriteLine($"x = {x}: K_b = {k:G14}, |K_b − K_b(EP)| = {residual:E2}");
        // Smoothness: K_b should be within ~1e-3 of EP value
        Assert.True(residual < 1e-3, $"K_b not stable near EP: {residual:E2}");
    }

    [Fact]
    public void XPeakPrecise_MatchesC2HwhmRatioRoundedConstant()
    {
        // The stored 6-place constant is the rounding of XPeakPrecise.
        double rounded = Math.Round(C2BareDoubledPtfClosedForm.XPeakPrecise, 6);
        Assert.Equal(C2HwhmRatio.BareDoubledPtfXPeak, rounded, precision: 6);
    }

    [Fact]
    public void HwhmLeftOverXPeakPrecise_MatchesC2HwhmRatioConstantToFiveDecimals()
    {
        // C2HwhmRatio.BareDoubledPtfHwhmRatio = 0.671535 (6 places, truncated, not
        // rounded — the precise value 0.671535518 would round to 0.671536). The two
        // agree to 5 decimal places (= 0.67154 vs 0.67154 after that precision).
        Assert.Equal(C2HwhmRatio.BareDoubledPtfHwhmRatio,
            C2BareDoubledPtfClosedForm.HwhmLeftOverXPeakPrecise, precision: 5);
    }

    [Fact]
    public void EvaluateKb_NonPositiveX_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => C2BareDoubledPtfClosedForm.EvaluateKb(0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => C2BareDoubledPtfClosedForm.EvaluateKb(-1.0));
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, new C2BareDoubledPtfClosedForm().Tier);
    }
}
