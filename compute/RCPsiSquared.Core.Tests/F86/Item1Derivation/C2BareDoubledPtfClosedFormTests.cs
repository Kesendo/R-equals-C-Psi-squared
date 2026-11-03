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

    [Theory]
    [InlineData(0.3, 0.4)]
    [InlineData(0.469, 1.9)]   // near the envelope's x_half and its τ*
    [InlineData(0.999, 3.0)]   // series branch, pre-EP
    [InlineData(1.001, 0.7)]   // series branch, post-EP
    [InlineData(1.714, 1.3)]   // near the envelope peak
    [InlineData(2.5, 2.4)]
    [InlineData(5.0, 0.2)]
    public void EvaluateKb_AnyTau_MatchesBruteDuhamel(double x, double tau)
    {
        // Two routes: the closed form / series against the eigendecomposition Duhamel at t = τ/4.
        // Their difference is the eigensolver's rounding, gated at the same 1e-12 as τ = 1.
        double closed = C2BareDoubledPtfClosedForm.EvaluateKb(x, tau);
        double brute = BruteKb(x, tPeak: tau / 4.0);
        _out.WriteLine($"x = {x}, τ = {tau}: closed = {closed:G17}, brute = {brute:G17}");
        Assert.True(Math.Abs(closed - brute) < BitExactTolerance, $"x = {x}, τ = {tau}: residual {Math.Abs(closed - brute):E2}");
    }

    [Theory]
    [InlineData(0.5)]
    [InlineData(1.0)]
    [InlineData(2.0)]
    public void EvaluateKb_AtTheEp_IsTheLeadingSeriesTerm(double tau)
    {
        // At the EP (s = 0) only c_2 survives: K_b(1, τ) = −e^(−2τ)·(τ³/3 + τ⁴/12).
        double expected = -Math.Exp(-2.0 * tau) * (tau * tau * tau / 3.0 + tau * tau * tau * tau / 12.0);
        Assert.Equal(expected, C2BareDoubledPtfClosedForm.EvaluateKb(1.0, tau), 15);
    }

    [Theory]
    [InlineData(0.5)]
    [InlineData(1.0)]
    [InlineData(2.196910)]
    public void EvaluateKb_AtTauOne_IsTheFixedTimeForm(double x)
    {
        Assert.Equal(C2BareDoubledPtfClosedForm.EvaluateKb(x), C2BareDoubledPtfClosedForm.EvaluateKb(x, 1.0), 15);
    }

    [Fact]
    public void EnvelopePeak_MatchesThePythonSecondRoute()
    {
        // simulations/f86_envelope_second_route.py --bare (expm Duhamel, τ scan + Brent, x scan +
        // Brent, x_half by brentq to 1e-15): x_peak 1.7141131714205489, x_half 0.4688283327303858,
        // ratio 0.7264892770517302. x_peak sits on a maximum flat to second order and is resolved to about √ε relative on either side, so it and the
        // ratio are gated at 1e-7; x_half sits on a flank and is gated at 1e-12.
        var (xPeak, kMax, xHalf, ratio) = C2BareDoubledPtfClosedForm.EnvelopePeak();
        _out.WriteLine($"x_peak = {xPeak:G17}, K_max = {kMax:G17}, x_half = {xHalf:G17}, ratio = {ratio:G17}");
        Assert.Equal(1.7141131714205489, xPeak, 1e-7);
        Assert.Equal(0.4688283327303858, xHalf, 1e-12);
        Assert.Equal(0.7264892770517302, ratio, 1e-7);
        Assert.Equal(0.09030018438431675, kMax, 1e-14);
        Assert.Equal(C2BareDoubledPtfClosedForm.EnvelopeXPeak, xPeak, 1e-7);
        Assert.Equal(C2BareDoubledPtfClosedForm.EnvelopeHwhmLeftOverXPeak, ratio, 1e-7);
    }

    [Fact]
    public void EnvelopeKb_IsNeverBelowTheFixedTimeSlice()
    {
        foreach (double x in new[] { 0.3, 0.72, 1.0, 1.714, 2.197, 4.0 })
            Assert.True(C2BareDoubledPtfClosedForm.EnvelopeKb(x).KMax >= Math.Abs(C2BareDoubledPtfClosedForm.EvaluateKb(x)));
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, new C2BareDoubledPtfClosedForm().Tier);
    }
}
