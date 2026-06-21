using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class ExceptionalPointClockTests
{
    [Fact]
    public void QEp_IsTwoOverGEff()
    {
        Assert.Equal(1.5, ExceptionalPointClock.QEp(4.0 / 3.0), 12);
        Assert.Equal(2.5, ExceptionalPointClock.QEp(0.8), 12);
    }

    [Fact]
    public void QPeak_IsXPeakTimesQEp()
    {
        double qEp = ExceptionalPointClock.QEp(4.0 / 3.0);
        Assert.Equal(2.196910329331 * qEp, ExceptionalPointClock.QPeak(4.0 / 3.0), 9);
    }

    [Fact]
    public void X_IsQOverQEp()
    {
        // x = Q/Q_EP = Q·g_eff/2
        Assert.Equal(1.0, ExceptionalPointClock.X(1.5, 4.0 / 3.0), 12);   // Q = Q_EP
        Assert.Equal(0.5, ExceptionalPointClock.X(0.75, 4.0 / 3.0), 12);
    }

    [Fact]
    public void Decay_RisesToFourGammaThenPins()
    {
        const double g0 = 1.0, gEff = 4.0 / 3.0;   // Q_EP = 1.5
        Assert.Equal(2.0 * g0 * (2.0 - System.Math.Sqrt(0.75)), ExceptionalPointClock.Decay(g0, 0.75, gEff), 9); // x=0.5 → 2.27γ₀
        Assert.Equal(4.0 * g0, ExceptionalPointClock.Decay(g0, 1.5, gEff), 9);   // at the EP, pinned at 4γ₀
        Assert.Equal(4.0 * g0, ExceptionalPointClock.Decay(g0, 3.0, gEff), 9);   // above, still pinned
    }

    [Fact]
    public void Omega_IsZeroBelow_GrowsAbove()
    {
        const double g0 = 1.0, gEff = 4.0 / 3.0;
        Assert.Equal(0.0, ExceptionalPointClock.Omega(g0, 0.75, gEff), 9);       // below the EP: no oscillation
        Assert.Equal(0.0, ExceptionalPointClock.Omega(g0, 1.5, gEff), 9);        // at the EP: still 0
        Assert.Equal(2.0 * g0 * System.Math.Sqrt(1.25), ExceptionalPointClock.Omega(g0, 2.25, gEff), 9); // x=1.5 → 2.236
    }

    [Fact]
    public void RotationAngle_ZeroBelow_LiftsOffAbove()
    {
        const double g0 = 1.0, gEff = 4.0 / 3.0;
        Assert.Equal(0.0, ExceptionalPointClock.RotationAngleDegrees(g0, 1.0, gEff), 9);   // below the EP
        Assert.Equal(0.0, ExceptionalPointClock.RotationAngleDegrees(g0, 1.5, gEff), 9);   // at the EP, θ still 0
        Assert.Equal(29.2, ExceptionalPointClock.RotationAngleDegrees(g0, 2.25, gEff), 1); // above: ≈29.2°
    }

    [Fact]
    public void EigenvectorOverlap_IsMinXOneOverX_PeaksAtTheToyEp()
    {
        // The toy 2×2 rate-channel reduction is genuinely defective at its EP (overlap→1); these numerics
        // are TRUE OF THE TOY. They are NOT a claim about the physical (n,n+1) chain block, which is
        // non-normal on the real axis but has NO real-axis defective EP (F86a-retraction; LocalGlobalEpLink).
        const double gEff = 4.0 / 3.0;   // Q_EP = 1.5
        Assert.Equal(1.0, ExceptionalPointClock.EigenvectorOverlap(1.5, gEff), 9);   // x=1, the toy EP: defective (parallel)
        Assert.Equal(0.5, ExceptionalPointClock.EigenvectorOverlap(0.75, gEff), 9);  // x=0.5
        Assert.Equal(1.0 / 1.5, ExceptionalPointClock.EigenvectorOverlap(2.25, gEff), 9); // x=1.5 → 1/x = 0.667
    }

    [Fact]
    public void Eigenvalues_MatchEpAlgebra()
    {
        // The helper reuses EpAlgebra; this pins the slow-mode decay/omega to the canonical source.
        const double g0 = 1.0, gEff = 4.0 / 3.0;
        foreach (double Q in new[] { 0.75, 1.5, 2.25 })
        {
            var (lamP, lamM) = RCPsiSquared.Core.Resonance.EpAlgebra.SlowestPairEigenvaluesComplex(g0, Q * g0, gEff);
            double slowDecay = System.Math.Min(-lamP.Real, -lamM.Real);
            Assert.Equal(slowDecay, ExceptionalPointClock.Decay(g0, Q, gEff), 9);
        }
    }

    [Fact]
    public void Kb_IsTheResonance_PeaksNearXPeak()
    {
        // K_b reuses the C2 closed form; magnitude is larger near x_peak than far below.
        double gEff = 4.0 / 3.0;
        double atPeak = System.Math.Abs(ExceptionalPointClock.Kb(ExceptionalPointClock.QPeak(gEff), gEff));
        double farBelow = System.Math.Abs(ExceptionalPointClock.Kb(0.3 * ExceptionalPointClock.QEp(gEff), gEff));
        Assert.True(atPeak > farBelow, $"K_b should peak near Q_peak; peak {atPeak} vs farBelow {farBelow}");
    }

    [Fact]
    public void Clock_IsGamma0Invariant_AnglePure_DecayLinear()
    {
        // The clock face is γ₀-free: λ_± = γ₀·(−4 ± √(4 − Q²g_eff²)) at fixed Q, so the Rotation
        // hand θ = arctan(ω/decay) cancels γ₀ exactly (θ = arctan(√(x²−1)/2) above the EP, x = Q/Q_EP),
        // while the Takt hand carries the unit (decay scales linearly with γ₀). Shape lives in Q alone;
        // γ₀ is the tick. (The journey work's "bit-identical at γ₀ = 1.0 / 0.05", now test-pinned,
        // 2026-06-10.)
        const double gEff = 4.0 / 3.0;
        foreach (double q in new[] { 0.75, 1.5, 1.6, 2.25, ExceptionalPointClock.QPeak(gEff), 20.0 })
        {
            Assert.Equal(
                ExceptionalPointClock.RotationAngleDegrees(1.0, q, gEff),
                ExceptionalPointClock.RotationAngleDegrees(0.05, q, gEff), 12);
            Assert.Equal(20.0,
                ExceptionalPointClock.Decay(1.0, q, gEff) / ExceptionalPointClock.Decay(0.05, q, gEff), 9);
        }
        // The angle at the K_b peak: θ_peak = arctan(ξ_peak/2) with ξ_peak = √(x_peak²−1) = 1.956122...,
        // = 44.3646°. Close to but provably NOT 45° (45° would need x = √5 = 2.23607, x_peak = 2.19691).
        Assert.Equal(44.3645555612,
            ExceptionalPointClock.RotationAngleDegrees(1.0, ExceptionalPointClock.QPeak(gEff), gEff), 9);
    }
}
