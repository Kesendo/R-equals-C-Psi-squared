using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The Symphony seam movement (movement 4): the calibration topology. Gate-first —
/// the γ- and J-anchors recover (γ₀, J) from the dimensionful spectrum and the over-determination
/// gate fires when the spectrum leaves the (XY, Q ≥ Q*(N)) domain. The test configs are pinned to the
/// closed-form coherence horizon Q*(N) = {1, √2, 1.879, 2.374} (CoherenceHorizonClaim): N=3,Q=1.5 is
/// protected (gate PASS); N=4,Q=1.5 is below Q*(4)=1.879 (gate FIRES).</summary>
public class SeamMovementTests
{
    // The protected config: N=3, J=0.075, γ=0.05 ⟹ Q=1.5 > Q*(3)=√2 (the ClockHandLadder default).
    private static Symphony Protected3() =>
        new(n: 3, j: 0.075, gamma: 0.05, calibrate: true);

    // Below the horizon: N=4, J=0.075, γ=0.05 ⟹ Q=1.5 < Q*(4)=1.879 (an overdamped mode takes the gap).
    private static Symphony Below4() =>
        new(n: 4, j: 0.075, gamma: 0.05, calibrate: true);

    [Fact]
    public void GammaAnchor_RecoversGamma_InProtectedRegime()
    {
        var seam = Protected3().Seam!;
        Assert.True(seam.Protected);
        Assert.Equal(0.05, seam.GammaRecovered, 9);   // γ₀_rec = gap/2 = γ₀ to ~1e-9
    }

    [Fact]
    public void GammaAnchor_UnderRecovers_BelowTheHorizon()
    {
        // R1: the γ-anchor is NOT regime-free. Below Q*(N) an overdamped q_min mode is slower
        // than 2γ₀, so gap/2 < γ₀ and the regime flag is false.
        var seam = Below4().Seam!;
        Assert.False(seam.Protected);
        Assert.True(seam.GammaRecovered < 0.9 * 0.05);   // under-recovers (ground truth ≈ 0.0204)
    }

    [Fact]
    public void SeamMovement_DoesNotReEvolve_TheTrajectoryIsBuiltOnce()
    {
        var s = Protected3();
        var seam = s.Seam!;
        _ = seam.GammaRecovered;
        _ = seam.Protected;
        _ = ((IInspectable)s).Children.Select(c => c.Summary).ToList();   // touch every lens
        Assert.Equal(1, s.EvolveCount);   // spectrum-only: ONE evolution, the seam reuses it
    }

    [Fact]
    public void Seam_IsAbsent_WhenCalibrateFlagIsOff()
    {
        var s = new Symphony(n: 3, j: 0.075, gamma: 0.05);   // no calibrate
        Assert.Null(s.Seam);
        Assert.DoesNotContain(((IInspectable)s).Children, c => c.DisplayName == "movement: the seam");
    }

    [Fact]
    public void Seam_AppearsAsAChild_WhenCalibrateIsOn()
    {
        var s = Protected3();
        Assert.Contains(((IInspectable)s).Children, c => c.DisplayName == "movement: the seam");
    }

    [Fact]
    public void JAnchor_RecoversJ_InProtectedRegime()
    {
        var seam = Protected3().Seam!;
        Assert.True(seam.XyOk);
        Assert.Equal(0.075, seam.JRecovered, 9);   // J_rec = ω_mem / (2 cos(π/4)) = J to ~1e-9
    }

    [Fact]
    public void JAnchor_IsNotApplicable_UnderNonXyNormalization()
    {
        // R2/Test 7: the band edge 2J·cos(π/(N+1)) is XY-specific; Heisenberg breaks it ~4×.
        var seam = new Symphony(n: 3, j: 0.075, gamma: 0.05,
            hType: HamiltonianType.Heisenberg, calibrate: true).Seam!;
        Assert.False(seam.XyOk);
        Assert.Equal(0.0, seam.JRecovered);   // guarded: no recovery off the XY normalization
        var coh = ((IInspectable)seam).Children.Single(c => c.DisplayName == "seam: coherence-hand");
        Assert.Contains("N/A", coh.Summary);
    }

    [Fact]
    public void Gate_Passes_InTheProtectedXyRegime()
    {
        var seam = Protected3().Seam!;
        Assert.True(seam.GatePass);
        Assert.Equal(1.5, seam.Ratio, 6);           // J_rec/γ₀_rec = Q
        Assert.True(seam.GateResidual < 1e-6);
    }

    [Fact]
    public void Gate_Fires_BelowTheCoherenceHorizon()
    {
        // Test 3 (the test that can fail): N=4, Q=1.5 < Q*(4)=1.879. The gap mode is overdamped,
        // ω_mem → 0, J_rec → 0, ratio → 0 ≠ Q, the gate FAILS.
        var seam = Below4().Seam!;
        Assert.False(seam.Protected);
        Assert.Equal(0.0, seam.JRecovered);
        Assert.Equal(0.0, seam.Ratio);
        Assert.True(seam.GateResidual > 1.0);       // |0 − 1.5| = 1.5
        Assert.False(seam.GatePass);
        var gate = ((IInspectable)seam).Children.Single(c => c.DisplayName == "the gate");
        Assert.Contains("FIRES", gate.Summary);
    }

    [Fact]
    public void Gate_Fires_OnNonXyNormalization()
    {
        var seam = new Symphony(n: 3, j: 0.075, gamma: 0.05,
            hType: HamiltonianType.Heisenberg, calibrate: true).Seam!;
        Assert.False(seam.GatePass);
    }

    [Fact]
    public void ChainCollapse_PredictsTau_AndRoundTripsAgainstTheModel()
    {
        var seam = Protected3().Seam!;
        // τ = 1/(2γ₀) = 1/gap; pin via the takt reading and round-trip against the model's own clock.
        Assert.True(seam.Protected);
        Assert.Equal(1.0 / (2.0 * 0.05), seam.TauPredicted, 6);   // 1/(2·0.05) = 10
        var collapse = ((IInspectable)seam).Children.Single(c => c.DisplayName == "the chain collapse");
        Assert.Contains("τ", collapse.Summary);
    }

    [Fact]
    public void ChainCollapse_PredictsTPeak_AtN2BellPlus()
    {
        // t_peak = 1/(4γ₀) is the N=2 Bell+ closed form (state-class-specific).
        var seam = new Symphony(n: 2, j: 0.1, gamma: 0.05, calibrate: true).Seam!;   // Q=2 > Q*(2)=1
        Assert.True(seam.Protected);
        Assert.Equal(1.0 / (4.0 * 0.05), seam.TPeakPredicted, 6);   // 1/(4·0.05) = 5
    }

    [Fact]
    public void Summary_SurfacesTheGateVerdict()
    {
        Assert.Contains("PASS", Protected3().Seam!.Summary);
        Assert.Contains("FIRES", Below4().Seam!.Summary);
    }
}
