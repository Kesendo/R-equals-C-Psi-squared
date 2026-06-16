using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The Symphony seam movement (movement 4): the calibration topology. Gate-first —
/// the γ- and J-anchors recover (γ₀, J) from the dimensionful spectrum and the over-determination
/// gate fires when the spectrum leaves the (XY, Q ≥ Q*(N)) domain. Ground truth:
/// simulations/_seam_movement_review.py.</summary>
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
}
