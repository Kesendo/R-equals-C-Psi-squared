using System;
using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The painters' movement inside the Symphony: the PTF read off ONE evolution, the piece
/// played a second (defected) and a guard third (δJ/2) time. The assertions guard the contract — the
/// closure Σ ln α stays in the ±0.05 window for canonical XY N=4 δJ=0.02; the K₁ chiral mirror is
/// EXACT (≤1e-10) on the defected run; each movement trajectory is built exactly once; a non-XY chain
/// declines honestly; and absent --defect-bond the Symphony is unchanged (the movement is absent).</summary>
public class PtfMovementTests
{
    private static System.Collections.Generic.List<IInspectable> Children(Symphony s) =>
        ((IInspectable)s).Children.ToList();

    private static PaintersMovement Movement(Symphony s) =>
        Children(s).OfType<PaintersMovement>().Single();

    // The canonical PTF protocol: XY chain, Z-dephasing, |δJ| ≤ 0.1, and a BONDING-MODE
    // initial state (the delocalized k=1 sine mode, F67). The bonding-state class shows up in
    // the CLOSURE (Σ ln α inside the ±0.05 window), NOT in the per-site reliability guard: Bell+
    // also rescales per-site (α ≈ 1, all reliable) but does not close — see
    // StateClass_PerSiteGuardDoesNotRefuseBellPair below.
    private static Symphony Canonical(int n = 4) =>
        new Symphony(n: n, j: 1.0, gamma: 0.05, hType: HamiltonianType.XY,
                     initialState: InitialStateKind.BondingMode, defectBond: 1, deltaJ: 0.02);

    [Fact]
    public void Closure_WithinWindow_ForCanonicalXY_N5()
    {
        var m = Movement(Canonical(5));
        Assert.True(m.HasLenses);
        Assert.Equal(5, m.Reliable.Count(r => r));
        Assert.True(m.ClosureInWindow,
            $"closure Σ ln α (reliable) = {m.ClosureSum} should be in ±{PaintersMovement.ClosureWindow}");
        Assert.Equal(-0.0444, m.ClosureSum, 3);
    }

    [Fact]
    public void Closure_FiniteSize_N4_JustOutsideWindow_ReportedHonestly()
    {
        // At N=4 the same canonical protocol lands at −0.067: all sites reliable, the number
        // finite and mirror-symmetric, but outside ±0.05. The lens must say so, not round it in.
        var m = Movement(Canonical(4));
        Assert.Equal(4, m.Reliable.Count(r => r));
        Assert.False(m.ClosureInWindow);
        Assert.Equal(-0.0671, m.ClosureSum, 3);
    }

    [Fact]
    public void StateClass_PerSiteGuardDoesNotRefuseBellPair()
    {
        // CORRECTED 2026-06-15 (arc ptf_bonding_class_guard). The prior assertion "the guard refuses
        // every Bell+ site" was a MINIMIZER-TRAP artifact, not physics: a bare golden-section trapped
        // on Bell+'s featureless trajectories (α ≈ 4–6, junk), and the |f| ≤ FMax sanity check caught
        // the junk → 0 reliable, which LOOKED like a state-class veto. With the grid-seed global fit
        // (FitAlpha), Bell+ rescales just like the bonding mode at the per-site level: α ≈ 1.016 on
        // every site, all RELIABLE (brute-grid + Python perspectives_panel agree; the global min is
        // α ≈ 1.016, MSE ~100–300× below the trapped α ≈ 4–6). So the per-site reliability guard does
        // NOT enforce the bonding-state class. The class distinction lives in the CLOSURE (Σ ln α):
        // Bell+ stays OUT of the ±0.05 window where the canonical bonding mode closes (N=5: −0.0444).
        var s = new Symphony(n: 4, j: 1.0, gamma: 0.05, hType: HamiltonianType.XY,
                             initialState: InitialStateKind.BellPair, defectBond: 1, deltaJ: 0.02);
        var m = Movement(s);
        Assert.Equal(4, m.Reliable.Count(r => r));            // the guard does NOT refuse Bell+
        Assert.All(m.Alphas, x => Assert.InRange(x, 0.9, 1.1));  // it rescales per-site, like bonding
        Assert.False(m.ClosureInWindow);                      // but Bell+ does not close (the discriminator)
    }

    [Fact]
    public void ChiralMirror_Exact_OnDefectedRun()
    {
        var m = Movement(Canonical());
        Assert.True(m.ChiralMirrorExact,
            $"K₁ mirror deviation {m.ChiralMirrorDeviation} should be ≤ {PaintersMovement.ChiralExactTol}");
        Assert.True(m.ChiralMirrorDeviation <= 1e-10);
    }

    [Fact]
    public void EachTrajectory_BuiltExactlyOnce()
    {
        var m = Movement(Canonical());
        // Touch everything that reads a trajectory repeatedly.
        _ = m.Summary;
        _ = m.Children.Select(c => c.Summary).ToList();
        _ = m.Alphas;
        _ = m.Reliable;
        _ = m.ChiralMirrorDeviation;
        _ = m.ClosureSum;
        // Four runs: clean P_A, defected P_B, guard P_guard, and the K₁ defected run — each once.
        Assert.Equal(4, m.BuildCount);
    }

    [Fact]
    public void ReadBack_DecodesOwnProfile_MatchesCanonicalCase()
    {
        // The closed-loop self-test: the painters' movement calibrates a decoder for its own (N, J, γ)
        // at the FIXED canonical δJ_cal (= 0.02) and decodes its OWN α-profile. For the canonical case
        // (N=4, bond 1, δJ=0.02) the movement's defect IS the calibration point, so the decoder returns
        // bond 1 with an exact match — the existing contract stays green.
        var m = Movement(Canonical(4));
        var readBack = m.Children.Single(c => c.DisplayName == "the decoder reads back");
        Assert.Contains("bond 1", readBack.Summary);
        Assert.Contains("match", readBack.Summary);
    }

    [Fact]
    public void ReadBack_OffCalibration_GivesHonestNonzeroNumbers()
    {
        // Fix A: the read-back must NOT calibrate at the movement's own δJ (that is circular — the
        // dictionary entry for this bond would be the very profile being decoded, residual exactly 0.0,
        // strength error exactly 0%). It calibrates at the fixed canonical δJ_cal = 0.02. An
        // off-calibration movement (bond 0, δJ ≠ δJ_cal) must read back to bond 0, with a strength error
        // strictly in (0, 25] % and a residual strictly > 0 — honest numbers, not the self-fulfilling 0/0.
        //
        // δJ = 0.025 (a small step off the 0.02 calibration point) is used rather than 0.05: the C#
        // painters' EDGE-bond (bond 0) α-response is steeply nonlinear — its f-profile peaks at the far
        // site, so δJ=0.05 already sits outside the linear window for this bond (the projected δĴ runs to
        // ~0.09, an 80% error, and the read even becomes ambiguous). At δJ=0.025 the bond is still well
        // inside the window: the point of the test is the no-longer-circular contract, which only needs
        // δJ ≠ δJ_cal, so it is exercised in the regime where the linear dictionary is honest.
        const double offDeltaJ = 0.025;
        var s = new Symphony(n: 5, j: 1.0, gamma: 0.05, hType: HamiltonianType.XY,
                             initialState: InitialStateKind.BondingMode, defectBond: 0, deltaJ: offDeltaJ);
        var m = Movement(s);
        Assert.True(m.HasLenses);

        var dec = DefectDecoder.Calibrate(5, 1.0, 0.05, DefectDecoder.DefaultDeltaJCal);
        var r = dec.Decode(m.Alphas);
        Assert.Equal(0, r.Bond);
        Assert.True(r.Residual > 0.0, $"off-calibration residual should be strictly positive, got {r.Residual}");
        double strengthErr = Math.Abs(r.DeltaJ - offDeltaJ) / offDeltaJ;
        Assert.True(strengthErr > 0.0 && strengthErr <= 0.25,
            $"off-calibration strength error should be in (0, 25] %, got {strengthErr:P2} (δĴ={r.DeltaJ})");

        // And the live read-back node text is honest, not the self-fulfilling 0.0 / 0%.
        var readBack = m.Children.Single(c => c.DisplayName == "the decoder reads back");
        Assert.Contains("bond 0", readBack.Summary);
        Assert.DoesNotContain("residual 0.00E+000", readBack.Summary);
    }

    [Fact]
    public void NonXY_DeclinesHonestly_NoLenses()
    {
        var s = new Symphony(n: 4, hType: HamiltonianType.Heisenberg, defectBond: 1);
        var m = Movement(s);
        Assert.False(m.HasLenses);
        Assert.Contains("needs the XY chain", m.Summary);
        Assert.Empty(m.Children);   // silent: no lenses
    }

    [Fact]
    public void NoDefectFlag_MovementAbsent_SymphonyUnchanged()
    {
        var s = new Symphony(n: 4, hType: HamiltonianType.XY);   // no defectBond
        Assert.Null(s.DefectBond);
        Assert.DoesNotContain(Children(s), c => c is PaintersMovement);
        Assert.DoesNotContain(Children(s), c => c.DisplayName == "movement: painters");
        // The first movement's children stay exactly the clean-run lenses (no painters).
        var labels = Children(s).Select(c => c.DisplayName).ToList();
        Assert.Equal(new[] { "score", "lens: palindrome", "lens: quarter (CΨ)",
                             "lens: quarter (local CΨ)", "lens: dose (K)", "lens: light", "clock", "events" }, labels);
    }

    [Fact]
    public void Movement_AppearsBeforeEvents_WhenDefectGiven()
    {
        var labels = Children(Canonical()).Select(c => c.DisplayName).ToList();
        Assert.Contains("movement: painters", labels);
        Assert.True(labels.IndexOf("movement: painters") < labels.IndexOf("events"),
            "the painters' movement should precede the events axis");
    }

    [Fact]
    public void Clock_TaktGap_IsTwoGamma_ForDephasingChain()
    {
        // The Takt floor is the slowest nonzero decay rate = 2γ for a Z-dephasing chain.
        var m = Movement(Canonical());
        var clock = m.Children.Single(c => c.DisplayName == "clock");
        Assert.Contains("Takt gap", clock.Summary);
        // 2γ = 0.1 at γ=0.05; read it off the typed scalar child.
        var gap = clock.Children
            .Select(c => c.Payload).OfType<InspectablePayload.Real>()
            .Single(r => r.Label == "Takt gap").Value;
        Assert.Equal(0.1, gap, 6);
    }

    [Fact]
    public void ClosureEvent_MergedIntoEventsAxis_AtWindowEnd()
    {
        var s = Canonical();
        var events = Children(s).Single(c => c.DisplayName == "events");
        Assert.Contains(events.Children, c => c.Summary.Contains("[painters]") && c.Summary.Contains("closure"));
    }

    [Fact]
    public void ChiralMirrorEvent_AbsentWhenExact()
    {
        // The chiral-mirror deviation is added to the events axis only if BROKEN; canonical is EXACT.
        var s = Canonical();
        var events = Children(s).Single(c => c.DisplayName == "events");
        Assert.DoesNotContain(events.Children, c => c.Summary.Contains("chiral mirror BROKEN"));
    }

    [Fact]
    public void FirstMovement_StillBuiltOnce_WithMovementPresent()
    {
        // The painters' movement must NOT trigger a second evolution of the first movement's ρ(t).
        var s = Canonical();
        _ = s.Summary;
        _ = Children(s).Select(c => c.Summary).ToList();
        var m = Movement(s);
        _ = m.Summary;
        Assert.Equal(1, s.EvolveCount);
    }

    [Fact]
    public void Alphas_MatchPythonTwin_AndAreMirrorSymmetric()
    {
        // Cross-language anchor: Python perspectives_panel on the identical case gives
        // α = [0.98124, 0.98551, 0.98551, 0.98124] (N=4 bonding, bond 1, δJ=0.02, γ=0.05).
        // The α mirror symmetry (site i ↔ N−1−i) is the K₁ chiral law showing up in the
        // painter rates for the reflection-symmetric middle-bond defect.
        var m = Movement(Canonical(4));
        var a = m.Alphas;
        Assert.Equal(4, a.Count);
        Assert.Equal(0.9812, a[0], 3);
        Assert.Equal(0.9855, a[1], 3);
        Assert.Equal(a[0], a[3], 6);
        Assert.Equal(a[1], a[2], 6);
        Assert.All(a, x => Assert.InRange(x, 0.9, 1.1));
    }

    [Fact]
    public void DeviationProfile_SignedFiniteRightLength_AndResponseIsProfileOverDeltaJ()
    {
        var m = Movement(Canonical(5));
        Assert.True(m.HasLenses);

        var g = m.DeviationProfile;
        Assert.Equal(5, g.Count);
        Assert.All(g, x => Assert.True(double.IsFinite(x), $"deviation entry {x} must be finite"));
        Assert.Contains(g, x => x != 0.0);   // signed, non-trivial

        var resp = m.DeviationResponse;
        Assert.Equal(5, resp.Count);
        for (int i = 0; i < 5; i++)
            Assert.Equal(g[i] / m.DeltaJ, resp[i], 12);   // response == profile / δJ, site-wise
    }

    [Fact]
    public void DeviationProfile_EmptyWhenMovementSilent()
    {
        // A non-XY chain declines (no lenses); the deviation accessors mirror F/Alphas and return empty.
        var s = new Symphony(n: 4, hType: HamiltonianType.Heisenberg, defectBond: 1);
        var m = Movement(s);
        Assert.False(m.HasLenses);
        Assert.Empty(m.DeviationProfile);
        Assert.Empty(m.DeviationResponse);
    }
}
