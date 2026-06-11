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
    // initial state (the delocalized k=1 sine mode, F67). Localized (|10…0⟩) or multi-sector
    // (Bell+) states break the rescaling picture and the guard rightly flags every site
    // unreliable; that discovery is pinned in StateClass_Matters_GuardRefusesBellPair below.
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
    public void StateClass_Matters_GuardRefusesBellPair()
    {
        // Bell+ is a two-sector state, outside the PTF class: the rescaling picture breaks,
        // α is δJ-independent junk, and the guard must refuse every site (f_guard ≈ 2f).
        var s = new Symphony(n: 4, j: 1.0, gamma: 0.05, hType: HamiltonianType.XY,
                             initialState: InitialStateKind.BellPair, defectBond: 1, deltaJ: 0.02);
        var m = Movement(s);
        Assert.Equal(0, m.Reliable.Count(r => r));
        Assert.False(m.ClosureInWindow);
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
        // The first movement's children stay exactly the original six.
        var labels = Children(s).Select(c => c.DisplayName).ToList();
        Assert.Equal(new[] { "score", "lens: palindrome", "lens: quarter (CΨ)",
                             "lens: dose (K)", "lens: light", "events" }, labels);
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
}
