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

    private static Symphony Canonical() =>
        new Symphony(n: 4, j: 1.0, gamma: 0.1, hType: HamiltonianType.XY,
                     initialState: InitialStateKind.BellPair, defectBond: 1, deltaJ: 0.02);

    [Fact]
    public void Closure_WithinWindow_ForCanonicalXY()
    {
        var m = Movement(Canonical());
        Assert.True(m.HasLenses);
        Assert.True(m.ClosureInWindow,
            $"closure Σ ln α (reliable) = {m.ClosureSum} should be in ±{PaintersMovement.ClosureWindow}");
        Assert.True(Math.Abs(m.ClosureSum) <= PaintersMovement.ClosureWindow);
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
        // 2γ = 0.2 at γ=0.1; read it off the typed scalar child.
        var gap = clock.Children
            .Select(c => c.Payload).OfType<InspectablePayload.Real>()
            .Single(r => r.Label == "Takt gap").Value;
        Assert.Equal(0.2, gap, 6);
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
    public void Alphas_MatchPythonGlobalMinimum_ForReliableSites()
    {
        // Cross-validation anchor (see simulations/_ptf_symphony_crossval.py): the C# golden-section
        // fit lands the genuine global MSE minima. Sites 0, 1, 3 agree with the Python twin to ~1e-2;
        // site 2 is a featureless trajectory whose global min is α≈1.016 (the Python reference's Brent
        // got trapped at α≈3.15, MSE 920× worse). Pin the three well-conditioned sites here.
        var m = Movement(Canonical());
        var a = m.Alphas;
        Assert.Equal(4, a.Count);
        Assert.Equal(8.878, a[0], 2);
        Assert.Equal(3.679, a[1], 2);
        Assert.Equal(6.556, a[3], 2);
    }
}
