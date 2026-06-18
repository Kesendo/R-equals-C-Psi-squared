using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class TrichotomyWitnessTests
{
    private static IReadOnlyList<double> Uniform(int n, double gamma) =>
        Enumerable.Repeat(gamma, n).ToList();

    [Fact]
    public void Ctor_RejectsOutOfRangeN() =>
        Assert.Throws<ArgumentOutOfRangeException>(() => new TrichotomyWitness(n: 9));

    [Fact] // the two reads MUST stay on DIFFERENT conventions — guards against re-merging them onto one
    public void TwoReads_StayOnDifferentConventions()  // (the Round-2 defect this whole feature fixed)
    {
        // At the same nominal q=1.5 the two reads see DIFFERENT survivors, because they are different
        // physical sweeps: the CARBON un-freeze read (ClassifyUnfreeze, J/γ=Q) sees the frozen (p,p)
        // interior below Q*(5)≈2.374; the ABSOLUTE read (SurvivorSector, what ClassifySeam uses) sees the
        // (0,1) band edge (J/γ=3, above Q*). If a future change re-merges them onto one convention, one of
        // these Δn flips and this fails — re-introducing the chain mislabel the two-read split removed.
        Assert.Equal(0, TrichotomyWitness.ClassifyUnfreeze(TopologyKind.Chain, 5, 1.5).Dn);   // carbon: interior
        var (pc, pr, _) = TrichotomyWitness.SurvivorSector(TopologyKind.Chain, 5, 1.5, Uniform(5, 0.5));
        Assert.Equal(1, Math.Abs(pc - pr));                                                  // absolute: band edge
    }

    [Fact]
    public void SurvivorRate_MatchesPostEpFlowField_AtN5_VacBlock()
    {
        // The convention pin (R2): the (0,1) survivor rate on the absolute SectorSlowest scale equals
        // the full-L PostEpFlowField slowest rate to 9 digits. Same (q, profile) convention as the
        // existing SectorReductionWitnessTests pin (VacBlock_SlowestRate_EqualsPostEpFlowField_AtN5):
        // q in {1.5, 1000.0}, tauGrid {0.0, 1.0}, ReadAssembly(q).SlowestRate as ground truth.
        const int n = 5; const double q = 1.5;
        var profile = Uniform(n, 0.5);
        var (pc, pr, rate) = TrichotomyWitness.SurvivorSector(TopologyKind.Chain, n, q, profile);
        var expected = new PostEpFlowField(n, new[] { 1.5, 1000.0 }, new[] { 0.0, 1.0 }, profile)
            .ReadAssembly(q).SlowestRate;
        Assert.Equal(expected, rate, 9);
    }

    [Fact]
    public void Survivor_Star5_IsTheOneOneCommutant()
    {
        var (pc, pr, _) = TrichotomyWitness.SurvivorSector(TopologyKind.Star, 5, 1.5, Uniform(5, 0.5));
        Assert.Equal((1, 1), (pc, pr));
    }

    [Fact]
    public void Deviation_CanalN5_IsGlobal_NotSectorPinned()
    {
        // The EXACT "canal N=5" profile from simulations/birth_canal_junction_nature.py line 281
        // (_stage2_seam cases). The global slowest mode stays the (0,1) edge-survivor (odd-drift,
        // no sector switch) — its rate-drift is 0.085, NOT the sector-pinned 0.304 (R2).
        var canal = new[] { 0.25, 1.5, 1.5, 1.5, 0.25 };
        Assert.Equal(0.0850, TrichotomyWitness.Deviation(TopologyKind.Chain, 5, canal), 3);
    }

    [Fact]
    public void Deviation_DeepEdgeN6_IsTheJunction()
    {
        // The EXACT "deep-edge N=6" profile from birth_canal_junction_nature.py line 280:
        // deep_edge(6, edge=0.25) = [0.25] + [(6-0.5)/4]*4 + [0.25] = [0.25, 1.375x4, 0.25].
        // Here the global slowest SWITCHES sectors (the junction): the interior (2,2) wins at low Q.
        var deep = new[] { 0.25, 1.375, 1.375, 1.375, 1.375, 0.25 };
        Assert.Equal(0.4079, TrichotomyWitness.Deviation(TopologyKind.Chain, 6, deep), 3);
    }

    // ====================== TASK 5: the two-read Classify ============================
    // CARBON read (Q = J/γ, uniform γ = 1/Q): ClassifyUnfreeze drives the un-freeze trichotomy.

    [Fact] // THE R1 GATE — now correct under the carbon read
    public void ClassifyUnfreeze_ChainBelowQStar_FrozenInterior_UnfreezingSeEp()
    {
        var r = TrichotomyWitness.ClassifyUnfreeze(TopologyKind.Chain, 5, 1.5); // Q=1.5 < Q*(5)≈2.374 (carbon)
        Assert.Equal(0, r.Dn);                 // frozen (p,p) interior
        Assert.True(r.ImMax < 1e-6, $"frozen below the horizon; got |Im|={r.ImMax}");
        Assert.Equal(TrichotomyWitness.Route.UnfreezingSeEp, r.Route);
    }

    [Fact]
    public void ClassifyUnfreeze_ChainAboveQStar_BandEdgeOscillates()
    {
        var r = TrichotomyWitness.ClassifyUnfreeze(TopologyKind.Chain, 5, 5.0);
        Assert.Equal(1, r.Dn);
        Assert.True(r.ImMax > 1e-2, $"band edge oscillates above Q*; got |Im|={r.ImMax}");
        Assert.Equal(TrichotomyWitness.Route.UnfreezingSeEp, r.Route);
    }

    [Fact]
    public void ClassifyUnfreeze_Star5_FrozenCommutant()
    {
        var r = TrichotomyWitness.ClassifyUnfreeze(TopologyKind.Star, 5, 8.0);
        Assert.Equal(TrichotomyWitness.Route.FrozenCommutant, r.Route);
        Assert.True(r.ImMax < 1e-6, $"star (1,1) commutant is frozen at every Q; got |Im|={r.ImMax}");
    }

    [Fact] // the star is FrozenCommutant at LOW Q too, not just at the high-Q ceiling (the tool-test catch)
    public void ClassifyUnfreeze_StarLowQ_IsFrozenCommutant_NotLevelCrossing()
    {
        // The star's frozen survivor is the [H,A]=0 commutant at EVERY Q; below the ceiling asymptote it must
        // NOT fall to the ring's FrozenLevelCrossing label (the bug the rendered tree exposed at Q≤6).
        var r = TrichotomyWitness.ClassifyUnfreeze(TopologyKind.Star, 6, 1.5);
        Assert.Equal(0, r.Dn);
        Assert.True(r.ImMax < 1e-6);
        Assert.Equal(TrichotomyWitness.Route.FrozenCommutant, r.Route);  // NOT FrozenLevelCrossing
    }

    [Fact]
    public void ClassifyUnfreeze_Ring5_FrozenLevelCrossing()
    {
        // GATE-FIRST CORRECTION (2026-06-18): the spec's Q=1.5 is ABOVE the ring N=5 interior handover —
        // at Q=1.5 the ring survivor is the (0,1) band edge (Dn=1), gap=2γ=1.333 (gate-measured). The ring
        // interior (the frozen V-Effect level crossing) only wins below Q≈1.3 (handover grows ~0.29N). At
        // Q=1.0 the survivor is the frozen (1,1) interior (slowIm=0), whose rate sits 76% below its commutant
        // ceiling (a level crossing, not a commutant ceiling) -> FrozenLevelCrossing. Diagnosed, not loosened.
        var r = TrichotomyWitness.ClassifyUnfreeze(TopologyKind.Ring, 5, 1.0);
        Assert.Equal(0, r.Dn);
        Assert.True(r.ImMax < 1e-6, $"ring interior is a frozen level crossing; got |Im|={r.ImMax}");
        Assert.Equal(TrichotomyWitness.Route.FrozenLevelCrossing, r.Route);
    }

    [Fact] // the carbon-mapping gate: the carbon block rate == the canonical Survivor.Gap
    public void CarbonBlock_RateMatchesSurvivor()
    {
        const int n = 5; const double Q = 1.5;
        var (gap, pc, pr, _) = IncompletenessSurvivorWitness.Survivor(n, Q, TopologyKind.Chain);
        // The carbon |Im|/rigidity block is built the SAME way Survivor builds its block (Qh=0.5,
        // uniform γ=1/Q). Its slowest non-kernel rate must equal Survivor.Gap to ~9 digits, else
        // the |Im| we read belongs to a different block. CarbonSlowestRate exposes that rate.
        double carbonRate = TrichotomyWitness.CarbonSlowestRate(TopologyKind.Chain, n, Q, pc, pr);
        Assert.Equal(gap, carbonRate, 9);
    }

    // ABSOLUTE read (fixed γ-profile): ClassifySeam drives the Δn-seam, reproducing the script.

    [Fact] // the seam read — absolute, matching birth_canal_junction_nature.py
    public void ClassifySeam_Uniform_Sterile_Canal_OddDrift_DeepEdge_Junction()
    {
        // uniform N=5: probe window above Q*(5), the (0,1) band edge is Q-flat -> sterile.
        Assert.Equal(TrichotomyWitness.SeamKind.Sterile,
            TrichotomyWitness.ClassifySeam(TopologyKind.Chain, 5, Uniform(5, 0.5)).Kind);
        // canal N=5: edge-protected (0,1) survivor drifts, no Δn switch -> odd-drift.
        var canal = new[] { 0.25, 1.5, 1.5, 1.5, 0.25 }; // exact from the script (line 281)
        Assert.Equal(TrichotomyWitness.SeamKind.OddDrift,
            TrichotomyWitness.ClassifySeam(TopologyKind.Chain, 5, canal).Kind);
        // deep-edge N=6: deep_edge(6,0.25)=[0.25,1.375x4,0.25]; interior (2,2) overtakes at low Q,
        // the survivor's Δn flips 0->1 -> junction.
        var deep = new[] { 0.25, 1.375, 1.375, 1.375, 1.375, 0.25 };
        Assert.Equal(TrichotomyWitness.SeamKind.Junction,
            TrichotomyWitness.ClassifySeam(TopologyKind.Chain, 6, deep).Kind);
    }

    // ====================== TASK 6: the IInspectable tree (four slices) ============================

    [Fact]
    public void Witness_RendersAndGoesPastN5()
    {
        var w = new TrichotomyWitness(n: 6, q: 1.5);
        var kids = ((IInspectable)w).Children.ToList();
        Assert.Equal(4, kids.Count);                       // RouteSweep, ThresholdLadder, ΔnSeam, Vocabulary
        foreach (var k in kids) _ = k.Children.ToList();   // force lazy enumeration; must not throw
        Assert.NotEmpty(w.DisplayName);
    }
}
