using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the walk-time step (experiments/COUPLING_DEFECT_WALK_TIME_STEP.md, 2026-07-12/13):
// a single defect bond J' = J(1+delta) edits exactly one summand of the front's walk-time. The arrival-time
// delay profile is a step: zero upstream of the bond, a flat plateau downstream near -delta/(2J) (first
// order; the broadband seed undershoots, the committed plateau numbers are pinned here). The dose the front
// pays is amplitude, not schedule: at gamma > 0 the timing stays ballistic and the step survives.
// Arrival time: first t with P_i(t) >= theta * max_t P_i(t) on the defect-free run at the SAME gamma.
public class WalkTimeTests
{
    static readonly World W = new();

    // one clean N=60 reference, shared by the gamma=0 tests (the expensive fixture).
    static readonly Lazy<double[][]> Clean60 = new(() => WalkTime.Reference(W, 60, 1.0, 0.0, seedSite: 0, dt: 0.02, tMax: 36.0));

    // P1, delta = +0.10: a stronger bond advances the front by the bond's walk-time, every site beyond it.
    // Committed plateau (sites 40-55): -0.0424 +/- 0.0004; upstream side stays at ~0.0012 or below.
    [Fact]
    public void Stronger_Bond_Advances_The_Front_By_One_Step()
    {
        var dt = WalkTime.Profile(W, 60, 1.0, 0.0, defect: (29, 30), delta: +0.10,
                                  clean: Clean60.Value, seedSite: 0, dt: 0.02, tMax: 36.0);
        double plateau = WalkTime.Mean(dt, 40, 55);
        Assert.InRange(plateau, -0.0464, -0.0384);            // the committed -0.0424, grid tolerance
        for (int i = 5; i <= 25; i++)                          // upstream of bond (29,30): no step
            Assert.True(Math.Abs(dt[i]) < 0.005, $"upstream site {i} moved: {dt[i]}");
    }

    // P1, delta = -0.10: a weaker bond delays the front (the sign of the metric step).
    // Committed plateau (sites 40-55): +0.0633 +/- 0.0017.
    [Fact]
    public void Weaker_Bond_Delays_The_Front()
    {
        var dt = WalkTime.Profile(W, 60, 1.0, 0.0, defect: (29, 30), delta: -0.10,
                                  clean: Clean60.Value, seedSite: 0, dt: 0.02, tMax: 36.0);
        double plateau = WalkTime.Mean(dt, 40, 55);
        Assert.InRange(plateau, 0.0563, 0.0703);              // the committed +0.0633, grid tolerance
    }

    // P4', gamma = 0.05: the watching eats the front's amplitude, not its schedule. The near-field step
    // survives dephasing: N=20, defect (4,5), sites 6-9 read -0.0488 to -0.0491 vs -0.05 (committed).
    [Fact]
    public void The_Watching_Pays_Amplitude_Not_Schedule()
    {
        var clean = WalkTime.Reference(W, 20, 1.0, 0.05, seedSite: 0, dt: 0.02, tMax: 12.0);
        var dt = WalkTime.Profile(W, 20, 1.0, 0.05, defect: (4, 5), delta: +0.10,
                                  clean: clean, seedSite: 0, dt: 0.02, tMax: 12.0);
        for (int i = 6; i <= 9; i++)
            Assert.InRange(dt[i], -0.0510, -0.0468);          // the committed -0.0488..-0.0491
    }

    // the knob itself: SetBond overwrites one bond's coupling, symmetrically, and only that bond.
    [Fact]
    public void SetBond_Edits_One_Bond_Only()
    {
        var cone = new Cone(W, 6, 1.0, 0.0);
        cone.SetBond(2, 3, 1.10);
        cone.Seed(0);
        cone.Step(0.05);
        Assert.Equal(1.0, cone.Structure, 9);                 // still trace-preserving with the defect
    }
}
