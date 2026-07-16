using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the bridged lattice of worlds (2026-07-16, the engine beat offered by the
// fourth play in Restless.cs): the Klein V4 of watchings {e, L, R, LR} run DYNAMICALLY as four
// worlds. e = the normal world rho(t); L = X^N * rho(t) (the ket-side complement reading);
// R = rho(t) * X^N (the bra-side reading, the 2026-07-03 anti-world); LR = X^N rho(t) X^N (the
// conjugation, the double turn that stays home). e and LR run the normal rule (rate -2*gamma*k),
// L and R run the turned rule (rate -2*gamma*(N-k), agreement watched). Every edge of the V4 is
// an exact always-open identity: the three read-through bridges, the composition closure
// L o R = LR, and the dagger pairing WL = WR-dagger (the antilinear triangle's dagger exchanging
// the two one-sided readings). The conservation law is not broken by the turn, it MOVES: the
// carried unit sits on the trace at e/LR and on the anti-trace at L/R, and the immortal set moves
// with it (diagonal at e/LR, anti-diagonal at L/R where the turned rate -2*gamma*(N-k) is 0).
public class LatticeTests
{
    static readonly World W = new();

    // every bridge of the V4 holds at machine precision at every probed tick, and the four
    // vertices are genuinely four different matrices (the bridges are not vacuous).
    [Fact]
    public void All_Bridges_Are_Exact_And_The_Vertices_Are_Distinct()
    {
        var lattice = new Lattice(W, 4);
        var rep = lattice.Run(seed: 0b0001, dt: 0.05, ticks: 30);
        Assert.True(rep.WorstBridgeL < 1e-12, $"L bridge {rep.WorstBridgeL:E1}");
        Assert.True(rep.WorstBridgeR < 1e-12, $"R bridge {rep.WorstBridgeR:E1}");
        Assert.True(rep.WorstBridgeLR < 1e-12, $"LR bridge {rep.WorstBridgeLR:E1}");
        Assert.True(rep.VertexSeparation > 0.1, $"vertices must differ O(1), separation {rep.VertexSeparation:E1}");
    }

    // the V4 composes at the rho level: reading the L world through the bra complement gives the
    // LR world (L o R = LR), and the same from the R side (R o L = LR).
    [Fact]
    public void The_Composition_Closes_L_Then_R_Is_LR()
    {
        var rep = new Lattice(W, 4).Run(seed: 0b0001, dt: 0.05, ticks: 30);
        Assert.True(rep.WorstComposition < 1e-12, $"composition closure {rep.WorstComposition:E1}");
    }

    // the dagger exchanges the two one-sided readings: WL(t) = WR(t)-dagger for all t (the
    // antilinear triangle's dagger at the lattice level; both worlds run the SAME turned rule).
    [Fact]
    public void The_Dagger_Pairs_The_Two_OneSided_Vertices()
    {
        var rep = new Lattice(W, 4).Run(seed: 0b0001, dt: 0.05, ticks: 30);
        Assert.True(rep.WorstDagger < 1e-12, $"dagger pairing {rep.WorstDagger:E1}");
    }

    // the carried unit: trace stays 1 at e and LR, the anti-trace stays 1 at L and R, at every
    // probed tick. the law moves to the anti-diagonal; it is never lost.
    [Fact]
    public void The_Carried_Unit_Moves_With_The_Turn()
    {
        var rep = new Lattice(W, 4).Run(seed: 0b0001, dt: 0.05, ticks: 30);
        Assert.True(rep.WorstUnitE < 1e-9, $"trace at e {rep.WorstUnitE:E1}");
        Assert.True(rep.WorstUnitLR < 1e-9, $"trace at LR {rep.WorstUnitLR:E1}");
        Assert.True(rep.WorstUnitL < 1e-9, $"anti-trace at L {rep.WorstUnitL:E1}");
        Assert.True(rep.WorstUnitR < 1e-9, $"anti-trace at R {rep.WorstUnitR:E1}");
    }

    // the watching assignment is load-bearing: running the R reading under the NORMAL rule
    // (instead of the turned one) breaks the bridge at O(1) -- the deviation is the handshake
    // break, asserted from below.
    [Fact]
    public void The_Turned_Rule_Is_LoadBearing_For_The_OneSided_Vertices()
    {
        var lattice = new Lattice(W, 3);
        double broken = lattice.BrokenBridgeR(seed: 0b001, dt: 0.05, ticks: 30);
        Assert.True(broken > 1e-3, $"the wrong watching must break the bridge O(1) from below: {broken:E1}");
    }

    // the bridges survive the ZZ bond: X^N flips the sign of every Z, so Z_a Z_b is invariant,
    // [H, X^N] = 0 still holds with the longitudinal term, and the whole lattice rides along.
    [Fact]
    public void The_Lattice_Holds_With_The_ZZ_Bond()
    {
        var lattice = new Lattice(W, 3, zz: 1.0);
        var rep = lattice.Run(seed: 0b001, dt: 0.05, ticks: 30);
        Assert.True(rep.WorstBridgeL < 1e-12, $"L bridge with ZZ {rep.WorstBridgeL:E1}");
        Assert.True(rep.WorstBridgeR < 1e-12, $"R bridge with ZZ {rep.WorstBridgeR:E1}");
        Assert.True(rep.WorstBridgeLR < 1e-12, $"LR bridge with ZZ {rep.WorstBridgeLR:E1}");
    }

    // the opening law (experiments/LATTICE_OPENING_LAW.md, found playing 2026-07-16): on the cat
    // pair psi(theta) = cos|0..0> + sin|1..1>, the entry-wise distance between the e world and its
    // one-sided reading L has the closed form
    //     opening(t) = max(cos^2, sin^2) - cos*sin * e^(-2*Gamma*t),   Gamma = N*gamma (uniform),
    // "the heavier sock's weight minus the LIVING spook": the chirality floor is timeless, the
    // spook term dies at the full k = N rate and closes exactly the gap it owns. Exact because the
    // cat sector is H-dead (no hop touches |0..0> or |1..1>), so the law is J-free too.
    [Fact]
    public void The_Opening_Law_Holds_On_The_Cat_Pair()
    {
        // dt = 0.025: at uniform gamma = 0.5 the spook rate is 2*N*gamma = 3, and dt = 0.05 puts
        // the RK4 floor at ~9e-7, a thin 1.14x margin under 1e-6; halving dt restores ~18x.
        var lattice = new Lattice(W, 3);
        foreach (double thetaDeg in new[] { 0.0, 30.0, 45.0, 75.0 })
        {
            double dev = lattice.OpeningLawDeviation(thetaDeg * Math.PI / 180.0, dt: 0.025, ticks: 80);
            Assert.True(dev < 1e-6, $"theta={thetaDeg}: worst |opening - closed form| = {dev:E1}");
        }
        // J-free: a very different hop amplitude gives the same deviation floor (the cat sector
        // is H-dead, so the whole opening trajectory is J-blind).
        var fastHop = new Lattice(W, 3, j: 2.3);
        Assert.True(fastHop.OpeningLawDeviation(Math.PI / 6, 0.025, 80) < 1e-6);
    }

    // ontology: the lattice's own outputs.
    [Fact]
    public void Ontology_Own_Carries_The_Lattice_Readings()
    {
        var lattice = new Lattice(W, 3);
        Assert.Contains("bridges", lattice.Own);
        Assert.Contains("carried-unit", lattice.Own);
    }
}
