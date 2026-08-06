using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the bridged lattice of worlds (2026-07-16, the engine beat offered by the
// fourth play in Restless.cs): the Klein V4 of watchings {e, L, R, LR} run DYNAMICALLY as four
// worlds. e = the normal world rho(t); L = X^N * rho(t) (the ket-side complement reading);
// R = rho(t) * X^N (the bra-side reading, the 2026-07-03 anti-world); LR = X^N rho(t) X^N (the
// conjugation, the double turn that stays home). e and LR run the normal rule (rate -2*gamma*k),
// L and R run the turned rule (rate -2*gamma*(N-k), agreement watched). Every edge of the V4 is
// an exact always-open identity IN EXACT ARITHMETIC: the three read-through bridges, the composition closure
// L o R = LR, and the dagger pairing WL = WR-dagger (the antilinear triangle's dagger exchanging
// the two one-sided readings). The conservation law is not broken by the turn, it MOVES: the
// carried unit sits on the trace at e/LR and on the anti-trace at L/R, and the immortal set moves
// with it (diagonal at e/LR, anti-diagonal at L/R where the turned rate -2*gamma*(N-k) is 0).
// In FLOAT the four worlds are integrated separately and the edges are exact only where the
// contraction has too few terms to reorder; see the BridgeBound note below, which is where the
// distinction is worked out and measured.
public class LatticeTests
{
    static readonly World W = new();

    const double Eps = 2.220446049250313e-16;

    // Where the V4 edges are exact, and what the bound below is NOT.
    //
    // In exact arithmetic every edge is an identity: [H, X^N] = 0 and k(~i,j) = N - k(i,j). The four
    // vertices are integrated SEPARATELY, though, and the complement relabelling reverses the order
    // of the contraction Restless.Rhs accumulates over ascending m (for the L bridge it is the hx
    // half that reverses, for R the xh half, for LR both; read out of BuildHandshake, not measured,
    // since h is private to Restless). Float addition is not associative, so a rounding can survive
    // wherever the row sums have enough nonzero terms to reorder.
    //
    // What is solid, and it is the only thing asserted as exact: where the row sums are too sparse
    // to reorder, a single-excitation seed at zz = 0, every edge is exactly 0.0. That has held at
    // every (gamma, J, N) tried, wide significands included, and it is compared with == 0.0.
    //
    // BridgeBound is a REGRESSION PIN at the settings the tests below actually run, and nothing
    // more. It is not a law and it has no derivation. Three different explanations for how the
    // residual behaves in time were written here and measured false in turn (dissipation contracting
    // the worlds together; an accumulation window closing at ln(1/eps)/(gap*dt); a growth exponent
    // above one half). What is measured is only this: at the settings below the worst edge is 2.0
    // eps, widening dt to 0.1 and ticks to 400 reaches 3.0 eps, and the pin sits at 32 eps. Move off
    // those settings and it does not hold: at N=4, seed 0b0011, zz=1, J=1 the worst edge reaches
    // 54.9 eps at gamma=0.3 with dt=0.05 and 16000 ticks, 61.3 eps at gamma=0.5 with dt=0.04, and
    // 130.9 eps at dt=0.075. So do NOT read 32 eps as a bound on the object; read it as "these runs
    // produce this, and a change that moves it is worth looking at". The exact route that would end
    // the whole question -- making the contraction order-canonical under the complement, or
    // contracting the turned world against the complemented index directly -- is an engine change
    // to Restless.Rhs and a PREDICTION, not a measurement. It is recorded in
    // compute/MirrorWorld/README.md so it does not live only in a test file.
    const double BridgeBound = 32.0 * Eps;

    // every bridge of the V4 is exact at every probed tick, and the four vertices are four
    // different matrices (the bridges are not vacuous). Exact rather than bounded because the
    // seed is a single excitation at zz = 0, the sparse regime where the contraction cannot
    // reorder; that is a property of THIS configuration, not of the identity in general.
    //
    // Read the separation guard for what it is. At a single-excitation seed the four worlds are
    // seeded into DISJOINT joint-popcount blocks and Restless.Alive keeps them there, so at every
    // cell at most one of any two worlds is nonzero and the pairwise max reduces to the larger of
    // the two peak magnitudes; the min over pairs is then the second-smallest of the four peaks.
    // The bridges asserted above make all four peaks EQUAL (l, r and lr are exact entry
    // relabellings of e), so here VertexSeparation is just max|e|, and this guard and the
    // liveness gate at the end of the file are two bounds on that one scalar.
    //
    // So it is a measurement of the run, but of the amplitude's DECAY and not of distinctness.
    // What keeps it off 0 is conservation rather than the seeding: trace(e) = 1 over a block of
    // C(N,1) non-negative diagonal cells forces max|e| >= 1/N, which is >= 1/6 at every N used in
    // this file and cannot approach the 0.1 below.
    [Fact]
    public void All_Bridges_Are_Exact_And_The_Vertices_Are_Distinct()
    {
        var lattice = new Lattice(W, 4);
        var rep = lattice.Run(seed: 0b0001, dt: 0.05, ticks: 30);
        Assert.True(rep.WorstBridgeL == 0.0, $"L bridge {rep.WorstBridgeL:E3}");
        Assert.True(rep.WorstBridgeR == 0.0, $"R bridge {rep.WorstBridgeR:E3}");
        Assert.True(rep.WorstBridgeLR == 0.0, $"LR bridge {rep.WorstBridgeLR:E3}");
        Assert.True(rep.VertexSeparation > 0.1, $"vertices must differ O(1), separation {rep.VertexSeparation:E1}");
    }

    // the V4 composes at the rho level: reading the L world through the bra complement gives the
    // LR world (L o R = LR), and the same from the R side (R o L = LR).
    [Fact]
    public void The_Composition_Closes_L_Then_R_Is_LR()
    {
        var rep = new Lattice(W, 4).Run(seed: 0b0001, dt: 0.05, ticks: 30);
        Assert.True(rep.WorstComposition == 0.0, $"composition closure {rep.WorstComposition:E3}");
    }

    // the dagger exchanges the two one-sided readings: WL(t) = WR(t)-dagger for all t (the
    // antilinear triangle's dagger at the lattice level; both worlds run the SAME turned rule).
    [Fact]
    public void The_Dagger_Pairs_The_Two_OneSided_Vertices()
    {
        var rep = new Lattice(W, 4).Run(seed: 0b0001, dt: 0.05, ticks: 30);
        Assert.True(rep.WorstDagger == 0.0, $"dagger pairing {rep.WorstDagger:E3}");
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

    // The bridges survive the ZZ bond: X^N flips the sign of every Z, so Z_a Z_b is invariant and
    // [H, X^N] = 0 still holds with the longitudinal term.
    //
    // The zz axis is swept from 0 up, and the seed density with it, because the two interact: the
    // longitudinal term reorders nothing by itself (h stays complement-symmetric), it adds a
    // diagonal entry to each row of the contraction, which is one more term the reordering can bite
    // on. So a single excitation is exact at zz = 0 and rounds at zz != 0, while at N = 5 it stays
    // exact at every zz tried. Holding the seed fixed here would have hidden that.
    //
    // dt is tied to the ZZ energy scale rather than left at 0.05, and the reason is a measurement
    // worth keeping: at N = 4, seed 0b011, zz = 16, dt = 0.05 the worst edge is 4e-7, about 1.8e9
    // eps, and at N = 5 the same cell reaches 1e17. That is not a broken bridge, it is RK4 outside
    // its stability region. The commutator's eigenvalues are imaginary and span the H spread of the
    // seed's popcount block, and RK4 is stable on the imaginary axis only up to dt*spread = 2*sqrt2;
    // at N=4, p=2, zz=16 the spread is 64.2, so dt=0.05 gives 3.21 and it diverges. The same cell at
    // dt = 0.01 gives 0.035 eps. TWO cells out of 30 did it, which is why the single-excitation
    // version of this test never saw it: p=1 keeps the spread under the limit. The bound below is a
    // statement about rounding and means nothing once the integrator is unstable, so the step is
    // chosen to keep dt*zz*(N-1) under 0.5, comfortably inside the RK4 limit.
    [Fact]
    public void The_Lattice_Holds_With_The_ZZ_Bond()
    {
        foreach (int n in new[] { 3, 4, 5 })
            foreach (int seed in new[] { 0b001, 0b011 })
                foreach (double zz in new[] { 0.0, 0.25, 1.0, 4.0, 16.0 })
                {
                    double dt = Math.Min(0.05, 0.5 / (zz * (n - 1) + 1.0));
                    var rep = new Lattice(W, n, zz: zz).Run(seed, dt, ticks: 30);
                    double w = Math.Max(Math.Max(rep.WorstBridgeL, rep.WorstBridgeR),
                        Math.Max(rep.WorstBridgeLR, Math.Max(rep.WorstComposition, rep.WorstDagger)));
                    Assert.True(w <= BridgeBound, $"ZZ bridge {w / Eps:F2} eps at N={n}, seed={seed}, zz={zz}, dt={dt}");
                }
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

    // The axis that governs is the SEED, i.e. how many nonzero terms the contraction's row sums
    // hold. gamma and J were swept too and never turned the effect on or off in anything measured,
    // but that is an observation and not a proof: gamma scales only the dephasing mask while J
    // scales only the hopping part of h, so the two knobs do move the relative sizes inside the sum
    // that gets rounded, and neither is a power of two at the values used here.
    [Fact]
    public void The_Sparse_Case_Is_Exact_And_The_Dense_Case_Rounds()
    {
        // the sparse end: exactly 0.0, so the residual below is a reordering artifact and not a
        // break in the identity. This is the break-input the bound needs to mean anything.
        foreach (int n in new[] { 3, 4, 5, 6 })
        {
            var s = new Lattice(W, n).Run(seed: 0b1, dt: 0.05, ticks: 30);
            double worstSparse = Math.Max(Math.Max(s.WorstBridgeL, s.WorstBridgeR),
                Math.Max(s.WorstBridgeLR, Math.Max(s.WorstComposition, s.WorstDagger)));
            Assert.True(worstSparse == 0.0,
                $"a single-excitation seed at zz=0 has too few terms to reorder, so every edge must "
                + $"be exact; worst edge {worstSparse:E3} at N={n}");
        }

        // and the dense end really does round, so the sparse rows above are not the whole story
        // Note for whoever takes the exact route named in the BridgeBound block above: this assertion pins that a rounding
        // error EXISTS, so making the contraction order-canonical will fail it by design. That is
        // the intended signal, not a regression; delete this row together with the bound when the
        // engine change lands, and the sparse rows above become the whole test. It has a second
        // failure mode with nothing to do with that: a target where the JIT contracts the multiply
        // and add into an FMA could make this exactly 0.0 with no code change at all.
        var dense = new Lattice(W, 4).Run(seed: 0b0011, dt: 0.05, ticks: 30);
        Assert.True(dense.WorstBridgeL > 0.0,
            "two excitations were measured to round and did not; the break-input is gone");

        // the grid: 320 cells, 304 distinct (the N=3 seed list repeats 7). Worst here 2.0 eps;
        // widening dt to 0.1 and ticks to 400 outside this test reaches 3.0 eps, the largest seen
        // inside the bound's domain and what the 32 keeps headroom over.
        //
        // 146 of these 320 cells return exactly 0.0, and 80 of them are the full seed (1<<n)-1,
        // which is structurally incapable of a residual: its occupied block is the single cell
        // (N,N), no hop can leave it, and the mask is 0 on both sides, so those cells integrate a
        // 1x1 constant. They are kept because the identity should hold there too and the cost is
        // small, but they are not what the bound is measured on; the residual comes from the
        // half-filled seeds.
        //
        // No VertexSeparation guard in this loop, and that is a fact about the lattice rather than
        // an omission, though the reason is RELAXATION and not the seed alone. Half filling is
        // necessary, since the seed and its complement then carry the same popcount and the four
        // worlds land in the SAME joint-popcount block instead of four disjoint ones, so their
        // supports overlap and the separation can fall. It is not sufficient: N=4, seed 0b0011,
        // zz=0, gamma=0.5, J=1 separates by 0.43 at 30 ticks and only falls to 7.6e-5 by 200, as
        // both worlds approach a common steady state. 17 of the 320 cells end under 0.1, the
        // smallest being seed 0b0101 at 200 ticks (2.5e-5). Away from half filling the supports
        // stay disjoint, the separation reduces to max|e| (see the note on the first test), and
        // conservation floors that at 1/N, so the guard could not fail there whatever the run
        // does. It is never the thing keeping a row honest. What keeps THIS test non-vacuous is the exactly-0.0
        // sparse rows at its top and The_Turned_Rule_Is_LoadBearing_For_The_OneSided_Vertices,
        // which breaks a bridge at O(1).
        double worst = 0.0;
        string where = "nothing ran";
        foreach (int n in new[] { 3, 4, 5, 6 })
            foreach (int seed in new[] { 1, 3, 5, 7, (1 << n) - 1 })
                foreach (double zz in new[] { 0.0, 1.0 })
                    foreach (var (g, j) in new[] { (0.5, 1.0), (0.05, 0.075), (0.07, 2.3), (1.0 / 3.0, 0.5) })
                        foreach (int ticks in new[] { 30, 200 })
                        {
                            var rep = new Lattice(W, n, j: j, gamma: g, zz: zz).Run(seed, 0.05, ticks);
                            double w = Math.Max(Math.Max(rep.WorstBridgeL, rep.WorstBridgeR),
                                Math.Max(rep.WorstBridgeLR, Math.Max(rep.WorstComposition, rep.WorstDagger)));
                            if (w > worst) { worst = w; where = $"N={n} seed={seed} zz={zz} gamma={g} J={j} ticks={ticks}"; }
                        }

        Assert.True(worst <= BridgeBound, $"worst V4 edge residual {worst / Eps:F2} eps at {where}");
    }

    // The one time-axis statement that survived measurement: without dephasing the residual keeps
    // growing, so the small numbers above are not a property of the identity but of short runs at a
    // dissipative point. Strict growth rather than a threshold, and robust across N, seed, zz, dt
    // and J (the only exception is the full seed, which is structurally zero, see the sweep above).
    //
    // Deliberately NOT asserted here: that the growth stops at a physical gamma. It does at
    // gamma=0.5 with dt=0.05, where 4000 and 16000 ticks agree bit for bit, and it does not at
    // gamma=0.3 or at dt=0.04, where the same cell reaches 55 and 61 eps. Whatever controls that is
    // not understood, and an assertion built on the one cell where it looked clean would be
    // measuring the cell.
    [Fact]
    public void Without_Dephasing_The_Rounding_Keeps_Accumulating()
    {
        double Worst(int ticks)
        {
            var r = new Lattice(W, 4, j: 1.0, gamma: 0.0, zz: 1.0).Run(0b0011, 0.05, ticks);
            return Math.Max(Math.Max(r.WorstBridgeL, r.WorstBridgeR),
                   Math.Max(r.WorstBridgeLR, Math.Max(r.WorstComposition, r.WorstDagger)));
        }
        double shortRun = Worst(1000), longRun = Worst(16000);
        Assert.True(longRun > shortRun,
            $"at gamma=0 nothing damps the accumulated rounding, so the running maximum must keep "
            + $"growing: {shortRun / Eps:F2} eps at 1000 ticks, {longRun / Eps:F2} eps at 16000");
    }

    // Liveness. A review stubbed Restless.Step to a no-op and SIX of the ten tests then in this
    // file still passed, because most of what they assert is a relabelling identity between four
    // matrices and a relabelling holds just as well when nothing moves. The VertexSeparation guard
    // in the first test does not catch it either: the four worlds are seeded into four DISJOINT
    // joint-popcount blocks (Restless.Alive keeps each one there forever), so no two of them ever
    // share a nonzero cell, and the separation is already 1.0 before a single step.
    //
    // The frozen reading is exact, which is why this compares with ==, and distinct seed CELLS are
    // all that takes; the disjoint blocks above are stronger than needed for it.
    // The claim is NOT that no arithmetic has run at ticks = 0 -- the bridge, unit and separation
    // loops all run at tick 0. It is that no STEP has run (Lattice.Run breaks before the first
    // Step when tick == ticks), so every cell is either 0 or the seeded 1.0, and every pairwise
    // difference is (1,0) - (0,0), whose Magnitude is returned unmodified. Nothing rounds.
    //
    // What this gates, measured by stubbing rather than argued: a lattice that stops moving as a
    // whole. Under a no-op Step the live separation is exactly 1.0 too, so the discriminator is
    // the exact `< 1.0` and no threshold is needed. What it does NOT gate, also measured, with
    // who does gate it, so nobody reads this as the file's liveness guarantee:
    //   - the dissipator. Dropping the mask term from Restless.Rhs leaves this test green (all
    //     four seed cells sit on the zero-rate set: e and LR are diagonal at k=0, L and R are at
    //     k=N where the turned rate -2*gamma*(N-k) vanishes), so the drop below is driven by the
    //     hopping alone. The_Opening_Law_Holds_On_The_Cat_Pair and
    //     The_Turned_Rule_Is_LoadBearing_For_The_OneSided_Vertices both fail on that mutation.
    //   - one or two of the four worlds freezing. A frozen world keeps the LARGEST peak, so the
    //     min-over-pairs lands on a live pair and this test stays green. Freezing L and R fails
    //     four of the bridge tests instead, which is where a broken vertex belongs.
    // For the record and not asserted: the live separation here is 0.438077.
    [Fact]
    public void The_Integrator_Actually_Moves_The_Worlds()
    {
        var lattice = new Lattice(W, 4);
        var frozen = lattice.Run(seed: 0b0001, dt: 0.05, ticks: 0);
        var live = lattice.Run(seed: 0b0001, dt: 0.05, ticks: 30);

        Assert.True(frozen.VertexSeparation == 1.0,
            $"before any step the four vertices are four distinct basis cells of amplitude 1 in "
            + $"disjoint blocks, so the separation is exactly 1.0; got {frozen.VertexSeparation:E17}");
        Assert.True(live.VertexSeparation < 1.0,
            $"30 ticks must move the worlds, and a stopped integrator leaves the separation at "
            + $"exactly 1.0; got {live.VertexSeparation:E17}");
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
