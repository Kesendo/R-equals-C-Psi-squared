using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the walk-time step (experiments/COUPLING_DEFECT_WALK_TIME_STEP.md, 2026-07-12/13):
// a single defect bond J' = J(1+delta) edits exactly one summand of the front's walk-time. The arrival-time
// delay profile is a step: zero upstream of the bond, a flat plateau downstream near -delta/(2J) (first
// order; the broadband seed undershoots, the committed plateau numbers are pinned here). The dose the front
// pays at the tested gamma=0.05, J=1, N=60 point is mainly amplitude: the timing step survives there.
// Beyond the arrival-dose fence the front crosses into a diffusive regime.
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

    // Two disconnected pairs let the edited pair move faster while the untouched pair keeps J=1.
    [Fact]
    public void SetBond_Edits_One_Bond_Only()
    {
        static Cone Run(int seed)
        {
            var cone = new Cone(W, 4, 1.0, 0.0, bonds: new[] { (0, 1), (2, 3) });
            cone.SetBond(2, 3, 1.10);
            cone.Seed(seed);
            cone.Step(0.01);
            return cone;
        }

        var untouched = Run(0);
        var edited = Run(2);
        Assert.Equal(Math.Pow(Math.Sin(0.01), 2), untouched.Population(1), 8);
        Assert.Equal(Math.Pow(Math.Sin(0.011), 2), edited.Population(3), 8);
        Assert.Equal(1.0, edited.Structure, 9);
    }

    // ---- the inverse (experiments/THE_DECODER_AND_THE_FRONT.md, run mode `walk decode`) ----
    // Seed site 0, J = 1, gamma = 0.05 (Q = 20, the DefectDecoder's protocol point), letters by the central
    // difference at +-0.02, dt = 0.005, tMax = 0.75 N + 3.
    static double TMax(int n) => 0.75 * n + 3.0;
    static readonly Dictionary<int, double[][]> LettersAt = new();
    static double[][] Letters(int n)
    {
        lock (LettersAt)
        {
            if (!LettersAt.TryGetValue(n, out var l))
                LettersAt[n] = l = WalkTime.Letters(W, n, 1.0, 0.05, seedSite: 0, dt: 0.005, tMax: TMax(n));
            return l;
        }
    }

    // The arrival letters are nested steps of full rank. A pin of the committed run at theta = 0.2, N = 3..10: the
    // worst signed cosine reads +0.824 down to -0.167 (the DefectDecoder's painted letters on the bonding mode read
    // -0.91..-0.99 at N = 3, 5..12), and the normalised smallest singular value 0.419 down to 0.158. Past N = 10 the
    // worst cosine keeps drifting down (about -0.54 at N = 20, simulations/decoder_carrier_dictionary.py), so this
    // gate is the window N <= 10, not a law. The N = 5 letters are pinned to the committed run (grid tolerance
    // 0.002): about -1/(2J) beyond the bond, 0 before it, the bump on the site just before the bond.
    [Fact]
    public void The_Front_Letters_Have_No_AntiCollinear_Pair()
    {
        for (int n = 3; n <= 10; n++)
        {
            Assert.True(WalkTime.WorstSignedCos(Letters(n)).Cos > -0.25, $"N={n}: an anti-collinear pair");
            Assert.True(WalkTime.Conditioning(Letters(n)) > 0.15, $"N={n}: a nearly degenerate dictionary");
        }
        double[][] committed =
        {
            new[] { -0.347, -0.458, -0.545, -0.528 },
            new[] { +0.013, -0.460, -0.554, -0.547 },
            new[] { +0.000, +0.062, -0.554, -0.547 },
            new[] { +0.000, -0.002, +0.163, -0.528 },
        };
        var l5 = Letters(5);
        for (int b = 0; b < 4; b++)
            for (int i = 1; i < 5; i++)
                Assert.InRange(l5[b][i], committed[b][i - 1] - 0.002, committed[b][i - 1] + 0.002);
    }

    // The DefectDecoder's confuser (N = 5: bond 3 weakened reads like bond 1 strengthened, the alpha path's residual
    // ratio 1.45) is read cleanly from the front, bond and sign, by the decoder's own single-letter rule.
    [Fact]
    public void The_Front_Reads_The_Decoders_Confuser()
    {
        var clean = WalkTime.Reference(W, 5, 1.0, 0.05, seedSite: 0, dt: 0.005, tMax: TMax(5));
        foreach (var (bond, delta) in new[] { (3, -0.02), (1, +0.02) })
        {
            var deltas = new double[4]; deltas[bond] = delta;
            var r = WalkTime.DecodeOne(Letters(5), WalkTime.Profile(W, 5, 1.0, 0.05, deltas, clean, 0, 0.005, TMax(5)));
            Assert.Equal(bond, r.Bond);
            Assert.InRange(r.Strength / delta, 0.95, 1.05);
            Assert.True(r.Ratio > 1000, $"bond {bond}: residual ratio {r.Ratio}");
        }
    }

    // Every bond at once, by least squares over all letters (the dictionary has full rank N - 1). The error is
    // the second order of the step law, so the gate is that law:
    // halving every delta quarters the largest error (ratio near 1/4), at N = 5 and N = 8.
    [Fact]
    public void The_Front_Reads_Every_Bond_At_Once()
    {
        var truths = new Dictionary<int, double[]>
        {
            [5] = new[] { +0.0076, -0.0108, -0.0035, -0.0129 },
            [8] = new[] { +0.0102, -0.0070, +0.0048, -0.0095, -0.0136, -0.0165, +0.0181 },
        };
        foreach (var (n, truth) in truths)
        {
            var clean = WalkTime.Reference(W, n, 1.0, 0.05, seedSite: 0, dt: 0.005, tMax: TMax(n));
            double Err(double scale)
            {
                var t = truth.Select(v => v * scale).ToArray();
                var x = WalkTime.Decode(Letters(n), WalkTime.Profile(W, n, 1.0, 0.05, t, clean, 0, 0.005, TMax(n)));
                return t.Zip(x, (a, b) => Math.Abs(a - b)).Max();
            }
            double full = Err(1.0), half = Err(0.5);
            Assert.True(full < 0.05 * truth.Max(Math.Abs), $"N={n}: error {full}");
            Assert.InRange(half / full, 0.18, 0.32);
        }
    }

    // With the observable held fixed (the whole population trajectory, every site, t = 0.25..20) the carrier decides.
    // Pins of the committed run at N = 5..9: the eigenmodes read anti-collinear (psi1 -0.882..-0.974, psi2
    // -0.897..-0.990) with a normalised smallest singular value of 0.06..0.11; the site seed's dictionary is at least
    // 3.5 times better conditioned than the bonding mode's (0.42..0.65 against 0.07..0.11); the end mixture, which
    // the mirror leaves fixed, reads well conditioned too (0.36..0.63), while (psi1 + psi2)/sqrt2, which it does not,
    // reads anti-collinear (-0.748..-0.799): mirror parity is not what separates them.
    [Fact]
    public void The_Carrier_Decides_On_The_Trajectory()
    {
        for (int n = 5; n <= 9; n++)
        {
            double[] Site(int s) { var v = new double[n]; v[s] = 1.0; return v; }
            var p12 = WalkTime.Mode(n, 1).Zip(WalkTime.Mode(n, 2), (a, b) => (a + b) / Math.Sqrt(2)).ToArray();
            double[][] Read(params (double, double[])[] rho) =>
                WalkTime.TrajectoryLetters(W, n, 1.0, 0.05, rho, dt: 0.01, sample: 0.25, tMax: 20.0);
            var bonding = Read((1.0, WalkTime.Mode(n, 1)));
            var seed = Read((1.0, Site(0)));
            Assert.True(WalkTime.WorstSignedCos(bonding).Cos < -0.85, $"N={n}: bonding mode");
            Assert.True(WalkTime.WorstSignedCos(Read((1.0, WalkTime.Mode(n, 2)))).Cos < -0.85, $"N={n}: psi2");
            Assert.True(WalkTime.Conditioning(seed) > 3.5 * WalkTime.Conditioning(bonding), $"N={n}: site seed against bonding mode");
            Assert.True(WalkTime.Conditioning(Read((0.5, Site(0)), (0.5, Site(n - 1)))) > 0.3, $"N={n}: end mixture");
            Assert.True(WalkTime.WorstSignedCos(Read((1.0, p12))).Cos < -0.7, $"N={n}: psi1 + psi2");
        }
    }

    // The conditioning reading itself, against two cases with a known answer: orthogonal letters read 1, and a
    // letter that is the sum of two others makes the dictionary singular (0 to rounding).
    [Fact]
    public void Conditioning_Reads_Orthogonal_As_One_And_Dependent_As_Zero()
    {
        Assert.Equal(1.0, WalkTime.Conditioning(new[] { new[] { 1.0, 0, 0 }, new[] { 0, 2.0, 0 }, new[] { 0, 0, 3.0 } }), 12);
        Assert.True(WalkTime.Conditioning(new[] { new[] { 1.0, 0, 1 }, new[] { 0, 1.0, 1 }, new[] { 1.0, 1, 2 } }) < 1e-7);
        Assert.Throws<ArgumentException>(() => WalkTime.WorstSignedCos(new[] { new[] { 1.0, double.NaN }, new[] { 0.0, 1.0 } }));
    }
}
