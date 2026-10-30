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

    // At equal cost (one shot = one position at one time, equal shots per time): the per-shot Fisher matrix of
    // every bond's delta. Handshake with the bond/seat study (experiments/THE_BOND_AND_THE_SEAT_WITH_FINITE_SHOTS.md):
    // N = 7, J = 0.075, gamma = 0.05 (Q = 1.5), the bonding mode on that page's baseline schedule
    // [2, 4, 6, 8, 10, 12, 15, 18, 24, 30] with 1e5 shots per time, in additive units (divided by J^2): its smallest
    // eigenvalue is 190,443 there, 190443.2487 on the exact route of simulations/decoder_carrier_dictionary.py. The
    // cone's RK4 meets it to its step error, and the gate is that error's law: halving dt cuts the deviation from
    // the exact value by about 2^4 = 16.
    [Fact]
    public void Fisher_Meets_The_Bond_Seat_Study_To_The_Step_Law()
    {
        const double exact = 190443.2487, jh = 0.075;
        double[] sched = { 2, 4, 6, 8, 10, 12, 15, 18, 24, 30 };
        double Lambda(double dt)
        {
            var f = WalkTime.PositionFisher(W, 7, jh, 0.05, new[] { (1.0, WalkTime.Mode(7, 1)) }, sched, dt);
            var total = new double[6, 6];
            for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) total[a, b] = f[a, b] * sched.Length * 1e5 / (jh * jh);
            return WalkTime.MinEigenvalue(total);
        }
        double e1 = Math.Abs(Lambda(0.5) - exact), e2 = Math.Abs(Lambda(0.25) - exact);
        Assert.True(e2 < 1.0, $"dt = 0.25 reads {exact + e2} against {exact}");
        Assert.InRange(e1 / e2, 10.0, 22.0);
    }

    // At Q = 20 on the decoder's window (60 equal-split times on (0, 1/gamma]), the localized carriers hold far more
    // information in the worst direction than the bonding mode: site 0 108, 119, 121 times and the centre site 93,
    // 76, 97 times at N = 5, 7, 9 on the exact route (part of it the uniform rescaling of every bond, the bonding
    // mode's nearly blind direction, a shrinking share as N grows). The cone's ratios are pinned to the exact route's to 1 % (the RK4 step at dt = 1/30
    // moves the eigenvalues by far less).
    [Fact]
    public void Localized_Carriers_Hold_More_Fisher_Information_At_Q20()
    {
        foreach (int n in new[] { 5, 7, 9 })
        {
            double[] Site(int s) { var v = new double[n]; v[s] = 1.0; return v; }
            var times = Enumerable.Range(1, 60).Select(k => 20.0 * k / 60).ToArray();
            double L(double[] wave) => WalkTime.MinEigenvalue(WalkTime.PositionFisher(W, n, 1.0, 0.05, new[] { (1.0, wave) }, times, dt: 1.0 / 30));
            double bonding = L(WalkTime.Mode(n, 1));
            var exact = new Dictionary<int, (double Seed, double Centre)> { [5] = (108.40, 93.08), [7] = (118.51, 75.58), [9] = (121.33, 97.47) }[n];
            Assert.InRange(L(Site(0)) / bonding / exact.Seed, 0.99, 1.01);
            Assert.InRange(L(Site(n / 2)) / bonding / exact.Centre, 0.99, 1.01);
        }
    }

    // F126's never-caught term is the carrier's clean evolution (times the universal e^{-Gamma t}). For the bonding
    // mode it is blind to the uniform rescaling of every bond: psi1 is an eigenstate of every rescaled hopping
    // matrix, so with every bond at J(1 + 0.1) its clean populations stay psi1(i)^2. The cone at gamma = 0 reads the
    // clean evolution; the residual is rounding (h psi1 = E psi1 to about 1e-16 per entry, carried through 500 RK4
    // steps of a commutator that vanishes in exact arithmetic: about 500 * 1e-16), gated at 1e-12. This is the linear
    // algebra of the never-caught term, not the F126 split itself. The control: the site seed's clean populations
    // do move under the same rescaling.
    [Fact]
    public void The_Bonding_Modes_Never_Caught_Term_Is_Blind_To_The_Uniform_Rescaling()
    {
        const int n = 7;
        double[] Run(double[] wave, double scale)
        {
            var cone = new Cone(W, n, 1.0, 0.0);
            for (int b = 0; b < n - 1; b++) cone.SetBond(b, b + 1, scale);
            cone.SeedPure(wave);
            for (int k = 0; k < 500; k++) cone.Step(0.01);
            return Enumerable.Range(0, n).Select(cone.Population).ToArray();
        }
        var psi1 = WalkTime.Mode(n, 1);
        var scaled = Run(psi1, 1.1);
        for (int i = 0; i < n; i++)
            Assert.True(Math.Abs(scaled[i] - psi1[i] * psi1[i]) < 1e-12, $"site {i}: {scaled[i]} against {psi1[i] * psi1[i]}");
        var seed = new double[n]; seed[0] = 1.0;
        var a = Run(seed, 1.0); var b2 = Run(seed, 1.1);
        Assert.True(a.Zip(b2, (x, y) => Math.Abs(x - y)).Max() > 0.01, "the site seed must feel the rescaling");
    }

    // The letters of every bond add up to the letter of the uniform rescaling of every bond, because h is the sum of
    // its bonds (experiments/THE_DECODER_AND_THE_FRONT.md §4). At gamma = 0 that rescaling is a rescaling of time, and
    // a carrier at rest under H does not move under it: the bonding mode's letters add up to nothing. What is left of
    // the sum is the central difference's own error, which goes as delta^2, so halving delta quarters it; the RK4 step
    // cannot break the rule, since its propagator keeps psi1 an eigenvector of every rescaled h. A first-order sum that
    // did not vanish would hold the ratio near 1. The control: the site seed moves, and its letters add up to its own
    // motion, t dp/dt, far from nothing (0.777 of their summed norms on the exact route).
    [Fact]
    public void The_Letters_Of_A_Carrier_At_Rest_Add_Up_To_Nothing()
    {
        const int n = 5;
        double Cancellation(double[] wave, double delta)
        {
            var l = WalkTime.TrajectoryLetters(W, n, 1.0, 0.0, new[] { (1.0, wave) }, dt: 0.01, sample: 0.25, tMax: 20.0, delta: delta);
            double sum2 = 0.0;
            for (int k = 0; k < l[0].Length; k++)
            {
                double s = 0.0;
                for (int b = 0; b < n - 1; b++) s += l[b][k];
                sum2 += s * s;
            }
            return Math.Sqrt(sum2) / l.Sum(row => Math.Sqrt(row.Sum(x => x * x)));
        }
        var psi1 = WalkTime.Mode(n, 1);
        double coarse = Cancellation(psi1, 0.02), fine = Cancellation(psi1, 0.01);
        Assert.InRange(coarse / fine, 3.9, 4.1);
        var seed = new double[n]; seed[0] = 1.0;
        Assert.True(Cancellation(seed, 0.01) > 0.5, "the site seed's letters add up to its motion");
    }

    // The bonding mode's weakest direction is the sum of its letters (experiments/THE_DECODER_AND_THE_FRONT.md §4):
    // along the norm vector (|l_0|, .., |l_N-2|) its normalised letters add up to the uniform rescaling's letter, of
    // which it has little (only the halo's), so at Q = 20 the norm vector's Rayleigh quotient on the cosine matrix sits
    // just above the smallest eigenvalue (it cannot lie below it): 1.20526709 times it at N = 5 on the exact route of
    // simulations/decoder_carrier_dictionary.py, letters by central differences at delta = 1e-4 (2e-8 below the
    // delta -> 0 value, far inside the gate; decoder_letters_sum_and_shots.py prints it as 1.205, and 1.206 and 1.196
    // at N = 7, 9). The cone's reading meets that value to
    // the central difference's law: doubling delta from 0.01 to 0.02 multiplies the deviation by about 4. The control:
    // the site seed's weakest direction lies elsewhere, its norm vector at 11.3 times the smallest eigenvalue there.
    [Fact]
    public void The_Bonding_Modes_Weakest_Direction_Is_The_Sum_Of_Its_Letters()
    {
        const int n = 5, m = n - 1;
        const double exact = 1.20526709;
        double Ratio(double[] wave, double delta)
        {
            var l = WalkTime.TrajectoryLetters(W, n, 1.0, 0.05, new[] { (1.0, wave) }, dt: 0.01, sample: 0.25, tMax: 20.0, delta: delta);
            var norm = l.Select(r => Math.Sqrt(r.Sum(x => x * x))).ToArray();
            var c = new double[m, m];
            for (int a = 0; a < m; a++)
                for (int b = 0; b < m; b++) c[a, b] = l[a].Zip(l[b], (x, y) => x * y).Sum() / (norm[a] * norm[b]);
            double num = 0.0, den = 0.0;
            for (int a = 0; a < m; a++)
            {
                den += norm[a] * norm[a];
                for (int b = 0; b < m; b++) num += norm[a] * c[a, b] * norm[b];
            }
            return num / den / WalkTime.MinEigenvalue(c);
        }
        var psi1 = WalkTime.Mode(n, 1);
        double e1 = Math.Abs(Ratio(psi1, 0.02) - exact), e2 = Math.Abs(Ratio(psi1, 0.01) - exact);
        var seed = new double[n]; seed[0] = 1.0;
        double site = Ratio(seed, 0.01);
        Assert.True(e2 < 1e-3, $"dt = 0.01, delta = 0.01 reads {exact + e2} against {exact}");
        Assert.InRange(e1 / e2, 3.5, 4.5);
        Assert.True(site > 5.0, $"the site seed's norm vector reads {site} times the smallest eigenvalue");
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
