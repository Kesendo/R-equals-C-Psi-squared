using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the first mirror in the world of mirrors (the fold lattice, adopted 2026-07-03
// from PROOF_CODIM1_BY_ADDITIVITY section 7). Every leg is an exact entry-wise identity between two
// block matrices -- a rearrangement, not an approximation -- so the pins demand machine zero, and the
// trajectory fold (the partner block running backward at the price 2*N*gamma) is pinned at RK4
// tolerance from two independent runs.
public class MirrorTests
{
    const double J = 1.0, G = 0.5;
    static readonly World W = new();

    // every leg holds on every block of the lattice, entry by entry, machine zero (N=4 and N=5).
    [Fact]
    public void Legs_Are_Exact_On_Every_Block()
    {
        foreach (int n in new[] { 4, 5 })
        {
            var mirror = new Mirror(W, n, J, G);
            for (int p = 0; p <= n; p++)
                for (int q = 0; q <= n; q++)
                {
                    Assert.True(mirror.TransposeResidual(p, q) < 1e-12, $"t leg broke at N={n} ({p},{q})");
                    Assert.True(mirror.BraFoldResidual(p, q) < 1e-12, $"f_P leg broke at N={n} ({p},{q})");
                    Assert.True(mirror.KetFoldResidual(p, q) < 1e-12, $"f_Q leg broke at N={n} ({p},{q})");
                    Assert.True(mirror.KleinResidual(p, q) < 1e-12, $"Klein leg broke at N={n} ({p},{q})");
                }
        }
    }

    // the fold's price is the deepest rate in the world: 2*N*gamma.
    [Fact]
    public void The_Price_Is_The_Deepest_Rate()
    {
        Assert.Equal(5.0, new Mirror(W, 5, J, G).Price, 12);
        Assert.Equal(2.0 * 8 * G, new Mirror(W, 8, J, G).Price, 12);
    }

    // the orbit of the interior core (1,1) at N=5: the four blocks the exclusion propagation walked
    // (RemainderR4InteriorExclusionTests in the main repo), two kept and two paid.
    [Fact]
    public void Orbit_Of_The_Interior_Core_Is_The_Propagation_Set()
    {
        var images = Mirror.OrbitImages(5, 1, 1);
        var kept = images.Where(x => x.FoldParity == 0).Select(x => (x.P, x.Q)).Distinct().OrderBy(x => x).ToArray();
        var paid = images.Where(x => x.FoldParity == 1).Select(x => (x.P, x.Q)).Distinct().OrderBy(x => x).ToArray();
        Assert.Equal(new[] { (1, 1), (4, 4) }, kept);
        Assert.Equal(new[] { (1, 4), (4, 1) }, paid);
    }

    // the lattice compresses to a fundamental domain: 36 blocks fall into 6 orbits at N=5.
    [Fact]
    public void The_Lattice_Folds_To_Six_Orbits_At_N5()
        => Assert.Equal(6, new Mirror(W, 5, J, G).OrbitCount());

    // a self-folded block (even N, q = N/2) is its own mirror image and pays the price out of its own
    // trace: trace L = -(price/2)*dim, exactly, for every p.
    [Fact]
    public void SelfFolded_Blocks_Pay_The_Price_From_Their_Own_Trace()
    {
        var mirror = new Mirror(W, 6, J, G);
        for (int p = 0; p <= 6; p++)
        {
            var (trace, law) = mirror.SelfFoldedTrace(p);
            Assert.Equal(law, trace, 9);
        }
    }

    // the trajectory fold: x forward under L(1,2), w backward under the partner -L(1,N-2); the mirror
    // predicts w(t) = exp(price*t) * fold(x(t)). Two independent RK4 runs agree tick by tick.
    [Fact]
    public void The_Mirror_Runs_The_Partner_Backward_At_The_Price()
    {
        var mirror = new Mirror(W, 4, J, G);
        var (_, nx, nw, worst) = mirror.TrajectoryFold(1, 2, dt: 0.005, ticks: 200);
        Assert.True(worst < 1e-6, $"the trajectory fold drifted: {worst:E2}");
        double price = mirror.Price;                               // |w|/|x| must be exp(price*t) exactly
        Assert.Equal(Math.Exp(price * 200 * 0.005), nw[200] / nx[200], 4);
    }

    // the rules turned around (the mirror's rho-level face): the anti-watched world (rate -2g(N-k))
    // is the normal world read through the bra complement, entry for entry, at every time.
    [Fact]
    public void The_AntiWatched_World_Is_The_World_Read_Through_The_Complement()
    {
        const int n = 3; const double dt = 0.05;
        int s = 1, sbar = (1 << n) - 1 - s, dim = 1 << n;
        var normal = new Restless(W, n, J, G);
        normal.Seed(s, 0.5); normal.Seed(sbar, 0.5);
        var anti = new Restless(W, n, J, G, antiWatching: true);
        anti.SeedCoherence(s, sbar, 0.5);
        for (int t = 0; t < 40; t++) { normal.Step(dt); anti.Step(dt); }
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                Assert.True((anti[i, j] - normal[i, dim - 1 - j]).Magnitude < 1e-12,
                    $"read-through broke at ({i},{j})");
    }

    // the conservation law is not taken, it moves: the anti-world's trace dies while its ANTI-trace
    // holds at 1, the exact twin of the normal world's conserved trace.
    [Fact]
    public void The_AntiWorld_Conserves_The_AntiTrace_Instead_Of_The_Trace()
    {
        const int n = 3;
        int s = 1, sbar = (1 << n) - 1 - s;
        var anti = new Restless(W, n, J, G, antiWatching: true);
        anti.SeedCoherence(s, sbar, 0.5);
        Assert.Equal(1.0, anti.AntiStructure, 12);
        for (int t = 0; t < 40; t++) anti.Step(0.05);
        Assert.Equal(1.0, anti.AntiStructure, 8);                  // the moved law holds
        Assert.True(anti.Structure < 0.2, $"the anti-world's trace must die; got {anti.Structure:0.000}");
    }

    // the two worlds' disagreement histograms are each other read backward (k <-> N-k).
    [Fact]
    public void The_Histograms_Mirror_Each_Other()
    {
        const int n = 3; const double dt = 0.05;
        int s = 1, sbar = (1 << n) - 1 - s;
        var normal = new Restless(W, n, J, G);
        normal.Seed(s, 0.5); normal.Seed(sbar, 0.5);
        var anti = new Restless(W, n, J, G, antiWatching: true);
        anti.SeedCoherence(s, sbar, 0.5);
        for (int t = 0; t < 40; t++) { normal.Step(dt); anti.Step(dt); }
        var hn = normal.WeightByDisagreement();
        var ha = anti.WeightByDisagreement();
        for (int k = 0; k <= n; k++)
            Assert.Equal(hn[k], ha[n - k], 10);
    }
}
