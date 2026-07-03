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
}
