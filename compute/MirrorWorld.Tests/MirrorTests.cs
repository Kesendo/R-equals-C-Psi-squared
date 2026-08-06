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
                    // The four legs do NOT all deserve the same gate, and the split is the
                    // content: the legs that COST NOTHING are exact, the legs that PAY THE
                    // PRICE are not. t and Klein only permute indices and flip signs, so no
                    // cell is ever recomputed by a second arithmetic route and there is nothing
                    // for a rounding error to enter through: exactly 0.0, at every gamma and J.
                    // f_P and f_Q additionally subtract the price on the diagonal, and paying
                    // is an arithmetic act, so they carry one rounding of a quantity the size
                    // of the price. See Folds_Are_Exact_Only_When_Gamma_Is_Dyadic below, which
                    // which is where that is measured and where the error model is stated.
                    Assert.True(mirror.TransposeResidual(p, q) == 0.0, $"t leg broke at N={n} ({p},{q})");
                    Assert.True(mirror.KleinResidual(p, q) == 0.0, $"Klein leg broke at N={n} ({p},{q})");
                    Assert.True(mirror.BraFoldResidual(p, q) <= FoldBound(mirror),
                        $"f_P leg broke at N={n} ({p},{q})");
                    Assert.True(mirror.KetFoldResidual(p, q) <= FoldBound(mirror),
                        $"f_Q leg broke at N={n} ({p},{q})");
                }
        }
    }

    // Machine epsilon, the unit the fold's error model is stated in.
    const double Eps = 2.220446049250313e-16;

    /// <summary>The fold leg's error bound, DERIVED rather than fitted. The identity being
    /// checked is -2γ(N−k) = −(−2γk) − 2Nγ, and the float route rounds four times: forming
    /// 2γk, forming 2γ(N−k), forming the price 2Nγ, and the subtraction itself. With u = eps/2
    /// the unit roundoff, each contributes at most u times a quantity bounded by the price, and
    /// the four sum to 3u·price = 1.5·eps·price.
    ///
    /// <para>Do not read the constant 1.5 as a measured margin. A random sweep over γ ∈ [1,2)
    /// (the ratio is invariant under γ → 2γ, so that binade is the whole space) reaches 0.99995
    /// of eps·price at N=6 and stays just under 1 out to N=300, rising with N. So a gate at
    /// 1.0·eps·price would pass with five decimal places of margin and no headroom at all; the
    /// 1.5 comes from the rounding count, and the fact that the measurement sits below it is a
    /// check on the derivation rather than the source of it.</para></summary>
    static double FoldBound(Mirror m) => 1.5 * Eps * m.Price;

    /// <summary>The exactness this world claims is REAL for two legs and CONDITIONAL for two, and
    /// the suite ran three months at gamma = 0.5 without seeing the difference.
    ///
    /// <para>t and Klein only permute indices and flip signs, so no cell is ever recomputed by a
    /// second arithmetic route: exactly 0.0 at every gamma AND every J tried. f_P and f_Q also
    /// subtract the price, and -2*gamma*k - 2*N*gamma is a different float route to
    /// -2*gamma*(N-k), so they round.</para>
    ///
    /// <para>The condition for the folds to stay bit-exact is NOT "gamma is dyadic". That label
    /// is wrong twice over: every finite double already IS a dyadic rational (0.07 is
    /// 1261007895663739/2^54), and gamma = 0.07, 0.123 and 0.001 are exact at N=2..8 while
    /// 1/7 is exact to N=4 and rounds from N=5. Exactness belongs to the PAIR (gamma, N). A
    /// sufficient condition, and the one gamma=0.5 satisfies so comfortably that it hid the
    /// question: write 2*gamma = m*2^e with m odd; if m*N &lt; 2^53 then every 2*gamma*k is
    /// representable and the subtraction is exact. gamma = 0.5 has m = 1, the narrowest possible
    /// significand; the project's canonical gamma = 0.05 has m = 3602879701896397, the widest,
    /// which is why it rounds from N=3. The condition is sufficient, not necessary: roundings can
    /// still cancel above it, so this test states which (gamma, N) it MEASURED and claims no
    /// characterisation.</para></summary>
    [Fact]
    public void The_Legs_That_Pay_The_Price_Are_The_Legs_That_Round()
    {
        // measured exact at N=2..8: a narrow significand of 2*gamma, m*N well under 2^53
        double[] narrow = { 0.25, 0.5, 1.0, 0.07, 0.001 };
        // measured to round from N=3 up: a full-width significand
        double[] wide = { 0.05, 0.1, 0.3, 1.0 / 3.0 };

        foreach (int n in new[] { 4, 5, 6 })
        {
            foreach (double g in narrow)
            {
                var m = new Mirror(W, n, J, g);
                for (int p = 0; p <= n; p++)
                    for (int q = 0; q <= n; q++)
                    {
                        Assert.True(m.BraFoldResidual(p, q) == 0.0,
                            $"f_P was measured exact at gamma={g}, N={n} ({p},{q}) and is not");
                        Assert.True(m.KetFoldResidual(p, q) == 0.0,
                            $"f_Q was measured exact at gamma={g}, N={n} ({p},{q}) and is not");
                    }
            }

            foreach (double g in wide)
            {
                var m = new Mirror(W, n, J, g);
                double worstFold = 0.0, worstFree = 0.0;
                for (int p = 0; p <= n; p++)
                    for (int q = 0; q <= n; q++)
                    {
                        worstFold = Math.Max(worstFold, Math.Max(m.BraFoldResidual(p, q), m.KetFoldResidual(p, q)));
                        worstFree = Math.Max(worstFree, Math.Max(m.TransposeResidual(p, q), m.KleinResidual(p, q)));
                    }

                // the discriminating half: the legs that cost nothing stay exact even here
                Assert.True(worstFree == 0.0,
                    $"t/Klein must stay exact at gamma={g}, N={n}, got {worstFree:E3}");
                // and the legs that pay obey the derived bound, not a chosen number
                Assert.True(worstFold <= FoldBound(m),
                    $"fold residual {worstFold:E3} exceeds the derived bound {FoldBound(m):E3} at gamma={g}, N={n}");
                // and this gamma really does break bit-exactness, so the row above measures something
                Assert.True(worstFold > 0.0,
                    $"gamma={g} was measured to round at N>=3 and did not at N={n}; the break-input " +
                    "is gone and the exactness half of this test now proves nothing");
            }
        }
    }

    /// <summary>The legs that cost nothing are exact along the OTHER axis too. The
    /// no-second-arithmetic-route argument says nothing about gamma in particular, so a J sweep
    /// is the check that can break it and the suite ran only at J = 1.</summary>
    [Fact]
    public void The_Free_Legs_Are_Exact_At_Every_Coupling_Too()
    {
        foreach (int n in new[] { 4, 5 })
            foreach (double j in new[] { 0.0, 0.1, 1.0 / 3.0, 2.71828, 1e-9 })
            {
                var m = new Mirror(W, n, j, 0.05);
                for (int p = 0; p <= n; p++)
                    for (int q = 0; q <= n; q++)
                    {
                        Assert.True(m.TransposeResidual(p, q) == 0.0, $"t leg broke at J={j}, N={n} ({p},{q})");
                        Assert.True(m.KleinResidual(p, q) == 0.0, $"Klein leg broke at J={j}, N={n} ({p},{q})");
                    }
            }
    }
    /// <summary>Where a RATE reading is still the whole truth. The answer is not "only in the
    /// empty world", which is what this test claimed on its first draft: it is whole on 4N of
    /// the (N+1)^2 blocks, at EVERY coupling.
    ///
    /// <para>Pin one side and every cell of the block carries the same disagreement, so the
    /// diagonal is a constant c times the identity; the surviving side's hops enter as -iJ/+iJ,
    /// and the block is c*Id + i*J*A with A real symmetric, whose eigenvalues are c + i*J*mu
    /// with mu real. The real part is exactly c however hard J is driven. The Hamiltonian is
    /// not absent there, it only turns the phase.</para>
    ///
    /// <para>The certificate is entry-wise, which is why no eigensolver appears: a constant
    /// diagonal, a purely imaginary off-diagonal, and that off-diagonal symmetric. Those three
    /// together ARE the statement that the real part cannot move, and each is a cell comparison.
    /// (Cross-checked in numpy while this was written, at N=4,5,6 and J=1,3: 16 / 20 / 24 blocks,
    /// exactly 4N, real parts matching the diagonal to 1e-15.)</para>
    ///
    /// <para>Relation to the main repo's pin, which is the weaker statement: there
    /// F1SpectrumStatisticsTests.MaxF1RatePairingDistance_IsBlind_ToABreakThatMovesOnlyAFrequency
    /// needs a constructed spectrum to exhibit a break the rate projection cannot see. Here the
    /// blindness is LOCATED rather than merely exhibited, and it does not say a rate reading is
    /// worthless: on 4N blocks it is exactly right, and on the other (N-1)^2 it is not.</para></summary>
    [Fact]
    public void The_Rate_Reading_Is_Whole_On_Exactly_Four_N_Blocks()
    {
        foreach (int n in new[] { 4, 5, 6 })
            foreach (double j in new[] { 0.25, 1.0, 3.0 })
            {
                var m = new Mirror(W, n, j, G);
                int whole = 0;
                for (int p = 0; p <= n; p++)
                    for (int q = 0; q <= n; q++)
                    {
                        var (spread, realMass, asym) = m.RateReadingSplit(p, q);
                        bool onePinned = p == 0 || p == n || q == 0 || q == n;

                        // the off-diagonal is ALWAYS purely imaginary and symmetric: the hops are
                        // the only thing there and they enter as -iJ / +iJ. This half is what makes
                        // the constant-diagonal test sufficient, and it holds on every block.
                        Assert.True(realMass == 0.0,
                            $"the off-diagonal grew a real part at N={n} J={j} ({p},{q}): {realMass:E3}");
                        Assert.True(asym == 0.0,
                            $"the off-diagonal is not symmetric at N={n} J={j} ({p},{q}): {asym:E3}");

                        if (onePinned)
                        {
                            Assert.True(spread == 0.0,
                                $"a pinned side must give a constant diagonal at N={n} ({p},{q}), spread {spread:E3}");
                            whole++;
                        }
                        else
                        {
                            // both sides free: the diagonal spreads, and that spread is the room
                            // J then has to move the real parts across
                            Assert.True(spread > 0.0,
                                $"both sides can hop at N={n} ({p},{q}); if the diagonal is constant " +
                                "there, the count below is wrong and the mechanism is not what this test says");
                        }
                    }

                Assert.Equal(4 * n, whole);
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

    // past the wall: the fold partner of the memory cut is again a memory cut (block (1,N-1), one hole
    // against one excitation, N^2 in site labels), and the fold leg holds cell by cell at N=60 --
    // where the full spectrum died at N=8. The wall was a property of the spectrum, not of the mirror.
    [Fact]
    public void The_Mirror_Walks_Past_The_Wall()
    {
        var mirror = new Mirror(W, 60, J, G);
        var (res, dim) = mirror.PastTheWallResidual();
        Assert.Equal(3600, dim);                                   // N^2, not 4^N
        // This IS a fold leg (BuildSiteBlock subtracts the same price), so it is owed the same
        // derived bound rather than the 1e-12 it carried until 2026-08-06, which was ~75x looser
        // than the arithmetic allows and would have passed on a construction that had started to
        // approximate. At N=60, gamma=0.5 the bound is 1.5*eps*120 = 4.0e-14.
        Assert.True(res <= FoldBound(mirror),
            $"the fold leg broke past the wall: {res:E2} against the derived bound {FoldBound(mirror):E2}");
    }

    // the trajectory fold past the wall: two independent RK4 runs at N=40, related by exp(price*t).
    [Fact]
    public void The_Trajectory_Fold_Holds_Past_The_Wall()
    {
        var mirror = new Mirror(W, 40, J, G);
        var (_, nx, nw, worst) = mirror.PastTheWallTrajectory(dt: 0.002, ticks: 50);
        Assert.True(worst < 1e-6, $"the past-the-wall trajectory fold drifted: {worst:E2}");
        Assert.Equal(Math.Exp(mirror.Price * 50 * 0.002), nw[50] / nx[50], 3);
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
