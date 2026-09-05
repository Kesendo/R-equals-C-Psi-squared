using MirrorWorld;

namespace MirrorWorldTests;

// Spooky action as a named composition (Field at N=2 + two Marginals): the k=2 sighting, the middle of
// the triptych the double slit (k=1) and the cat (k=N) are the ends of. The sock drawer (the 00/11
// records) is the immortal diagonal; the spook |00><11| pays -2*gamma*2 = -2(gamma_A + gamma_B); the
// local pages are blind to it (F70 kinematics: a k=2 disagreement cannot reach a one-site page, and the
// Bell skeleton has no k=1 content). The one-sided-watching pin runs the living world (Restless) with
// only one site watched: the carrier pays exactly the watched site's rate, while page A stays I/2
// EXACTLY, as it does under every Z-watching profile (the pages carry no k=1 content and H keeps the
// poles), so what the one-sided profile changes is the carrier's price and nothing on the page.
// Meaning in docs/quantum/SPOOKY_ACTION_TRANSLATED.md, not here.
public class SpookyActionTests
{
    const double G = 0.05;
    static readonly World W = new();

    [Fact]
    public void The_Sock_Drawer_Is_Immortal_The_Spook_Pays_Twice_The_Slit()
    {
        var sa = new SpookyAction(W, G);
        Assert.True(sa.SockDrawer == 2.0);                    // |00><00| + |11><11|: the two records (literal 1.0s)
        Assert.True(sa.Spook == 2.0);                         // the one spook |00><11|, with its mirror twin
        Assert.True(sa.SpookRate == -2.0 * G * 2);            // k=2: the product Pair forms, exact
        Assert.True(sa.SpookRate == 2 * new DoubleSlit(W, G).BetweenRate);   // twice the slit's k=1 (an exact doubling)
        Assert.True(sa.SpookRate == new Cat(W, 2, G).CoherenceRate);         // = the N=2 cat: same law, same doubles
        double sock0 = sa.SockDrawer, spook0 = sa.Spook;
        for (int s = 0; s < 200; s++) sa.Watch(0.05);         // t = 10, two 1/e times
        Assert.True(sa.SockDrawer == sock0);                  // the records never move (k=0, never stepped)
        Assert.True(sa.Spook < spook0);                       // the spook pays
    }

    [Fact]
    public void The_Spook_Decays_On_The_Exp_Law()
    {
        // Field.Step is forward Euler, so the spook has the EXACT closed form (1 + rate*dt)^n under the scheme;
        // compare exactly (the argument is in DoubleSlitTests), then read the exp law as the first-order
        // limit: the deviation from e^{-4 gamma t} halves when dt halves (t = 5, one 1/e time;
        // e^{-1} * rate^2 t/2 * dt ~ 3.7e-5 at dt = 1e-3).
        double ratio1 = WatchToT(new SpookyAction(W, G), dt: 1e-3, t: 5.0, out double euler1);
        Assert.True(ratio1 == euler1, $"the Euler product is exact: {ratio1:R} vs {euler1:R}");
        double exact = Math.Exp(-4.0 * G * 5.0);
        double ratio2 = WatchToT(new SpookyAction(W, G), dt: 5e-4, t: 5.0, out _);
        double dev1 = Math.Abs(ratio1 - exact), dev2 = Math.Abs(ratio2 - exact);
        double pred1 = exact * (4.0 * G) * (4.0 * G) * 5.0 / 2.0 * 1e-3;   // e^{rt} r^2 t dt / 2
        Assert.InRange(dev1 / pred1, 0.995, 1.005);  // the closed form, to the O(r dt) it leaves out
        Assert.InRange(dev1 / dev2, 1.95, 2.05);
    }

    static double WatchToT(SpookyAction sa, double dt, double t, out double eulerProduct)
    {
        int n = (int)Math.Round(t / dt);
        double s0 = sa.Spook, factor = 1.0 + sa.SpookRate * dt;
        eulerProduct = 1.0;
        for (int s = 0; s < n; s++) { sa.Watch(dt); eulerProduct *= factor; }
        return sa.Spook / s0;
    }

    [Fact]
    public void The_Pages_Hold_While_The_Spook_Dies()
    {
        // F70 kinematics from below: the k=2 spook cannot reach a one-site page (and the skeleton has
        // no k=1 content), so both pages hold their immortal record -- off-diagonal exactly zero,
        // diagonal constant -- at every sampled time, while the spook between them visibly decays.
        var sa = new SpookyAction(W, G);
        double spook0 = sa.Spook;
        for (int tick = 0; tick <= 200; tick++)
        {
            foreach (var pg in new[] { sa.PageA, sa.PageB })
            {
                // exactly zero, not small: the page entry (0,1) sums cells that differ in ONE bit, and the
                // skeleton never seeded a k=1 cell, so every term of the partial trace is the literal 0.0.
                Assert.True(pg.Novelty == 0.0, $"the page must be structurally blind, read {pg.Novelty:R}");
                Assert.True(pg[0, 0].Real == 1.0);            // the record (bare weights: the page is I),
                Assert.True(pg[1, 1].Real == 1.0);            // untouched because k=0 cells are never stepped
            }
            if (tick < 200) sa.Watch(0.05);
        }
        Assert.True(sa.Spook < 0.4 * spook0, "the spook must visibly decay while the pages hold");

        // the same door with the answer nonzero: a k=1 cell |00><01| (the two indices differ at site 0) IS
        // seen by site 0's page and is invisible to site 1's page; so a zero above means blindness, not a
        // page that reads nothing.
        var one = new Field(W, 2, G);
        one[0b00, 0b01] = 1.0;
        Assert.True(new Marginal(one, new[] { 0 }).Novelty == 2.0);
        Assert.True(new Marginal(one, new[] { 1 }).Novelty == 0.0);
    }

    [Fact]
    public void One_Sided_Watching_The_Carrier_Pays_Bobs_Rate_And_Page_A_Holds_Under_Every_Profile()
    {
        // the normalized Bell pair in the LIVING world, H on (J=1, ZZ on too). Page A is I/2 EXACTLY under
        // EVERY Z-watching profile: the flip-flop's matrix element between |00> and |11> is the literal 0.0,
        // both poles carry the same ZZ energy, and Z-dephasing never touches populations, so every RK4
        // stage adds exact zeros to the populations and to the page's off-diagonal cells. What the profile
        // changes is the carrier's price alone: -2*gamma_B with only site 1 watched, -2*(gamma_A+gamma_B)
        // with both. So the page assertions are the same under both profiles (that is the point: the
        // adjective "unwatched" does no work on the page) and the carrier's rate is what tells them apart.
        const double dt = 0.01;
        foreach (var (profile, rate) in new[] { (new[] { 0.0, G }, -2.0 * G), (new[] { G, G }, -4.0 * G) })
        {
            var w = new Restless(W, 2, 1.0, G, siteGammas: profile, zz: 1.0);
            w.Seed(0b00, 0.5);
            w.Seed(0b11, 0.5);
            w.SeedCoherence(0b00, 0b11, 0.5);                 // the true Bell density matrix, trace 1
            var pageA = new Marginal(w, new[] { 0 });
            double carrier0 = w[0b00, 0b11].Magnitude;
            for (int tick = 0; tick <= 1000; tick++)
            {
                if (tick % 250 == 0)
                {
                    Assert.True(pageA[0, 0].Real == 0.5, $"the locked page, read {pageA[0, 0].Real:R}");
                    Assert.True(pageA[1, 1].Real == 0.5);
                    Assert.True(pageA[0, 1].Magnitude == 0.0);
                    // the carrier: on this entry the generator is the scalar `rate` (the H commutator is the exact
                    // zero above), so RK4 runs the scalar recursion rho <- rho*(1 - x + x^2/2 - x^3/6 + x^4/24),
                    // x = |rate|*dt. The bound is a ROUNDING BUDGET, 16 eps per step (the truncation x^5/120 per
                    // step is 0.2 percent of it), linear in the step count; measured deviations sit ~500x under it,
                    // while a rate off by 1 percent sits ~4e-3 away at t = 10, nine orders past the bound.
                    double x = -rate * dt;
                    double tol = tick * (Math.Pow(x, 5) / 120.0 + 16.0 * 2.220446049250313e-16) + 2.220446049250313e-16;
                    double dev = Math.Abs(w[0b00, 0b11].Magnitude / carrier0 - Math.Exp(rate * w.T));
                    Assert.True(dev <= tol, $"tick {tick}, rate {rate}: |carrier/c0 - e^(rate t)| = {dev:R} exceeds the budget {tol:R}");
                }
                if (tick < 1000) w.Step(dt);                  // to t = 10
            }
            Assert.True(w[0b00, 0b11].Magnitude < 0.5 * carrier0, "the spook must actually decay");
        }
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(5)]
    [InlineData(8)]
    public void The_Distance_Never_Enters_The_Price(int n)
    {
        // "spooky action AT A DISTANCE" -- where is the distance for us? (Tom, 2026-07-12). Only in
        // the WAY: the bonds H walks to build the pair (correlation walks, it does not jump --
        // experiments/RELAY_PROTOCOL.md). The finished between keeps no column for distance: put the
        // pair at the two ENDS of an N-chain (middle sites definite |0>) and the carrier still has
        // disagreement k=2 and price -2*gamma*2, for every N -- while the cat's carrier over the SAME
        // register pays -2*gamma*N. The watching reads disagreement, never distance.
        var world = new World();
        int far = 1 | (1 << (n - 1));                    // |10..01>: the two partners at the ends
        var pair = new Pair(world, 0, far, G);
        Assert.Equal(2, pair.Disagreement);              // k = 2, however long the chain between
        Assert.True(pair.Rate == -2.0 * G * 2);          // the price never sees the distance (the same product Pair forms)...
        Assert.True(new Pair(world, 0, (1 << n) - 1, G).Rate == -2.0 * G * n);   // ...the cat's does

        // and the end pages still lock, while the end-to-end page holds the whole far spook
        var cloud = new Field(world, n, G);
        cloud[0, 0] = 1.0;
        cloud[far, far] = 1.0;
        cloud[0, far] = 1.0;
        var endA = new Marginal(cloud, new[] { 0 });
        var ends = new Marginal(cloud, new[] { 0, n - 1 });
        for (int s = 0; s < 100; s++) cloud.Step(0.05);
        // exact, all three: the end page's off-diagonal sums cells one bit apart, none seeded (0.0 exactly);
        // the two-end page's (00,11) entry has exactly one nonzero term in its trace sum, the far spook itself.
        Assert.True(endA.Novelty == 0.0);                             // one end alone: blind (F70)
        Assert.True(ends[0b00, 0b11].Real == cloud[0, far]);          // both ends together: the whole spook
        Assert.True(ends.Novelty == 2.0 * Math.Abs(cloud[0, far]));
        Assert.True(cloud[0, far] > 0.0 && cloud[0, far] < 1.0, "the spook must have decayed, and not to zero");
    }

    [Fact]
    public void Own_Reads_The_Phenomenon_And_The_Pages_Inherit_Through_The_Cloud()
    {
        var sa = new SpookyAction(W, G);
        Assert.Equal(new[] { "sockDrawer", "spook", "pages" }, sa.Own);
        Assert.Equal(new[] { "x", "y", "z" }, sa.Inherited);
        Assert.Equal(new[] { "page", "structure", "novelty" }, sa.PageA.Own);
        Assert.Equal(new[] { "structure", "novelty", "x", "y", "z" }, sa.PageA.Inherited);
    }
}
