using MirrorWorld;

namespace MirrorWorldTests;

// Spooky action as a named composition (Field at N=2 + two Marginals): the k=2 sighting, the middle of
// the triptych the double slit (k=1) and the cat (k=N) are the ends of. The sock drawer (the 00/11
// records) is the immortal diagonal; the spook |00><11| pays -2*gamma*2 = -2(gamma_A + gamma_B); the
// local pages are blind to it (F70 kinematics: a k=2 disagreement cannot reach a one-site page, and the
// Bell skeleton has no k=1 content). The one-sided-watching pin runs the living world (Restless) with
// only one site watched: the carrier pays exactly the watched site's rate while the unwatched page
// stays I/2 to machine precision. Meaning in docs/quantum/SPOOKY_ACTION_TRANSLATED.md, not here.
public class SpookyActionTests
{
    const double G = 0.05;
    static readonly World W = new();

    [Fact]
    public void The_Sock_Drawer_Is_Immortal_The_Spook_Pays_Twice_The_Slit()
    {
        var sa = new SpookyAction(W, G);
        Assert.Equal(2.0, sa.SockDrawer, 12);                 // |00><00| + |11><11|: the two records
        Assert.Equal(2.0, sa.Spook, 12);                      // the one spook |00><11|, with its mirror twin
        Assert.Equal(-4.0 * G, sa.SpookRate, 12);             // k=2: -2*gamma*2 = -2(gamma+gamma)
        Assert.Equal(2 * new DoubleSlit(W, G).BetweenRate, sa.SpookRate, 12);   // twice the slit's k=1
        Assert.Equal(new Cat(W, 2, G).CoherenceRate, sa.SpookRate, 12);         // = the N=2 cat: same law
        double sock0 = sa.SockDrawer, spook0 = sa.Spook;
        for (int s = 0; s < 200; s++) sa.Watch(0.05);         // t = 10, two 1/e times
        Assert.Equal(sock0, sa.SockDrawer, 12);               // the records never move (k=0, immortal)
        Assert.True(sa.Spook < spook0);                       // the spook pays
    }

    [Fact]
    public void The_Spook_Decays_On_The_Exp_Law()
    {
        var sa = new SpookyAction(W, G);
        double spook0 = sa.Spook;
        const double dt = 1e-3; const int n = 5000;           // t = 5, the 1/e time 1/(4*gamma)
        for (int s = 0; s < n; s++) sa.Watch(dt);
        Assert.Equal(Math.Exp(-4.0 * G * 5.0), sa.Spook / spook0, 2);   // -> e^{-4*gamma*t} = 1/e
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
                Assert.Equal(0.0, pg.Novelty, 15);            // machine zero: structurally absent, not small
                Assert.Equal(1.0, pg[0, 0].Real, 15);         // the record (bare weights: the page is I)
                Assert.Equal(1.0, pg[1, 1].Real, 15);
            }
            if (tick < 200) sa.Watch(0.05);
        }
        Assert.True(sa.Spook < 0.4 * spook0, "the spook must visibly decay while the pages hold");
    }

    [Fact]
    public void One_Sided_Watching_The_Carrier_Pays_Bobs_Rate_And_The_Unwatched_Page_Never_Moves()
    {
        // the normalized Bell pair in the LIVING world, H on (J=1, ZZ on too), only site 1 watched:
        // siteGammas = {0, G}. The poles |00>,|11> are H-eigenstates of EQUAL energy (the flip-flop
        // annihilates both; the ZZ bond gives both the same diagonal), so the carrier rotates nothing
        // and pays exactly the watched site's -2*gamma_B (RK4-accurate) -- while the unwatched page
        // stays I/2 to machine precision at all times. The watching on one side never reaches the
        // other side's page; only the between pays.
        var w = new Restless(W, 2, 1.0, G, siteGammas: new[] { 0.0, G }, zz: 1.0);
        w.Seed(0b00, 0.5);
        w.Seed(0b11, 0.5);
        w.SeedCoherence(0b00, 0b11, 0.5);                     // the true Bell density matrix, trace 1
        var pageA = new Marginal(w, new[] { 0 });
        double carrier0 = w[0b00, 0b11].Magnitude;
        for (int tick = 0; tick <= 1000; tick++)
        {
            if (tick % 250 == 0)
            {
                Assert.Equal(0.5, pageA[0, 0].Real, 12);      // I/2, machine-exact: the locked page
                Assert.Equal(0.5, pageA[1, 1].Real, 12);
                Assert.Equal(0.0, pageA[0, 1].Magnitude, 12);
                Assert.Equal(Math.Exp(-2.0 * G * w.T), w[0b00, 0b11].Magnitude / carrier0, 6);
            }
            if (tick < 1000) w.Step(0.01);                    // to t = 10, one 1/e time of -2*gamma_B
        }
        Assert.True(w[0b00, 0b11].Magnitude < 0.5 * carrier0, "the spook must actually decay");
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
        Assert.Equal(-4.0 * G, pair.Rate, 12);           // the price never sees the distance...
        Assert.Equal(-2.0 * G * n, new Pair(world, 0, (1 << n) - 1, G).Rate, 12);   // ...the cat's does

        // and the end pages still lock, while the end-to-end page holds the whole far spook
        var cloud = new Field(world, n, G);
        cloud[0, 0] = 1.0;
        cloud[far, far] = 1.0;
        cloud[0, far] = 1.0;
        var endA = new Marginal(cloud, new[] { 0 });
        var ends = new Marginal(cloud, new[] { 0, n - 1 });
        for (int s = 0; s < 100; s++) cloud.Step(0.05);
        Assert.Equal(0.0, endA.Novelty, 15);                          // one end alone: blind (F70)
        Assert.Equal(cloud[0, far], ends[0b00, 0b11].Real, 12);       // both ends together: the whole spook
        Assert.Equal(2.0 * Math.Abs(cloud[0, far]), ends.Novelty, 12);
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
