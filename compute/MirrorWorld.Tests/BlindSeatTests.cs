using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the blind seat (adopted 2026-08-24 from experiments/THE_SEAT_THAT_CUTS.md +
// experiments/THE_BLIND_SITE.md; no F-number yet, the GammaFold-per-site-turn precedent). Everything here
// is a count or a rank by elimination over GF(p) at two primes: no eigensolver anywhere, the Seed genre.
// The literals are the committed run's numbers (simulations/results/seat_cut_blindness/
// seat_cut_blindness_run.txt), asserted BESIDE the closed form so a gate built only from the formula
// cannot pass by construction.
public class BlindSeatTests
{
    static readonly World W = new();

    static BlindSeat Uniform(int n, bool zz) =>
        new(W, n, Enumerable.Repeat(1L, n - 1).ToArray(), heisenberg: zz);

    static BlindSeat Chain(long[] bonds, bool zz) =>
        new(W, bonds.Length + 1, bonds, heisenberg: zz);

    // the committed uniform ZZ law, every seat, count against closed form (both asserted).
    [Theory]
    [InlineData(3)] [InlineData(5)] [InlineData(7)] [InlineData(9)]
    [InlineData(11)] [InlineData(12)] [InlineData(13)]
    public void Uniform_ZzBook_Count_Is_The_Committed_Law(int n)
    {
        var bs = Uniform(n, zz: true);
        for (int s = 0; s < n; s++)
        {
            Assert.Equal(bs.UniformLaw(s), bs.Blind(s));
            Assert.Equal((Cyclotomy.Gcd(2 * s + 1, n) - 1) / 2, bs.Blind(s));
        }
    }

    // the committed XY law, every seat; the N = 11 row is the run's own list, literal.
    [Theory]
    [InlineData(6)] [InlineData(7)] [InlineData(9)] [InlineData(11)]
    [InlineData(12)] [InlineData(13)]
    public void Uniform_XyBook_Count_Is_The_Committed_Law(int n)
    {
        var bs = Uniform(n, zz: false);
        for (int s = 0; s < n; s++)
        {
            Assert.Equal(bs.UniformLaw(s), bs.Blind(s));
            Assert.Equal(Cyclotomy.Gcd(s + 1, n + 1) - 1, bs.Blind(s));
        }
    }

    [Fact]
    public void Uniform_Xy_N11_Is_The_Committed_List()
    {
        var bs = Uniform(11, zz: false);
        var expected = new[] { 0, 1, 2, 3, 0, 5, 0, 3, 2, 1, 0 };
        Assert.Equal(expected, Enumerable.Range(0, 11).Select(s => bs.Blind(s)).ToArray());
    }

    // the four committed profile rows of the criterion table (run's `criterion` part), literal.
    [Fact]
    public void The_Committed_Profile_Rows_Reproduce()
    {
        Assert.Equal(new[] { 0, 1, 0, 3, 0, 1, 0 },
            Row(Chain(new long[] { 1, 2, 1, 1, 2, 1 }, zz: true)));
        Assert.Equal(new[] { 0, 0, 0, 0, 4, 0, 0, 0, 0 },
            Row(Chain(new long[] { 1, 2, 3, 4, 4, 3, 2, 1 }, zz: true)));
        Assert.Equal(new[] { 0, 0, 2, 1, 0 },
            Row(Chain(new long[] { 1, 4, 2, 2 }, zz: true)));
        Assert.Equal(new[] { 0, 1, 0, 0 },
            Row(Chain(new long[] { 1, 3, 2 }, zz: true)));
    }

    static int[] Row(BlindSeat bs) => Enumerable.Range(0, bs.N).Select(s => bs.Blind(s)).ToArray();

    // a zero bond: the count survives (the fence is the two-halves PHRASING's, not the count's).
    // [1,0,1] carries blind 2 at every seat on both books, past the two-halves ceiling min(j, N-1-j).
    [Theory]
    [InlineData(true)] [InlineData(false)]
    public void Zero_Bond_Count_Is_The_Committed_List(bool zzBook)
    {
        var bs = Chain(new long[] { 1, 0, 1 }, zzBook);
        Assert.Equal(new[] { 2, 2, 2, 2 }, Row(bs));
    }

    [Fact]
    public void Zero_Bond_N6_Is_The_Committed_List()
    {
        var bs = Chain(new long[] { 1, 1, 0, 1, 1 }, zz: true);
        Assert.Equal(new[] { 3, 4, 3, 3, 4, 3 }, Row(bs));
    }

    // the parity forcing: XY, odd chain, odd seat -- blind at EVERY zero-free profile, the signs free.
    [Theory]
    [InlineData(5)] [InlineData(7)] [InlineData(9)]
    public void Xy_Odd_Seats_Of_Odd_Chains_Are_Parity_Forced(int n)
    {
        var profiles = new[]
        {
            new long[] { 2, -5, 7, 3, -1, 4, 6, -2 },
            new long[] { 1, 1, 1, 1, 1, 1, 1, 1 },
            new long[] { 9, 2, -9, 5, 3, -7, 1, 8 },
        };
        foreach (var p in profiles)
        {
            var bs = Chain(p.Take(n - 1).ToArray(), zz: false);
            for (int s = 1; s < n; s += 2)
                Assert.True(bs.Blind(s) >= 1, $"N={n} seat {s} must be blind");
        }
    }

    // the Heisenberg control: no parity forcing, and the irregular signed profile above has no blind seat.
    [Theory]
    [InlineData(5)] [InlineData(7)] [InlineData(9)]
    public void ZzBook_Irregular_Profile_Has_No_Blind_Seat(int n)
    {
        var bs = Chain(new long[] { 2, -5, 7, 3, -1, 4, 6, -2 }.Take(n - 1).ToArray(), zz: true);
        for (int s = 0; s < n; s++)
            Assert.Equal(0, bs.Blind(s));
    }

    // the count does not depend on the coupling's scale or sign, only on the profile's shape.
    [Fact]
    public void The_Count_Is_Scale_And_Sign_Free()
    {
        var a = Chain(new long[] { 1, 4, 2, 2 }, zz: true);
        var b = Chain(new long[] { 3, 12, 6, 6 }, zz: true);
        var c = Chain(new long[] { -1, -4, -2, -2 }, zz: true);
        Assert.Equal(Row(a), Row(b));
        Assert.Equal(Row(a), Row(c));
    }

    // the span: 1 + blind on the zero-free chain (uniform and irregular). The identity can break only at
    // a degenerate spectrum and the zero bond is NOT the discriminator: [1,1,0,1,1] gives span 7 against
    // 1+blind = 5 at seats 1 and 4, while [1,0,1] holds at every seat (asserted below).
    [Theory]
    [InlineData(3)] [InlineData(5)] [InlineData(7)]
    public void Span_Is_One_Plus_Blind_On_The_ZeroFree_Chain(int n)
    {
        foreach (var bs in new[] { Uniform(n, true), Chain(new long[] { 1, 4, 2, 2, 3, 1 }.Take(n - 1).ToArray(), true) })
            for (int s = 0; s < n; s++)
                Assert.Equal(1 + bs.Blind(s), bs.Span(s));
    }

    [Fact]
    public void Span_Exceeds_One_Plus_Blind_At_The_Committed_Seats()
    {
        var bs = Chain(new long[] { 1, 1, 0, 1, 1 }, zz: true);
        Assert.Equal(7, bs.Span(1));
        Assert.Equal(7, bs.Span(4));
        Assert.Equal(5, 1 + bs.Blind(1));
        foreach (int s in new[] { 0, 2, 3, 5 })
            Assert.Equal(1 + bs.Blind(s), bs.Span(s));
    }

    // the two-halves ceiling holds where the fence holds: blind <= min(j, N-1-j) on the zero-free chain.
    [Fact]
    public void The_Ceiling_Holds_ZeroFree_And_Falls_At_A_Zero_Bond()
    {
        var zf = Chain(new long[] { 1, 2, 1, 1, 2, 1 }, zz: true);
        for (int s = 0; s < zf.N; s++)
            Assert.True(zf.Blind(s) <= Math.Min(s, zf.N - 1 - s));
        var zb = Chain(new long[] { 1, 0, 1 }, zz: true);
        Assert.True(zb.Blind(0) > 0);   // past min(0, 3) = 0
    }

    // it walks past the wall: the laws hold at an N no eigendecomposition reaches, both books.
    [Theory]
    [InlineData(60)] [InlineData(200)]
    public void It_Walks_Past_The_Wall(int n)
    {
        foreach (bool book in new[] { true, false })
        {
            var bs = Uniform(n, book);
            for (int s = 0; s < n; s++)
                Assert.Equal(bs.UniformLaw(s), bs.Blind(s));
        }
    }

    // the count predicts the running world: under a one-seat light on the Cone, a blind mode is
    // stationary and a sighted one decays. XY N = 7, seat 3 (blind 3): mode m = 2 has a node at the seat
    // (sin(pi*2*(3+1)/8) = 0), mode m = 1 does not.
    [Fact]
    public void A_Blind_Mode_Survives_The_Cone_And_A_Sighted_One_Decays()
    {
        const int n = 7; const int seat = 3; const double dt = 0.02;
        var g = new double[n]; g[seat] = 0.5;

        double[] Mode(int m) => Enumerable.Range(0, n)
            .Select(i => Math.Sqrt(2.0 / (n + 1)) * Math.Sin(Math.PI * m * (i + 1) / (n + 1))).ToArray();

        var blind = new Cone(W, n, 1.0, 0.0, siteGammas: g);
        blind.SeedPure(Mode(2));
        double before = blind.Novelty;
        for (int t = 0; t < 200; t++) blind.Step(dt);
        Assert.True(Math.Abs(blind.Novelty - before) < 1e-13,   // stationary: the light never touches it
            $"blind-mode drift {Math.Abs(blind.Novelty - before):E2}");

        var sighted = new Cone(W, n, 1.0, 0.0, siteGammas: g);
        sighted.SeedPure(Mode(1));
        double sBefore = sighted.Novelty;
        for (int t = 0; t < 200; t++) sighted.Step(dt);
        Assert.True(sighted.Novelty < 0.8 * sBefore, "a sighted mode must lose coherence");
    }

    // the per-site Cone stays faithful to Restless's site mask where both run.
    [Fact]
    public void PerSite_Cone_Agrees_With_Restless_Single_Excitation()
    {
        const int n = 4; const double j = 1.0, dt = 0.05;
        var g = new double[] { 0.0, 0.5, 0.0, 0.0 };
        var cone = new Cone(W, n, j, 0.0, siteGammas: g);
        var rest = new Restless(W, n, j, 0.0, siteGammas: g);
        cone.Seed(0); rest.Seed(1);
        for (int t = 0; t < 30; t++) { cone.Step(dt); rest.Step(dt); }
        for (int a = 0; a < n; a++)
            Assert.Equal(rest[1 << a, 1 << a].Real, cone.Population(a), 9);
    }

    // the parity predicate itself, all four ways it can be off: the book, the chain's parity, the seat's.
    [Fact]
    public void ParityForced_Names_Exactly_The_Odd_Odd_Xy_Cell()
    {
        Assert.True(Uniform(7, zz: false).ParityForced(3));    // XY, odd chain, odd seat
        Assert.False(Uniform(7, zz: false).ParityForced(2));   // even seat
        Assert.False(Uniform(8, zz: false).ParityForced(3));   // even chain
        Assert.False(Uniform(7, zz: true).ParityForced(3));    // the ZZ book has no parity forcing
    }

    // the two buckets stay pure: the seat owns its count, law and span; the sector is the Cone's.
    [Fact]
    public void Ontology_The_Seat_Owns_Its_Count_And_Hangs_On_The_Frame()
    {
        var bs = new BlindSeat(W, 5, new long[] { 1, 1, 1, 1 });
        Assert.Equal(new[] { "blind", "law", "span" }, bs.Own);
        Assert.DoesNotContain("structure", bs.Own);       // the Cone's own stays the Cone's
        Assert.IsType<World>(bs.Parent);                  // the sector is rebuilt, not borrowed
        Assert.Contains("x", bs.Inherited);               // only the frame arrives as inheritance
    }

    // and [1,0,1] HOLDS the span identity at every seat -- the zero bond is not the discriminator.
    [Fact]
    public void The_Zero_Bond_Path_Holds_The_Span_Identity()
    {
        var bs = Chain(new long[] { 1, 0, 1 }, zz: true);
        for (int s = 0; s < 4; s++)
            Assert.Equal(1 + bs.Blind(s), bs.Span(s));
    }
}
