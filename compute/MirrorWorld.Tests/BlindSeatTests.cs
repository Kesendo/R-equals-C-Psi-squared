using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the blind seat (adopted 2026-08-24 from experiments/THE_SEAT_THAT_CUTS.md +
// experiments/THE_BLIND_SITE.md; adopted on the genre gate the GammaFold-per-site-turn way, F157 since the same evening). Everything here
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

    // the law is a UNIFORM-chain law, and it is a function of (N, book, seat) with the bond list nowhere
    // in it: off the uniform chain it parts company with the count in BOTH directions, so `law` is a
    // reading of the chain the object was told about and never a prediction for this one. The rows are
    // the committed run's `scope` part, which is where the reflection-symmetry reading died.
    [Fact]
    public void UniformLaw_Is_Not_The_Count_Off_The_Uniform_Chain()
    {
        var a = Chain(new long[] { 1, 4, 2, 2 }, zz: true);            // N = 5: MORE blindness than the law
        Assert.Equal(new[] { 0, 0, 2, 0, 0 }, Enumerable.Range(0, 5).Select(a.UniformLaw).ToArray());
        Assert.Equal(new[] { 0, 0, 2, 1, 0 }, Row(a));

        var b = Chain(new long[] { 1, 2, 1, 1, 2, 1 }, zz: true);      // N = 7: more again, at seats 1 and 5
        Assert.Equal(new[] { 0, 0, 0, 3, 0, 0, 0 }, Enumerable.Range(0, 7).Select(b.UniformLaw).ToArray());
        Assert.Equal(new[] { 0, 1, 0, 3, 0, 1, 0 }, Row(b));

        var c = Chain(new long[] { 1, 2, 3, 4, 4, 3, 2, 1 }, zz: true); // N = 9: LESS, the other direction
        Assert.Equal(new[] { 0, 1, 0, 0, 4, 0, 0, 1, 0 }, Enumerable.Range(0, 9).Select(c.UniformLaw).ToArray());
        Assert.Equal(new[] { 0, 0, 0, 0, 4, 0, 0, 0, 0 }, Row(c));
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

    // the span: 1 + blind on the zero-free chain (uniform and irregular), which holds because a zero-free
    // chain makes H simple on the Krylov complement (Corollary C). The rows are LITERALS, an exact
    // rational nullity taken off this object's route, and they are what makes this a gate at all: Span and
    // Blind read the same H at the same two primes, so a bad reduction inflates BOTH and the identity
    // alone survives it unbroken. The proof's own G6 takes both sides over the rationals for exactly that
    // reason. The identity is asserted BESIDE the literal here, never instead of it.
    [Theory]
    [InlineData(3, new[] { 1, 2, 1 },             new[] { 1, 1, 1 })]
    [InlineData(5, new[] { 1, 1, 3, 1, 1 },       new[] { 1, 1, 3, 2, 1 })]
    [InlineData(7, new[] { 1, 1, 1, 4, 1, 1, 1 }, new[] { 1, 1, 1, 2, 1, 1, 1 })]
    public void Span_Is_One_Plus_Blind_On_The_ZeroFree_Chain(int n, int[] uniformSpan, int[] irregularSpan)
    {
        var cases = new[]
        {
            (bs: Uniform(n, true), expected: uniformSpan),
            (bs: Chain(new long[] { 1, 4, 2, 2, 3, 1 }.Take(n - 1).ToArray(), true), expected: irregularSpan),
        };
        foreach (var (bs, expected) in cases)
            for (int s = 0; s < n; s++)
            {
                Assert.Equal(expected[s], bs.Span(s));         // the exact rational value, off this route
                Assert.Equal(1 + bs.Blind(s), bs.Span(s));     // and the identity, beside it
            }
    }

    // and where it breaks, the whole row is the literal: H is DEGENERATE on the Krylov complement at
    // seats 1 and 4 (the cut chain's two halves repeat each other), which is Corollary C failing, and
    // not the zero bond as such -- [1,0,1] carries a zero bond and holds at every seat.
    [Fact]
    public void Span_Exceeds_One_Plus_Blind_At_The_Committed_Seats()
    {
        var bs = Chain(new long[] { 1, 1, 0, 1, 1 }, zz: true);
        Assert.Equal(new[] { 4, 7, 4, 4, 7, 4 }, Enumerable.Range(0, 6).Select(s => bs.Span(s)).ToArray());
        Assert.Equal(5, 1 + bs.Blind(1));
        Assert.Equal(5, 1 + bs.Blind(4));
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

    // the count predicts the running world, and the guard must be able to watch the prediction fail.
    // A blind mode that is ALSO an eigenvector is a fixed point of the whole generator (||rho_dot|| =
    // 1.1e-16), so asserting that it does not move tests RK4 rather than blindness -- and the content of
    // the words "largest H-INVARIANT subspace" is exactly the blind states that DO move. So the seeds
    // here run over every blind mode AND a superposition of two of them, which rotates, and the
    // observables are the two the physics fixes: the seat's amplitude stays exactly zero (algebraically:
    // (h psi)[seat] = 0 for a blind psi, and it survives RK4 because every stage keeps the seat row of
    // rho zero) and the state stays PURE. Novelty is NOT that observable: it is not conserved by the
    // unitary part, and a blind superposition moves it by O(1) while staying perfectly blind.
    // The blind modes are counted by an INTEGER node criterion, (N+1) | m*(seat+1), and that count is
    // asserted against Blind() -- the earlier version of this test never called Blind() at all.
    [Fact]
    public void The_Cone_Cannot_Touch_A_Blind_Mode_That_Moves_And_Does_Touch_A_Sighted_One()
    {
        const int n = 7; const int seat = 3; const double dt = 0.02; const int steps = 200;
        var g = new double[n]; g[seat] = 0.5;
        var bs = Uniform(n, zz: false);                  // the Cone runs the XY book, the only one it has

        // The DST-I mode, built through its own exact symmetry rather than site by site. Two roundings
        // are avoided deliberately, and both matter: the analytic node is written as an exact 0 instead
        // of a sine that returns 1.2e-16, and the far half is MIRRORED from the near half, because
        // psi_m(n-1-i) = (-1)^(m+1) psi_m(i) holds exactly while Sin(3*pi/4) and Sin(5*pi/4) differ in
        // the last bit. Without the mirror the seat amplitude starts at 9e-18 rather than at zero and no
        // exact assertion below is available; with it, (h psi)[seat] cancels bit for bit.
        double[] Mode(int m)
        {
            var v = new double[n];
            double sign = m % 2 == 0 ? -1.0 : 1.0;          // (-1)^(m+1)
            for (int i = 0; i <= n / 2; i++)
            {
                v[i] = (m * (i + 1)) % (n + 1) == 0
                    ? 0.0
                    : Math.Sqrt(2.0 / (n + 1)) * Math.Sin(Math.PI * m * (i + 1) / (n + 1));
                v[n - 1 - i] = sign * v[i];
            }
            return v;
        }

        var blindModes = Enumerable.Range(1, n).Where(m => (m * (seat + 1)) % (n + 1) == 0).ToArray();
        Assert.Equal(bs.Blind(seat), blindModes.Length);   // the node count IS the GF(p) rank's answer

        double Purity(Cone c)
        {
            double p = 0;
            for (int a = 0; a < c.N; a++)
                for (int b = 0; b < c.N; b++) p += (c[a, b] * c[b, a]).Real;
            return p;
        }

        Cone Run(double[] psi)
        {
            var c = new Cone(W, n, 1.0, 0.0, siteGammas: g);
            c.SeedPure(psi);
            for (int t = 0; t < steps; t++) c.Step(dt);
            return c;
        }

        // every blind mode, and a superposition of two that is blind WITHOUT being an eigenvector
        var seeds = blindModes.Select(m => (name: $"mode {m}", psi: Mode(m))).ToList();
        var sup = Enumerable.Range(0, n)
            .Select(i => (Mode(blindModes[0])[i] + Mode(blindModes[1])[i]) / Math.Sqrt(2.0)).ToArray();
        seeds.Add(("a blind superposition, not an eigenvector", sup));

        double worstBlindLoss = 0.0;
        foreach (var (name, psi) in seeds)
        {
            Assert.Equal(0.0, psi[seat]);                       // blind: no amplitude at the seat, exactly
            var c = Run(psi);
            Assert.Equal(0.0, c.Population(seat));              // and it stays exactly zero, not merely small
            worstBlindLoss = Math.Max(worstBlindLoss, 1.0 - Purity(c));
        }

        // the superposition must actually MOVE, or this test has quietly become the old fixed-point one
        var moving = new Cone(W, n, 1.0, 0.0, siteGammas: g);
        moving.SeedPure(sup);
        double n0 = moving.Novelty;
        for (int t = 0; t < steps; t++) moving.Step(dt);
        Assert.True(Math.Abs(moving.Novelty - n0) > 0.1,
            $"the blind superposition must be non-stationary; Novelty moved {Math.Abs(moving.Novelty - n0):E2}");

        // a sighted mode decoheres, and the two losses are not the same kind of number: the blind one is
        // RK4 truncation, the sighted one is physics. The gate is their ratio, not a threshold on either.
        double sightedLoss = 1.0 - Purity(Run(Mode(1)));
        Assert.True(sightedLoss > 1e5 * worstBlindLoss,
            $"sighted loss {sightedLoss:E2} must dwarf the blind residual {worstBlindLoss:E2}");
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

    // the guards, which are the one-sidedness's premises rather than decoration: an H built past
    // MaxCoupling wraps int64 and the count then comes out TOO SMALL, the one direction the file's whole
    // safety argument forbids. Measured before the guard existed: N = 6 on the ZZ book at
    // |J| = 4378862956477877167 reported blind 0 at seats 1 and 4 where the count is 1.
    [Fact]
    public void A_Coupling_That_Would_Wrap_Int64_Is_Refused_And_So_Is_A_Seat_Off_The_Chain()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            new BlindSeat(W, 6, Enumerable.Repeat(4378862956477877167L, 5).ToArray()));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            new BlindSeat(W, 3, new[] { BlindSeat.MaxCoupling + 1, 1L }));
        var ok = new BlindSeat(W, 3, new[] { BlindSeat.MaxCoupling, -BlindSeat.MaxCoupling });
        Assert.Equal(0, ok.Blind(0));                       // the boundary itself still computes

        var bs = Uniform(5, zz: true);
        Assert.Throws<ArgumentOutOfRangeException>(() => bs.Blind(-1));
        Assert.Throws<ArgumentOutOfRangeException>(() => bs.Span(5));
        Assert.Throws<ArgumentOutOfRangeException>(() => bs.UniformLaw(-1));   // once returned -1 silently
        Assert.Throws<ArgumentOutOfRangeException>(() => bs.ParityForced(5));
        Assert.Throws<ArgumentException>(() => new BlindSeat(W, 1, Array.Empty<long>()));
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
        Assert.Equal(new[] { 3, 3, 3, 3 }, Enumerable.Range(0, 4).Select(s => bs.Span(s)).ToArray());
        for (int s = 0; s < 4; s++)
            Assert.Equal(1 + bs.Blind(s), bs.Span(s));
    }
}
