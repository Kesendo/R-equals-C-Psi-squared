using System.Numerics;
using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the crack (built 2026-09-01 beside experiments/THE_CRACKED_BELL.md section E,
// whose gate is simulations/cracked_bell_gate.py stage E; registry F160). The pins are exact where an
// exact route exists (integer polynomial identities, Descartes counts, a vanishing at the band edge
// compared with == 0) and windows around committed numbers where the object is a reading (a bisected
// root). No eigensolver anywhere, Seed's and Divisor's genre.
public class CrackTests
{
    static readonly Cyclotomy C = new();

    static Crack Make(int n, long uNum, long uDen = 1) => new(C, n, uNum, uDen);

    // ---- the identity, exactly: the matrix's characteristic polynomial IS the road polynomial ----

    [Theory]
    [InlineData(3, 1, 2)] [InlineData(4, 1, 2)] [InlineData(5, 3, 4)] [InlineData(6, 9, 10)]
    [InlineData(7, 0, 1)] [InlineData(8, 1, 1)] [InlineData(9, 5, 4)] [InlineData(11, 6, 5)]
    [InlineData(12, 7, 5)] [InlineData(13, 4, 1)] [InlineData(16, 9999, 10000)] [InlineData(21, 11, 10)]
    public void The_Characteristic_Polynomial_Is_The_Road_Exactly(int n, long p, long q)
    {
        var residual = Make(n, p, q).IdentityResidual();
        Assert.All(residual, c => Assert.True(c.IsZero, $"coefficient residual {c}"));
    }

    // the two ends of the road are the parent's two combs, as polynomials from two DIFFERENT recursions:
    // u = 0 (chain) is U_N(x/2), F2b's polynomial; u = 1 (ring) is 2 T_N(x/2) - 2. Both met by the matrix.
    [Theory] [InlineData(3)] [InlineData(4)] [InlineData(7)] [InlineData(10)] [InlineData(15)]
    public void The_Ends_Of_The_Road_Are_The_Two_Combs(int n)
    {
        Assert.Equal(Crack.ChebyshevSecondKind(n), Make(n, 0).CharacteristicPolynomial());
        Assert.Equal(Crack.RingPolynomial(n), Make(n, 1).CharacteristicPolynomial());
        Assert.True(Make(n, 0).IdentityHolds() && Make(n, 1).IdentityHolds() && Make(n, 3, 5).IdentityHolds());
    }

    // the two printed lines of G are one curve: a typo in either would part them. What parts them in
    // double is the argument reduction of sin((N+1)k) against sin(Nk)cos k + cos(Nk)sin k, an error
    // LINEAR in N (worst |diff| / (eps (1+u^2)) on one libm: 4.3, 20, 39, 380, 1816 at N = 4, 12, 25,
    // 200, 1001; another libm gave 4.6, 20, 48, 376, 1799), so the gate is the model N (1+u^2) eps and
    // that the ratio to it stays inside one band across three decades of N. The band [0.25, 4] is wide
    // against the measured spread (1.1 to 1.9) and is there to discriminate the EXPONENT: an N^0 or an
    // N^2 law would miss it by ~250x at the far end.
    [Fact]
    public void The_Two_Written_Forms_Of_The_Curve_Agree_To_A_Law_Linear_In_N()
    {
        const double eps = 2.220446049250313e-16;
        foreach (int n in new[] { 4, 12, 25, 200, 1001 })
        {
            double worstRatio = 0;
            foreach (double u in new[] { 0.0, 0.3, 0.999, 1.0, 1.5, 4.0 })
                for (int i = 1; i < 400; i++)
                {
                    double k = Math.PI * i / 400;
                    double ratio = Math.Abs(Crack.G(k, n, u) - Crack.GSecondForm(k, n, u)) / (n * (1 + u * u) * eps);
                    worstRatio = Math.Max(worstRatio, ratio);
                }
            Assert.InRange(worstRatio, 0.25, 4.0);            // the model is linear in N, and tight both ways
        }
    }

    // ---- the departure count, exactly ----

    [Theory]
    [InlineData(5, 1, 2)] [InlineData(8, 1, 2)] [InlineData(9, 999, 1000)] [InlineData(12, 0, 1)] [InlineData(31, 1, 1)]
    public void Nothing_Leaves_The_Band_While_The_Bond_Is_Not_Strengthened(int n, long p, long q)
    {
        var c = Make(n, p, q);
        Assert.Equal(0, c.Departures);
        Assert.Equal(0, c.LawCount);
    }

    [Theory]
    [InlineData(8, 1000000001, 1000000000)] [InlineData(8, 1001, 1000)] [InlineData(10, 3, 2)] [InlineData(12, 3, 1)] [InlineData(16, 7, 5)]
    public void Even_Rings_Shed_Two_At_Every_Strengthening(int n, long p, long q)
    {
        var c = Make(n, p, q);
        Assert.Equal(2, c.Departures);
        Assert.Equal(1, c.DeparturesAbove());
        Assert.Equal(1, c.DeparturesBelow());
        Assert.Equal(2, c.LawCount);
    }

    // odd N: one leaves at once, the second only PAST u = (N+1)/(N-1); AT the threshold the level sits
    // on the band edge, P(-2) = 0 exactly, and is not a departure. Straddled by 1/1000 in u, and met at
    // the exact rational threshold itself.
    [Theory] [InlineData(5)] [InlineData(7)] [InlineData(9)] [InlineData(11)] [InlineData(15)] [InlineData(25)] [InlineData(41)]
    public void Odd_Rings_Shed_The_Second_Only_Past_The_Threshold(int n)
    {
        var (tn, td) = Crack.OddThreshold(n);
        Assert.Equal(1, Make(n, 1000000001, 1000000000).Departures);
        Assert.Equal(1, Make(n, 1000 * tn - td, 1000 * td).Departures);      // u = thr - 1/1000
        var at = Make(n, tn, td);
        Assert.Equal(1, at.Departures);
        Assert.Equal((false, true), at.OnTheEdge());                         // P(-2) = 0, exactly
        Assert.Equal(2, Make(n, 1000 * tn + td, 1000 * td).Departures);      // u = thr + 1/1000
        Assert.Equal(2, Make(n, 4).Departures);
        // the sibling's convention u = 1 + delta reads delta = (tn - td)/td = 2/(N-1), cross-multiplied
        Assert.Equal(2 * td, (tn - td) * (n - 1));
    }

    // Descartes against the closed law, over a grid of rationals on both sides of every threshold
    [Fact]
    public void The_Descartes_Count_Is_The_Law_Everywhere()
    {
        foreach (int n in Enumerable.Range(3, 18))
            foreach (long p in new long[] { 0, 1, 3, 5, 7, 9, 10, 11, 12, 13, 15, 20, 30, 41 })
            {
                var c = Make(n, p, 10);
                Assert.Equal(c.LawCount, c.Departures);
            }
    }

    // the band edge: P(2) = 0 at u = 1 for every N, P(-2) = 0 at u = 1 for even N only
    [Theory] [InlineData(6)] [InlineData(7)] [InlineData(12)] [InlineData(13)]
    public void The_Perfect_Ring_Sits_On_The_Edge_By_Parity(int n)
    {
        Assert.Equal((true, n % 2 == 0), Make(n, 1).OnTheEdge());
        Assert.Equal((false, false), Make(n, 9, 10).OnTheEdge());
    }

    // ---- the identity past the wall: mod two primes at N+1 points, O(N) per point ----

    [Theory] [InlineData(8, 9, 10)] [InlineData(21, 11, 10)] [InlineData(200, 1, 2)] [InlineData(200, 3, 2)] [InlineData(1001, 7, 5)]
    public void The_Identity_Holds_Mod_Two_Primes_Past_The_Wall(int n, long p, long q)
        => Assert.Equal(0, Make(n, p, q).IdentityMismatchesModTwoPrimes());

    // the modular route is not vacuous: fed the wrong ring's road (u swapped for 1/u), it mismatches at
    // every one of the N+1 points, at both primes; a road of the wrong degree is refused
    [Fact]
    public void The_Modular_Route_Sees_A_Wrong_Ring()
    {
        var ring = Make(40, 3, 2);
        Assert.Equal(0, ring.IdentityMismatchesModTwoPrimes());
        var wrongU = Make(40, 2, 3).RoadPolynomial();
        foreach (long p in new[] { 2147483647L, 999999937L })
            Assert.Equal(41, ring.MismatchesAgainst(wrongU, p));
        var wrongN = Make(42, 3, 2).RoadPolynomial();
        Assert.Throws<ArgumentException>(() => ring.MismatchesAgainst(wrongN, 999999937L));
    }

    // ---- the reading: roots of the curve, tied to the exact count ----

    [Theory]
    [InlineData(4, 999, 1000)] [InlineData(5, 98, 100)] [InlineData(8, 85, 100)] [InlineData(11, 7, 10)]
    [InlineData(12, 3, 10)] [InlineData(16, 1, 1000)] [InlineData(25, 1, 2)] [InlineData(40, 999, 1000)]
    public void Exactly_N_Levels_Come_Off_The_Curve(int n, long p, long q)
    {
        // the check is that Levels() RETURNS: InBandMomenta throws when its scan count misses the exact
        // count, so the length assert below is the absence of that throw, written out
        var levels = Make(n, p, q).Levels();
        Assert.Equal(n, levels.Length);
        Assert.True(levels.All(e => Math.Abs(e) < 2.0), "a weakened bond keeps every level strictly inside the band (Perron-Frobenius)");
    }

    [Theory] [InlineData(8, 3, 2)] [InlineData(9, 3, 2)] [InlineData(9, 2, 1)] [InlineData(12, 4, 1)]
    public void The_Departed_Levels_Are_On_The_Same_Curve_Continued(int n, long p, long q)
    {
        var c = Make(n, p, q);
        var levels = c.Levels();
        Assert.Equal(n, levels.Length);
        Assert.Equal(c.Departures, levels.Count(e => Math.Abs(e) > 2.0));
        // Gershgorin: no level beyond 1 + u
        Assert.True(levels.All(e => Math.Abs(e) <= 1.0 + c.U));
        // the judge is the exact polynomial: a Newton step of q^N P at each departed level, evaluated in
        // double (q <= 2, N <= 12, so the coefficients are tame). The step is the distance to the exact
        // root, and the model for it is the bisection's own resolution: the root is found in kappa to
        // adjacent doubles (1 ulp of kappa) and E = 2 cosh kappa carries that as dE ~ kappa * eps * E,
        // so the step is a small multiple of eps * kappa * |E| (kappa ~ 1 here, ~ln u at large u); measured
        // 3.0e-16 to 9.0e-16 on these four rows (the low end moves by an ulp between libms). Gated at 8 eps max(1, kappa) |E|, and the parity mutations
        // of GBelow land at 3e-2 and 9e-3.
        var road = c.RoadPolynomial();
        foreach (double e in levels.Where(x => Math.Abs(x) > 2.0))
            Assert.True(NewtonStep(road, e) <= 8 * 2.220446049250313e-16 * Math.Max(1.0, Math.Acosh(Math.Abs(e) / 2)) * Math.Abs(e),
                $"N={n}, u={p}/{q}: departed level {e} is {NewtonStep(road, e)} off the exact polynomial's root");
    }

    // a scan too coarse for a near-double pair THROWS: the exact count is the oracle the scan is held to
    [Fact]
    public void A_Coarser_Scan_Must_Not_Silently_Undercount()
        => Assert.Throws<InvalidOperationException>(() => Make(12, 9999, 10000).InBandMomenta(points: 500));

    // the pair reading is fenced by the pair's own geometry (a straddling, isolated pair): a crack too
    // deep for it THROWS instead of returning the nearest two levels (N = 25, u = 1/2 used to return
    // 0.04859 for a true 0.05266), while the 0.971754 row and the decade-law rows pass the fence
    [Fact]
    public void Too_Deep_A_Crack_Refuses_The_Pair_Reading()
    {
        Assert.Throws<InvalidOperationException>(() => Make(25, 1, 2).SplitOfPair(1));
        Assert.Throws<InvalidOperationException>(() => Make(12, 1, 100).SplitOfPair(1));
        Assert.InRange(Make(12, 9, 10).SplitOfPair(1), 0.0336, 0.0337);
    }

    // the continued curve in its grouped form reads the departed level at LARGE u to the same model as at
    // small u: the exact polynomial's Newton step at the returned level, a small multiple of
    // eps * kappa * |E| with kappa ~ ln u. The first form of the curve returned 1e5.0107 at N = 12,
    // u = 1e5 for an exact 1e5.00001 (step 1.1e-2 against this gate's 2e-10) and a test row certified it
    // as fine; the first row here is that row, judged. The strong bond is a dimer: the top level is
    // u + 1/u + O(u^-3) (1e5.00001 and 1000.001 are that), so it is pinned to that asymptote as well,
    // to the same kappa * eps * u resolution.
    [Theory] [InlineData(12, 100000)] [InlineData(20, 10000000000)] [InlineData(30, 1000)]
    public void Large_U_Reads_The_Departed_Level_To_The_Same_Model(int n, long u)
    {
        const double eps = 2.220446049250313e-16;
        var c = Make(n, u, 1);
        var road = c.RoadPolynomial();
        var departed = c.DepartedLevels();
        Assert.Equal(c.Departures, departed.Length);
        foreach (double e in departed)
            Assert.True(NewtonStep(road, e) <= 8 * eps * Math.Acosh(Math.Abs(e) / 2) * Math.Abs(e), $"N={n}, u={u}: departed level {e} is {NewtonStep(road, e)} off the exact polynomial's root");
        double top = departed.Max(), kappa = Math.Log(u);
        Assert.True(Math.Abs(top - (u + 1.0 / u)) <= 8 * eps * kappa * u + 4.0 / (u * (double)u * u), $"top {top} against the dimer asymptote u + 1/u");
    }

    // a Newton step of the integer polynomial at e, evaluated in double by Horner
    static double NewtonStep(BigInteger[] road, double e)
    {
        double pv = 0, dv = 0;
        for (int i = road.Length - 1; i >= 0; i--) { dv = dv * e + pv; pv = pv * e + (double)road[i]; }
        return Math.Abs(pv / dv);
    }

    // just past u = 1 the distance from the edge falls under one ulp of 2 (E - 2 ~ (u-1)^2): the count
    // says two levels left, the reading refuses rather than returning 2.0 for them
    [Fact]
    public void Just_Past_U_1_The_Reading_Refuses_And_The_Count_Stands()
    {
        var c = Make(12, 1000000001, 1000000000);          // u - 1 = 1e-9
        Assert.Equal(2, c.Departures);
        Assert.Throws<InvalidOperationException>(() => c.DepartedLevels());
        Assert.True(Make(12, 1000001, 1000000).DepartedLevels().All(e => Math.Abs(e) > 2.0));   // u - 1 = 1e-6 reads
    }

    // and at N = 200, u = 1000, where the first form overflowed, the reading simply returns
    [Fact]
    public void Past_The_Wall_The_Continued_Curve_Still_Reads()
    {
        var departed = Make(200, 1000, 1).DepartedLevels();
        Assert.Equal(2, departed.Length);
        Assert.InRange(departed.Max(), 1000.0, 1001.0);
        Assert.InRange(departed.Min(), -1001.0, -1000.0);
    }

    // the proof's committed detuned Perron root: 1.9999500012500 at N = 4, delta = 1e-4 (PROOF_RING_GAP_DOMINANCE
    // section Scope prints 1.9999500; THE_CRACKED_BELL section The three closed forms carries the rest)
    [Fact]
    public void The_Perron_Root_Reproduces_The_Proofs_Detuned_Floor()
    {
        double top = Make(4, 9999, 10000).Levels()[0];
        Assert.InRange(top, 1.99995000124, 1.99995000126);                    // committed 1.9999500012500
        // at N = 4 the polynomial is biquadratic, x^4 - (3+u^2) x^2 + (1-u)^2, so the top root has a closed
        // form; two double evaluations of one number, each within a few ulp: 4 ulp of 2 is the model
        double u = 0.9999, a = 3 + u * u, b = (1 - u) * (1 - u);
        double closed = Math.Sqrt((a + Math.Sqrt(a * a - 4 * b)) / 2);
        Assert.True(Math.Abs(top - closed) <= 4 * 2.220446049250313e-16 * 2.0, $"{top} vs closed form {closed}");
        // the series over the Perron shift, top = 2 - delta/2 + delta^2/8 - delta^4/128 + ... (the delta^3
        // term vanishes, 60-digit check), read as a law: (top - 2 + delta/2)/delta^2 - 1/8 falls TWO decades
        // per decade of delta. A first draft of this gate wrote one decade and the run said 100.03. Only ONE
        // decade pair is gated, and the reason is the reading's own floor: the root is bisected to ~2e-16,
        // which over delta^2 is 2e-10 at delta = 1e-3 (2.5% of the 7.8e-9 residual there) and 2e-8 at
        // delta = 1e-4, where it swamps the 7.8e-11 residual entirely (the run said 37 for that pair).
        double Coefficient(long q) { double d = 1.0 / q; return (Make(4, q - 1, q).Levels()[0] - 2 + d / 2) / (d * d); }
        double r1 = Math.Abs(Coefficient(100) - 0.125), r2 = Math.Abs(Coefficient(1000) - 0.125);
        Assert.InRange(r1 / r2, 85.0, 120.0);
    }

    // at u = 1 the in-band momenta are the INTERIOR of the ring comb: the edge levels (k = 0 at every N,
    // k = pi at even N) belong to OnTheEdge, and the tie in-band + departures + edge = N still holds
    [Theory] [InlineData(6)] [InlineData(7)]
    public void The_Ring_Ends_Interior_Comb_Excludes_The_Edge(int n)
    {
        var c = Make(n, 1);
        var ks = c.InBandMomenta();
        Assert.Equal(n - c.EdgeLevels(), ks.Length);
        Assert.All(ks, k => Assert.True(k > 0 && k < Math.PI));
        Assert.Equal(n, c.Levels().Length);
    }

    // the cross-dock with Warble: the ratio the (1,1) clock keeps (0.9705 propagated, WarbleTests) is the
    // exact split's, 0.971754 read off the curve alone (gate E6b); no matrix on this side
    [Fact]
    public void The_Split_Ratio_Is_The_Curves_And_Meets_The_Warbles_Pin()
    {
        var c = Make(12, 9, 10);
        double s1 = c.SplitOfPair(1), s2 = c.SplitOfPair(2);
        Assert.InRange(s1 / s2, 0.97170, 0.97180);                            // committed 0.971754
        Assert.InRange(Math.PI / (2 * s1), 46.69, 46.71);                     // T_zero from the law, committed 46.697
        Assert.InRange(Math.PI / (2 * s2), 45.37, 45.39);                     // committed 45.378
    }

    // the split's next order as a LAW, not a number: the residual to (4 delta/N)[1 + delta c_m] falls one
    // decade per decade of delta (gate E5a's band 8.5..12), at an even and an odd N
    [Theory] [InlineData(12, 1)] [InlineData(25, 3)] [InlineData(40, 1)]
    public void The_Split_Correction_Is_A_Decade_Law(int n, int m)
    {
        double Residual(double delta)
        {
            long q = 100000, p = (long)Math.Round(q * (1 - delta));
            double split = Make(n, p, q).SplitOfPair(m);
            // gate E5a's c_measured minus c_closed: the measured (split/(4 delta/N) - 1)/delta against c_m
            return Math.Abs((split / (4 * delta / n) - 1) / delta - Formulas.F160_SplitCorrection(n, m));
        }
        double ratio = Residual(1e-2) / Residual(1e-3);
        Assert.InRange(ratio, 8.5, 12.0);
    }

    // c_m changes sign at N sin^2 k_m = 2: for m = 1 between N = 19 and 20 (root 19.03)
    [Fact]
    public void The_Correction_Changes_Sign_Between_Nineteen_And_Twenty()
    {
        Assert.True(Formulas.F160_SplitCorrection(19, 1) > 0);
        Assert.True(Formulas.F160_SplitCorrection(20, 1) < 0);
        Assert.InRange(Formulas.F160_SplitCorrection(19, 1), 7.8e-4, 8.0e-4);   // committed +7.9e-4
        Assert.True(Formulas.F160_SplitCorrection(12, 1) > 0);
        Assert.True(Formulas.F160_SplitCorrection(25, 1) < 0);
    }

    // ---- the ontology ----

    [Fact]
    public void The_Crack_Owns_The_Road_And_Inherits_The_Two_Combs_And_Nothing_Else()
    {
        var c = Make(12, 17, 20);
        Assert.Equal(new[] { "road", "departures" }, c.Own);
        Assert.IsType<Cyclotomy>(c.Parent);
        Assert.Equal(2, c.Inherited.Count);
        Assert.All(c.Inherited, o => Assert.Contains("turn fractions", o));
        // the frame does NOT arrive: Cyclotomy has no parent (it inherits nothing itself, since 2026-08-03),
        // so this is the first object whose inherited bucket is non-empty and has no frame in it; the
        // ontology reporting where the crack hangs, not an omission.
        foreach (var frame in new[] { "x", "y", "z" })
            Assert.DoesNotContain(frame, c.Inherited);
        // u is reduced on entry
        Assert.Equal((17L, 20L), (c.UNum, c.UDen));
        Assert.Equal((1L, 1L), (Make(5, 30, 30).UNum, Make(5, 30, 30).UDen));
    }
}
