using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>From-below pins for <see cref="CrackedRingExactCurveClaim"/> (F160, proof
/// <c>docs/proofs/PROOF_CRACKED_RING_EXACT_CURVE.md</c>). Every closed form is met by a SECOND route: the
/// polynomial against a fraction-free Bareiss determinant of the integer matrix, the factor forms against the
/// unfactored polynomial at x = +-2, the departure count against Descartes' rule on the shifted integer polynomial
/// (exact for a real-rooted polynomial, and it sees a root ON the edge as a vanishing constant term), c_m against
/// roots of the curve at two crack depths, and the chain-end velocity against a finite difference of the
/// polynomial's roots. The Bareiss, first-kind-recursion, edge-factor and Descartes rows are exact integer
/// comparisons; the float-route rows are six: four carry an error MODEL, a rounding bound or a decade law, with no
/// bare threshold (the curve at 2 cos k, c_m, the velocity, the F65 magnitudes), and two are counts or signs with
/// no tolerance at all (the N real roots found by bisection, the sign of c_1 at N = 19 and 20). A first version
/// classified roots against the band edge with a literal 1e-9, where Descartes was the exact route.</summary>
public class CrackedRingExactCurveClaimTests
{
    private static CrackedRingExactCurveClaim BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var qubitAnchor = new QubitDimensionalAnchorClaim();
        var f66 = new F66PoleModesPi2Inheritance(ladder, qubitAnchor);
        var f65 = new F65XxChainSpectrumPi2Inheritance(ladder, f66);
        var f2b = new F2bXyChainSpectrumPi2Inheritance(ladder, f65);
        var absorption = new AbsorptionTheoremClaim(ladder);
        var discriminant = new PolynomialDiscriminantAnchorClaim(new PolynomialFoundationClaim(), qubitAnchor, ladder);
        var carrier = new UniversalCarrierClaim(absorption, ladder, discriminant);
        var clock = new ClockHandLadderClaim(f2b, absorption, carrier);
        var bandEdge = new TopologyBandEdgeClaim(clock, absorption);
        return new CrackedRingExactCurveClaim(bandEdge, f2b);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void Anchor_NamesTheProofTheTwoPagesAndTheRegistry()
    {
        var anchor = BuildClaim().Anchor;
        Assert.Contains("docs/proofs/PROOF_CRACKED_RING_EXACT_CURVE.md", anchor);
        Assert.Contains("experiments/THE_CRACKED_BELL.md", anchor);
        Assert.Contains("experiments/COUPLING_DEFECT_WALK_TIME_STEP.md", anchor);
        Assert.Contains("experiments/THE_COMB_ON_THE_ROAD.md", anchor);
        Assert.Contains("docs/ANALYTICAL_FORMULAS.md", anchor);
        Assert.Contains("compute/MirrorWorld/Crack.cs", anchor);
    }

    [Fact]
    public void Parents_AreHeldAsTypedEdges_AndNullThrows()
    {
        var claim = BuildClaim();
        Assert.NotNull(claim.BandEdge);
        Assert.NotNull(claim.ChainEnd);
        Assert.Throws<ArgumentNullException>(() => new CrackedRingExactCurveClaim(null!, claim.ChainEnd));
        Assert.Throws<ArgumentNullException>(() => new CrackedRingExactCurveClaim(claim.BandEdge, null!));
    }

    [Fact]
    public void TheStatement_KeepsItsFences()
    {
        var name = BuildClaim().Name;
        Assert.Contains("boundary-condition parameter", name);
        Assert.Contains("XY adjacency book only", name);
        Assert.Contains("gamma absent from G", name);
        Assert.DoesNotContain("topolog", name);       // u is not a topology knob; the word is spent twice in the repo
        Assert.DoesNotContain("Heisenberg", name);
    }

    // ---------------------------------------------------------------- the polynomial against an independent determinant

    /// <summary>Fraction-free Bareiss elimination over the integers: det of an integer matrix, exactly.</summary>
    private static BigInteger BareissDeterminant(BigInteger[,] a)
    {
        int n = a.GetLength(0);
        var m = (BigInteger[,])a.Clone();
        BigInteger sign = BigInteger.One, prev = BigInteger.One;
        for (int k = 0; k < n - 1; k++)
        {
            if (m[k, k].IsZero)
            {
                int swap = -1;
                for (int i = k + 1; i < n; i++) if (!m[i, k].IsZero) { swap = i; break; }
                if (swap < 0) return BigInteger.Zero;
                for (int j = 0; j < n; j++) (m[k, j], m[swap, j]) = (m[swap, j], m[k, j]);
                sign = -sign;
            }
            for (int i = k + 1; i < n; i++)
                for (int j = k + 1; j < n; j++)
                    m[i, j] = (m[i, j] * m[k, k] - m[i, k] * m[k, j]) / prev;   // exact division, Bareiss
            prev = m[k, k];
        }
        return sign * m[n - 1, n - 1];
    }

    /// <summary>The integer matrix uDen * (x I - H) for the cracked ring: uDen x on the diagonal, -uDen on the
    /// N-1 chain bonds, -uNum on the two corners.</summary>
    private static BigInteger[,] ScaledPencil(int n, long uNum, long uDen, long x)
    {
        var a = new BigInteger[n, n];
        for (int i = 0; i < n; i++) a[i, i] = (BigInteger)uDen * x;
        for (int i = 0; i < n - 1; i++) { a[i, i + 1] = -uDen; a[i + 1, i] = -uDen; }
        a[0, n - 1] = -uNum;
        a[n - 1, 0] = -uNum;
        return a;
    }

    private static BigInteger Evaluate(BigInteger[] ascending, long x)
    {
        BigInteger acc = BigInteger.Zero;
        for (int i = ascending.Length - 1; i >= 0; i--) acc = acc * x + ascending[i];
        return acc;
    }

    // det(uDen (x I - H)) = uDen^N P(x); the claim's polynomial is uDen^2 P(x). Two routes, compared as integers.
    [Theory]
    [InlineData(3, 1, 2)]
    [InlineData(4, 2, 1)]
    [InlineData(5, 1, 1)]
    [InlineData(6, 3, 4)]
    [InlineData(7, 1, 3)]
    [InlineData(8, 5, 1)]
    [InlineData(9, 0, 1)]
    [InlineData(11, 7, 5)]
    public void RoadPolynomial_MeetsTheBareissDeterminant_AtIntegerPoints(int n, long uNum, long uDen)
    {
        var p = CrackedRingExactCurveClaim.RoadPolynomialScaled(n, uNum, uDen);
        Assert.Equal(n + 1, p.Length);
        Assert.Equal(((BigInteger)uDen) * uDen, p[n]);                 // leading coefficient uDen^2
        BigInteger denPow = BigInteger.Pow(uDen, n - 2);
        foreach (long x in new long[] { -5, -3, -2, -1, 0, 1, 2, 3, 7 })
        {
            BigInteger fromClaim = denPow * Evaluate(p, x);
            BigInteger fromMatrix = BareissDeterminant(ScaledPencil(n, uNum, uDen, x));
            Assert.True(fromClaim == fromMatrix, $"N={n}, u={uNum}/{uDen}, x={x}: polynomial {fromClaim} vs determinant {fromMatrix}");
        }
    }

    // ---------------------------------------------------------------- the two ends

    private static BigInteger[] TwoChebyshevFirstKind(int n)
    {
        // 2 T_n(x/2): q_0 = 2, q_1 = x, q_n = x q_{n-1} - q_{n-2}
        var prev = new BigInteger[] { 2 };
        var cur = new BigInteger[] { 0, 1 };
        for (int m = 2; m <= n; m++)
        {
            var next = new BigInteger[m + 1];
            for (int i = 0; i < cur.Length; i++) next[i + 1] += cur[i];
            for (int i = 0; i < prev.Length; i++) next[i] -= prev[i];
            prev = cur; cur = next;
        }
        return cur;
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(7)]
    [InlineData(12)]
    public void TheRingEnd_IsTwoT_N_MinusTwo_AndTheChainEnd_IsThePathDeterminant(int n)
    {
        var ring = CrackedRingExactCurveClaim.RoadPolynomialScaled(n, 1, 1);
        var twoT = TwoChebyshevFirstKind(n);
        twoT[0] -= 2;
        Assert.Equal(twoT, ring);
        // The chain end against the Bareiss determinant of the PATH pencil (no corners) at integer x: an
        // independent route. (A first version compared it with ChebyshevSecondKindMonic, which is the same code.)
        var chain = CrackedRingExactCurveClaim.RoadPolynomialScaled(n, 0, 1);
        foreach (long x in new long[] { -3, -1, 0, 2, 5 })
            Assert.True(Evaluate(chain, x) == BareissDeterminant(ScaledPencil(n, 0, 1, x)), $"N={n}, x={x}");
    }

    [Fact]
    public void Chebyshev_SmallCases_AreTheTextbookOnes()
    {
        Assert.Equal(new BigInteger[] { 1 }, CrackedRingExactCurveClaim.ChebyshevSecondKindMonic(0));
        Assert.Equal(new BigInteger[] { 0, 1 }, CrackedRingExactCurveClaim.ChebyshevSecondKindMonic(1));
        Assert.Equal(new BigInteger[] { -1, 0, 1 }, CrackedRingExactCurveClaim.ChebyshevSecondKindMonic(2));
        Assert.Equal(new BigInteger[] { 0, -2, 0, 1 }, CrackedRingExactCurveClaim.ChebyshevSecondKindMonic(3));
    }

    // ---------------------------------------------------------------- the band-edge factor forms

    [Theory]
    [InlineData(3, 1, 1)]
    [InlineData(4, 3, 2)]
    [InlineData(5, 3, 2)]
    [InlineData(6, 1, 4)]
    [InlineData(7, 4, 3)]
    [InlineData(9, 5, 4)]
    [InlineData(10, 11, 10)]
    [InlineData(13, 7, 6)]
    public void FactorForms_MeetTheUnfactoredPolynomial_AtTheBandEdges(int n, long uNum, long uDen)
    {
        var p = CrackedRingExactCurveClaim.RoadPolynomialScaled(n, uNum, uDen);
        Assert.True(Evaluate(p, 2) == CrackedRingExactCurveClaim.BandTopValueScaled(n, uNum, uDen));
        Assert.True(Evaluate(p, -2) == CrackedRingExactCurveClaim.BandBottomValueScaled(n, uNum, uDen));
    }

    // ---------------------------------------------------------------- the departure count, against the roots themselves

    /// <summary>All real roots of the scaled polynomial by sign changes on a fine grid plus bisection (the
    /// polynomial is real-rooted with N distinct roots for u != 1, so this route is complete when it finds N).</summary>
    private static double[] RootsByBisection(BigInteger[] p, double lo, double hi, int cells)
    {
        double Eval(double x)
        {
            double acc = 0;
            for (int i = p.Length - 1; i >= 0; i--) acc = acc * x + (double)p[i];
            return acc;
        }
        var roots = new List<double>();
        double h = (hi - lo) / cells, a = lo, fa = Eval(a);
        for (int c = 0; c < cells; c++)
        {
            double b = a + h, fb = Eval(b);
            if (fa == 0) roots.Add(a);
            else if (fa * fb < 0)
            {
                double x0 = a, x1 = b;
                for (int it = 0; it < 200; it++)
                {
                    double mid = 0.5 * (x0 + x1), fm = Eval(mid);
                    if (fm == 0) { x0 = x1 = mid; break; }
                    if (Eval(x0) * fm < 0) x1 = mid; else x0 = mid;
                }
                roots.Add(0.5 * (x0 + x1));
            }
            a = b; fa = fb;
        }
        return roots.ToArray();
    }

    /// <summary>Descartes' rule on q(y) = P(y + 2) (roots above the band) or q(y) = P(-y - 2) (below): the number
    /// of sign changes among the nonzero coefficients. For a real-rooted P that number IS the number of positive
    /// roots of q, and a root exactly at the edge is q(0) = 0, a vanishing constant term, returned separately.
    /// Exact integer arithmetic; shares no code with the claim's factor forms.</summary>
    private static (int Beyond, bool OnEdge) DescartesBeyondEdge(BigInteger[] p, bool above)
    {
        int n = p.Length - 1;
        // q(y) = sum_i p_i (s y + t)^i with (s, t) = (1, 2) above and (-1, -2) below
        long sgn = above ? 1 : -1, shift = above ? 2 : -2;
        var q = new BigInteger[n + 1];
        for (int i = 0; i <= n; i++)
        {
            // (s y + t)^i = sum_j C(i, j) s^j t^(i-j) y^j
            BigInteger binom = BigInteger.One;
            for (int j = 0; j <= i; j++)
            {
                if (j > 0) binom = binom * (i - j + 1) / j;
                BigInteger term = p[i] * binom * BigInteger.Pow(sgn, j) * BigInteger.Pow(shift, i - j);
                q[j] += term;
            }
        }
        bool onEdge = q[0].IsZero;
        int changes = 0, last = 0;
        for (int j = 0; j <= n; j++)
        {
            if (q[j].IsZero) continue;
            int sg = q[j].Sign;
            if (last != 0 && sg != last) changes++;
            last = sg;
        }
        return (changes, onEdge);
    }

    // The count from the factor forms against Descartes' rule on the shifted integer polynomial, over N = 3..14 and
    // u on both sides of 1 and of the odd threshold; three rows hit their own threshold u = (N+1)/(N-1) (N = 3 at
    // 2, N = 5 at 3/2, N = 9 at 5/4) and Descartes sees the edge root as a vanishing constant term. Two exact
    // routes, compared as integers; nothing about a root's position is read in floating point here.
    [Fact]
    public void DepartureCount_MatchesDescartesOnTheShiftedPolynomial()
    {
        for (int n = 3; n <= 14; n++)
        {
            foreach (var (uNum, uDen) in new (long, long)[] { (1, 2), (9, 10), (11, 10), (5, 4), (3, 2), (2, 1), (5, 1) })
            {
                var p = CrackedRingExactCurveClaim.RoadPolynomialScaled(n, uNum, uDen);
                var (above, topEdge) = DescartesBeyondEdge(p, above: true);
                var (below, bottomEdge) = DescartesBeyondEdge(p, above: false);
                Assert.Equal(above, CrackedRingExactCurveClaim.DeparturesAbove(n, uNum, uDen));
                Assert.Equal(below, CrackedRingExactCurveClaim.DeparturesBelow(n, uNum, uDen));
                Assert.Equal((topEdge, bottomEdge), CrackedRingExactCurveClaim.OnTheEdge(n, uNum, uDen));
                Assert.True(above + below <= 2, $"N={n}, u={uNum}/{uDen}: {above + below} departures");
            }
        }
    }

    // And the roots themselves, as a reading beside the exact count: bisection finds N real roots at every row.
    [Fact]
    public void TheRoots_AreAllReal_AndCountN()
    {
        for (int n = 3; n <= 14; n++)
        {
            foreach (var (uNum, uDen) in new (long, long)[] { (1, 2), (11, 10), (2, 1), (5, 1) })
            {
                var p = CrackedRingExactCurveClaim.RoadPolynomialScaled(n, uNum, uDen);
                double u = (double)uNum / uDen;
                var roots = RootsByBisection(p, -(u + 3), u + 3, 200000);
                Assert.True(roots.Length == n, $"N={n}, u={uNum}/{uDen}: found {roots.Length} roots");
            }
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    [InlineData(9)]
    [InlineData(15)]
    public void OddN_AtTheThreshold_TheBottomLevelSitsOnTheEdge_AndIsNotCountedAsDeparted(int n)
    {
        var (num, den) = CrackedRingExactCurveClaim.OddThreshold(n);
        Assert.Equal((long)(n + 1) * den, (long)(n - 1) * num);  // num/den == (N+1)/(N-1)
        Assert.Equal(1, (int)System.Numerics.BigInteger.GreatestCommonDivisor(num, den));   // in lowest terms
        Assert.True(CrackedRingExactCurveClaim.BandBottomValueScaled(n, num, den).IsZero);
        Assert.Equal((false, true), CrackedRingExactCurveClaim.OnTheEdge(n, num, den));
        Assert.Equal(0, CrackedRingExactCurveClaim.DeparturesBelow(n, num, den));
        Assert.Equal(1, CrackedRingExactCurveClaim.DeparturesAbove(n, num, den));
        // one step past: gone; one step before: still inside
        Assert.Equal(1, CrackedRingExactCurveClaim.DeparturesBelow(n, num * 100 + 1, den * 100));
        Assert.Equal(0, CrackedRingExactCurveClaim.DeparturesBelow(n, num * 100 - 1, den * 100));
    }

    [Fact]
    public void EvenN_HasNoThreshold_TheBottomLevelLeavesAtOnce_AndBothEdgesAreMetAtTheRing()
    {
        Assert.Throws<ArgumentException>(() => CrackedRingExactCurveClaim.OddThreshold(8));
        Assert.Equal(1, CrackedRingExactCurveClaim.DeparturesBelow(8, 1000001, 1000000));
        Assert.Equal(0, CrackedRingExactCurveClaim.DeparturesBelow(8, 999999, 1000000));
        Assert.Equal((true, true), CrackedRingExactCurveClaim.OnTheEdge(8, 1, 1));
        Assert.Equal((true, false), CrackedRingExactCurveClaim.OnTheEdge(7, 1, 1));
        Assert.Equal(0, CrackedRingExactCurveClaim.DeparturesAbove(7, 1, 1) + CrackedRingExactCurveClaim.DeparturesBelow(7, 1, 1));
    }

    // ---------------------------------------------------------------- the curve at x = 2 cos k

    [Theory]
    [InlineData(6, 0.7, 0.9)]
    [InlineData(9, 0.25, 2.3)]
    [InlineData(12, 0.99, 1.1)]
    public void TheCurve_IsThePolynomialAtTwoCosK_OverSinK(int n, double u, double k)
    {
        // The identity is exact (gate P3 of the proof script has it symbolically); in doubles the only thing to
        // tolerate is rounding, so the bound is a rounding MODEL: Horner on N + 1 terms accumulates at most
        // (N + 1) roundings of size eps times the magnitude sum of the terms, and the curve side a handful more.
        long den = 100, num = (long)Math.Round(u * 100);
        var p = CrackedRingExactCurveClaim.RoadPolynomialScaled(n, num, den);
        double x = 2 * Math.Cos(k), acc = 0, magnitude = 0;
        for (int i = p.Length - 1; i >= 0; i--) { acc = acc * x + (double)p[i]; magnitude += Math.Abs((double)p[i]) * Math.Pow(Math.Abs(x), i); }
        double lhs = acc / (den * (double)den);                          // P(2 cos k)
        double rhs = CrackedRingExactCurveClaim.Curve(k, n, (double)num / den) / Math.Sin(k);
        double eps = Math.Pow(2, -52);
        double uu = (double)num / den;
        double gTerms = (Math.Abs(1 - uu * uu) + (1 + uu * uu) + 2 * uu) / Math.Abs(Math.Sin(k));   // the magnitude sum of G's own terms over sin k
        double bound = eps * (2 * (n + 1) * magnitude / (den * (double)den) + 8 * gTerms);
        Assert.True(Math.Abs(lhs - rhs) <= bound, $"{lhs} vs {rhs}: residual {Math.Abs(lhs - rhs)} against the rounding model {bound}");
    }

    // ---------------------------------------------------------------- c_m, against roots of the curve at two depths

    private static double PairSplit(int n, double delta, int m)
    {
        double u = 1 - delta, km = 2 * Math.PI * m / n;
        // the two roots of G nearest k_m, by bisection in a window that contains both and no other (delta small)
        double w = 6.0 * delta / n + 1e-9;
        var ks = new List<double>();
        int cells = 4000;
        double a = km - w, fa = CrackedRingExactCurveClaim.Curve(a, n, u), h = 2 * w / cells;
        for (int c = 0; c < cells; c++)
        {
            double b = a + h, fb = CrackedRingExactCurveClaim.Curve(b, n, u);
            if (fa * fb < 0)
            {
                double x0 = a, x1 = b;
                for (int it = 0; it < 100; it++)
                {
                    double mid = 0.5 * (x0 + x1);
                    if (CrackedRingExactCurveClaim.Curve(x0, n, u) * CrackedRingExactCurveClaim.Curve(mid, n, u) < 0) x1 = mid; else x0 = mid;
                }
                ks.Add(0.5 * (x0 + x1));
            }
            a = b; fa = fb;
        }
        Assert.True(ks.Count == 2, $"N={n}, m={m}, delta={delta}: {ks.Count} roots in the pair window");
        return 2 * Math.Cos(ks[0]) - 2 * Math.Cos(ks[1]);
    }

    // (Delta E/(4 delta/N) - 1)/delta -> c_m, with the residual falling one decade per decade of delta. The
    // depths are 1e-2 and 1e-3: at 1e-4 the double-precision root of the curve is no longer good enough to read
    // a residual of c_3 delta^2 ~ 1e-7 (measured: the residual stopped falling there), so that depth would test
    // the arithmetic and not the law. The window (8, 12.5) is symmetric about 10 on a log scale (a factor 1.25
    // either way).
    [Theory]
    [InlineData(8, 1)]
    [InlineData(12, 1)]
    [InlineData(12, 2)]
    [InlineData(16, 2)]
    [InlineData(20, 3)]
    public void SplitCorrection_IsTheSecondOrderOfTheSplit_ByADecadeLaw(int n, int m)
    {
        double c = CrackedRingExactCurveClaim.SplitCorrection(n, m);
        double Read(double delta) => (Math.Abs(PairSplit(n, delta, m)) / (4 * delta / n) - 1) / delta;
        double r1 = Math.Abs(Read(1e-2) - c), r2 = Math.Abs(Read(1e-3) - c);
        double ratio = r1 / r2;
        Assert.True(ratio > 8 && ratio < 12.5, $"decade ratio {ratio} (r1={r1}, r2={r2}, c={c})");
    }

    // N sin^2(2 pi/N) falls with N from N = 5 on and crosses 2 at N = 19.03, so c_1 is positive (band centre side) at N = 19 and
    // negative (band edge side) at N = 20. A first version of this test had the two signs the other way round.
    [Fact]
    public void SplitCorrection_ChangesSign_ForMEqualsOne_BetweenN19AndN20()
    {
        Assert.True(CrackedRingExactCurveClaim.SplitCorrection(19, 1) > 0);
        Assert.True(CrackedRingExactCurveClaim.SplitCorrection(20, 1) < 0);
        Assert.True(CrackedRingExactCurveClaim.SplitCorrection(6, 1) > 0);        // band centre: c = 1/2 - 1/(6 sin^2 60deg) = 1/2 - 2/9
        Assert.Throws<ArgumentOutOfRangeException>(() => CrackedRingExactCurveClaim.SplitCorrection(8, 4));
    }

    // ---------------------------------------------------------------- the chain-end velocity, against the roots' finite difference

    // The one-sided difference (x_k(u) - x_k(0))/u against the closed form: its residual is the second-order
    // term, (u/2) d^2x_k/du^2 + O(u^2), so it falls one decade when u does, or two where that coefficient
    // vanishes (the middle level of an odd chain, k = (N+1)/2, does that). The depths are 1e-3 and 1e-4: at
    // 1e-5 the middle level's difference quotient, whose residual is O(u^2), sits on the bisection's own rounding
    // floor (1e-16 in x times 1e5; the other levels' O(u) residuals are still above it) and a first version of
    // this test read that floor as a failed law.
    [Theory]
    [InlineData(5)]
    [InlineData(8)]
    [InlineData(11)]
    public void ChainEndVelocity_MatchesTheFiniteDifferenceOfTheRoots_ByADecadeLaw(int n)
    {
        double[] Residuals(long den)
        {
            var plus = RootsByBisection(CrackedRingExactCurveClaim.RoadPolynomialScaled(n, 1, den), -3, 3, 200000);
            Assert.Equal(n, plus.Length);
            Array.Sort(plus); Array.Reverse(plus);
            var r = new double[n];
            for (int k = 1; k <= n; k++)
            {
                double e0 = 2 * Math.Cos(k * Math.PI / (n + 1));
                r[k - 1] = (plus[k - 1] - e0) * den - CrackedRingExactCurveClaim.ChainEndVelocity(n, k);
            }
            return r;
        }
        var coarse = Residuals(1000);       // u = 1e-3
        var fine = Residuals(10000);        // u = 1e-4
        for (int k = 0; k < n; k++)
        {
            double ratio = Math.Abs(coarse[k]) / Math.Abs(fine[k]);
            bool one = Math.Abs(ratio / 10.0 - 1.0) <= 0.15, two = Math.Abs(ratio / 100.0 - 1.0) <= 0.15;
            Assert.True(one || two, $"N={n}, k={k + 1}: residual ratio {ratio} (coarse {coarse[k]}, fine {fine[k]}) is neither one nor two decades per decade");
        }
    }

    [Fact]
    public void ChainEndVelocity_MagnitudesAreF65sRates_AtN5()
    {
        // F65's verified row at N = 5: alpha/gamma_0 = 1/6, 1/2, 2/3, 1/2, 1/6, with the alternating sign.
        // The only error is the rounding of Math.Sin, one square and one division: a handful of ulps of the value.
        double[] expected = { 1.0 / 6, -1.0 / 2, 2.0 / 3, -1.0 / 2, 1.0 / 6 };
        for (int k = 1; k <= 5; k++)
            Assert.True(Math.Abs(CrackedRingExactCurveClaim.ChainEndVelocity(5, k) - expected[k - 1]) <= 8 * Math.Pow(2, -52) * Math.Abs(expected[k - 1]));
    }

    [Fact]
    public void Guards_RefuseWhatTheObjectDoesNotSpeakAbout()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => CrackedRingExactCurveClaim.RoadPolynomialScaled(2, 1, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => CrackedRingExactCurveClaim.RoadPolynomialScaled(5, 1, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => CrackedRingExactCurveClaim.RoadPolynomialScaled(5, -1, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => CrackedRingExactCurveClaim.ChainEndVelocity(5, 6));
    }
}
