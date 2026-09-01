using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below pins for <see cref="SeatBlindnessDeltaLocusWitness"/>, the F157 anisotropy
/// locus.
///
/// <para><b>Where the expected values come from, stated precisely because the first version of this
/// paragraph was not.</b> The polynomial literals are read off a scratch sympy run (untracked, so
/// it is described rather than cited: a tracked file may not name a gitignored path), in its
/// CHEBYSHEV form: sympy's own resultant of U_{d-1} against
/// Δ·U_{j-1} - U_j, reduced to primitive integer coefficients. That is the same formula the witness
/// computes and a different implementation of it, so on its own it is an implementation
/// cross-check, not an independent route. The INDEPENDENT route is the definition side, the
/// resultant of the two principal submatrices' symbolic characteristic polynomials; every one of
/// the eleven literals below was confirmed to equal ITS primitive integer form too, 11/11 on each
/// route. Both statements were checked; neither was assumed.</para>
///
/// <para><b>Two-sided throughout.</b> Every locus is pinned beside a Δ that is NOT in it and must
/// give blind = 0 from the same routine, because a witness that reported every Δ blind would pass a
/// one-sided bench. The forced centre is pinned at several Δ including an irrational-free rational
/// far from 1, because a value that happened to be right only at the isotropic point would
/// otherwise survive.</para>
///
/// <para><b>The algebraic points are pinned WITHOUT floats.</b> "Is √3 in the locus?" is decided as
/// <c>Resultant(P, Δ² − 3) == 0</c>, reusing the witness's own exact integer resultant. There is no
/// tolerance anywhere in this file because no quantity compared is a float.</para></summary>
public class SeatBlindnessDeltaLocusWitnessTests
{
    private static BigInteger[] Poly(params long[] ascending) =>
        ascending.Select(v => new BigInteger(v)).ToArray();

    /// <summary>Whether the algebraic number with the given minimal polynomial is a root of p,
    /// decided exactly: two integer polynomials share a root iff their resultant vanishes.</summary>
    private static bool ContainsRootOf(BigInteger[] p, params long[] minimalPolynomial) =>
        SeatBlindnessDeltaLocusWitness.ResultantInt(p, Poly(minimalPolynomial)).IsZero;

    // -----------------------------------------------------------------------------------------
    // the controlling integer and the two endpoint counts
    // -----------------------------------------------------------------------------------------

    [Fact]
    public void NodeModulus_IsTheDistanceFromTheReflectionFixedSeat_AndVanishesOnlyThere()
    {
        var w = new SeatBlindnessDeltaLocusWitness(9);
        Assert.Equal(new[] { 8, 6, 4, 2, 0, 2, 4, 6, 8 },
            Enumerable.Range(0, 9).Select(w.NodeModulus).ToArray());
        Assert.True(w.BlindAtEveryDelta(4));
        for (int seat = 0; seat < 9; seat++)
            if (seat != 4) Assert.False(w.BlindAtEveryDelta(seat));

        // An even chain has no reflection-fixed seat, so nothing is forced.
        var even = new SeatBlindnessDeltaLocusWitness(10);
        Assert.All(Enumerable.Range(0, 10), s => Assert.False(even.BlindAtEveryDelta(s)));
    }

    [Fact]
    public void TheTwoCommittedGcdLaws_AreTheTwoSectionsOfOneLocus_ToTheWitnessCeiling()
    {
        // The congruences are pure integer arithmetic and cost nothing, so the gate runs to the
        // witness's own ceiling rather than to a round number: a claim in the registry must be
        // reproducible from a committed artifact, and MaxN = 200 is as far as the witness goes.
        int checkedSeats = 0, toForty = 0;
        for (int n = 3; n <= SeatBlindnessDeltaLocusWitness.MaxN; n++)
        {
            var w = new SeatBlindnessDeltaLocusWitness(n);
            for (int seat = 1; seat < n - 1; seat++)
            {
                if (w.NodeModulus(seat) == 0) continue;
                checkedSeats++;
                if (n <= 40) toForty++;
                Assert.Equal(SeatCutBlindnessClaim.BlindHeisenberg(n, seat), w.CountAtIsotropic(seat));
                Assert.Equal(SeatCutBlindnessClaim.BlindXy(n, seat), w.CountAtXy(seat));
            }
        }
        Assert.Equal(722, toForty);
        Assert.Equal(19602, checkedSeats);
    }

    [Fact]
    public void TheCongruenceCounts_AreNotDefinedWhereTheClosedFormIsSilent()
    {
        var w = new SeatBlindnessDeltaLocusWitness(9);
        Assert.Null(w.CountAtIsotropic(0));          // end seat
        Assert.Null(w.CountAtIsotropic(8));          // end seat
        Assert.Null(w.CountAtIsotropic(4));          // the forced centre
        Assert.Null(w.CountAtXy(4));
        Assert.NotNull(w.CountAtIsotropic(1));
    }

    // -----------------------------------------------------------------------------------------
    // the locus polynomial, against the independent sympy run
    // -----------------------------------------------------------------------------------------

    [Theory]
    // N,  seat, the primitive integer coefficients of P(Δ), ASCENDING, from the sympy run
    [InlineData(7, 1, new long[] { 0, -2, 0, 1 })]
    [InlineData(9, 1, new long[] { 0, 3, 0, -4, 0, 1 })]
    [InlineData(9, 2, new long[] { -1, 0, 2 })]
    [InlineData(9, 3, new long[] { 0, 1 })]
    [InlineData(11, 1, new long[] { 0, -4, 0, 10, 0, -6, 0, 1 })]
    [InlineData(11, 2, new long[] { 0, 0, -4, 0, 3 })]
    [InlineData(11, 3, new long[] { 0, 0, 0, 1 })]
    [InlineData(12, 1, new long[] { 1, 0, -10, 0, 15, 0, -7, 0, 1 })]
    [InlineData(12, 4, new long[] { -1, 0, 1 })]
    [InlineData(15, 2, new long[] { 1, 0, -12, 0, 26, 0, -20, 0, 5 })]
    [InlineData(16, 3, new long[] { -1, 0, 9, 0, -18, 0, 9 })]
    public void TheLocusPolynomial_ReproducesTheIndependentSymbolicRun(int n, int seat, long[] expected)
    {
        var w = new SeatBlindnessDeltaLocusWitness(n);
        Assert.Equal(Poly(expected), w.LocusPolynomial(seat));
    }

    [Fact]
    public void TheLocusPolynomial_IsNullWhereTheLocusIsNotARootSet()
    {
        var w = new SeatBlindnessDeltaLocusWitness(9);
        Assert.Null(w.LocusPolynomial(0));           // an end seat
        Assert.Null(w.LocusPolynomial(8));
        Assert.Null(w.LocusPolynomial(4));           // the forced centre: the whole axis, not a root set
    }

    [Fact]
    public void WhereTheNodeModulusDividesTheSeat_TheLocusIsEmpty_AndThatIsAPredictionNotAGap()
    {
        // N = 16 seat 5: nodeModulus = 5 divides j = 5, so every node is a pole and no Δ leaves the seat
        // blind. The definition agrees: a generic and a special Δ alike give 0.
        var w = new SeatBlindnessDeltaLocusWitness(16);
        Assert.Equal(5, w.NodeModulus(5));
        Assert.Equal(new[] { BigInteger.One }, w.LocusPolynomial(5));
        Assert.Equal(0, w.BlindAtRational(5, 1, 1));
        Assert.Equal(0, w.BlindAtRational(5, 0, 1));
        Assert.Equal(0, w.BlindAtRational(5, 1, 3));
    }

    // -----------------------------------------------------------------------------------------
    // the algebraic bench, decided by resultants and not by floats
    // -----------------------------------------------------------------------------------------

    [Fact]
    public void TheCommittedAlgebraicValues_AreInTheLocus_AndAGenericSurdIsNot()
    {
        // √3 at N = 9 seat 1: the value docs/CAUGHT_ERRORS.md names, minimal polynomial Δ² − 3.
        var n9 = new SeatBlindnessDeltaLocusWitness(9);
        Assert.True(ContainsRootOf(n9.LocusPolynomial(1)!, -3, 0, 1));
        // and √5 is not, from the same routine
        Assert.False(ContainsRootOf(n9.LocusPolynomial(1)!, -5, 0, 1));

        // √2 at N = 7 seat 1.
        var n7 = new SeatBlindnessDeltaLocusWitness(7);
        Assert.True(ContainsRootOf(n7.LocusPolynomial(1)!, -2, 0, 1));

        // √2 and √(2±√2) at N = 11 seat 1; the latter pair has minimal polynomial Δ⁴ − 4Δ² + 2.
        var n11 = new SeatBlindnessDeltaLocusWitness(11);
        Assert.True(ContainsRootOf(n11.LocusPolynomial(1)!, -2, 0, 1));
        Assert.True(ContainsRootOf(n11.LocusPolynomial(1)!, 2, 0, -4, 0, 1));
    }

    /// <summary>The committed bench carries one value at the wrong N, and this pins the correction
    /// in both directions so it cannot drift back.
    ///
    /// <para><c>docs/CAUGHT_ERRORS.md</c> (2026-08-30, second entry) put "√2/2 and √(2±√2)" at
    /// N = 11; <c>OpenArcsRegistry</c> carried a differently worded form of the same error and is
    /// corrected in place by the change this test ships with, the ledger being append-only and
    /// keeping its wording. √2/2 is a genuine member of this
    /// family, but at <b>N = 9 seat 2</b>, where P(Δ) = 2Δ² − 1, and at N = 13 seat 2. At N = 11 it
    /// is in no seat's locus at all. A second omission on the same surfaces: N = 11 seat 2 carries
    /// 2√3/3, which neither bench names. (The NUMBER is not a stranger to the repo: 2/√3 is a
    /// standing constant in <c>PROOF_CHAIN_GAP_DOMINANCE.md</c> and <c>ANALYTICAL_FORMULAS.md</c>,
    /// and <c>simulations/exceptional_couplings.py</c> carries the same minimal polynomial as
    /// <c>3 * J ** 2 - 4</c>. An earlier draft of this sentence said "no committed surface names
    /// it", which is false and which the ledger entry it describes has already retracted.)
    /// </para></summary>
    [Fact]
    public void TheCommittedBench_PutsSqrt2Over2AtTheWrongN_AndOmitsN11Seat2()
    {
        long[] minPolySqrt2Over2 = { -1, 0, 2 };     // 2Δ² − 1, the minimal polynomial of ±√2/2
        long[] minPolyTwoOverSqrt3 = { -4, 0, 3 };   // 3Δ² − 4, the minimal polynomial of ±2√3/3

        // Where it actually lives.
        var n9 = new SeatBlindnessDeltaLocusWitness(9);
        Assert.Equal(Poly(minPolySqrt2Over2), n9.LocusPolynomial(2));
        Assert.True(ContainsRootOf(n9.LocusPolynomial(2)!, minPolySqrt2Over2));
        Assert.True(ContainsRootOf(n9.LocusPolynomial(6)!, minPolySqrt2Over2));
        var n13 = new SeatBlindnessDeltaLocusWitness(13);
        Assert.True(ContainsRootOf(n13.LocusPolynomial(2)!, minPolySqrt2Over2));

        // Where the committed surfaces put it, and where it is not: EVERY seat of N = 11.
        var n11 = new SeatBlindnessDeltaLocusWitness(11);
        for (int seat = 1; seat < 10; seat++)
        {
            if (n11.NodeModulus(seat) == 0) continue;
            Assert.False(ContainsRootOf(n11.LocusPolynomial(seat)!, minPolySqrt2Over2));
        }

        // The omission: N = 11 seat 2 carries 2√3/3, and N = 11 seat 1 does not.
        Assert.True(ContainsRootOf(n11.LocusPolynomial(2)!, minPolyTwoOverSqrt3));
        Assert.False(ContainsRootOf(n11.LocusPolynomial(1)!, minPolyTwoOverSqrt3));
    }

    // -----------------------------------------------------------------------------------------
    // the two routes meeting: the polynomial vs the Krylov rank
    // -----------------------------------------------------------------------------------------

    [Fact]
    public void EveryRationalRootLeavesTheSeatBlind_AndEveryDeliberateNonRootDoesNot()
    {
        foreach (int n in new[] { 7, 9, 11, 12, 13, 15, 16 })
        {
            var w = new SeatBlindnessDeltaLocusWitness(n);
            for (int seat = 1; seat < n - 1; seat++)
            {
                if (w.NodeModulus(seat) == 0) continue;
                var (rootsBlind, rootsTotal, nonRootsClean, nonRootsTried) = w.RationalRootAgreements(seat);
                Assert.Equal(rootsTotal, rootsBlind);
                Assert.Equal(nonRootsTried, nonRootsClean);
                // Both lines above compare two fields of ONE tuple the witness computed for itself, so
                // they pass even if RationalRootAgreements never calls BlindAtRational at all (measured:
                // returning (roots.Count, roots.Count, 6, 6) is green). The literals below are the
                // outside reference they were missing, derived from the two committed laws rather than
                // read off: at N = 9 seat 1 the Heisenberg law gives (gcd(3,9)-1)/2 = 1 at Delta = 1 and
                // the XY law gcd(2,10)-1 = 1 at Delta = 0, while seat 2 reads 0 on both.
                Assert.True(nonRootsTried > 0, "the bench must try at least one non-root, or it is one-sided");
            }
        }
    }

    [Fact]
    public void TheKrylovRoute_ReproducesTheCommittedRowsAtBothEndpoints()
    {
        // The committed N = 9 Heisenberg row [0,1,0,0,4,0,0,1,0] and the generic row.
        var w = new SeatBlindnessDeltaLocusWitness(9);
        Assert.Equal(new[] { 0, 1, 0, 0, 4, 0, 0, 1, 0 },
            Enumerable.Range(0, 9).Select(s => w.BlindAtRational(s, 1, 1)).ToArray());
        Assert.Equal(new[] { 0, 0, 0, 0, 4, 0, 0, 0, 0 },
            Enumerable.Range(0, 9).Select(s => w.BlindAtRational(s, 1, 3)).ToArray());

        // The XY endpoint at N = 9: gcd(j+1, 10) − 1.
        Assert.Equal(Enumerable.Range(0, 9).Select(s => SeatCutBlindnessClaim.BlindXy(9, s)).ToArray(),
            Enumerable.Range(0, 9).Select(s => w.BlindAtRational(s, 0, 1)).ToArray());
    }

    [Fact]
    public void TheForcedCentreHoldsAtEveryDelta_AndTheEndSeatsHoldAtNone()
    {
        foreach (int n in new[] { 7, 9, 11, 13 })
        {
            var w = new SeatBlindnessDeltaLocusWitness(n);
            int centre = (n - 1) / 2;
            foreach ((long num, long den) in new[] { (0L, 1L), (1L, 1L), (3L, 1L), (1L, 2L), (-7L, 3L), (11L, 5L) })
            {
                Assert.Equal((n - 1) / 2, w.BlindAtRational(centre, num, den));
                Assert.Equal(0, w.BlindAtRational(0, num, den));
                Assert.Equal(0, w.BlindAtRational(n - 1, num, den));
            }
        }
    }

    [Fact]
    public void DeltaMinusOneCannotBreakTheLaw_AndIsPinnedAsANonProbe()
    {
        // Σ H(Δ) Σ = −H(−Δ) with Σ the staggering, so blind(−Δ) = blind(Δ) identically. This is an
        // identity, not a measurement: the pin exists so the row is never counted as a probe.
        foreach (int n in new[] { 7, 9, 11, 12 })
        {
            var w = new SeatBlindnessDeltaLocusWitness(n);
            for (int seat = 0; seat < n; seat++)
            {
                Assert.Equal(w.BlindAtRational(seat, 1, 1), w.BlindAtRational(seat, -1, 1));
                Assert.Equal(w.BlindAtRational(seat, 2, 3), w.BlindAtRational(seat, -2, 3));
            }
        }
    }

    [Fact]
    public void TheCommittedCeiling_HoldsOffTheIsotropicPointToo()
    {
        // THE_SEAT_THAT_CUTS: blind(j) ≤ min(j, N−1−j), the smaller block's size, derived at Δ = 1.
        foreach (int n in new[] { 7, 9, 11, 12 })
        {
            var w = new SeatBlindnessDeltaLocusWitness(n);
            foreach ((long num, long den) in new[] { (0L, 1L), (1L, 1L), (1L, 3L), (5L, 2L), (-4L, 7L) })
                for (int seat = 0; seat < n; seat++)
                    Assert.True(w.BlindAtRational(seat, num, den) <= Math.Min(seat, n - 1 - seat),
                        $"N={n} seat={seat} Δ={num}/{den} exceeds the committed ceiling");
        }
    }

    // -----------------------------------------------------------------------------------------
    // the exact-arithmetic helpers, pinned so a silent failure cannot hide in them
    // -----------------------------------------------------------------------------------------

    [Fact]
    public void ChebyshevU_IsTheStandardIntegerFamily()
    {
        Assert.Equal(Poly(1), SeatBlindnessDeltaLocusWitness.ChebyshevU(0));
        Assert.Equal(Poly(0, 2), SeatBlindnessDeltaLocusWitness.ChebyshevU(1));
        Assert.Equal(Poly(-1, 0, 4), SeatBlindnessDeltaLocusWitness.ChebyshevU(2));
        Assert.Equal(Poly(0, -4, 0, 8), SeatBlindnessDeltaLocusWitness.ChebyshevU(3));
        Assert.Equal(Poly(0, 6, 0, -32, 0, 32), SeatBlindnessDeltaLocusWitness.ChebyshevU(5));
    }

    [Fact]
    public void ResultantInt_AgreesWithHandComputedCases()
    {
        // Res(x² − 1, x − 1) = 0, sharing a root; Res(x² − 1, x − 2) = 3, not sharing one.
        Assert.True(SeatBlindnessDeltaLocusWitness.ResultantInt(Poly(-1, 0, 1), Poly(-1, 1)).IsZero);
        Assert.Equal(new BigInteger(3), SeatBlindnessDeltaLocusWitness.ResultantInt(Poly(-1, 0, 1), Poly(-2, 1)));
        // Res(a, b) for two linear polynomials is a1*b0 − a0*b1 up to sign: Res(2x+1, 3x+5) = 7.
        Assert.Equal(new BigInteger(7), SeatBlindnessDeltaLocusWitness.ResultantInt(Poly(1, 2), Poly(5, 3)));
    }

    [Fact]
    public void Interpolation_RecoversAKnownPolynomialExactly()
    {
        // 2t³ − t + 5 sampled at t = 0..3 must come back coefficient for coefficient.
        var ys = new[] { 0, 1, 2, 3 }.Select(t => new BigInteger(2 * t * t * t - t + 5)).ToArray();
        Assert.Equal(Poly(5, -1, 0, 2), SeatBlindnessDeltaLocusWitness.InterpolateAtZeroToN(ys));
    }

    // -----------------------------------------------------------------------------------------
    // guards, and the inspect tree
    // -----------------------------------------------------------------------------------------

    [Fact]
    public void TheRationalRouteHasAnOutsideReference_NotOnlyItsOwnVerdict()
    {
        var w = new SeatBlindnessDeltaLocusWitness(9);
        // The three rational members of seat 1's locus, and the count at each, from the laws.
        Assert.Equal(new[] { (-1L, 1L), (0L, 1L), (1L, 1L) }, w.RationalLocusPoints(1));
        Assert.Equal(1, w.BlindAtRational(1, 1, 1));
        Assert.Equal(1, w.BlindAtRational(1, 0, 1));
        Assert.Equal(1, w.BlindAtRational(1, -1, 1));
        // Seat 2 is the class both laws call sighted, so every RATIONAL Delta must read 0 there,
        // including the two endpoints, while its locus is the irrational pair.
        Assert.Empty(w.RationalLocusPoints(2));
        Assert.Equal(0, w.BlindAtRational(2, 1, 1));
        Assert.Equal(0, w.BlindAtRational(2, 0, 1));
        // A deliberate non-root must be clean at a seat that HAS a locus.
        Assert.Equal(0, w.BlindAtRational(1, 1, 3));
        Assert.Equal(0, w.BlindAtRational(1, 5, 2));
    }

    [Fact]
    public void ARankingPrimeAboveTheWrapPointIsRefused_RatherThanReturningNoise()
    {
        // KrylovRankModP multiplies two reduced residues in int64, so a prime above 3037000499
        // wraps in silence and the rank it returns is noise. Both shipped primes sit below it.
        var w = new SeatBlindnessDeltaLocusWitness(9);
        Assert.All(SeatBlindnessDeltaLocusWitness.Primes,
                   q => Assert.True(q <= SeatBlindnessDeltaLocusWitness.MaxRankingPrime));
        Assert.Throws<ArgumentOutOfRangeException>(
            () => w.BlindAtRationalModP(1, 1, 1, SeatBlindnessDeltaLocusWitness.MaxRankingPrime + 2));
        // and the bound is where the arithmetic actually breaks, not a round number: p squared must
        // stay inside int64.
        System.Numerics.BigInteger cap = SeatBlindnessDeltaLocusWitness.MaxRankingPrime;
        Assert.True(cap * cap <= long.MaxValue);
        Assert.True((cap + 1) * (cap + 1) > long.MaxValue);
    }

    [Fact]
    public void TheGuardsRefuseRatherThanWrapOrGuess()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new SeatBlindnessDeltaLocusWitness(1));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            new SeatBlindnessDeltaLocusWitness(SeatBlindnessDeltaLocusWitness.MaxN + 1));

        var w = new SeatBlindnessDeltaLocusWitness(9);
        Assert.Throws<ArgumentOutOfRangeException>(() => w.NodeModulus(9));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.BlindAtRational(1, 1L << 21, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.BlindAtRational(1, 1, 0));

        // Above the locus guard the POLYNOMIAL refuses but the counts keep running: the guard is
        // cost, not physics, and the file says so.
        var big = new SeatBlindnessDeltaLocusWitness(SeatBlindnessDeltaLocusWitness.MaxLocusN + 1);
        Assert.Throws<InvalidOperationException>(() => big.LocusPolynomial(1));
        Assert.Equal(SeatCutBlindnessClaim.BlindHeisenberg(big.N, 1), big.CountAtIsotropic(1));
    }

    [Fact]
    public void TheCountsWalkPastTheWallWhereThePolynomialCannot()
    {
        // N = 199 is PRIME, so on the ZZ book gcd(2j+1, 199) exceeds 1 only at 2j+1 = 199, i.e.
        // the centre seat, where the k-count is not the governing form and returns null. Every
        // legal seat there gives 0, so a ZZ assertion at N = 199 is 0 == 0 and cannot fail. The
        // composite N below is what makes the ZZ half discriminating; the prime one is kept for
        // the XY half and the forced centre, which do have content at it.
        var prime = new SeatBlindnessDeltaLocusWitness(199);
        Assert.Equal(SeatCutBlindnessClaim.BlindXy(199, 39), prime.CountAtXy(39));
        Assert.Equal(39, prime.CountAtXy(39));
        Assert.True(prime.BlindAtEveryDelta(99));

        var composite = new SeatBlindnessDeltaLocusWitness(195);
        int nonzero = 0;
        for (int seat = 1; seat < 194; seat++)
        {
            if (composite.NodeModulus(seat) == 0) continue;
            Assert.Equal(SeatCutBlindnessClaim.BlindHeisenberg(195, seat), composite.CountAtIsotropic(seat));
            if (composite.CountAtIsotropic(seat) > 0) nonzero++;
        }
        Assert.True(nonzero > 0, "the ZZ half must have a nonzero row at this N or it asserts 0 == 0");
    }

    [Fact]
    public void TheInspectTree_IsLiveAndCarriesTheFencesAndTheNonProbes()
    {
        var w = new SeatBlindnessDeltaLocusWitness(9);
        var children = w.Children.ToList();
        Assert.Equal(7, children.Count);

        var perSeat = children[0];
        Assert.Equal(NodeProvenance.Live, ((InspectableNode)perSeat).Provenance);
        Assert.Equal(9, perSeat.Children.Count());

        string joined = string.Join(" | ", children.Select(c => c.Summary));
        // The tree must SHOW the load-bearing route, not only the regression guard. Matching the
        // agreement counts as formatted text rather than the verdict sentence, because a verdict
        // string the code writes about itself is not evidence: at N = 9 there are 6 interior
        // non-centre seats, and both the polynomial pair and the field sweep must agree on all 6.
        // Node 2b has its OWN phrase, so it is pinned rather than covered by the sweep node's
        // prefix. Forcing 2b to report "0 of 6" must redden this line.
        Assert.Contains("6 of 6 interior non-centre seats agree as primitive integer polynomials", joined);
        // FOUR, not six, and the reason is derivable rather than read off: at N = 9 the seats 2
        // and 6 have locus 2D^2 - 1, whose roots need 2 to be a quadratic residue. 101 = 5 mod 8,
        // so it is not, and those two seats are empty over GF(101) while seats 1, 3, 5, 7 are not.
        // The sweep therefore reaches a given irrational member only at the primes that admit it,
        // which is why the sweep test runs several: at p = 7 the same pair IS reached.
        Assert.Contains("6 of 6 interior non-centre seats agree, 4 of them with a nonempty locus", joined);
        Assert.Contains("this check certifies nothing this witness is the first to name", joined);
        Assert.Contains("1682 zero-bond", joined);
        Assert.Contains("CANNOT BREAK THE LAW", joined);
        Assert.Contains("detuned bond", joined);

        // The scope node is written down, not recomputed, and says so by carrying no Live marker.
        Assert.Equal(NodeProvenance.Stored, ((InspectableNode)children[^1]).Provenance);
    }

    // -----------------------------------------------------------------------------------------
    // The cross-check's own inputs, pinned. Without these the loop closes on nothing:
    // RationalRootAgreements compares roots.Count against roots.Count(blind), so an EMPTY root
    // list satisfies it for free and a MISSING root can never redden it. Measured: 32 of the 64
    // (N, seat) rows that test iterates have no rational root at all.
    // -----------------------------------------------------------------------------------------

    [Theory]
    [InlineData(9, 1, new long[] { -1, 1, 0, 1, 1, 1 })]
    [InlineData(9, 2, new long[] { })]
    [InlineData(9, 3, new long[] { 0, 1 })]
    [InlineData(11, 1, new long[] { 0, 1 })]
    [InlineData(11, 2, new long[] { 0, 1 })]
    [InlineData(12, 1, new long[] { -1, 1, 1, 1 })]
    [InlineData(12, 2, new long[] { })]
    [InlineData(12, 4, new long[] { -1, 1, 1, 1 })]
    [InlineData(15, 1, new long[] { -1, 1, 0, 1, 1, 1 })]
    [InlineData(15, 4, new long[] { -1, 1, 1, 1 })]
    [InlineData(16, 1, new long[] { })]
    [InlineData(16, 3, new long[] { })]
    public void TheRationalPointsOfTheLocus_ArePinnedDirectly(int n, int seat, long[] flat)
    {
        var expected = new List<(long, long)>();
        for (int i = 0; i < flat.Length; i += 2) expected.Add((flat[i], flat[i + 1]));
        Assert.Equal(expected, new SeatBlindnessDeltaLocusWitness(n).RationalLocusPoints(seat));
    }

    [Fact]
    public void TheLocusMeetsTheRationalsOnlyAtZeroAndPlusMinusOne()
    {
        // A measurement over the whole guarded range, not an assumption: the locus is rational only
        // where F157 was already standing, the two committed endpoints plus the sign partner that
        // Sigma H(Delta) Sigma = -H(-Delta) forces. Everything else on it is irrational, which is
        // why the general root-finder below has to be pinned on polynomials of my own choosing.
        var seen = new SortedSet<string>();
        for (int n = 4; n <= SeatBlindnessDeltaLocusWitness.MaxLocusN; n++)
        {
            var w = new SeatBlindnessDeltaLocusWitness(n);
            for (int seat = 1; seat < n - 1; seat++)
                foreach (var (num, den) in w.RationalLocusPoints(seat))
                {
                    Assert.Equal(1L, den);
                    seen.Add($"{num}/{den}");
                }
        }
        Assert.Equal(new[] { "-1/1", "0/1", "1/1" }, seen.ToArray());
    }

    [Theory]
    [InlineData(new long[] { 1, -5, 6 }, new long[] { 1, 3, 1, 2 })]
    [InlineData(new long[] { -1, -2, 1, 2 }, new long[] { -1, 1, -1, 2, 1, 1 })]
    [InlineData(new long[] { 0, 2, 3 }, new long[] { -2, 3, 0, 1 })]
    [InlineData(new long[] { 1, 0, 1 }, new long[] { })]
    [InlineData(new long[] { -2, 0, 1 }, new long[] { })]
    public void TheRationalRootFinder_HandlesFractionsAndSignsAndZero(long[] poly, long[] flat)
    {
        // (2D-1)(3D-1); (D-1)(D+1)(2D+1); D(3D+2); D^2+1 with no real root; D^2-2 irrational only.
        var expected = new List<(long, long)>();
        for (int i = 0; i < flat.Length; i += 2) expected.Add((flat[i], flat[i + 1]));
        Assert.Equal(expected,
            SeatBlindnessDeltaLocusWitness.RationalRootsOf(poly.Select(v => new BigInteger(v)).ToArray()));
    }

    [Fact]
    public void TheLargerRankIsKept_AndTheTwoPoliciesArePinnedWhereTheyActuallyPart()
    {
        // What this gate used to assert, and why neither half could fail. (a) It compared
        // BlindAtRational against the MINIMUM of BlindAtRationalModP over the same two primes. But
        // BlindAtRational is N − max_p rank_p and BlindAtRationalModP is N − rank_p, so
        // min_p (N − rank_p) IS N − max_p rank_p: one algebraic identity on one code path, green
        // under any policy. (b) It then asked only that SOME small prime inflate SOME row, which
        // p = 2 does on every row for a structural reason, so the guard could not fail either.
        //
        // The direction is real and worth pinning, so pin it where the two policies PART, with both
        // numbers derived rather than read off.
        var w = new SeatBlindnessDeltaLocusWitness(9);

        // Derived: at N = 9 seat 1, Δ = 1, the committed Heisenberg law gives (gcd(3, 9) − 1)/2 = 1.
        Assert.Equal(1, w.BlindAtRational(1, 1, 1));

        // Derived: every off-diagonal of the cleared matrix is 2·den, hence even, so mod 2 the
        // matrix is DIAGONAL, the Krylov space spanned by e_seat is one-dimensional, the rank is 1
        // and the count is N − 1 = 8. This holds at every seat and every Δ, which is exactly why
        // "some small prime inflates somewhere" was not a test.
        Assert.Equal(8, w.BlindAtRationalModP(1, 1, 1, 2));
        foreach (int seat in Enumerable.Range(0, 9))
            Assert.Equal(8, w.BlindAtRationalModP(seat, 1, 1, 2));

        // So a prime set containing 2 separates the two policies by 7 on this row. Keeping the
        // LARGER rank (the smaller count) is what makes a bad prime harmless; keeping the smaller
        // rank would report 8 where the truth is 1.
        // and the policy is reached where it LIVES, not restated beside it: hand BlindAtRational a
        // prime set containing 2. Keeping the larger rank makes the bad prime harmless and the
        // answer is 1; keeping the smaller rank would return 8. Flipping Math.Max to Math.Min in
        // the witness reddens this line, and reddens nothing else in the file (measured).
        long[] withBad = { 2L, SeatBlindnessDeltaLocusWitness.Primes.Max() };
        Assert.Equal(1, w.BlindAtRational(1, 1, 1, withBad));
        Assert.Equal(8, withBad.Max(q => w.BlindAtRationalModP(1, 1, 1, q)));   // the inverse policy

        // And the honest disclosure the old comment made and the assertions then buried: on the two
        // REAL primes the maximum never has to do any work. Pinned as a count, so a future prime
        // change that starts making it work is visible rather than silent.
        long[] realPrimes = SeatBlindnessDeltaLocusWitness.Primes.ToArray();
        int rows = 0, disagreements = 0;
        foreach (int n in new[] { 7, 9, 11 })
        {
            var v = new SeatBlindnessDeltaLocusWitness(n);
            for (int seat = 0; seat < n; seat++)
                foreach ((long num, long den) in new[] { (1L, 1L), (0L, 1L), (1L, 3L), (3L, 1L) })
                {
                    rows++;
                    if (v.BlindAtRationalModP(seat, num, den, realPrimes[0])
                        != v.BlindAtRationalModP(seat, num, den, realPrimes[1])) disagreements++;
                }
        }
        Assert.Equal(108, rows);           // (7 + 9 + 11) seats times four Δ
        Assert.Equal(0, disagreements);
    }

    [Fact]
    public void TheGuards_RefuseOnBothSidesAndOnTheNegativeBranch()
    {
        var w = new SeatBlindnessDeltaLocusWitness(9);
        long cap = SeatBlindnessDeltaLocusWitness.MaxDeltaTerm;
        Assert.Throws<ArgumentOutOfRangeException>(() => w.BlindAtRational(1, cap + 1, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.BlindAtRational(1, 1, cap + 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.BlindAtRational(1, -(cap + 1), 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.BlindAtRationalModP(1, 1, 1, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => w.BlindAtRationalModP(9, 1, 1, 7));

        // The cap's stated reason is exact and checkable rather than a threshold that merely
        // passes: it holds every matrix entry strictly below both ranking primes, which is what
        // makes the sibling witness's fake-blindness mode impossible here rather than defended.
        long worstEntry = Math.Max(2 * cap, (SeatBlindnessDeltaLocusWitness.MaxN - 1) * cap);
        Assert.True(worstEntry < SeatBlindnessDeltaLocusWitness.Primes.Min(),
            $"worst entry {worstEntry} is not below the smallest ranking prime "
            + SeatBlindnessDeltaLocusWitness.Primes.Min());
    }

    // -----------------------------------------------------------------------------------------
    // The headline number, which nothing tested and which was therefore wrong at every N.
    // Summary computed nodeModulus % j instead of j % nodeModulus -- the reverse divisibility. The tree below it
    // printed the truth the whole time, so only a pin on the SUMMARY ITSELF could catch it.
    // -----------------------------------------------------------------------------------------

    [Theory]
    [InlineData(7, new[] { 1, 5 })]
    [InlineData(9, new[] { 1, 2, 3, 5, 6, 7 })]
    [InlineData(11, new[] { 1, 2, 3, 7, 8, 9 })]
    [InlineData(12, new[] { 1, 2, 3, 4, 7, 8, 9, 10 })]
    [InlineData(15, new[] { 1, 2, 3, 4, 5, 9, 10, 11, 12, 13 })]
    [InlineData(16, new[] { 1, 2, 3, 4, 11, 12, 13, 14 })]
    public void TheSeatsWithANonemptyLocus_ArePinnedAsASetAndNotOnlyAsACount(int n, int[] expected)
    {
        var w = new SeatBlindnessDeltaLocusWitness(n);
        Assert.Equal(expected, Enumerable.Range(0, n).Where(w.HasNonemptyLocus).ToArray());

        // and the predicate must agree with the polynomial it is a shortcut for, seat by seat
        for (int seat = 1; seat < n - 1; seat++)
        {
            if (w.NodeModulus(seat) == 0) continue;
            bool nonconstant = w.LocusPolynomial(seat)!.Length > 1;
            Assert.Equal(nonconstant, w.HasNonemptyLocus(seat));
        }
    }

    [Fact]
    public void TheSummaryReportsTheSameCountTheTreeDoes()
    {
        // The bug this pins: a headline derived from its own predicate rather than from the
        // objects it summarises. Both sides are recomputed here, and they are computed from
        // different things -- the polynomial's degree, and the predicate.
        foreach (int n in new[] { 7, 9, 11, 12, 15, 16 })
        {
            var w = new SeatBlindnessDeltaLocusWitness(n);
            int fromPolynomials = Enumerable.Range(1, n - 2)
                .Count(seat => w.NodeModulus(seat) != 0 && w.LocusPolynomial(seat)!.Length > 1);
            // ".": without the terminator, "locus: 10" also matches a summary reporting 100.
            Assert.Contains($"nonempty finite locus: {fromPolynomials}.", w.Summary);
            // The forced list is the other half of the same sentence and had no gate: emptying it
            // passed everything. At odd N the centre seat is forced and nothing else is; at even N
            // there is none. Derived from N_node = 0 <=> 2j = N-1.
            if (n % 2 == 1) Assert.Contains($"Forced (blind at every Delta): {n / 2}", w.Summary);
            else Assert.Contains("Forced (blind at every Delta): none", w.Summary);
        }
    }

    [Fact]
    public void SeatOneAlwaysHasANonemptyLocus_WhichIsTheCaseTheInvertedTestAlwaysDropped()
    {
        // nodeModulus % 1 == 0 for every nodeModulus, so the reverse divisibility excluded seat 1 at every N without
        // exception. It is the seat both worked examples in the registry entry are about.
        for (int n = 5; n <= SeatBlindnessDeltaLocusWitness.MaxLocusN; n++)
        {
            var w = new SeatBlindnessDeltaLocusWitness(n);
            if (w.NodeModulus(1) == 0) continue;
            Assert.True(w.HasNonemptyLocus(1), $"N={n}: seat 1 should carry a nonempty locus");
            Assert.True(w.LocusPolynomial(1)!.Length > 1);
        }
    }

    // -----------------------------------------------------------------------------------------
    // The cross-check that makes the closed form a CLAIM rather than an arithmetic coincidence.
    // LocusPolynomial builds Res_x(U_{nodeModulus-1}, D*U_{j-1} - U_j), which presupposes that the shared
    // eigenvalues sit at the roots of U_{nodeModulus-1} -- the very thing being asserted. DefinitionPolynomial
    // forms Res_lambda(chi_L, chi_R) from the two submatrices instead, by the Jacobi three-term
    // recursion, using no Chebyshev identity at all. Agreement is the closed form being CORRECT.
    // -----------------------------------------------------------------------------------------

    [Fact]
    public void TheChebyshevRouteAndTheDefinitionRouteAgree_AsPrimitiveIntegerPolynomials()
    {
        int sharp = 0, empty = 0;
        for (int n = 4; n <= 16; n++)
        {
            var w = new SeatBlindnessDeltaLocusWitness(n);
            for (int seat = 1; seat < n - 1; seat++)
            {
                if (w.NodeModulus(seat) == 0)
                {
                    Assert.Null(w.LocusPolynomial(seat));
                    Assert.Null(w.DefinitionPolynomial(seat));
                    continue;
                }
                var chebyshev = w.LocusPolynomial(seat)!;
                var definition = w.DefinitionPolynomial(seat)!;
                Assert.Equal(chebyshev, definition);
                if (chebyshev.Length > 1) sharp++; else empty++;
            }
        }
        // Reported rather than averaged away: an agreement on {1} == {1} is a real prediction but
        // a weaker one, so the count that carries the claim is the SHARP one.
        Assert.Equal(70, sharp);
        Assert.Equal(28, empty);
    }

    [Fact]
    public void TheDefinitionRouteIsNotTheChebyshevRouteInDisguise()
    {
        // What this gate used to be, and why it was replaced: it compared the two routes at
        // DIFFERENT SEATS, so it only ever showed that seat 1 and seat 2 have different
        // polynomials. Replacing DefinitionPolynomial's whole body with
        //     return LocusPolynomial(seat);
        // passed it verbatim, which is the one mutation it was named for. Measured, not argued.
        //
        // The real question is whether the definition route reads the CHAIN at all. So feed it a
        // physically wrong end shift (N-4 where the ZZ term gives N-3) and require the routes to
        // part; the disguised version cannot part, because it never looks at the shift.
        // WHICH wrong shift, and why this one. N-5 is the interior diagonal itself, so feeding it as
        // the END shift removes the defect and leaves both halves with a CONSTANT diagonal. Two paths
        // of sizes p and q with equal diagonals share an eigenvalue whenever gcd(p+1, q+1) > 1, at
        // every Δ, so the resultant collapses to the zero polynomial there. That makes N-5 the
        // sharpest probe available and the only one that separates at EVERY seat, including the seat
        // whose locus is {0}. Verified independently in sympy at N = 9 over shifts 3..12: seat 1 and
        // seat 3 give the zero polynomial at shift 4 = N-5 and nothing else does.
        //
        // (An earlier version of this gate used N-4 only, and its comment claimed seat 3 gives Δ "for
        // EVERY shift, measured" on the strength of that single shift. It does not: at N-5 it gives
        // the zero polynomial. The claim was a measurement asserted from one sample.)
        var w = new SeatBlindnessDeltaLocusWitness(9);
        foreach (int seat in new[] { 1, 2, 3 })
        {
            Assert.Equal(w.LocusPolynomial(seat), w.DefinitionPolynomial(seat));
            Assert.NotEqual(w.LocusPolynomial(seat), w.DefinitionPolynomial(seat, w.N - 5));
        }
        // The degenerate shift really does collapse the route, rather than merely moving it: at the
        // two seats whose halves have gcd(p+1, q+1) > 1 the answer is the zero polynomial.
        Assert.Equal(new[] { System.Numerics.BigInteger.Zero }, w.DefinitionPolynomial(1, w.N - 5));
        Assert.Equal(new[] { System.Numerics.BigInteger.Zero }, w.DefinitionPolynomial(3, w.N - 5));
        // A shift that is wrong but NOT degenerate still separates, at the seats with room to move.
        foreach (int seat in new[] { 1, 2 })
        {
            Assert.NotEqual(w.LocusPolynomial(seat), w.DefinitionPolynomial(seat, w.N - 4));
            Assert.NotEqual(w.LocusPolynomial(seat), w.DefinitionPolynomial(seat, w.N - 2));
        }
    }

    // ---------------------------------------------------------------------------------------
    // The route the polynomial never touches: a rank sweep over the whole field.
    // ---------------------------------------------------------------------------------------

    [Fact]
    public void TheRankSweepReachesAnIrrationalMemberThroughItsImageModP()
    {
        // N = 9 seat 2 is the headline class: blind at NEITHER committed endpoint, locus ±√2/2.
        // Mod 7 that locus is 2Δ² − 1 = 0, i.e. Δ² = 2⁻¹ = 4, i.e. Δ = ±2. Derived by hand here,
        // not read off the witness, so the literal survives any change to the polynomial code.
        var w = new SeatBlindnessDeltaLocusWitness(9);
        var (byRank, byPoly) = w.LocusOverFp(2, 7)!.Value;
        Assert.Equal(new long[] { 2, 5 }, byRank);          // the RANK, with no polynomial in sight
        Assert.Equal(new long[] { 2, 5 }, byPoly);
    }

    [Fact]
    public void TheRankSweepAndThePolynomialNameTheSameSet_AtEveryInteriorSeatAndSeveralPrimes()
    {
        int rows = 0, nonEmpty = 0;
        foreach (int n in new[] { 7, 9, 11, 12 })
        {
            var w = new SeatBlindnessDeltaLocusWitness(n);
            foreach (long prime in new long[] { 7, 11, 13, 101, 211 })
                for (int seat = 1; seat < n - 1; seat++)
                {
                    var got = w.LocusOverFp(seat, prime);
                    if (got is null) continue;              // end seat or the forced centre
                    var (byRank, byPoly) = got.Value;
                    Assert.Equal(byPoly, byRank);
                    rows++;
                    if (byPoly.Length > 0) nonEmpty++;
                }
        }
        // Without these the test could pass on rows that are all empty on both sides.
        Assert.Equal(140, rows);   // 4+6+8+10 interior non-centre seats, times 5 primes
        // MEASURED, not derived: how many of the 140 rows have a locus at all mod p. Without it
        // the equality above could pass on rows that are empty on both sides.
        Assert.Equal(94, nonEmpty);
    }

    [Fact]
    public void TheRankSweepAgreesWithAnEliminationWrittenHere_NotWithTheWitnessArithmetic()
    {
        // The gate this replaces compared the sweep's output against a wrong polynomial evaluated by
        // the TEST's own loop, and then against a set the test had already pinned twelve lines above,
        // so the only surviving content was NotEqual({3,4}, {2,5}). It asserted nothing about the
        // witness. This one rebuilds the rank side from the PHYSICS, in this file, and compares.
        //
        // The single-excitation matrix on a uniform open chain: hop 2 on every bond, and the ZZ term
        // pays -Delta at each of a bond's two endpoints and +Delta at every other site, so a site
        // touching k of the N-1 bonds carries Delta*(N-1-2k). That is 6*Delta at the two chain ends
        // and 4*Delta inside, at N = 9. Written from that sentence, not from BuildCleared.
        const int N = 9, prime = 7;
        static int BlindHere(int n, int seat, long delta, long p)
        {
            var h = new long[n, n];
            for (int b = 0; b < n - 1; b++) { h[b, b + 1] += 2; h[b + 1, b] += 2; }
            for (int site = 0; site < n; site++)
            {
                int touches = (site > 0 ? 1 : 0) + (site < n - 1 ? 1 : 0);
                h[site, site] = ((n - 1 - 2 * touches) * delta % p + p) % p;
            }
            // Krylov rows from e_seat, then a plain elimination over GF(p).
            var rows = new List<long[]>();
            var v = new long[n]; v[seat] = 1;
            for (int k = 0; k <= n; k++)
            {
                rows.Add((long[])v.Clone());
                var next = new long[n];
                for (int a = 0; a < n; a++)
                {
                    long acc = 0;
                    for (int b = 0; b < n; b++) acc = (acc + ((h[a, b] % p + p) % p) * v[b]) % p;
                    next[a] = acc;
                }
                v = next;
            }
            int rank = 0;
            for (int col = 0; col < n && rank < rows.Count; col++)
            {
                int pivot = -1;
                for (int r = rank; r < rows.Count; r++) if (rows[r][col] % p != 0) { pivot = r; break; }
                if (pivot < 0) continue;
                (rows[rank], rows[pivot]) = (rows[pivot], rows[rank]);
                long inv = 1, bse = ((rows[rank][col] % p) + p) % p, e = p - 2;
                while (e > 0) { if ((e & 1) == 1) inv = inv * bse % p; bse = bse * bse % p; e >>= 1; }
                for (int r = 0; r < rows.Count; r++)
                {
                    if (r == rank || rows[r][col] % p == 0) continue;
                    long f = rows[r][col] * inv % p;
                    for (int c = 0; c < n; c++) rows[r][c] = ((rows[r][c] - f * rows[rank][c]) % p + p) % p;
                }
                rank++;
            }
            return n - rank;
        }

        var w = new SeatBlindnessDeltaLocusWitness(N);
        foreach (int seat in new[] { 1, 2, 3, 5, 6, 7 })
        {
            var mine = Enumerable.Range(0, prime).Where(d => BlindHere(N, seat, d, prime) > 0)
                                 .Select(d => (long)d).ToArray();
            var (byRank, byPoly) = w.LocusOverFp(seat, prime)!.Value;
            Assert.Equal(mine, byRank);      // two eliminations, written independently
            Assert.Equal(mine, byPoly);      // and the polynomial names the same set
        }
        // And the row that carries the point: seat 2's locus is +-sqrt(2)/2, whose image mod 7 is
        // +-2 because 2^-1 = 4 and 2^2 = 4. Derived, then found by an elimination that never saw a
        // Chebyshev polynomial.
        Assert.Equal(new long[] { 2, 5 },
            Enumerable.Range(0, prime).Where(d => BlindHere(N, 2, d, prime) > 0).Select(d => (long)d));
    }

    [Fact]
    public void TheSweepGuardsRefuseRatherThanReturnAnEmptyComparison()
    {
        var w = new SeatBlindnessDeltaLocusWitness(9);
        Assert.Null(w.LocusOverFp(0, 7));                                    // end seat
        Assert.Null(w.LocusOverFp(8, 7));                                    // the other end
        Assert.Null(w.LocusOverFp(4, 7));                                    // the forced centre
        Assert.Throws<ArgumentOutOfRangeException>(() => w.LocusOverFp(1, 2));
        Assert.Throws<ArgumentOutOfRangeException>(
            () => w.LocusOverFp(1, SeatBlindnessDeltaLocusWitness.MaxSweepPrime + 2));
    }

    [Fact]
    public void ThePerSeatBodiesAreAsserted_NotOnlyTheSevenTopLevelSummaries()
    {
        // The tree test builds its string from the seven TOP-LEVEL summaries, so every per-seat body
        // was unasserted prose: the formatted polynomial, the printed counts, the rational-point list,
        // and the end-seat and centre-seat wordings. Format has real branching (unit coefficients
        // suppressed, sign placement, D^k) and a bug in it was invisible.
        var w = new SeatBlindnessDeltaLocusWitness(9);
        var seats = w.Children.First().Children.ToList();
        Assert.Equal(9, seats.Count);

        Assert.Contains("an END seat", seats[0].Summary);
        Assert.Contains("blind = 0 at Delta = 1 and 0 at Delta = 0", seats[0].Summary);

        // Seat 1: N_node = 6, polynomial D^5 - 4D^3 + 3D. The formatting is the assertion here:
        // the leading coefficient 1 is suppressed, the interior signs are separated, and the linear
        // term carries no exponent.
        Assert.Contains("N_node = 6; locus polynomial P(D) = D^5 - 4D^3 + 3D;", seats[1].Summary);
        Assert.Contains("its RATIONAL points: -1, 0, 1;", seats[1].Summary);
        Assert.Contains("1 at Delta = 1, 1 at Delta = 0", seats[1].Summary);

        // Seat 2: the class both laws call sighted, so the counts must be 0 and the rational list empty
        // while the polynomial is nonconstant. A leading 2 IS printed, unlike seat 1's leading 1.
        Assert.Contains("N_node = 4; locus polynomial P(D) = 2D^2 - 1;", seats[2].Summary);
        Assert.Contains("its RATIONAL points: none;", seats[2].Summary);
        Assert.Contains("0 at Delta = 1, 0 at Delta = 0", seats[2].Summary);

        Assert.Contains("the reflection-fixed CENTRE seat, N_node = 0", seats[4].Summary);
        Assert.Contains("the count is (N-1)/2 = 4", seats[4].Summary);
    }

    [Fact]
    public void TheKrylovRouteMeetsTheCommittedLaws_AtMoreThanOneChainLength()
    {
        // The rank was checked against the two committed closed forms at N = 9 alone; everywhere else
        // it was checked against itself (the ceiling, the centre, the +-Delta identity) or against the
        // polynomial. This walks the DEFINITION against both laws over several N and both endpoints.
        int rows = 0, nonzero = 0;
        foreach (int n in new[] { 5, 7, 8, 9, 11, 12 })
        {
            var w = new SeatBlindnessDeltaLocusWitness(n);
            for (int seat = 0; seat < n; seat++)
            {
                int h = w.BlindAtRational(seat, 1, 1), x = w.BlindAtRational(seat, 0, 1);
                Assert.Equal(SeatCutBlindnessClaim.BlindHeisenberg(n, seat), h);
                Assert.Equal(SeatCutBlindnessClaim.BlindXy(n, seat), x);
                rows++;
                if (h > 0 || x > 0) nonzero++;
            }
        }
        Assert.Equal(52, rows);                 // 5+7+8+9+11+12
        Assert.True(nonzero >= 20, $"only {nonzero} of {rows} rows are nonzero, so this is mostly 0 == 0");
    }

    // ---------------------------------------------------------------------------------------
    // The rows at IRRATIONAL Delta, which every other route here cannot reach.
    // ---------------------------------------------------------------------------------------

    [Fact]
    public void TheCommittedIrrationalRows_AreComputedHere_NotQuotedFromAScout()
    {
        // These two rows are what the entry's new content rests on, and until this test they were
        // produced by no committed code: the count routes decide RATIONAL Delta and the field sweep
        // decides membership, so both rows came from a local gitignored scout. They are now counted,
        // by the multiplicity of the minimal polynomial in P_j.
        var w = new SeatBlindnessDeltaLocusWitness(9);

        // Delta = sqrt(3), minimal polynomial D^2 - 3. The committed sentence is that this row
        // "reproduces the isotropic row exactly", so it must equal the Delta = 1 row as well.
        var atSqrt3 = w.BlindRowAtAlgebraic(new BigInteger[] { -3, 0, 1 });
        Assert.Equal(new[] { 0, 1, 0, 0, 4, 0, 0, 1, 0 }, atSqrt3);
        Assert.Equal(Enumerable.Range(0, 9).Select(s => w.BlindAtRational(s, 1, 1)).ToArray(), atSqrt3);

        // Delta = sqrt(2)/2, minimal polynomial 2D^2 - 1. This is the row for the class that is
        // genuinely new: seats 2 and 6 are blind here and BOTH committed gcd laws call them sighted.
        var atHalfSqrt2 = w.BlindRowAtAlgebraic(new BigInteger[] { -1, 0, 2 });
        Assert.Equal(new[] { 0, 0, 1, 0, 4, 0, 1, 0, 0 }, atHalfSqrt2);
        Assert.Equal(0, SeatCutBlindnessClaim.BlindHeisenberg(9, 2));
        Assert.Equal(0, SeatCutBlindnessClaim.BlindXy(9, 2));

        // Two-sided: an irrational Delta OFF the locus must leave only the forced centre.
        Assert.Equal(new[] { 0, 0, 0, 0, 4, 0, 0, 0, 0 }, w.BlindRowAtAlgebraic(new BigInteger[] { -5, 0, 1 }));
        Assert.Equal(new[] { 0, 0, 0, 0, 4, 0, 0, 0, 0 }, w.BlindRowAtAlgebraic(new BigInteger[] { -7, 0, 1 }));
    }

    [Fact]
    public void TheCountIsTheMULTIPLICITY_AndTheDistinctReadingWouldBreakTheCommittedLaw()
    {
        // The entry says this is load-bearing rather than bookkeeping, and nothing gated it: the
        // rational bench asks "blind > 0", never "blind == multiplicity", and its root list is a SET,
        // so the double root is flattened exactly where the distinction lives.
        var w = new SeatBlindnessDeltaLocusWitness(11);
        // P_2 = 3D^4 - 4D^2, so 0 is a DOUBLE root: k = 2 and k = 4 both give Delta = 0.
        Assert.Equal(new BigInteger[] { 0, 0, -4, 0, 3 }, w.LocusPolynomial(2));
        Assert.Equal(2, w.BlindAtAlgebraic(2, new BigInteger[] { 0, 1 }));
        // and that is exactly what the committed XY law needs, gcd(3, 12) - 1 = 2. The distinct-root
        // reading would give 1 and the law would fail.
        Assert.Equal(2, SeatCutBlindnessClaim.BlindXy(11, 2));
        Assert.Equal(2, w.BlindAtRational(2, 0, 1));
        // The rational route agrees at this Delta because 0 is rational; the point is that the
        // POLYNOMIAL side had to be read with multiplicity to match it.
        Assert.Single(w.RationalLocusPoints(2));
    }

    [Fact]
    public void ACompositeModulusIsRefused_BecauseItWouldAgreeInSilence()
    {
        // The elimination inverts by Fermat, so a composite modulus has no inverses and the whole
        // GF(p) argument is void. It does not fail loudly: measured at 9, 15, 25 and 2001, BOTH sides
        // of the sweep came back empty and the comparison reported agreement. That is a gate agreeing
        // with itself about nothing, which is the shape this repo has a name for.
        var w = new SeatBlindnessDeltaLocusWitness(9);
        foreach (long composite in new long[] { 9, 15, 25, 2001 })
            Assert.Throws<ArgumentOutOfRangeException>(() => w.LocusOverFp(1, composite));
        // and the primes on either side of one of them still work
        Assert.NotNull(w.LocusOverFp(1, 23));
        Assert.NotNull(w.LocusOverFp(1, 29));
    }

    [Fact]
    public void WhereTheSweepCarriesNothing_ItSaysSo_RatherThanClaimingTheOpposite()
    {
        // The sweep node used to print "0 of them with a nonempty locus so the comparison is not
        // 0 == 0", which asserts its own refutation: zero rows with a locus is exactly a comparison
        // of 0 against 0. At N = 4 both interior seats have node modulus 1, so their locus is empty
        // by the d | j rule and there is nothing to compare at any prime.
        var small = new SeatBlindnessDeltaLocusWitness(4);
        Assert.Equal(1, small.NodeModulus(1));
        Assert.Equal(1, small.NodeModulus(2));
        Assert.Empty(small.LocusOverFp(1, 101)!.Value.PolynomialSaysBlind);
        string smallJoined = string.Join(" | ", small.Children.Select(c => c.Summary));
        Assert.Contains("the comparison here IS 0 == 0 and carries nothing", smallJoined);
        Assert.DoesNotContain("so the comparison is not 0 == 0", smallJoined);

        // and the other branch still reads the other way where the comparison does carry something
        var w = new SeatBlindnessDeltaLocusWitness(9);
        string joined = string.Join(" | ", w.Children.Select(c => c.Summary));
        Assert.Contains("so the comparison is not 0 == 0", joined);
        Assert.DoesNotContain("carries nothing", joined);
    }

    [Theory]
    // The nine-row bench the local number-field scout measured on 2026-08-30 and that lived only in
    // a gitignored file: (N, seat, minimal polynomial ascending, expected blind). Both directions on
    // every value, so a route that answered "blind" everywhere would fail half the rows. Row six is
    // degree FOUR, which is also the row that exercises the irreducibility guard's quadratic search.
    [InlineData(9,  1, new long[] { -3, 0, 1 },        1)]   // sqrt(3)
    [InlineData(9,  1, new long[] { -5, 0, 1 },        0)]   // sqrt(5), control
    [InlineData(9,  2, new long[] { -1, 0, 2 },        1)]   // sqrt(2)/2
    [InlineData(11, 2, new long[] { -1, 0, 2 },        0)]   // sqrt(2)/2 is NOT at N = 11, control
    [InlineData(7,  1, new long[] { -2, 0, 1 },        1)]   // sqrt(2)
    [InlineData(11, 1, new long[] { 2, 0, -4, 0, 1 },  1)]   // sqrt(2 +- sqrt 2), degree 4
    [InlineData(11, 2, new long[] { -4, 0, 3 },        1)]   // 2*sqrt(3)/3
    [InlineData(11, 1, new long[] { -4, 0, 3 },        0)]   // and not at seat 1, control
    [InlineData(13, 2, new long[] { -1, 0, 2 },        1)]   // sqrt(2)/2's second home
    public void TheNumberFieldBench_IsReproducedByTheMultiplicityRoute(
        int n, int seat, long[] minimalPolynomial, int expected)
    {
        var w = new SeatBlindnessDeltaLocusWitness(n);
        Assert.Equal(expected, w.BlindAtAlgebraic(seat, minimalPolynomial.Select(c => (BigInteger)c).ToArray()));
    }

    [Fact]
    public void AReducibleMinimalPolynomialIsRefused_BecauseItNamesSeveralDeltaAtOnce()
    {
        // The prototype this route ports raises on a reducible input and tells its porter to guard
        // rather than return a wrong rank. The count is the multiplicity of ONE Delta's minimal
        // polynomial; a product names several, and the number it would return answers no question.
        var w = new SeatBlindnessDeltaLocusWitness(9);

        // The realistic mistake: handing over the whole locus polynomial. At N = 9 seat 1 that is
        // D^5 - 4D^3 + 3D = D(D-1)(D+1)(D^2-3), caught by its rational root 0.
        Assert.Throws<ArgumentException>(() => w.BlindAtAlgebraic(1, w.LocusPolynomial(1)!));
        // A reducible quartic with NO rational root, which only the quadratic search can catch:
        // (D^2 - 2)(D^2 - 3) = D^4 - 5D^2 + 6.
        Assert.Throws<ArgumentException>(() => w.BlindAtAlgebraic(1, new BigInteger[] { 6, 0, -5, 0, 1 }));
        // Not squarefree: (D^2 - 3)^2.
        Assert.Throws<ArgumentException>(() => w.BlindAtAlgebraic(1, new BigInteger[] { 9, 0, -6, 0, 1 }));
        // Above degree 5 the guard refuses rather than guessing, and says so.
        Assert.Throws<ArgumentException>(
            () => w.BlindAtAlgebraic(1, new BigInteger[] { -1, 0, 0, 0, 0, 0, 1 }));
        // A constant after trimming would make the counting loop run forever; refused.
        Assert.Throws<ArgumentException>(() => w.BlindAtAlgebraic(1, new BigInteger[] { 5, 0 }));

        // and the irreducible ones the bench uses still pass, including the quartic
        Assert.Equal(1, w.BlindAtAlgebraic(1, new BigInteger[] { -3, 0, 1 }));
        Assert.Equal(1, new SeatBlindnessDeltaLocusWitness(11)
            .BlindAtAlgebraic(1, new BigInteger[] { 2, 0, -4, 0, 1 }));
    }
}
