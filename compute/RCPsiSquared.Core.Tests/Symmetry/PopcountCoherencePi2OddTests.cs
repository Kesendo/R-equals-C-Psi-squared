using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class PopcountCoherencePi2OddTests
{
    // ─────────────── α anchors (general, all popcount pairs) ───────────────

    [Theory]
    // Adjacent (n_q = n_p + 1): preserves the original three-anchor structure
    [InlineData(3, 1, 2, 0.0)]               // adjacent inter popcount-mirror N=3 central
    [InlineData(5, 2, 3, 0.0)]               // adjacent inter popcount-mirror N=5 central
    [InlineData(7, 3, 4, 0.0)]               // adjacent inter popcount-mirror N=7 central
    [InlineData(9, 4, 5, 0.0)]               // adjacent inter popcount-mirror N=9 central
    [InlineData(4, 1, 2, 3.0 / 10)]          // adjacent K-intermediate (near-mirror half) N=4
    [InlineData(6, 2, 3, 2.0 / 7)]           // adjacent K-intermediate N=6
    [InlineData(8, 3, 4, 5.0 / 18)]          // adjacent K-intermediate N=8
    [InlineData(10, 4, 5, 3.0 / 11)]         // adjacent K-intermediate N=10
    [InlineData(5, 0, 1, 0.5)]               // adjacent generic boundary
    [InlineData(5, 1, 2, 0.5)]               // adjacent generic interior
    [InlineData(7, 2, 3, 0.5)]               // adjacent generic
    // Non-adjacent: covered only by the generalised formula
    [InlineData(2, 0, 2, 0.0)]               // popcount-mirror N=2 (= GHZ_2 popcount)
    [InlineData(3, 0, 3, 0.0)]               // popcount-mirror N=3 (GHZ_3 popcount)
    [InlineData(5, 1, 4, 0.0)]               // non-adjacent inter-mirror N=5
    [InlineData(6, 2, 4, 0.0)]               // non-adjacent inter-mirror N=6
    [InlineData(4, 2, 2, 0.0)]               // intra-mirror n_p = n_q = N/2 = 2
    [InlineData(6, 3, 3, 0.0)]               // intra-mirror n_p = n_q = N/2 = 3
    [InlineData(4, 0, 2, 3.0 / 7)]           // non-adjacent K-intermediate N=4 (n_q = N/2)
    [InlineData(6, 0, 3, 10.0 / 21)]         // non-adjacent K-intermediate N=6 (n_q = N/2)
    [InlineData(6, 1, 3, 5.0 / 13)]          // non-adjacent K-intermediate N=6 (n_q = N/2)
    [InlineData(6, 3, 5, 5.0 / 13)]          // K-intermediate, mirror image of (1, 3)
    [InlineData(3, 1, 1, 0.5)]               // intra-non-mirror generic
    [InlineData(5, 1, 1, 0.5)]               // intra-non-mirror generic
    [InlineData(5, 0, 2, 0.5)]               // inter-non-mirror non-K generic
    public void AlphaAnchor_MatchesExpectedAcrossAllPopcountPairs(int N, int np, int nq, double expected)
    {
        Assert.Equal(expected, PopcountCoherencePi2Odd.AlphaAnchor(N, np, nq), 12);
    }

    [Theory]
    // K-intermediate closed form α = C(N, N/2) / (2·(C(N, n_other) + C(N, N/2)))
    [InlineData(4, 0, 2, 3.0 / 7)]
    [InlineData(4, 1, 2, 3.0 / 10)]          // adjacent near-mirror half: recovers (N+2)/(4(N+1)) = 6/20
    [InlineData(4, 2, 3, 3.0 / 10)]
    [InlineData(4, 2, 4, 3.0 / 7)]           // mirror image of (0, 2)
    [InlineData(6, 0, 3, 10.0 / 21)]
    [InlineData(6, 1, 3, 5.0 / 13)]
    [InlineData(6, 2, 3, 2.0 / 7)]           // adjacent near-mirror half: 8/28
    [InlineData(6, 3, 4, 2.0 / 7)]
    [InlineData(6, 3, 5, 5.0 / 13)]
    [InlineData(8, 3, 4, 5.0 / 18)]          // adjacent near-mirror half: 10/36
    [InlineData(10, 4, 5, 3.0 / 11)]         // adjacent near-mirror half: 12/44
    public void AlphaKIntermediateClosed_MatchesFormula(int N, int np, int nq, double expected)
    {
        Assert.Equal(expected, PopcountCoherencePi2Odd.AlphaKIntermediateClosed(N, np, nq), 12);
    }

    [Theory]
    [InlineData(5, 2, 3)]                    // odd N: never K-intermediate
    [InlineData(4, 0, 1)]                    // neither n_p nor n_q is N/2
    [InlineData(4, 2, 2)]                    // intra-mirror, not K-intermediate
    public void AlphaKIntermediateClosed_ThrowsOnNonKIntermediate(int N, int np, int nq)
    {
        Assert.Throws<ArgumentException>(() => PopcountCoherencePi2Odd.AlphaKIntermediateClosed(N, np, nq));
    }

    [Fact]
    public void AlphaAnchor_AgreesWithKrawtchouk_AllCompatiblePairs_NUpTo10()
    {
        for (int N = 2; N <= 10; N++)
        {
            for (int np = 0; np <= N; np++)
            {
                for (int nq = np; nq <= N; nq++)
                {
                    double alphaK = PopcountCoherencePi2Odd.AlphaKrawtchouk(N, np, nq);
                    double alphaA = PopcountCoherencePi2Odd.AlphaAnchor(N, np, nq);
                    Assert.True(Math.Abs(alphaK - alphaA) < 1e-10,
                        $"N={N} (n_p={np}, n_q={nq}): Krawtchouk={alphaK}, anchor={alphaA}");
                }
            }
        }
    }

    // ─────────────── Pi2OddInMemory closed-form values ───────────────

    [Theory]
    [InlineData(3, 1, 2, 1, 3.0 / 5)]            // adjacent mirror, HD<N: 3/5
    [InlineData(5, 2, 3, 1, 10.0 / 19)]          // adjacent mirror, HD<N: 10/19
    [InlineData(5, 2, 3, 3, 10.0 / 19)]          // same, HD=3 still <N=5
    [InlineData(7, 3, 4, 1, 35.0 / 69)]          // adjacent mirror N=7
    [InlineData(5, 0, 1, 1, 0.5)]                // adjacent generic
    [InlineData(5, 1, 2, 1, 0.5)]                // adjacent generic interior
    [InlineData(2, 0, 2, 2, 0.0)]                // GHZ_2 popcount-(0,2) HD=N=2 → Π²-classical
    [InlineData(3, 0, 3, 3, 0.0)]                // GHZ_3 HD=N=3 → 0
    [InlineData(4, 0, 4, 4, 0.0)]                // GHZ_4 HD=N=4 → 0
    [InlineData(5, 0, 5, 5, 0.0)]                // GHZ_5 HD=N=5 → 0
    [InlineData(2, 1, 1, 2, 0.0)]                // Singlet/Triplet at N=2: HD=N=2 → 0
    [InlineData(4, 2, 2, 2, 3.0 / 5)]            // intra-mirror N=4 HD<N: 3/5
    [InlineData(4, 0, 2, 2, 9.0 / 17)]           // K-intermediate N=4 popcount-(0,2): 9/17
    [InlineData(6, 0, 3, 3, 30.0 / 59)]          // K-intermediate N=6 popcount-(0,3): 30/59
    [InlineData(5, 1, 4, 3, 5.0 / 9)]            // non-adjacent inter-mirror N=5 popcount-(1,4) HD=3
    [InlineData(6, 2, 4, 2, 15.0 / 29)]          // non-adjacent inter-mirror N=6 popcount-(2,4) HD=2
    public void Pi2OddInMemory_KnownClosedFormValues(int N, int np, int nq, int hd, double expected)
    {
        Assert.Equal(expected, PopcountCoherencePi2Odd.Pi2OddInMemory(N, np, nq, hd), 10);
    }

    // ─────────────── Krawtchouk known values ───────────────

    [Theory]
    [InlineData(1, 0, 3, 3)]                 // K_1(0; 3) = coef. of z in (1+z)^3 = 3
    [InlineData(1, 1, 3, 1)]                 // K_1(1; 3) = coef. of z in (1−z)(1+z)^2 = 1
    [InlineData(1, 2, 3, -1)]
    [InlineData(1, 3, 3, -3)]                // K_1(N; N) = −K_1(0; N) by reflection
    [InlineData(2, 0, 5, 10)]                // K_2(0; 5) = C(5, 2) = 10
    [InlineData(2, 1, 5, 2)]
    public void Krawtchouk_KnownValues(int n, int s, int N, long expected)
    {
        Assert.Equal(expected, PopcountCoherencePi2Odd.Krawtchouk(n, s, N));
    }

    [Theory]
    [InlineData(2, 1)]                       // K_1(1; 2) = 0
    [InlineData(4, 1)]                       // K_2(1; 4) = 0
    [InlineData(4, 3)]                       // K_2(3; 4) = 0
    [InlineData(6, 1)]                       // K_3(1; 6) = 0
    [InlineData(6, 3)]                       // K_3(3; 6) = 0
    [InlineData(6, 5)]                       // K_3(5; 6) = 0
    public void Krawtchouk_HalfPopcountVanishesOnOddS(int N, int sOdd)
    {
        Assert.Equal(0, PopcountCoherencePi2Odd.Krawtchouk(N / 2, sOdd, N));
    }

    [Fact]
    public void Krawtchouk_ReflectionOrthogonalityLemma_AllPairsAtNUpTo7()
    {
        // Σ_s (−1)^s · C(N, s) · K_n(s; N) · K_m(s; N) = 2^N · C(N, n) · [n + m = N].
        // This is the core lemma underlying the closed form for E − O at all
        // anchor categories in PopcountCoherencePi2Odd.
        for (int N = 2; N <= 7; N++)
        {
            for (int n = 0; n <= N; n++)
            {
                for (int m = 0; m <= N; m++)
                {
                    long lhs = 0;
                    for (int s = 0; s <= N; s++)
                    {
                        long term = PopcountCoherencePi2Odd.Binomial(N, s)
                                  * PopcountCoherencePi2Odd.Krawtchouk(n, s, N)
                                  * PopcountCoherencePi2Odd.Krawtchouk(m, s, N);
                        lhs += (s % 2 == 0) ? term : -term;
                    }
                    long rhs = (n + m == N) ? (1L << N) * PopcountCoherencePi2Odd.Binomial(N, n) : 0;
                    Assert.True(lhs == rhs,
                        $"N={N} n={n} m={m}: LHS={lhs}, RHS={rhs}");
                }
            }
        }
    }

    // ─────────────── Static fraction ───────────────

    [Theory]
    [InlineData(5, 0, 1)]
    [InlineData(5, 2, 3)]
    [InlineData(7, 3, 4)]
    public void StaticFraction_InterSector_MatchesDirectBinomial(int N, int np, int nq)
    {
        double expected = 0.25 / PopcountCoherencePi2Odd.Binomial(N, np)
                        + 0.25 / PopcountCoherencePi2Odd.Binomial(N, nq);
        Assert.Equal(expected, PopcountCoherencePi2Odd.StaticFraction(N, np, nq), 14);
    }

    [Theory]
    [InlineData(2, 1, 0.5)]                  // intra at N=2, n=1: 1/C(2,1) = 1/2
    [InlineData(3, 1, 1.0 / 3)]              // intra at N=3, n=1: 1/C(3,1) = 1/3
    [InlineData(4, 2, 1.0 / 6)]              // intra-mirror at N=4, n=N/2=2: 1/C(4,2) = 1/6
    public void StaticFraction_IntraSector_UsesSingleSectorFormula(int N, int n, double expected)
    {
        Assert.Equal(expected, PopcountCoherencePi2Odd.StaticFraction(N, n, n), 14);
    }

    // ─────────────── Anchor predicates ───────────────

    [Theory]
    [InlineData(3, 1, 2, true)]              // odd N central inter-mirror
    [InlineData(5, 2, 3, true)]
    [InlineData(7, 3, 4, true)]
    [InlineData(4, 2, 2, true)]              // intra-mirror n_p+n_q=N at n=N/2
    [InlineData(4, 1, 2, false)]
    [InlineData(5, 0, 1, false)]
    public void IsPopcountMirror_TrueIffSumEqualsN(int N, int np, int nq, bool expected)
    {
        Assert.Equal(expected, PopcountCoherencePi2Odd.IsPopcountMirror(N, np, nq));
    }

    [Theory]
    [InlineData(4, 1, 2, true)]              // adjacent near-mirror half
    [InlineData(4, 0, 2, true)]              // non-adjacent K-intermediate (n_q = N/2)
    [InlineData(6, 1, 3, true)]              // non-adjacent K-intermediate
    [InlineData(6, 3, 5, true)]              // K-intermediate (n_p = N/2)
    [InlineData(4, 2, 2, false)]             // intra-mirror is popcount-mirror, not K-intermediate
    [InlineData(5, 2, 3, false)]             // odd N never K-intermediate
    [InlineData(4, 0, 1, false)]             // neither n_p nor n_q is N/2=2
    public void IsKIntermediate_TrueIffOneOfPairIsHalfPopcountAtEvenNNoMirror(int N, int np, int nq, bool expected)
    {
        Assert.Equal(expected, PopcountCoherencePi2Odd.IsKIntermediate(N, np, nq));
    }

    [Theory]
    [InlineData(2, 2, true)]
    [InlineData(3, 3, true)]
    [InlineData(5, 5, true)]
    [InlineData(5, 3, false)]
    [InlineData(5, 1, false)]
    public void IsHdComplement_TrueIffHdEqualsN(int N, int hd, bool expected)
    {
        Assert.Equal(expected, PopcountCoherencePi2Odd.IsHdComplement(N, hd));
    }
}
