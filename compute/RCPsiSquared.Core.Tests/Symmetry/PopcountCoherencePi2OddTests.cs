using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class PopcountCoherencePi2OddTests
{
    [Theory]
    [InlineData(3, 1, 2, 0.0)]               // mirror N=3 central pair
    [InlineData(5, 2, 3, 0.0)]               // mirror N=5 central pair
    [InlineData(7, 3, 4, 0.0)]               // mirror N=7 central pair
    [InlineData(9, 4, 5, 0.0)]               // mirror N=9 central pair
    [InlineData(4, 1, 2, 3.0 / 10)]          // near-mirror half N=4
    [InlineData(6, 2, 3, 2.0 / 7)]           // near-mirror half N=6
    [InlineData(8, 3, 4, 5.0 / 18)]          // near-mirror half N=8
    [InlineData(10, 4, 5, 3.0 / 11)]         // near-mirror half N=10
    [InlineData(5, 0, 1, 0.5)]               // generic boundary
    [InlineData(5, 1, 2, 0.5)]               // generic interior, N=5
    [InlineData(7, 2, 3, 0.5)]               // generic, N=7
    [InlineData(8, 0, 1, 0.5)]               // generic boundary, N=8
    public void AlphaThreeAnchor_MatchesExpectedAtKnownAnchors(int N, int np, int nq, double expected)
    {
        Assert.Equal(expected, PopcountCoherencePi2Odd.AlphaThreeAnchor(N, np, nq), 12);
    }

    [Fact]
    public void ThreeAnchor_AgreesWithKrawtchoukClosedForm_AllPopcountCoherencePairs_NUpTo16()
    {
        for (int N = 3; N <= 16; N++)
        {
            for (int n = 0; n < N; n++)
            {
                double alphaK = PopcountCoherencePi2Odd.AlphaKrawtchouk(N, n, n + 1);
                double alphaA = PopcountCoherencePi2Odd.AlphaThreeAnchor(N, n, n + 1);
                Assert.True(Math.Abs(alphaK - alphaA) < 1e-10,
                    $"N={N} (n_p={n}, n_q={n + 1}): Krawtchouk={alphaK}, three-anchor={alphaA}");
            }
        }
    }

    [Theory]
    [InlineData(3, 1, 2, 3.0 / 5)]           // mirror N=3: 3/5 = 0.6
    [InlineData(5, 2, 3, 10.0 / 19)]         // mirror N=5: 10/19 ≈ 0.5263
    [InlineData(7, 3, 4, 35.0 / 69)]         // mirror N=7: 35/69 ≈ 0.5072
    [InlineData(9, 4, 5, 126.0 / 251)]       // mirror N=9: 126/251 ≈ 0.5020
    [InlineData(5, 0, 1, 0.5)]               // generic gives 0.5 exactly
    [InlineData(5, 1, 2, 0.5)]               // generic interior
    public void Pi2OddInMemory_KnownClosedFormValues(int N, int np, int nq, double expected)
    {
        Assert.Equal(expected, PopcountCoherencePi2Odd.Pi2OddInMemory(N, np, nq), 10);
    }

    [Theory]
    // Krawtchouk(n, s, N): coefficient of z^n in (1−z)^s (1+z)^(N−s).
    [InlineData(1, 0, 3, 3)]                 // K_1(0; 3) = coef. of z in (1+z)^3 = 3
    [InlineData(1, 1, 3, 1)]                 // K_1(1; 3) = coef. of z in (1−z)(1+z)^2 = 1
    [InlineData(1, 2, 3, -1)]                // K_1(2; 3) = -1
    [InlineData(1, 3, 3, -3)]                // K_1(3; 3) = -3 (reflection K_1(N; N) = -K_1(0; N))
    [InlineData(2, 0, 5, 10)]                // K_2(0; 5) = C(5, 2) = 10
    [InlineData(2, 1, 5, 2)]                 // K_2(1; 5) = 2
    public void Krawtchouk_KnownValues(int n, int s, int N, long expected)
    {
        Assert.Equal(expected, PopcountCoherencePi2Odd.Krawtchouk(n, s, N));
    }

    [Theory]
    [InlineData(5, 0, 1)]
    [InlineData(5, 2, 3)]
    [InlineData(7, 3, 4)]
    public void StaticFraction_MatchesDirectBinomial(int N, int np, int nq)
    {
        double expected = 0.25 / PopcountCoherencePi2Odd.Binomial(N, np)
                        + 0.25 / PopcountCoherencePi2Odd.Binomial(N, nq);
        Assert.Equal(expected, PopcountCoherencePi2Odd.StaticFraction(N, np, nq), 14);
    }

    [Theory]
    [InlineData(3, 1, 2, true)]   // odd N central
    [InlineData(5, 2, 3, true)]
    [InlineData(7, 3, 4, true)]
    [InlineData(4, 1, 2, false)]  // even N has no exact mirror in (n, n+1)
    [InlineData(5, 0, 1, false)]
    public void IsPopcountMirror_OnlyAtOddNCentralPair(int N, int np, int nq, bool expected)
    {
        Assert.Equal(expected, PopcountCoherencePi2Odd.IsPopcountMirror(N, np, nq));
    }

    [Theory]
    [InlineData(4, 1, 2, true)]
    [InlineData(4, 2, 3, true)]
    [InlineData(6, 2, 3, true)]
    [InlineData(6, 3, 4, true)]
    [InlineData(5, 2, 3, false)]  // odd N never near-mirror-half
    [InlineData(4, 0, 1, false)]  // not the half-pair
    public void IsNearMirrorHalf_OnlyAtEvenNCentralPairs(int N, int np, int nq, bool expected)
    {
        Assert.Equal(expected, PopcountCoherencePi2Odd.IsNearMirrorHalf(N, np, nq));
    }
}
