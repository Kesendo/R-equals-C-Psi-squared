using System;
using System.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class LiouvilleOperatorContentTests
{
    [Fact]
    public void VacToSingleExcitation_IsAllNDiff1()
    {
        // |00⟩⟨01| at N=2: a=00, b=01, a⊕b=01, popcount=1 → {1:1}, the |vac⟩⟨ψ_k| band coherence.
        int n = 2, d = 4;
        var vec = ComplexVector.Build.Dense(d * d);
        vec[0b00 * d + 0b01] = Complex.One;
        var (mean, hist) = LiouvilleOperatorContent.NDiffHistogram(vec, n);
        Assert.Equal(1.0, mean, 6);
        Assert.Equal(1.0, hist[1], 6);
    }

    [Fact]
    public void PopulationAndPair_IsHalfZeroHalfTwo()
    {
        // Equal weight on |00⟩⟨00| (n_diff 0) and |01⟩⟨10| (a⊕b=11, popcount 2) → {0:½,2:½}, mean 1:
        // the {0,2}-coherence (population/antisymmetric block) signature.
        int n = 2, d = 4;
        var vec = ComplexVector.Build.Dense(d * d);
        vec[0b00 * d + 0b00] = new Complex(Math.Sqrt(0.5), 0);
        vec[0b01 * d + 0b10] = new Complex(Math.Sqrt(0.5), 0);
        var (mean, hist) = LiouvilleOperatorContent.NDiffHistogram(vec, n);
        Assert.Equal(0.5, hist[0], 6);
        Assert.Equal(0.5, hist[2], 6);
        Assert.Equal(1.0, mean, 6);
    }

    [Fact]
    public void Population_HasBinomialPauliWeight_ButZeroLight()
    {
        // |0><0| is diagonal: Pauli content = the 2^N strings of {I,Z}^N, so the TOTAL-weight
        // histogram is the binomial C(N,w)/2^N (w = #Z letters), mean N/2 - while n_diff (XY-weight)
        // is pure {0:1}. Same vector, opposite gradings: the V-Effect identity refutation in miniature.
        foreach (int n in new[] { 2, 3, 4 })
        {
            int d = 1 << n;
            var vec = ComplexVector.Build.Dense(d * d);
            vec[0] = Complex.One;
            var (meanW, wh) = LiouvilleOperatorContent.PauliWeightHistogram(vec, n);
            var (_, ndh) = LiouvilleOperatorContent.NDiffHistogram(vec, n);
            Assert.Equal(1.0, ndh[0], 9);                                     // no light at all
            for (int w = 0; w <= n; w++)
                Assert.Equal(Binom(n, w) / Math.Pow(2, n), wh.GetValueOrDefault(w), 9);
            Assert.Equal(n / 2.0, meanW, 9);                                  // binomial mean N/2
        }
    }

    [Fact]
    public void SingleCoherence_TotalWeightIsLightPlusZShadow()
    {
        // |001><000| at N=3: XY-weight 1 fixed (the disagreeing site), Z-shadow on the 2 agreeing
        // sites => total weight 1 + Binom(2,k)/4 = {1:¼, 2:½, 3:¼}, while n_diff stays {1:1}. The two
        // histograms differ by exactly the Z-shadow popcount(a&b) - the V-Effect identity in miniature.
        int n = 3, d = 1 << n;
        var vec = ComplexVector.Build.Dense(d * d);
        vec[0b001 * d + 0b000] = Complex.One;
        var (_, wh) = LiouvilleOperatorContent.PauliWeightHistogram(vec, n);
        var (meanNd, ndh) = LiouvilleOperatorContent.NDiffHistogram(vec, n);
        Assert.Equal(1.0, ndh[1], 9);
        Assert.Equal(1.0, meanNd, 9);
        Assert.False(wh.ContainsKey(0));                                      // weight >= n_diff = 1
        Assert.Equal(0.25, wh.GetValueOrDefault(1), 9);
        Assert.Equal(0.50, wh.GetValueOrDefault(2), 9);
        Assert.Equal(0.25, wh.GetValueOrDefault(3), 9);
    }

    private static double Binom(int n, int k)
    {
        double r = 1;
        for (int i = 0; i < k; i++) r = r * (n - i) / (i + 1);
        return r;
    }
}
