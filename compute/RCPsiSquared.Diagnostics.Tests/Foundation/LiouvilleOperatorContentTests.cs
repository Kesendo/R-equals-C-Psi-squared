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
}
