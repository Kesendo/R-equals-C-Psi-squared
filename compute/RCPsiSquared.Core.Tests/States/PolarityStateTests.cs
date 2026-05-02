using RCPsiSquared.Core.States;

namespace RCPsiSquared.Core.Tests.States;

public class PolarityStateTests
{
    [Fact]
    public void UniformPlus_AllSitesEqualToInverseSqrtTwoPow_N()
    {
        // |+⟩^3 has every basis amplitude equal to 1/√(2^N) = 1/(2√2).
        int N = 3;
        var psi = PolarityState.Uniform(N, +1);
        double expected = 1.0 / Math.Sqrt(1 << N);
        Assert.Equal(1 << N, psi.Count);
        for (int i = 0; i < psi.Count; i++)
            Assert.Equal(expected, psi[i].Real, 12);
    }

    [Fact]
    public void Build_Normalized()
    {
        var psi = PolarityState.Build(N: 3, signs: new[] { +1, -1, +1 });
        double norm = psi.ConjugateDotProduct(psi).Real;
        Assert.Equal(1.0, norm, 12);
    }

    [Fact]
    public void Build_RejectsInvalidSigns()
    {
        Assert.Throws<ArgumentException>(() => PolarityState.Build(3, new[] { 1, 0, -1 }));
        Assert.Throws<ArgumentException>(() => PolarityState.Build(3, new[] { 1, -1 })); // wrong length
    }

    [Fact]
    public void UniformPlus_VsUniformMinus_HaveOppositeSignsInOddSites()
    {
        // |−⟩^3 = ⊗(|0⟩−|1⟩)/√2 — basis amplitudes alternate sign by Hamming weight.
        var plus = PolarityState.Uniform(3, +1);
        var minus = PolarityState.Uniform(3, -1);
        for (int i = 0; i < 8; i++)
        {
            int popcount = System.Numerics.BitOperations.PopCount((uint)i);
            double expectedSign = popcount % 2 == 0 ? +1 : -1;
            Assert.Equal(expectedSign * Math.Abs(plus[i].Real), minus[i].Real, 12);
        }
    }
}
