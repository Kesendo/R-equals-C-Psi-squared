using RCPsiSquared.Core.States;

namespace RCPsiSquared.Core.Tests.States;

public class BondingModeTests
{
    [Theory]
    [InlineData(3, 1)]
    [InlineData(4, 2)]
    [InlineData(5, 3)]
    public void Build_IsNormalized(int N, int k)
    {
        var psi = BondingMode.Build(N, k);
        double norm = psi.ConjugateDotProduct(psi).Real;
        Assert.Equal(1.0, norm, 10);
    }

    [Fact]
    public void Build_Amplitudes_FollowF65SineFormula()
    {
        // ψ_k(j) = √(2/(N+1)) · sin(π·k·(j+1)/(N+1)) at single-excitation index 2^(N-1-j).
        int N = 4;
        int k = 1;
        var psi = BondingMode.Build(N, k);
        double norm = Math.Sqrt(2.0 / (N + 1));
        for (int j = 0; j < N; j++)
        {
            double expected = norm * Math.Sin(Math.PI * k * (j + 1) / (N + 1));
            int idx = 1 << (N - 1 - j);
            Assert.Equal(expected, psi[idx].Real, 10);
        }
    }

    [Theory]
    [InlineData(3, 0)]
    [InlineData(3, 4)]
    public void Build_RejectsKOutOfRange(int N, int k)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BondingMode.Build(N, k));
    }

    [Fact]
    public void PairState_IsNormalized()
    {
        var psi = BondingMode.PairState(N: 4, k: 2);
        double norm = psi.ConjugateDotProduct(psi).Real;
        Assert.Equal(1.0, norm, 10);
    }
}
