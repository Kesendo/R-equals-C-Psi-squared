using RCPsiSquared.Core.Lindblad;

namespace RCPsiSquared.Core.Tests.Lindblad;

public class PalindromeResidualScalingTests
{
    [Theory]
    [InlineData(2, HamiltonianClass.Main, 1.0)]      // (2-1)·4^0 = 1
    [InlineData(3, HamiltonianClass.Main, 8.0)]      // (3-1)·4^1 = 8
    [InlineData(4, HamiltonianClass.Main, 48.0)]     // (4-1)·4^2 = 48
    [InlineData(5, HamiltonianClass.Main, 256.0)]    // (5-1)·4^3 = 256
    [InlineData(2, HamiltonianClass.SingleBody, 1.0)] // (2·2-3)·4^0 = 1
    [InlineData(3, HamiltonianClass.SingleBody, 12.0)] // (2·3-3)·4^1 = 12
    public void FactorChain_MatchesClosedForm(int N, HamiltonianClass cls, double expected)
    {
        Assert.Equal(expected, PalindromeResidualScaling.FactorChain(N, cls), 12);
    }

    [Theory]
    [InlineData(2, HamiltonianClass.Main, 8.0)]              // 4·2/1 = 8
    [InlineData(3, HamiltonianClass.Main, 6.0)]              // 4·3/2 = 6
    [InlineData(4, HamiltonianClass.Main, 16.0 / 3.0)]       // 4·4/3
    [InlineData(2, HamiltonianClass.SingleBody, 12.0)]       // 4·(2·2−1)/(2·2−3) = 12/1 = 12
    [InlineData(3, HamiltonianClass.SingleBody, 20.0 / 3.0)] // 4·(2·3−1)/(2·3−3) = 20/3
    [InlineData(4, HamiltonianClass.SingleBody, 28.0 / 5.0)] // 4·(2·4−1)/(2·4−3) = 28/5
    public void AdjacentRatio_MatchesFormula(int N, HamiltonianClass cls, double expected)
    {
        Assert.Equal(expected, PalindromeResidualScaling.AdjacentRatio(N, cls), 12);
    }

    [Fact]
    public void FactorFromGraph_Chain_AgreesWithFactorChain()
    {
        // For an N-site chain: B = N-1, D2 = 4N-6.
        for (int N = 3; N <= 6; N++)
        {
            int B = N - 1;
            int D2 = 4 * N - 6;
            Assert.Equal(
                PalindromeResidualScaling.FactorChain(N, HamiltonianClass.Main),
                PalindromeResidualScaling.FactorFromGraph(N, B, D2, HamiltonianClass.Main), 12);
            Assert.Equal(
                PalindromeResidualScaling.FactorChain(N, HamiltonianClass.SingleBody),
                PalindromeResidualScaling.FactorFromGraph(N, B, D2, HamiltonianClass.SingleBody), 12);
        }
    }
}
