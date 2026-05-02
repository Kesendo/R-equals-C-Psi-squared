using RCPsiSquared.Core.CoherenceBlocks;

namespace RCPsiSquared.Core.Tests.CoherenceBlocks;

public class ChromaticityTests
{
    // F74: c(n, N) = min(n, N-1-n) + 1
    [Theory]
    [InlineData(7, 0, 1)]
    [InlineData(7, 1, 2)]
    [InlineData(7, 2, 3)]
    [InlineData(7, 3, 4)]   // central, c_max for odd N
    [InlineData(7, 4, 3)]
    [InlineData(7, 6, 1)]
    [InlineData(8, 3, 4)]   // even N: two adjacent c_max blocks
    [InlineData(8, 4, 4)]
    [InlineData(5, 1, 2)]   // c=2 N=5 (used in Statement 2 c=2 verification)
    [InlineData(8, 1, 2)]   // c=2 N=8
    public void Compute_MatchesF74Formula(int N, int n, int expected)
    {
        Assert.Equal(expected, Chromaticity.Compute(N, n));
    }

    [Fact]
    public void HammingDistances_AreOddAndAscending()
    {
        var hds = Chromaticity.HammingDistances(N: 7, n: 3); // c = 4
        Assert.Equal(new[] { 1, 3, 5, 7 }, hds);
    }

    [Fact]
    public void Compute_RejectsInvalidN()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Chromaticity.Compute(0, 0));
    }

    [Theory]
    [InlineData(5, -1)]
    [InlineData(5, 5)]
    public void Compute_RejectsInvalidLowerPopcount(int N, int n)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Chromaticity.Compute(N, n));
    }
}
