using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F2bXyChainSpectrumPi2InheritanceTests
{
    private static F2bXyChainSpectrumPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var qubitAnchor = new QubitDimensionalAnchorClaim();
        var f66 = new F66PoleModesPi2Inheritance(ladder, qubitAnchor);
        var f65 = new F65XxChainSpectrumPi2Inheritance(ladder, f66);
        return new F2bXyChainSpectrumPi2Inheritance(ladder, f65);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void HoppingFactor_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().HoppingFactor, precision: 14);
    }

    [Theory]
    // F2b: E_k = 2J·cos(πk/(N+1)) for k=1..N
    // N=3, J=1, k=1: 2·cos(π/4) = √2 ≈ 1.4142
    // N=3, J=1, k=2: 2·cos(π/2) = 0
    // N=3, J=1, k=3: 2·cos(3π/4) = -√2
    // N=5, J=1, k=3: 2·cos(π/2) = 0 (middle mode is always zero for odd N+1=even-N... actually N=5 → N+1=6 → middle k=3, cos(3π/6)=cos(π/2)=0)
    [InlineData(3, 1.0, 1, 1.4142135623730951)]
    [InlineData(3, 1.0, 2, 0.0)]
    [InlineData(3, 1.0, 3, -1.4142135623730951)]
    [InlineData(5, 1.0, 3, 0.0)]
    [InlineData(3, 0.5, 1, 0.7071067811865476)]   // J scaling
    public void Eigenvalue_MatchesClosedForm(int N, double J, int k, double expected)
    {
        Assert.Equal(expected, BuildClaim().Eigenvalue(N, J, k), precision: 12);
    }

    [Theory]
    // ψ_k(i) = √(2/(N+1))·sin(πk(i+1)/(N+1))
    // N=3, k=1, i=0: √(2/4)·sin(π/4) = (1/√2)·(1/√2) = 1/2
    // N=3, k=1, i=1: √(2/4)·sin(π/2) = (1/√2)·1 = 1/√2 ≈ 0.7071
    // N=3, k=1, i=2: √(2/4)·sin(3π/4) = 1/2 (palindromic to i=0)
    [InlineData(3, 1, 0, 0.5)]
    [InlineData(3, 1, 1, 0.7071067811865476)]
    [InlineData(3, 1, 2, 0.5)]
    public void EigenvectorAmplitude_MatchesClosedForm(int N, int k, int site, double expected)
    {
        Assert.Equal(expected, BuildClaim().EigenvectorAmplitude(N, k, site), precision: 12);
    }

    [Theory]
    [InlineData(2, 2)]
    [InlineData(3, 3)]
    [InlineData(5, 5)]
    [InlineData(10, 10)]
    public void ModeCount_EqualsN(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().ModeCount(N));
    }

    [Theory]
    [InlineData(3, 1, 0)]
    [InlineData(3, 1, 1)]
    [InlineData(3, 1, 2)]
    [InlineData(5, 2, 3)]
    public void EigenvectorMatchesF65_HoldsForVariousNK(int N, int k, int site)
    {
        // F2b's |ψ_k(i)|² should match F65's BondingModePopulation exactly.
        Assert.True(BuildClaim().EigenvectorMatchesF65(N, k, site));
    }

    [Theory]
    [InlineData(3, 1)]
    [InlineData(3, 2)]
    [InlineData(5, 3)]
    public void ObcBoundaryConditionsHold(int N, int k)
    {
        Assert.True(BuildClaim().ObcBoundaryConditionsHold(N, k));
    }

    [Fact]
    public void Eigenvalue_PalindromicPairing()
    {
        // E_k + E_{N+1-k} = 0 (XY spectrum is symmetric around 0).
        var f = BuildClaim();
        for (int N = 2; N <= 7; N++)
            for (int k = 1; k <= N; k++)
                Assert.Equal(0.0, f.Eigenvalue(N, 1.0, k) + f.Eigenvalue(N, 1.0, N + 1 - k), precision: 12);
    }

    [Fact]
    public void Eigenvalue_KOutOfRange_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().Eigenvalue(N: 3, J: 1.0, k: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().Eigenvalue(N: 3, J: 1.0, k: 4));
    }

    [Fact]
    public void EigenvectorAmplitude_SiteOutOfRange_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().EigenvectorAmplitude(N: 3, k: 1, site: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().EigenvectorAmplitude(N: 3, k: 1, site: 3));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f66 = new F66PoleModesPi2Inheritance(ladder, new QubitDimensionalAnchorClaim());
        var f65 = new F65XxChainSpectrumPi2Inheritance(ladder, f66);
        Assert.Throws<ArgumentNullException>(() =>
            new F2bXyChainSpectrumPi2Inheritance(null!, f65));
    }

    [Fact]
    public void Constructor_NullF65_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F2bXyChainSpectrumPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }
}
