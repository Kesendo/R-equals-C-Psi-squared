using RCPsiSquared.Core.ChainSystems;

namespace RCPsiSquared.Core.Tests.ChainSystems;

public class BondTopologyTests
{
    [Fact]
    public void Chain_N5_HasFourBonds_LinearlyOrdered()
    {
        var sys = new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05);
        Assert.Equal(4, sys.Bonds.Count);
        Assert.Equal(new Bond(0, 1, 1.0), sys.Bonds[0]);
        Assert.Equal(new Bond(1, 2, 1.0), sys.Bonds[1]);
        Assert.Equal(new Bond(3, 4, 1.0), sys.Bonds[3]);
    }

    [Fact]
    public void Star_N5_HasFourSpokesFromCenter()
    {
        var sys = new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05, Topology: TopologyKind.Star);
        Assert.Equal(4, sys.Bonds.Count);
        Assert.All(sys.Bonds, b => Assert.Equal(0, b.Site1));
    }

    [Fact]
    public void Ring_N5_HasFiveBondsClosingTheLoop()
    {
        var sys = new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05, Topology: TopologyKind.Ring);
        Assert.Equal(5, sys.Bonds.Count);
        Assert.Equal(new Bond(4, 0, 1.0), sys.Bonds[4]);
    }

    [Fact]
    public void OverrideBonds_SupportsNonUniformJ()
    {
        var sys = new ChainSystem(N: 4, J: 1.0, GammaZero: 0.05);
        var custom = new[] { new Bond(0, 1, 0.5), new Bond(1, 2, 1.5), new Bond(2, 3, 0.7) };
        var withCustom = sys.WithBonds(custom);
        Assert.Equal(custom, withCustom.Bonds);
    }
}
