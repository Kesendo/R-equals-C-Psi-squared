using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The self-mirror fixed point as an OBJECT that inherits its x/y/z frame from the system.
/// Two-sided gate: the object sits at Re λ = −σ (inherited center), and its own delta (the self-paired
/// sector) is populated for even N and EMPTY for odd N , a sign/parity error fails the odd case.</summary>
public class SelfMirrorObjectTests
{
    private static SelfMirrorObject Build(int n, double gamma = 0.1)
    {
        var h = new ChainSystem(n, 1.0, gamma, HamiltonianType.XY, TopologyKind.Chain).BuildHamiltonian();
        var channels = Enumerable.Range(0, n).Select(l => new ChannelRate($"q{l}", gamma)).ToList();
        return new SelfMirrorObject(new MirrorSystem(n, h, channels));
    }

    [Fact]
    public void Object_InheritsFrame_AndSitsAtMinusSigma()
    {
        var obj = Build(4, gamma: 0.1);
        Assert.Equal(4, obj.N);              // N inherited from the system, not owned
        Assert.Equal(0.4, obj.Sigma, 10);    // σ = Σγ = Nγ inherited
        Assert.Equal(-0.4, obj.Center, 10);  // the object sits at Re λ = −σ
    }

    [Theory]
    [InlineData(2, true)]
    [InlineData(4, true)]
    [InlineData(3, false)]
    [InlineData(5, false)]
    public void FixedPointSector_PopulatedIffEvenN(int n, bool populated)
    {
        int count = Build(n).SelfPairedCount;
        if (populated)
            Assert.True(count > 0, $"N={n} (even): the k=N/2 self-mirror sector should be populated, got {count}");
        else
            Assert.Equal(0, count);          // odd N: half-integer w_XY = N/2, the sector is empty
    }
}
