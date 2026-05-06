using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Runtime.F1Family;

namespace RCPsiSquared.Runtime.Tests.F1Family;

public class ChainSystemPrimitiveTests
{
    [Fact]
    public void ChainSystemPrimitive_ExposesUnderlyingChain()
    {
        var primitive = new ChainSystemPrimitive(
            new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05,
                HType: HamiltonianType.XY, Topology: TopologyKind.Chain));

        Assert.Equal(5, primitive.System.N);
        Assert.Equal(1.0, primitive.System.J);
        Assert.Equal(0.05, primitive.System.GammaZero);
        Assert.Equal(Tier.Tier1Derived, primitive.Tier);
    }
}
