using RCPsiSquared.Core.F86.Item1Derivation;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class BondSubClassTests
{
    [Theory]
    [InlineData(5, 0, 2.517071370614458, BondSubClass.Endpoint)]
    [InlineData(5, 1, 1.4831232668717471, BondSubClass.Flanking)]
    [InlineData(6, 0, 2.580604998574177, BondSubClass.Endpoint)]
    [InlineData(6, 1, 1.6366804166405813, BondSubClass.Flanking)]
    [InlineData(6, 2, 1.4283108372718716, BondSubClass.CentralSelfPaired)]
    [InlineData(7, 2, 1.5010592087436854, BondSubClass.Mid)]
    [InlineData(7, 1, 7.255447920341961, BondSubClass.Orbit2Escape)]
    [InlineData(8, 3, 20.091690478319904, BondSubClass.CentralEscapeOrbit3)]
    public void Classify_ReturnsExpectedSubClass(int n, int b, double qPeak, BondSubClass expected)
    {
        Assert.Equal(expected, BondSubClassExtensions.Classify(n, b, qPeak));
    }
}
