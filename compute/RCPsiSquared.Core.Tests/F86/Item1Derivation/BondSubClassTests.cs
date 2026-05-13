using RCPsiSquared.Core.F86.Item1Derivation;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class BondSubClassTests
{
    [Theory]
    [InlineData(5, 0, 2.40, BondSubClass.Endpoint)]
    [InlineData(5, 1, 1.50, BondSubClass.Flanking)]
    [InlineData(6, 0, 2.52, BondSubClass.Endpoint)]
    [InlineData(6, 1, 1.65, BondSubClass.Flanking)]
    [InlineData(6, 2, 1.44, BondSubClass.CentralSelfPaired)]
    [InlineData(7, 2, 1.54, BondSubClass.Mid)]
    [InlineData(7, 1, 7.27, BondSubClass.Orbit2Escape)]
    [InlineData(8, 3, 16.79, BondSubClass.CentralEscapeOrbit3)]
    public void Classify_ReturnsExpectedSubClass(int n, int b, double qPeak, BondSubClass expected)
    {
        Assert.Equal(expected, BondSubClassExtensions.Classify(n, b, qPeak));
    }
}
