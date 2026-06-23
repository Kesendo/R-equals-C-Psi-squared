using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>The exact Z[i] value type the path-k charpoly/isolation is built on. The only
/// non-trivial operation is multiplication (the Gaussian product rule).</summary>
public class GaussianIntegerTests
{
    [Fact]
    public void Multiply_FollowsTheGaussianProductRule()
    {
        // (2 + 3i)(1 − i) = 2 − 2i + 3i − 3i² = 5 + i.
        var product = new GaussianInteger(2, 3) * new GaussianInteger(1, -1);
        Assert.Equal(new GaussianInteger(5, 1), product);
    }

    [Fact]
    public void Multiply_ISquared_IsMinusOne()
    {
        Assert.Equal(new GaussianInteger(-1, 0), GaussianInteger.I * GaussianInteger.I);
    }
}
