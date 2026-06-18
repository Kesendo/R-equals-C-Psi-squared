using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class TrichotomyWitnessTests
{
    [Fact]
    public void Ctor_RejectsOutOfRangeN() =>
        Assert.Throws<ArgumentOutOfRangeException>(() => new TrichotomyWitness(n: 9));
}
