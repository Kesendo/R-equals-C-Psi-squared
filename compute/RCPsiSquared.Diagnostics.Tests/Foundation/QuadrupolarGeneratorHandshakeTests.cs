using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Quadrupolar;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using M = RCPsiSquared.Core.Numerics.GaussianRationalMatrix;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The shared builder in Core against the chain witness's string builder: at zero coupling the
/// chain's generator is the Kronecker sum of two nuclei, entry by entry, ==. (The sodium witness builds
/// on the shared builder itself; its second route at 3/2 is the Pauli-string construction in Core,
/// <c>Spin32StringsTests</c>, and its gate on the Lamb shift reads the term the generator consumes.)</summary>
public class QuadrupolarGeneratorHandshakeTests
{
    private static GaussianRational[,] Shared(long j0, long j1, long j2, long q, long w) =>
        new QuadrupolarGenerator(SpinSpace.Of(3), new QuadrupolarRates(j0, j1, j2), SecularHamiltonian.Of(q, w)).Physical();

    [Theory]
    [InlineData(1, 2, 1, 2, 3)]
    [InlineData(0, 0, 1, 0, 2)]
    [InlineData(2, -1, 2, 3, 1)]
    public void ChainBuilder_AtZeroCoupling_IsTheKroneckerSumOfTwoNuclei(long q, long w, long j0, long j1, long j2)
    {
        var one = Shared(j0, j1, j2, q, w);
        var sum = GeneratorComposition.KroneckerSum(one, 4, one, 4);
        Assert.True(M.AreEqual(SodiumChainWitness.UncoupledGeneratorMatrix(q, w, j0, j1, j2), sum));
    }
}
