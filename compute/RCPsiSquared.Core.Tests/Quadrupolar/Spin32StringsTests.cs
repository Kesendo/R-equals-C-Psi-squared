using System;
using RCPsiSquared.Core.Quadrupolar;
using Xunit;
using M = RCPsiSquared.Core.Numerics.GaussianRationalMatrix;

namespace RCPsiSquared.Core.Tests.Quadrupolar;

/// <summary>The second route at 3/2: the generator built from the five Pauli-string jumps and H in
/// strings equals the spin-matrix builder's physical view entry by entry, including the Lamb shift
/// (−3K₂·ZI − 3K₁·IZ, read off T₂q†T₂q by hand).</summary>
public class Spin32StringsTests
{
    [Theory]
    [InlineData(1, 0, 0, 3, 6, 0, 0)]
    [InlineData(0, 1, 2, 2, 3, 0, 0)]
    [InlineData(1, 2, 3, 1, 1, 0, 0)]
    [InlineData(1, 4, 2, 3, 1, 0, 0)]
    [InlineData(1, 2, 4, 3, 1, 0, 0)]
    [InlineData(1, 4, 0, 0, 1, 0, 0)]
    [InlineData(0, -1, 2, 1, 2, 0, 0)]
    [InlineData(1, 2, 3, 1, 1, 1, 3)]
    [InlineData(2, 1, 5, -2, 4, -1, 2)]
    public void StringsGenerator_EqualsThePhysicalView(long j0, long j1, long j2, long q, long w, long k1, long k2)
    {
        var rates = new QuadrupolarRates(j0, j1, j2);
        var h = SecularHamiltonian.Of(q, w).WithLambShift(k1, k2);
        var physical = new QuadrupolarGenerator(SpinSpace.Of(3), rates, h).Physical();
        Assert.True(M.AreEqual(Spin32Strings.Generator(rates, h), physical));
    }

    [Fact]
    public void Strings_RefuseTheLock()
    {
        Assert.Throws<NotSupportedException>(() => Spin32Strings.Hamiltonian(SecularHamiltonian.Of(1, 1).WithLock(1)));
    }
}
