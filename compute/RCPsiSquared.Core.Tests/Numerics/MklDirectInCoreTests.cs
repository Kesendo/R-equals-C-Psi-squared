using System.Numerics;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

public class MklDirectInCoreTests
{
    [Fact]
    public void MklDirect_IsAccessibleFromCoreNamespace()
    {
        // 4x4 identity, column-major. Eigenvalues should all be 1.
        var a = new Complex[16];
        for (int i = 0; i < 4; i++) a[i * 4 + i] = Complex.One;
        var (values, _, _) = MklDirect.EigenvaluesLeftRightDirectRaw(a, 4);
        Assert.Equal(4, values.Length);
        foreach (var v in values)
            Assert.True((v - Complex.One).Magnitude < 1e-10, $"eigenvalue {v} should be 1");
    }

    [Fact]
    public void GcAllowVeryLargeObjects_PermitsComplexArrayOver2GB()
    {
        // 145M complex = 2.32 GB, over the default 2 GB single-object cap.
        // With gcAllowVeryLargeObjects this allocates; without it throws OutOfMemoryException.
        long count = 145_000_000L;
        Complex[] big = new Complex[count];
        Assert.Equal(count, big.LongLength);
        big[count - 1] = new Complex(1.0, 2.0);   // touch the far end
        Assert.Equal(new Complex(1.0, 2.0), big[count - 1]);
    }
}
