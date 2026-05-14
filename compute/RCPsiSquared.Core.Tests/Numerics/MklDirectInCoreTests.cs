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

    [Fact]
    public void LuFactorizeRaw_ThenLuSolveRaw_RecoversKnownSolution()
    {
        // Strongly diagonally dominant 4x4 complex A (column-major), non-symmetric.
        // Known x; b = A x; LU-factorize then solve must recover x.
        int n = 4;
        var A = new Complex[]
        {
            new(10, 1), new(1, 0), new(2, -1), new(0, 1),     // col 0
            new(1, 1), new(10, -2), new(1, 0), new(1, 1),     // col 1
            new(0, -1), new(2, 1), new(10, 0), new(1, -1),    // col 2
            new(1, 0), new(0, 1), new(1, 1), new(10, 2),      // col 3
        };
        var x = new Complex[] { new(1, 2), new(-1, 1), new(3, 0), new(0, -2) };

        // b = A x  (column-major: A[i,j] = A[j*n + i]).
        var b = new Complex[n];
        for (int i = 0; i < n; i++)
        {
            Complex s = Complex.Zero;
            for (int j = 0; j < n; j++) s += A[j * n + i] * x[j];
            b[i] = s;
        }

        var lu = (Complex[])A.Clone();
        var ipiv = new int[n];
        MklDirect.LuFactorizeRaw(lu, n, ipiv);
        MklDirect.LuSolveRaw(lu, n, ipiv, b, conjugateTranspose: false);

        for (int i = 0; i < n; i++)
            Assert.True((b[i] - x[i]).Magnitude < 1e-10,
                $"x[{i}]: expected {x[i]}, got {b[i]}");
    }

    [Fact]
    public void LuSolveRaw_ConjugateTranspose_SolvesAHermitianTransposeSystem()
    {
        int n = 4;
        var A = new Complex[]
        {
            new(10, 1), new(1, 0), new(2, -1), new(0, 1),     // col 0
            new(1, 1), new(10, -2), new(1, 0), new(1, 1),     // col 1
            new(0, -1), new(2, 1), new(10, 0), new(1, -1),    // col 2
            new(1, 0), new(0, 1), new(1, 1), new(10, 2),      // col 3
        };
        var x = new Complex[] { new(1, 2), new(-1, 1), new(3, 0), new(0, -2) };

        // b = A^H x.  A^H[i,j] = conj(A[j,i]); A[j,i] column-major is A[i*n + j].
        var b = new Complex[n];
        for (int i = 0; i < n; i++)
        {
            Complex s = Complex.Zero;
            for (int j = 0; j < n; j++) s += Complex.Conjugate(A[i * n + j]) * x[j];
            b[i] = s;
        }

        var lu = (Complex[])A.Clone();
        var ipiv = new int[n];
        MklDirect.LuFactorizeRaw(lu, n, ipiv);
        MklDirect.LuSolveRaw(lu, n, ipiv, b, conjugateTranspose: true);

        for (int i = 0; i < n; i++)
            Assert.True((b[i] - x[i]).Magnitude < 1e-10,
                $"x[{i}]: expected {x[i]}, got {b[i]}");
    }
}
