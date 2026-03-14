using System.Numerics;
using System.Runtime.InteropServices;

namespace RCPsiSquared.Compute;

/// <summary>
/// Direct MKL P/Invoke with pinned arrays to bypass the .NET 2GB marshaling limit.
/// </summary>
public static class MklDirect
{
    [DllImport("libMathNetNumericsMKL", CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)]
    private static extern unsafe int z_eigen(
        [MarshalAs(UnmanagedType.Bool)] bool isSymmetric,
        int n,
        Complex* a,
        Complex* vectors,
        Complex* values,
        Complex* d);

    /// <summary>
    /// Eigenvalues directly from column-major Complex[] array.
    /// No MathNet Matrix involved. Minimum memory path for N >= 7.
    /// WARNING: input array 'a' is destroyed by LAPACK.
    /// </summary>
    public static unsafe Complex[] EigenvaluesRaw(Complex[] a, int n)
    {
        var vectors = new Complex[(long)n * n];
        var values = new Complex[n];
        var d = new Complex[(long)n * n];

        fixed (Complex* pA = a)
        fixed (Complex* pVectors = vectors)
        fixed (Complex* pValues = values)
        fixed (Complex* pD = d)
        {
            int result = z_eigen(false, n, pA, pVectors, pValues, pD);
            if (result != 0)
                throw new InvalidOperationException($"MKL z_eigen failed with code {result}");
        }

        return values;
    }
}
