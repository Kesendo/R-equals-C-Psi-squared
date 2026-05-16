using System.Numerics;

namespace RCPsiSquared.Core.BlockSpectrum.JordanWigner;

/// <summary>Shared low-level complex-vector operations used by the Krylov-subspace
/// primitives <see cref="JwSlaterPairArnoldiEig"/> and <see cref="JwSlaterPairShiftInvertArnoldi"/>.
/// Hand-rolled for-loops (no LINQ) to avoid per-call delegate allocations on the BiCGStab
/// and MGS hot paths; thread-safe (no static state).</summary>
internal static class KrylovOps
{
    /// <summary>Σ |x_i|² — the squared L² norm of a complex vector.</summary>
    public static double NormSquared(Complex[] x)
    {
        double sum = 0.0;
        int n = x.Length;
        for (int i = 0; i < n; i++) sum += x[i].Real * x[i].Real + x[i].Imaginary * x[i].Imaginary;
        return sum;
    }

    /// <summary>Σ conj(a_i) · b_i — the Hermitian inner product ⟨a, b⟩.</summary>
    public static Complex ConjugateDot(Complex[] a, Complex[] b)
    {
        Complex sum = Complex.Zero;
        int n = a.Length;
        for (int i = 0; i < n; i++) sum += Complex.Conjugate(a[i]) * b[i];
        return sum;
    }

    /// <summary>y ← y + α·x (in-place AXPY).</summary>
    public static void AxpyInPlace(Complex[] y, Complex[] x, Complex alpha)
    {
        int n = y.Length;
        for (int i = 0; i < n; i++) y[i] += alpha * x[i];
    }

    /// <summary>Allocate a length-<paramref name="dim"/> complex vector with entries
    /// drawn from a real/imag-independent uniform distribution on [−0.5, +0.5) and
    /// renormalise to unit L² norm. Deterministic given <paramref name="seed"/>.</summary>
    public static Complex[] RandomNormalized(int dim, int seed)
    {
        var rng = new Random(seed);
        var v = new Complex[dim];
        double normSq = 0.0;
        for (int i = 0; i < dim; i++)
        {
            v[i] = new Complex(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
            normSq += v[i].Real * v[i].Real + v[i].Imaginary * v[i].Imaginary;
        }
        double inv = 1.0 / Math.Sqrt(normSq);
        for (int i = 0; i < dim; i++) v[i] *= inv;
        return v;
    }
}
