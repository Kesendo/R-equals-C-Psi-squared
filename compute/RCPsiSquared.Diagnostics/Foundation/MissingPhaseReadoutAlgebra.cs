using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

// Local exact arithmetic for the readout witness: a + b sqrt(2), a,b in Q(i).
// Keeping both coefficients avoids rounding the uniform frequency or its residue.
internal readonly record struct ReadoutScalar(GaussianRational A, GaussianRational B)
{
    internal static readonly ReadoutScalar Zero = new(GaussianRational.Zero, GaussianRational.Zero);
    internal static readonly ReadoutScalar One = new(GaussianRational.One, GaussianRational.Zero);
    internal static readonly ReadoutScalar I = new(GaussianRational.I, GaussianRational.Zero);
    internal static readonly ReadoutScalar RootTwo = new(GaussianRational.Zero, GaussianRational.One);
    internal bool IsZero => A.IsZero && B.IsZero;
    internal ReadoutScalar Conjugate => new(A.Conjugate, B.Conjugate);
    internal BigRational Rational => B.IsZero && A.Im.IsZero ? A.Re
        : throw new InvalidOperationException($"Expected a real rational, got {this}.");

    public static implicit operator ReadoutScalar(int value) => new(value, GaussianRational.Zero);
    public static implicit operator ReadoutScalar(BigRational value) => new(value, GaussianRational.Zero);
    public static implicit operator ReadoutScalar(GaussianRational value) => new(value, GaussianRational.Zero);
    public static ReadoutScalar operator +(ReadoutScalar x, ReadoutScalar y) => new(x.A + y.A, x.B + y.B);
    public static ReadoutScalar operator -(ReadoutScalar x, ReadoutScalar y) => new(x.A - y.A, x.B - y.B);
    public static ReadoutScalar operator -(ReadoutScalar x) => new(-x.A, -x.B);
    public static ReadoutScalar operator *(ReadoutScalar x, ReadoutScalar y) =>
        new(x.A * y.A + 2 * x.B * y.B, x.A * y.B + x.B * y.A);
    public static ReadoutScalar operator /(ReadoutScalar x, ReadoutScalar y)
    {
        var denominator = y.A * y.A - 2 * y.B * y.B;
        if (denominator.IsZero) throw new DivideByZeroException("Division by zero in Q(sqrt(2),i).");
        return new((x.A * y.A - 2 * x.B * y.B) / denominator,
                   (x.B * y.A - x.A * y.B) / denominator);
    }
    public override string ToString() => B.IsZero ? A.ToString()
        : A.IsZero ? $"({B})sqrt(2)" : $"{A} + ({B})sqrt(2)";
}

internal static class MissingPhaseReadoutAlgebra
{
    internal static ReadoutScalar[,] Zero(int rows, int columns)
    {
        var a = new ReadoutScalar[rows, columns];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++) a[i, j] = ReadoutScalar.Zero;
        return a;
    }

    internal static ReadoutScalar[,] Lift(GaussianRational[,] input)
    {
        var a = Zero(input.GetLength(0), input.GetLength(1));
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) a[i, j] = input[i, j];
        return a;
    }

    internal static ReadoutScalar[,] Add(ReadoutScalar[,] a, ReadoutScalar[,] b)
    {
        var c = Zero(a.GetLength(0), a.GetLength(1));
        for (int i = 0; i < c.GetLength(0); i++)
            for (int j = 0; j < c.GetLength(1); j++) c[i, j] = a[i, j] + b[i, j];
        return c;
    }

    internal static ReadoutScalar[,] Scale(ReadoutScalar value, ReadoutScalar[,] a)
    {
        var c = Zero(a.GetLength(0), a.GetLength(1));
        for (int i = 0; i < c.GetLength(0); i++)
            for (int j = 0; j < c.GetLength(1); j++) c[i, j] = value * a[i, j];
        return c;
    }

    internal static ReadoutScalar[,] Subtract(ReadoutScalar[,] a, ReadoutScalar[,] b) => Add(a, Scale(-1, b));

    internal static ReadoutScalar[,] Multiply(ReadoutScalar[,] a, ReadoutScalar[,] b)
    {
        if (a.GetLength(1) != b.GetLength(0)) throw new ArgumentException("Matrix dimensions disagree.");
        var c = Zero(a.GetLength(0), b.GetLength(1));
        for (int i = 0; i < c.GetLength(0); i++)
            for (int k = 0; k < a.GetLength(1); k++)
                if (!a[i, k].IsZero)
                    for (int j = 0; j < c.GetLength(1); j++)
                        if (!b[k, j].IsZero) c[i, j] += a[i, k] * b[k, j];
        return c;
    }

    internal static ReadoutScalar[,] Adjoint(ReadoutScalar[,] a)
    {
        var c = Zero(a.GetLength(1), a.GetLength(0));
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) c[j, i] = a[i, j].Conjugate;
        return c;
    }

    internal static ReadoutScalar Inner(ReadoutScalar[,] a, ReadoutScalar[,] b)
    {
        ReadoutScalar result = 0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) result += a[i, j].Conjugate * b[i, j];
        return result;
    }

    internal static ReadoutScalar Trace(ReadoutScalar[,] a)
    {
        ReadoutScalar trace = 0;
        for (int i = 0; i < a.GetLength(0); i++) trace += a[i, i];
        return trace;
    }

    internal static int NonzeroCount(ReadoutScalar[,] a)
    {
        int count = 0;
        foreach (var value in a) if (!value.IsZero) count++;
        return count;
    }

    internal static ReadoutScalar[,] Outer(ReadoutScalar[,] a, ReadoutScalar[,] b) => Multiply(a, Adjoint(b));

    // Exact Gaussian elimination with pivot search; no numerical tolerance or pseudoinverse.
    internal static ReadoutScalar[,] Solve(ReadoutScalar[,] matrix, ReadoutScalar[,] rhs)
    {
        int n = matrix.GetLength(0), columns = rhs.GetLength(1);
        if (matrix.GetLength(1) != n || rhs.GetLength(0) != n)
            throw new ArgumentException("A square system with matching right-hand side is required.");
        var a = Zero(n, n + columns);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++) a[i, j] = matrix[i, j];
            for (int j = 0; j < columns; j++) a[i, n + j] = rhs[i, j];
        }
        for (int k = 0; k < n; k++)
        {
            int pivot = k;
            while (pivot < n && a[pivot, k].IsZero) pivot++;
            if (pivot == n) throw new InvalidOperationException("Singular exact bordered system.");
            if (pivot != k)
                for (int j = 0; j < n + columns; j++) (a[k, j], a[pivot, j]) = (a[pivot, j], a[k, j]);
            var divisor = a[k, k];
            for (int j = k; j < n + columns; j++) a[k, j] /= divisor;
            for (int i = 0; i < n; i++)
            {
                if (i == k || a[i, k].IsZero) continue;
                var factor = a[i, k];
                for (int j = k; j < n + columns; j++) a[i, j] -= factor * a[k, j];
            }
        }
        var solution = Zero(n, columns);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < columns; j++) solution[i, j] = a[i, n + j];
        return solution;
    }
}
