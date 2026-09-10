using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Small exact F_101 arithmetic for the 90-dimensional F163 certificate.
/// Matrix sizes are below the characteristic, so Faddeev-LeVerrier divides only by units.
/// Rank and inverse use the same finite-field elimination; no floating rank decisions.</summary>
internal static class RouteBN6Modular
{
    internal const int Prime = 101;
    internal static int Mod(long value) => (int)((value % Prime + Prime) % Prime);
    private static int Inv(int value) => value == 0
        ? throw new InvalidOperationException("Singular F163 modular denominator.")
        : (int)CrossFormCertificate.Inv(value, Prime);

    internal static int[,] Identity(int n)
    {
        var result = new int[n, n];
        for (int i = 0; i < n; i++) result[i, i] = 1;
        return result;
    }

    internal static int[,] Add(int[,] a, int[,] b, int scaleB = 1)
    {
        var result = new int[a.GetLength(0), a.GetLength(1)];
        for (int i = 0; i < result.GetLength(0); i++)
        for (int j = 0; j < result.GetLength(1); j++)
            result[i, j] = Mod(a[i, j] + (long)scaleB * b[i, j]);
        return result;
    }

    internal static int[,] Scale(int[,] a, int scalar)
    {
        var result = new int[a.GetLength(0), a.GetLength(1)];
        for (int i = 0; i < result.GetLength(0); i++)
        for (int j = 0; j < result.GetLength(1); j++)
            result[i, j] = Mod((long)scalar * a[i, j]);
        return result;
    }

    internal static int[,] Multiply(int[,] a, int[,] b)
    {
        var result = new int[a.GetLength(0), b.GetLength(1)];
        for (int i = 0; i < result.GetLength(0); i++)
        for (int j = 0; j < result.GetLength(1); j++)
        {
            long sum = 0;
            for (int k = 0; k < a.GetLength(1); k++) sum += a[i, k] * b[k, j];
            result[i, j] = Mod(sum);
        }
        return result;
    }

    internal static bool Equal(int[,] a, int[,] b) =>
        a.GetLength(0) == b.GetLength(0) && a.GetLength(1) == b.GetLength(1)
        && a.Cast<int>().SequenceEqual(b.Cast<int>());
    internal static bool IsZero(int[,] a) => a.Cast<int>().All(c => c == 0);
    internal static int Trace(int[,] a) => Mod(Enumerable.Range(0, a.GetLength(0)).Sum(i => a[i, i]));
    internal static int Discriminant(int[,] a) => Mod(2 * Trace(Multiply(a, a)) - Trace(a) * Trace(a));

    internal static RouteBN6ProfileInvariants Invariants(int[,] a, int[,] b)
    {
        int alpha = Discriminant(a), gamma = Discriminant(b);
        int beta = Mod(4 * Trace(Multiply(a, b)) - 2 * Trace(a) * Trace(b));
        return new(alpha, beta, gamma, Mod(beta * beta - 4 * alpha * gamma));
    }

    private static int Reduce(int[,] a, int pivotColumns)
    {
        int rank = 0, rows = a.GetLength(0), columns = a.GetLength(1);
        for (int column = 0; column < pivotColumns && rank < rows; column++)
        {
            int pivot = rank;
            while (pivot < rows && a[pivot, column] == 0) pivot++;
            if (pivot == rows) continue;
            for (int j = 0; j < columns; j++) (a[rank, j], a[pivot, j]) = (a[pivot, j], a[rank, j]);
            int inverse = Inv(a[rank, column]);
            for (int j = 0; j < columns; j++) a[rank, j] = Mod(a[rank, j] * inverse);
            for (int i = 0; i < rows; i++)
            {
                if (i == rank) continue;
                int factor = a[i, column];
                if (factor == 0) continue;
                for (int j = 0; j < columns; j++) a[i, j] = Mod(a[i, j] - factor * a[rank, j]);
            }
            rank++;
        }
        return rank;
    }

    internal static int Rank(int[,] a) => Reduce((int[,])a.Clone(), a.GetLength(1));
    internal static int[,] Inverse(int[,] a)
    {
        int n = a.GetLength(0);
        var augmented = new int[n, 2 * n];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++) augmented[i, j] = a[i, j];
            augmented[i, n + i] = 1;
        }
        if (Reduce(augmented, n) != n) throw new InvalidOperationException("Singular F163 reduced resolvent.");
        var result = new int[n, n];
        for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) result[i, j] = augmented[i, n + j];
        return result;
    }

    internal static int[] Characteristic(int[,] a)
    {
        int n = a.GetLength(0);
        if (n >= Prime || a.GetLength(1) != n) throw new ArgumentException("F163 characteristic needs a square matrix smaller than 101.");
        var identity = Identity(n);
        var b = identity;
        var result = new int[n + 1];
        result[n] = 1;
        for (int k = 1; k <= n; k++)
        {
            b = Multiply(a, b);
            result[n - k] = Mod(-Trace(b) * Inv(k));
            b = Add(b, identity, result[n - k]);
        }
        return result;
    }

    internal static int[,] MatrixPolynomial(int[] coefficients, int[,] a)
    {
        var result = new int[a.GetLength(0), a.GetLength(0)];
        var identity = Identity(a.GetLength(0));
        for (int k = coefficients.Length - 1; k >= 0; k--)
            result = Add(Multiply(result, a), identity, coefficients[k]);
        return result;
    }

    internal static int Evaluate(int[] coefficients, int x)
    {
        int result = 0;
        for (int k = coefficients.Length - 1; k >= 0; k--) result = Mod(result * x + coefficients[k]);
        return result;
    }

    internal static int[] Derivative(int[] coefficients) => coefficients.Length <= 1 ? new[] { 0 }
        : Trim(Enumerable.Range(1, coefficients.Length - 1).Select(k => Mod(k * coefficients[k])).ToArray());

    private static int[] Trim(int[] a)
    {
        int length = a.Length;
        while (length > 1 && a[length - 1] == 0) length--;
        return a.Take(length).ToArray();
    }

    internal static (int[] Quotient, int[] Remainder) Divide(int[] dividend, int[] divisor)
    {
        divisor = Trim(divisor);
        if (divisor.All(c => c == 0)) throw new DivideByZeroException();
        var remainder = Trim((int[])dividend.Clone());
        var quotient = new int[Math.Max(1, remainder.Length - divisor.Length + 1)];
        int inverse = Inv(divisor[^1]);
        for (int k = remainder.Length - divisor.Length; k >= 0; k--)
        {
            int c = Mod(remainder[k + divisor.Length - 1] * inverse);
            quotient[k] = c;
            for (int j = 0; j < divisor.Length; j++) remainder[k + j] = Mod(remainder[k + j] - c * divisor[j]);
        }
        return (Trim(quotient), Trim(remainder));
    }

    internal static int[] Gcd(int[] a, int[] b)
    {
        while (b.Any(c => c != 0)) (a, b) = (b, Divide(a, b).Remainder);
        int inverse = Inv(a[^1]);
        return a.Select(c => Mod(c * inverse)).ToArray();
    }

    internal static int[] MultiplyPolynomials(int[] a, int[] b)
    {
        var result = new int[a.Length + b.Length - 1];
        for (int i = 0; i < a.Length; i++)
        for (int j = 0; j < b.Length; j++) result[i + j] = Mod(result[i + j] + a[i] * b[j]);
        return result;
    }
}
