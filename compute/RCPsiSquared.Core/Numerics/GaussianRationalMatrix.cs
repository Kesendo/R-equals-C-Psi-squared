using System;

namespace RCPsiSquared.Core.Numerics;

/// <summary>Dense exact matrices over ℚ(i), held as <c>GaussianRational[,]</c>: the operations an
/// exact generator needs, entry by entry. Every entry is set explicitly (a default
/// <see cref="BigRational"/> has denominator 0); nothing is rounded.</summary>
public static class GaussianRationalMatrix
{
    public static GaussianRational[,] Zero(int rows, int cols)
    {
        var m = new GaussianRational[rows, cols];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++) m[i, j] = GaussianRational.Zero;
        return m;
    }

    public static GaussianRational[,] Zero(int n) => Zero(n, n);

    public static GaussianRational[,] Identity(int n)
    {
        var m = Zero(n);
        for (int i = 0; i < n; i++) m[i, i] = GaussianRational.One;
        return m;
    }

    public static GaussianRational[,] Multiply(GaussianRational[,] a, GaussianRational[,] b)
    {
        int n = a.GetLength(0), k = a.GetLength(1), p = b.GetLength(1);
        if (b.GetLength(0) != k) throw new ArgumentException("inner dimensions differ.");
        var r = Zero(n, p);
        for (int i = 0; i < n; i++)
            for (int t = 0; t < k; t++)
            {
                if (a[i, t].IsZero) continue;
                for (int j = 0; j < p; j++)
                    if (!b[t, j].IsZero) r[i, j] += a[i, t] * b[t, j];
            }
        return r;
    }

    public static GaussianRational[,] Add(GaussianRational[,] a, GaussianRational[,] b) => Combine(a, b, (u, v) => u + v);

    public static GaussianRational[,] Subtract(GaussianRational[,] a, GaussianRational[,] b) => Combine(a, b, (u, v) => u - v);

    public static GaussianRational[,] Scale(GaussianRational[,] a, GaussianRational s)
    {
        int n = a.GetLength(0), p = a.GetLength(1);
        var r = Zero(n, p);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++) r[i, j] = s * a[i, j];
        return r;
    }

    public static GaussianRational[,] Transpose(GaussianRational[,] a)
    {
        int n = a.GetLength(0), p = a.GetLength(1);
        var r = Zero(p, n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++) r[j, i] = a[i, j];
        return r;
    }

    public static GaussianRational[,] Adjoint(GaussianRational[,] a)
    {
        var r = Transpose(a);
        for (int i = 0; i < r.GetLength(0); i++)
            for (int j = 0; j < r.GetLength(1); j++) r[i, j] = r[i, j].Conjugate;
        return r;
    }

    public static bool AreEqual(GaussianRational[,] a, GaussianRational[,] b)
    {
        if (a.GetLength(0) != b.GetLength(0) || a.GetLength(1) != b.GetLength(1)) return false;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
                if (a[i, j] != b[i, j]) return false;
        return true;
    }

    /// <summary>The superoperator of ρ ↦ A·ρ·B on the row-major vec of a d × d matrix (index d·i + j
    /// for |i⟩⟨j|): entry [d·x + y, d·c + e] = A[x, c]·B[e, y].</summary>
    public static GaussianRational[,] Sandwich(GaussianRational[,] a, GaussianRational[,] b)
    {
        int d = a.GetLength(0);
        var s = Zero(d * d);
        for (int x = 0; x < d; x++)
            for (int c = 0; c < d; c++)
            {
                if (a[x, c].IsZero) continue;
                for (int e = 0; e < d; e++)
                    for (int y = 0; y < d; y++)
                        if (!b[e, y].IsZero) s[d * x + y, d * c + e] += a[x, c] * b[e, y];
            }
        return s;
    }

    private static GaussianRational[,] Combine(GaussianRational[,] a, GaussianRational[,] b,
                                               Func<GaussianRational, GaussianRational, GaussianRational> f)
    {
        int n = a.GetLength(0), p = a.GetLength(1);
        if (b.GetLength(0) != n || b.GetLength(1) != p) throw new ArgumentException("shapes differ.");
        var r = Zero(n, p);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++) r[i, j] = f(a[i, j], b[i, j]);
        return r;
    }
}
