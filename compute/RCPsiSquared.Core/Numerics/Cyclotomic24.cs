using System;
using System.Collections.Generic;
using System.Linq;

namespace RCPsiSquared.Core.Numerics;

/// <summary>Exact arithmetic in ℚ(ζ₂₄), ζ = e^{iπ/12}: rational coefficient vectors of length 8 on the basis
/// {1, ζ, …, ζ⁷} that <see cref="CyclotomicRing"/> uses at m = 24 (Φ₂₄ = x⁸ − x⁴ + 1; its polynomial and its table
/// of root powers are taken from there), so vector equality is equality of the numbers. The field holds i = ζ⁶,
/// √2 = ζ³ + ζ⁻³, √3 = ζ² + ζ⁻² and e^{ikπ/12} for every integer k, so every phase that is a multiple of 15°, and
/// cos and sin of those phases. Complex conjugation is the automorphism ζ ↦ ζ⁻¹. The sodium page's pulses live here:
/// the spin-3/2 ladder carries √3, the 90° pulse a factor 1/√2, the phase cycle steps of 45°.</summary>
public readonly struct Cyclotomic24 : IEquatable<Cyclotomic24>
{
    public const int Order = 24;
    public const int Degree = 8;

    private static readonly long[] Phi = CyclotomicRing.Phi(Order);
    private readonly BigRational[]? _c;

    private Cyclotomic24(BigRational[] c) => _c = c;

    /// <summary>The coefficient of ζ^j, j = 0 … 7.</summary>
    public BigRational this[int j] => _c is null ? BigRational.Zero : _c[j];

    public static Cyclotomic24 Zero => default;
    public static Cyclotomic24 One => FromRational(BigRational.One);
    public static Cyclotomic24 I => RootPower(6);
    public static Cyclotomic24 Sqrt2 => RootPower(3) + RootPower(-3);
    public static Cyclotomic24 Sqrt3 => RootPower(2) + RootPower(-2);

    public static Cyclotomic24 FromRational(BigRational q)
    {
        var c = Empty();
        c[0] = q;
        return new Cyclotomic24(c);
    }

    public static Cyclotomic24 FromGaussian(GaussianRational z) => FromRational(z.Re) + I * FromRational(z.Im);

    private static readonly Cyclotomic24[] Powers = Enumerable.Range(0, Order).Select(e =>
    {
        var acc = CyclotomicRing.Zero(Order);
        CyclotomicRing.AddRootPower(acc, Order, e, 1);
        return new Cyclotomic24(acc.Select(v => new BigRational(v)).ToArray());
    }).ToArray();

    /// <summary>ζ^k = e^{ikπ/12}, any integer k (the 24 powers are read once from <see cref="CyclotomicRing"/>'s
    /// table; a value never changes after construction, so they are shared).</summary>
    public static Cyclotomic24 RootPower(long k) => Powers[(int)(((k % Order) + Order) % Order)];

    /// <summary>cos(kπ/12) = (ζ^k + ζ^{−k})/2.</summary>
    public static Cyclotomic24 Cos(long k) => (RootPower(k) + RootPower(-k)) * new BigRational(1, 2);

    /// <summary>sin(kπ/12) = (ζ^k − ζ^{−k})/(2i).</summary>
    public static Cyclotomic24 Sin(long k) => (RootPower(k) - RootPower(-k)) * (I * new BigRational(-1, 2));

    public bool IsZero => _c is null || _c.All(v => v.IsZero);

    /// <summary>Whether the number is rational (only the coefficient of 1 non-zero), and which.</summary>
    public bool TryRational(out BigRational q)
    {
        q = this[0];
        for (int j = 1; j < Degree; j++) if (!this[j].IsZero) return false;
        return true;
    }

    /// <summary>The complex conjugate, ζ^j ↦ ζ^{−j}.</summary>
    public Cyclotomic24 Conjugate
    {
        get
        {
            var r = Zero;
            for (int j = 0; j < Degree; j++)
                if (!this[j].IsZero) r += RootPower(-j) * this[j];
            return r;
        }
    }

    public static Cyclotomic24 operator +(Cyclotomic24 a, Cyclotomic24 b)
    {
        var c = Empty();
        for (int j = 0; j < Degree; j++) c[j] = a[j] + b[j];
        return new Cyclotomic24(c);
    }

    public static Cyclotomic24 operator -(Cyclotomic24 a, Cyclotomic24 b)
    {
        var c = Empty();
        for (int j = 0; j < Degree; j++) c[j] = a[j] - b[j];
        return new Cyclotomic24(c);
    }

    public static Cyclotomic24 operator -(Cyclotomic24 a) => Zero - a;

    public static Cyclotomic24 operator *(Cyclotomic24 a, BigRational s)
    {
        var c = Empty();
        for (int j = 0; j < Degree; j++) c[j] = a[j] * s;
        return new Cyclotomic24(c);
    }

    public static Cyclotomic24 operator *(Cyclotomic24 a, Cyclotomic24 b)
    {
        if (a.IsZero || b.IsZero) return Zero;
        var prod = new BigRational[2 * Degree - 1];
        for (int e = 0; e < prod.Length; e++) prod[e] = BigRational.Zero;
        for (int i = 0; i < Degree; i++)
        {
            if (a[i].IsZero) continue;
            for (int j = 0; j < Degree; j++)
                if (!b[j].IsZero) prod[i + j] += a[i] * b[j];
        }
        for (int e = prod.Length - 1; e >= Degree; e--)          // fold by the monic Φ₂₄
        {
            var lead = prod[e];
            if (lead.IsZero) continue;
            prod[e] = BigRational.Zero;
            for (int j = 0; j < Degree; j++)
                if (Phi[j] != 0) prod[e - Degree + j] -= lead * Phi[j];
        }
        var c = Empty();
        Array.Copy(prod, c, Degree);
        return new Cyclotomic24(c);
    }

    /// <summary>The 8 × 8 rational matrix of multiplication by this number on the basis {1, ζ, …, ζ⁷}: column j holds
    /// the coefficients of this·ζ^j.</summary>
    public BigRational[,] RegularRepresentation()
    {
        var m = new BigRational[Degree, Degree];
        for (int j = 0; j < Degree; j++)
        {
            var col = this * RootPower(j);
            for (int r = 0; r < Degree; r++) m[r, j] = col[r];
        }
        return m;
    }

    public static bool operator ==(Cyclotomic24 a, Cyclotomic24 b) => a.Equals(b);
    public static bool operator !=(Cyclotomic24 a, Cyclotomic24 b) => !a.Equals(b);

    public bool Equals(Cyclotomic24 other)
    {
        for (int j = 0; j < Degree; j++) if (this[j] != other[j]) return false;
        return true;
    }

    public override bool Equals(object? obj) => obj is Cyclotomic24 z && Equals(z);

    public override int GetHashCode()
    {
        var h = new HashCode();
        for (int j = 0; j < Degree; j++) h.Add(this[j]);
        return h.ToHashCode();
    }

    public override string ToString()
    {
        var terms = new List<string>();
        for (int j = 0; j < Degree; j++)
            if (!this[j].IsZero) terms.Add(j == 0 ? $"{this[j]}" : $"({this[j]})·ζ{(j == 1 ? "" : Superscript(j))}");
        return terms.Count == 0 ? "0" : string.Join(" + ", terms);
    }

    private static string Superscript(int j) => "⁰¹²³⁴⁵⁶⁷"[j].ToString();

    private static BigRational[] Empty()
    {
        var c = new BigRational[Degree];
        for (int j = 0; j < Degree; j++) c[j] = BigRational.Zero;
        return c;
    }
}

/// <summary>Matrices over ℚ(ζ₂₄) (<see cref="Cyclotomic24"/>), rows × columns, with the exact rank read through the
/// regular representation: an m × n matrix over the field of rank r is an 8m × 8n rational matrix of rank 8r.</summary>
public static class Cyclotomic24Matrix
{
    public static Cyclotomic24[,] Zero(int rows, int cols) => new Cyclotomic24[rows, cols];   // default is zero

    public static Cyclotomic24[,] Zero(int n) => Zero(n, n);

    public static Cyclotomic24[,] Identity(int n)
    {
        var m = Zero(n);
        for (int i = 0; i < n; i++) m[i, i] = Cyclotomic24.One;
        return m;
    }

    public static Cyclotomic24[,] FromGaussian(GaussianRational[,] a)
    {
        var m = Zero(a.GetLength(0), a.GetLength(1));
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) m[i, j] = Cyclotomic24.FromGaussian(a[i, j]);
        return m;
    }

    public static Cyclotomic24[,] Multiply(Cyclotomic24[,] a, Cyclotomic24[,] b)
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

    public static Cyclotomic24[,] Multiply(params Cyclotomic24[][,] factors) => factors.Aggregate(Multiply);

    public static Cyclotomic24[,] Add(Cyclotomic24[,] a, Cyclotomic24[,] b) => Combine(a, b, (u, v) => u + v);

    public static Cyclotomic24[,] Subtract(Cyclotomic24[,] a, Cyclotomic24[,] b) => Combine(a, b, (u, v) => u - v);

    public static Cyclotomic24[,] Scale(Cyclotomic24[,] a, Cyclotomic24 s) => Combine(a, a, (u, _) => s * u);

    public static Cyclotomic24[,] Scale(Cyclotomic24[,] a, BigRational s) => Combine(a, a, (u, _) => u * s);

    public static Cyclotomic24[,] Transpose(Cyclotomic24[,] a)
    {
        var r = Zero(a.GetLength(1), a.GetLength(0));
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) r[j, i] = a[i, j];
        return r;
    }

    public static Cyclotomic24[,] Adjoint(Cyclotomic24[,] a)
    {
        var r = Zero(a.GetLength(1), a.GetLength(0));
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) r[j, i] = a[i, j].Conjugate;
        return r;
    }

    /// <summary>A ⊗ B, the row index of A the high digit.</summary>
    public static Cyclotomic24[,] Kron(Cyclotomic24[,] a, Cyclotomic24[,] b)
    {
        int ar = a.GetLength(0), ac = a.GetLength(1), br = b.GetLength(0), bc = b.GetLength(1);
        var r = Zero(ar * br, ac * bc);
        for (int i = 0; i < ar; i++)
            for (int j = 0; j < ac; j++)
            {
                if (a[i, j].IsZero) continue;
                for (int k = 0; k < br; k++)
                    for (int l = 0; l < bc; l++) r[i * br + k, j * bc + l] = a[i, j] * b[k, l];
            }
        return r;
    }

    public static Cyclotomic24 Trace(Cyclotomic24[,] a)
    {
        var t = Cyclotomic24.Zero;
        for (int i = 0; i < a.GetLength(0); i++) t += a[i, i];
        return t;
    }

    public static bool AreEqual(Cyclotomic24[,] a, Cyclotomic24[,] b)
    {
        if (a.GetLength(0) != b.GetLength(0) || a.GetLength(1) != b.GetLength(1)) return false;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
                if (a[i, j] != b[i, j]) return false;
        return true;
    }

    public static bool IsZero(Cyclotomic24[,] a)
    {
        foreach (var z in a) if (!z.IsZero) return false;
        return true;
    }

    /// <summary>The exact rank over ℚ(ζ₂₄): the rational rank of the regular representation, divided by 8 (the
    /// division is checked).</summary>
    public static int Rank(Cyclotomic24[,] a)
    {
        int rows = a.GetLength(0), cols = a.GetLength(1), d = Cyclotomic24.Degree;
        var big = new BigRational[d * rows, d * cols];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
            {
                var block = a[i, j].RegularRepresentation();
                for (int r = 0; r < d; r++)
                    for (int c = 0; c < d; c++) big[d * i + r, d * j + c] = block[r, c];
            }
        int rank = d * cols - BigRationalLinearAlgebra.Nullspace(big).Count;
        if (rank % d != 0) throw new InvalidOperationException($"the regular representation has rank {rank}, not a multiple of {d}.");
        return rank / d;
    }

    private static Cyclotomic24[,] Combine(Cyclotomic24[,] a, Cyclotomic24[,] b, Func<Cyclotomic24, Cyclotomic24, Cyclotomic24> f)
    {
        if (a.GetLength(0) != b.GetLength(0) || a.GetLength(1) != b.GetLength(1)) throw new ArgumentException("shapes differ.");
        var r = Zero(a.GetLength(0), a.GetLength(1));
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) r[i, j] = f(a[i, j], b[i, j]);
        return r;
    }
}
