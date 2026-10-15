using System;
using System.Linq;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>A square matrix of Gaussian integers held as two long arrays (real part, imaginary part), for the
/// exact small-N checks of the F87 batteries (the golden router's split, the converse's lit chiral). Core's
/// <see cref="Core.Numerics.GaussianInteger"/> is BigInteger-backed; these checks stay far inside long range
/// (single-digit letters; the callers stay at N ≤ 5) and run a matrix-unit loop, so long arithmetic is used. Site 0 is the most
/// significant bit of a basis index, matching <see cref="Product"/>'s Kronecker order.</summary>
internal readonly record struct GMat(long[,] Re, long[,] Im)
{
    public int D => Re.GetLength(0);

    public static GMat Zero(int d) => new(new long[d, d], new long[d, d]);

    public static GMat Identity(int d)
    {
        var m = Zero(d);
        for (int i = 0; i < d; i++) m.Re[i, i] = 1;
        return m;
    }

    public static GMat Mul(GMat x, GMat y)
    {
        int d = x.D;
        var p = Zero(d);
        for (int i = 0; i < d; i++)
            for (int k = 0; k < d; k++)
            {
                long xr = x.Re[i, k], xi = x.Im[i, k];
                if (xr == 0 && xi == 0) continue;
                for (int j = 0; j < d; j++)
                {
                    long yr = y.Re[k, j], yi = y.Im[k, j];
                    if (yr == 0 && yi == 0) continue;
                    p.Re[i, j] += xr * yr - xi * yi;
                    p.Im[i, j] += xr * yi + xi * yr;
                }
            }
        return p;
    }

    public static GMat Kron(GMat x, GMat y)
    {
        int dx = x.D, dy = y.D;
        var p = Zero(dx * dy);
        for (int i = 0; i < dx; i++)
            for (int j = 0; j < dx; j++)
                for (int k = 0; k < dy; k++)
                    for (int l = 0; l < dy; l++)
                    {
                        p.Re[i * dy + k, j * dy + l] = x.Re[i, j] * y.Re[k, l] - x.Im[i, j] * y.Im[k, l];
                        p.Im[i * dy + k, j * dy + l] = x.Re[i, j] * y.Im[k, l] + x.Im[i, j] * y.Re[k, l];
                    }
        return p;
    }

    /// <summary>The letters I, X, Y, Z and the integer metallic directions at c = 3/2 (r = 2): a = 2X + Y and
    /// the perpendicular b = X − 2Y.</summary>
    public static GMat Letter(char c) => c switch
    {
        'I' => new(new long[,] { { 1, 0 }, { 0, 1 } }, new long[2, 2]),
        'X' => new(new long[,] { { 0, 1 }, { 1, 0 } }, new long[2, 2]),
        'Y' => new(new long[2, 2], new long[,] { { 0, -1 }, { 1, 0 } }),
        'Z' => new(new long[,] { { 1, 0 }, { 0, -1 } }, new long[2, 2]),
        'a' => new(new long[,] { { 0, 2 }, { 2, 0 } }, new long[,] { { 0, -1 }, { 1, 0 } }),   // 2X + Y
        'b' => new(new long[,] { { 0, 1 }, { 1, 0 } }, new long[,] { { 0, 2 }, { -2, 0 } }),   // X − 2Y
        _ => throw new ArgumentOutOfRangeException(nameof(c)),
    };

    public static GMat Product(string letters) =>
        letters.Skip(1).Aggregate(Letter(letters[0]), (m, c) => Kron(m, Letter(c)));

    public static GMat Combine(GMat x, long kx, GMat y, long ky)
    {
        int d = x.D;
        var s = Zero(d);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                s.Re[i, j] = kx * x.Re[i, j] + ky * y.Re[i, j];
                s.Im[i, j] = kx * x.Im[i, j] + ky * y.Im[i, j];
            }
        return s;
    }

    public static GMat Transpose(GMat x)
    {
        int d = x.D;
        var t = Zero(d);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                t.Re[j, i] = x.Re[i, j];
                t.Im[j, i] = x.Im[i, j];
            }
        return t;
    }

    /// <summary>−i·x.</summary>
    public static GMat TimesMinusI(GMat x) => new(Copy(x.Im), Negate(x.Re));

    public static bool Equal(GMat x, GMat y)
    {
        for (int i = 0; i < x.D; i++)
            for (int j = 0; j < x.D; j++)
                if (x.Re[i, j] != y.Re[i, j] || x.Im[i, j] != y.Im[i, j]) return false;
        return true;
    }

    public static (long Re, long Im) Trace(GMat x)
    {
        long re = 0, im = 0;
        for (int i = 0; i < x.D; i++) { re += x.Re[i, i]; im += x.Im[i, i]; }
        return (re, im);
    }

    private static long[,] Copy(long[,] a) => (long[,])a.Clone();

    private static long[,] Negate(long[,] a)
    {
        var n = (long[,])a.Clone();
        for (int i = 0; i < n.GetLength(0); i++)
            for (int j = 0; j < n.GetLength(1); j++) n[i, j] = -n[i, j];
        return n;
    }
}
