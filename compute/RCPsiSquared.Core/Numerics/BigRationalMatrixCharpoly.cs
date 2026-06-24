using System;

namespace RCPsiSquared.Core.Numerics;

/// <summary>The characteristic polynomial det(λI − A) of a square matrix over BigRational, via the
/// division-free Samuelson-Berkowitz algorithm (the same recursion as
/// <see cref="GaussianMatrixCharpoly"/>, over the field of rationals). Coefficients lowest-first,
/// monic. Used to read the charpoly of the rate-sector restriction M|W in the full-D AT
/// reconstruction.</summary>
public static class BigRationalMatrixCharpoly
{
    public static BigRational[] Characteristic(BigRational[,] matrix)
    {
        int n = matrix.GetLength(0);
        if (matrix.GetLength(1) != n) throw new ArgumentException("matrix must be square.");
        if (n == 0) return new[] { BigRational.One };

        var p = new[] { BigRational.One };
        for (int i = n; i >= 1; i--)
        {
            int m = n - i + 1, sub = m - 1;
            BigRational aii = matrix[i - 1, i - 1];

            var q = new BigRational[m + 1];
            q[0] = BigRational.One;
            q[1] = -aii;
            if (sub > 0)
            {
                var vec = new BigRational[sub];
                for (int row = 0; row < sub; row++) vec[row] = matrix[i + row, i - 1];
                for (int k = 0; k < sub; k++)
                {
                    BigRational term = BigRational.Zero;
                    for (int j = 0; j < sub; j++) term += matrix[i - 1, i + j] * vec[j];
                    q[k + 2] = -term;
                    if (k < sub - 1)
                    {
                        var next = new BigRational[sub];
                        for (int row = 0; row < sub; row++)
                        {
                            BigRational s = BigRational.Zero;
                            for (int col = 0; col < sub; col++) s += matrix[i + row, i + col] * vec[col];
                            next[row] = s;
                        }
                        vec = next;
                    }
                }
            }

            var pNew = new BigRational[m + 1];
            for (int r = 0; r <= m; r++)
            {
                BigRational s = BigRational.Zero;
                for (int c = 0; c < m; c++)
                {
                    int d = r - c;
                    if (d >= 0 && d <= m) s += q[d] * p[c];
                }
                pNew[r] = s;
            }
            p = pNew;
        }

        Array.Reverse(p);
        return p;
    }
}
