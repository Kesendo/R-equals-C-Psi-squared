using System;

namespace RCPsiSquared.Core.Numerics;

/// <summary>The characteristic polynomial det(λI − A) of a square matrix over Z[i], via the
/// division-free Samuelson-Berkowitz algorithm (valid over any commutative ring, so it keeps the
/// coefficients exact in Z[i] with no p&gt;n constraint — unlike Faddeev-LeVerrier, which divides
/// by 1..n). Coefficients are returned lowest-first (index = power), monic (leading = 1).
/// This is the F89 path-k charpoly engine feeding the AT-factor division and the DDF.</summary>
public static class GaussianMatrixCharpoly
{
    /// <summary>det(λI − A) as coefficients lowest-first (length n+1, leading coefficient 1).</summary>
    public static GaussianInteger[] Characteristic(GaussianInteger[,] matrix)
    {
        int n = matrix.GetLength(0);
        if (matrix.GetLength(1) != n) throw new ArgumentException("matrix must be square.");
        if (n == 0) return new[] { GaussianInteger.One };

        // p = C_1 · C_2 · … · C_n, each C_i a lower-triangular Toeplitz matrix built from the
        // bottom-right principal submatrix M_i. Built right-to-left; result is highest-first.
        var p = new[] { GaussianInteger.One };
        for (int i = n; i >= 1; i--)
        {
            int m = n - i + 1;                              // size of M_i (bottom-right block)
            int sub = m - 1;                                // size of M_{i+1} / R / S
            GaussianInteger aii = matrix[i - 1, i - 1];

            // q[0]=1, q[1]=−a_ii, q[t]=−(R · M_{i+1}^{t−2} · S) for t = 2..m
            var q = new GaussianInteger[m + 1];
            q[0] = GaussianInteger.One;
            q[1] = -aii;
            if (sub > 0)
            {
                var vec = new GaussianInteger[sub];         // M_{i+1}^k · S, starting at k=0 (= S)
                for (int row = 0; row < sub; row++) vec[row] = matrix[i + row, i - 1];
                for (int k = 0; k < sub; k++)
                {
                    GaussianInteger term = GaussianInteger.Zero;        // R · (M_{i+1}^k · S)
                    for (int j = 0; j < sub; j++) term += matrix[i - 1, i + j] * vec[j];
                    q[k + 2] = -term;
                    if (k < sub - 1)
                    {
                        var next = new GaussianInteger[sub];            // M_{i+1} · vec
                        for (int row = 0; row < sub; row++)
                        {
                            GaussianInteger s = GaussianInteger.Zero;
                            for (int col = 0; col < sub; col++)
                                s += matrix[i + row, i + col] * vec[col];
                            next[row] = s;
                        }
                        vec = next;
                    }
                }
            }

            // pNew = C_i · p, with C_i the (m+1)×m Toeplitz C_i[r,c] = q[r−c] (0 ≤ r−c ≤ m).
            var pNew = new GaussianInteger[m + 1];
            for (int r = 0; r <= m; r++)
            {
                GaussianInteger s = GaussianInteger.Zero;
                for (int c = 0; c < m; c++)
                {
                    int d = r - c;
                    if (d >= 0 && d <= m) s += q[d] * p[c];
                }
                pNew[r] = s;
            }
            p = pNew;
        }

        Array.Reverse(p);                                   // highest-first → lowest-first
        return p;
    }
}
