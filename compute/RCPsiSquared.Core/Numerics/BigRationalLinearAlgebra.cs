using System;
using System.Collections.Generic;

namespace RCPsiSquared.Core.Numerics;

/// <summary>Exact linear algebra over <see cref="BigRational"/>: matrix product/transpose, the right
/// nullspace (via reduced row echelon form), and square linear solve. The exact-rational foundation
/// for the F89 full-D AT reconstruction — the rate-confined invariant subspace W = nullspace of
/// [P_Uc·K^m], and the restriction M|W = solve((WᵀW)·R = WᵀMW).</summary>
public static class BigRationalLinearAlgebra
{
    /// <summary>Matrix product a·b.</summary>
    public static BigRational[,] Multiply(BigRational[,] a, BigRational[,] b)
    {
        int n = a.GetLength(0), m = a.GetLength(1), p = b.GetLength(1);
        if (b.GetLength(0) != m) throw new ArgumentException("matrix dimension mismatch.");
        var r = new BigRational[n, p];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++)
            {
                BigRational s = BigRational.Zero;
                for (int k = 0; k < m; k++) s += a[i, k] * b[k, j];
                r[i, j] = s;
            }
        return r;
    }

    /// <summary>Transpose.</summary>
    public static BigRational[,] Transpose(BigRational[,] a)
    {
        int n = a.GetLength(0), m = a.GetLength(1);
        var r = new BigRational[m, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) r[j, i] = a[i, j];
        return r;
    }

    /// <summary>A basis (each a column vector, length = #columns) for the right nullspace
    /// { x : a·x = 0 }, computed from the reduced row echelon form.</summary>
    public static List<BigRational[]> Nullspace(BigRational[,] a)
    {
        int rows = a.GetLength(0), cols = a.GetLength(1);
        var m = (BigRational[,])a.Clone();
        var pivotCol = new List<int>();
        int r = 0;
        for (int c = 0; c < cols && r < rows; c++)
        {
            int piv = -1;
            for (int i = r; i < rows; i++) if (!m[i, c].IsZero) { piv = i; break; }
            if (piv < 0) continue;
            if (piv != r)
                for (int j = 0; j < cols; j++) (m[r, j], m[piv, j]) = (m[piv, j], m[r, j]);
            var inv = BigRational.One / m[r, c];
            for (int j = 0; j < cols; j++) m[r, j] = m[r, j] * inv;
            for (int i = 0; i < rows; i++)
            {
                if (i == r || m[i, c].IsZero) continue;
                var f = m[i, c];
                for (int j = 0; j < cols; j++) m[i, j] = m[i, j] - f * m[r, j];
            }
            pivotCol.Add(c);
            r++;
        }

        var pivotSet = new HashSet<int>(pivotCol);
        var basis = new List<BigRational[]>();
        for (int free = 0; free < cols; free++)
        {
            if (pivotSet.Contains(free)) continue;
            var v = new BigRational[cols];
            for (int j = 0; j < cols; j++) v[j] = BigRational.Zero;
            v[free] = BigRational.One;
            for (int k = 0; k < pivotCol.Count; k++) v[pivotCol[k]] = -m[k, free];
            basis.Add(v);
        }
        return basis;
    }

    /// <summary>Solve a·x = b for a square invertible a (b may have several columns); Gauss-Jordan.</summary>
    public static BigRational[,] Solve(BigRational[,] a, BigRational[,] b)
    {
        int n = a.GetLength(0);
        if (a.GetLength(1) != n) throw new ArgumentException("matrix must be square.");
        int p = b.GetLength(1);
        var aug = new BigRational[n, n + p];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++) aug[i, j] = a[i, j];
            for (int j = 0; j < p; j++) aug[i, n + j] = b[i, j];
        }
        for (int c = 0; c < n; c++)
        {
            int piv = -1;
            for (int i = c; i < n; i++) if (!aug[i, c].IsZero) { piv = i; break; }
            if (piv < 0) throw new InvalidOperationException("matrix is singular.");
            if (piv != c)
                for (int j = 0; j < n + p; j++) (aug[c, j], aug[piv, j]) = (aug[piv, j], aug[c, j]);
            var inv = BigRational.One / aug[c, c];
            for (int j = 0; j < n + p; j++) aug[c, j] = aug[c, j] * inv;
            for (int i = 0; i < n; i++)
            {
                if (i == c || aug[i, c].IsZero) continue;
                var f = aug[i, c];
                for (int j = 0; j < n + p; j++) aug[i, j] = aug[i, j] - f * aug[c, j];
            }
        }
        var x = new BigRational[n, p];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++) x[i, j] = aug[i, n + j];
        return x;
    }
}
