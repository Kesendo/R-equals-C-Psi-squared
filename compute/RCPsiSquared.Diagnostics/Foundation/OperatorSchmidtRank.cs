using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The operator Schmidt rank of an operator held exactly on the Pauli strings, across any split of the
/// sites into two sets: the rank over ℚ of the coefficient matrix c[left string, right string]. The strings are
/// numbered in base 4 with site 0 the most significant digit and the letters 0 = I, 1 = X, 2 = Y, 3 = Z, the
/// numbering of <see cref="MirrorLocalityBoundWitness"/> and <see cref="TrianglePalindromeLocusWitness"/>, which both
/// read their ranks here. Rank one across every split is a product operator ⊗_l a_l.
///
/// <para>Why a rank decides locality: a map that factorises over the sites, linear or antilinear, after a frame
/// change per site, a transposition or a relabelling of the sites by conjugation, sends a product to a product; a map
/// of operator Schmidt rank r across a split sends a product to an operator of rank at most r across it. So a product
/// carried onto a non-product rules every such map out, the lemma of finding 1 of
/// <c>experiments/THE_CUBE_SITE_BY_SITE.md</c> and, applied to the identity, of F116 §4
/// (<c>docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md</c>).</para></summary>
public static class OperatorSchmidtRank
{
    /// <summary>The rank across the split of <paramref name="n"/> sites into <paramref name="left"/> and the rest,
    /// for an operator with coefficient vector <paramref name="v"/> of length 4^n.</summary>
    public static int Across(BigRational[] v, int n, IReadOnlyCollection<int> left)
    {
        if (v.Length != (int)Math.Pow(4, n)) throw new ArgumentException($"a coefficient vector on {n} sites has length 4^{n}; got {v.Length}.");
        var leftSites = left.Distinct().OrderBy(l => l).ToArray();
        if (leftSites.Any(l => l < 0 || l >= n)) throw new ArgumentException("a site of the split is outside the chain.");
        var rightSites = Enumerable.Range(0, n).Where(l => !leftSites.Contains(l)).ToArray();
        int rows = (int)Math.Pow(4, leftSites.Length), cols = (int)Math.Pow(4, rightSites.Length);
        var m = new BigRational[rows, cols];
        for (int i = 0; i < rows; i++) for (int k = 0; k < cols; k++) m[i, k] = BigRational.Zero;
        for (int s = 0; s < v.Length; s++)
        {
            if (v[s].IsZero) continue;
            int i = 0, k = 0;
            foreach (int l in leftSites) i = i * 4 + Digit(s, n, l);
            foreach (int l in rightSites) k = k * 4 + Digit(s, n, l);
            m[i, k] = v[s];
        }
        return RankOverQ(m);
    }

    /// <summary>The rank across the first <paramref name="cut"/> sites and the rest.</summary>
    public static int AcrossCut(BigRational[] v, int n, int cut) => Across(v, n, Enumerable.Range(0, cut).ToArray());

    private static int Digit(int s, int n, int site) => s / (int)Math.Pow(4, n - 1 - site) % 4;

    /// <summary>Exact rank over ℚ by Gaussian elimination in BigRational.</summary>
    public static int RankOverQ(BigRational[,] A)
    {
        int n = A.GetLength(0), m = A.GetLength(1);
        var rows = new BigRational[n][];
        for (int i = 0; i < n; i++) { rows[i] = new BigRational[m]; for (int k = 0; k < m; k++) rows[i][k] = A[i, k]; }
        int r = 0;
        for (int c = 0; c < m && r < n; c++)
        {
            int pr = -1;
            for (int i = r; i < n; i++) if (!rows[i][c].IsZero) { pr = i; break; }
            if (pr < 0) continue;
            (rows[r], rows[pr]) = (rows[pr], rows[r]);
            for (int i = r + 1; i < n; i++)
            {
                if (rows[i][c].IsZero) continue;
                var f = rows[i][c] / rows[r][c];
                for (int k = c; k < m; k++) rows[i][k] = rows[i][k] - f * rows[r][k];
            }
            r++;
        }
        return r;
    }
}
