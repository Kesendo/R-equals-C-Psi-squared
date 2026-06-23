using System.Collections.Generic;
using System.Numerics;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>The general path-k (SE, DE) S_2-symmetric Liouvillian sub-block at integer q0, built
/// EXACTLY over Z[i] via the integer mirror basis e_n + e_m (BᵀB diagonal with entries ∈ {1,2}).
/// The raw symmetric block M = (BᵀB)⁻¹ BᵀLB carries denominator-2 entries, so this returns the
/// ×2-cleared integer block 2M — a Gaussian-INTEGER matrix whose eigenvalues are 2·λ_k (twice the
/// physical H_B-mixed decay rates). The factor of 2 is tracked by the caller; the Galois group is
/// scale-invariant, so the S_d certificate is unaffected. Port of
/// simulations/f89_pathk_galois.py:build_pathk_sym_over_qi. path-k ⟹ nBlock = k + 1.</summary>
public static class F89PathKSeDeBlock
{
    /// <summary>The ×2-cleared Gaussian-integer S_2-sym (SE,DE) block at integer q0 (γ = 1),
    /// dim_sym × dim_sym. Eigenvalues are 2·λ_k.</summary>
    public static GaussianInteger[,] BuildTwoTimesSymBlock(int q0, int nBlock)
    {
        // DE pairs (j < k) in lexicographic order, and their index lookup.
        var dePairs = new List<(int J, int K)>();
        for (int a = 0; a < nBlock; a++)
            for (int b = a + 1; b < nBlock; b++)
                dePairs.Add((a, b));
        var pairIndex = new Dictionary<(int, int), int>();
        for (int pj = 0; pj < dePairs.Count; pj++) pairIndex[dePairs[pj]] = pj;

        // Basis (SE index i, DE pair index pj): i outer, pair inner — matches the Python order.
        var basis = new List<(int I, int Pair)>();
        for (int i = 0; i < nBlock; i++)
            for (int pj = 0; pj < dePairs.Count; pj++)
                basis.Add((i, pj));
        int nb = basis.Count;
        var idxOf = new Dictionary<(int, int), int>();
        for (int t = 0; t < nb; t++) idxOf[basis[t]] = t;

        GaussianInteger hop = new(0, 2 * q0);          // +i·2q (DE/bra hop)
        GaussianInteger negHop = new(0, -2 * q0);      // −i·2q (SE/ket hop)

        // Build L over Z[i] (nb × nb).
        var L = new GaussianInteger[nb, nb];
        for (int col = 0; col < nb; col++)
        {
            var (i, pj) = basis[col];
            var (j, k) = dePairs[pj];

            foreach (int i2 in new[] { i - 1, i + 1 })          // SE hop (ket): −i·2q
                if (i2 >= 0 && i2 < nBlock)
                    L[idxOf[(i2, pj)], col] += negHop;

            foreach (int nj in new[] { j - 1, j + 1 })          // DE hop on j (bra): +i·2q
                if (nj >= 0 && nj < nBlock && nj != k)
                {
                    var np = nj < k ? (nj, k) : (k, nj);
                    L[idxOf[(i, pairIndex[np])], col] += hop;
                }
            foreach (int nk in new[] { k - 1, k + 1 })          // DE hop on k (bra): +i·2q
                if (nk >= 0 && nk < nBlock && nk != j)
                {
                    var np = j < nk ? (j, nk) : (nk, j);
                    L[idxOf[(i, pairIndex[np])], col] += hop;
                }

            L[col, col] += (i == j || i == k)                   // γ = 1 diagonal: −2 (overlap) / −6
                ? new GaussianInteger(-2, 0)
                : new GaussianInteger(-6, 0);
        }

        // S_2-mirror symmetrization: orbit {t, t2} (or {t}) under i ↦ nBlock−1−i, pair mirrored.
        var orbits = new List<int[]>();
        var handled = new bool[nb];
        for (int t = 0; t < nb; t++)
        {
            if (handled[t]) continue;
            var (i, pj) = basis[t];
            var (j, k) = dePairs[pj];
            int mj = nBlock - 1 - j, mk = nBlock - 1 - k;
            var mp = mj < mk ? (mj, mk) : (mk, mj);
            int t2 = idxOf[(nBlock - 1 - i, pairIndex[mp])];
            handled[t] = true;
            if (t2 != t) { handled[t2] = true; orbits.Add(new[] { t, t2 }); }
            else orbits.Add(new[] { t });
        }

        // 2M[r,c] = (2/d_r) · (BᵀLB)[r,c], where (BᵀLB)[r,c] = Σ_{s∈orbit_r,u∈orbit_c} L[s,u].
        int dim = orbits.Count;
        var twoM = new GaussianInteger[dim, dim];
        for (int r = 0; r < dim; r++)
        {
            int scale = 2 / orbits[r].Length;               // d_r ∈ {1,2} ⟹ 2/d_r ∈ {2,1}
            for (int c = 0; c < dim; c++)
            {
                GaussianInteger sum = GaussianInteger.Zero;
                foreach (int s in orbits[r])
                    foreach (int u in orbits[c])
                        sum += L[s, u];
                twoM[r, c] = new GaussianInteger(sum.Re * scale, sum.Im * scale);
            }
        }
        return twoM;
    }
}
