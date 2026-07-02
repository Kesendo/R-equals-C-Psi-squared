using System;
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
    /// dim_sym × dim_sym. Eigenvalues are 2·λ_k.
    ///
    /// <para><b>METRIC WARNING (non-orthonormal basis).</b> The S₂-symmetric basis is a reflection-orbit-sum
    /// basis; orbits have size 1 (reflection-fixed) or 2, so the Hilbert-Schmidt metric on it is diag(1,2),
    /// NOT the identity. The block is therefore related to the physically HS-orthonormal representation by a
    /// non-unitary diagonal similarity S = diag(√orbit-size). SIMILARITY-INVARIANT quantities are UNAFFECTED:
    /// eigenvalues, the Galois/monodromy group, the residual-vs-AT split, and the diabolic-vs-defective
    /// character (M₂ = λI is basis-free) are all correct on this block. But any INNER-PRODUCT-SENSITIVE
    /// diagnostic (e.g. "is operator X scalar on a sub-plane", the Theorem-A twin-scalar / D-half test; a
    /// departure-from-normality restricted to a coalescing plane; a Petermann/phase-rigidity factor) MUST be
    /// computed in an HS-orthonormal basis (use the raw coherence block <see cref="WeightCoherenceBlock.Build"/>,
    /// or apply the congruence S before restricting), or it reports a spurious non-scalarity at ODD nBlock
    /// (N=5, N=7; even nBlock has a uniform metric and is unaffected). This cost one false twin-scalar reading
    /// before it was caught, 2026-07-02.</para></summary>
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

    // zz(c) = Σ_{bond (b,b+1)} ⟨c|Z_bZ_{b+1}|c⟩ = Σ_b (+1 if bits b,b+1 equal, −1 if differ).
    private static int Zz(int n, int c)
    {
        int s = 0;
        for (int b = 0; b < n - 1; b++)
            s += (((c >> b) & 1) == ((c >> (b + 1)) & 1)) ? 1 : -1;
        return s;
    }

    /// <summary>The per-orbit ZZ-frequency weights zzDiag = Σ_b⟨ket|Z_bZ_{b+1}|ket⟩ − Σ_b⟨bra|Z_bZ_{b+1}|bra⟩
    /// of the S_2-sym (SE,DE) block, in the SAME orbit ordering as <see cref="BuildTwoTimesSymBlock"/> (the basis
    /// + orbit construction is replicated verbatim so index r here is index r there). The XXZ anisotropy adds the
    /// diagonal frequency −i·qΔ·zzDiag (the ZZ term is Hermitian, so the absorption-theorem real rate is
    /// untouched). The generator is reflection-invariant (a chain's bond set is reversal-symmetric), so zzDiag is
    /// constant on each reflection orbit and any representative gives it. Real-integer valued. Consumed by
    /// PathKMonodromyScout's exact-residual XXZ port (the ×2-cleared generator is −2i·zzDiag, matching 2M).</summary>
    public static int[] BuildZzFrequencyDiag(int nBlock)
    {
        var dePairs = new List<(int J, int K)>();
        for (int a = 0; a < nBlock; a++)
            for (int b = a + 1; b < nBlock; b++)
                dePairs.Add((a, b));
        var pairIndex = new Dictionary<(int, int), int>();
        for (int pj = 0; pj < dePairs.Count; pj++) pairIndex[dePairs[pj]] = pj;

        var basis = new List<(int I, int Pair)>();
        for (int i = 0; i < nBlock; i++)
            for (int pj = 0; pj < dePairs.Count; pj++)
                basis.Add((i, pj));
        int nb = basis.Count;
        var idxOf = new Dictionary<(int, int), int>();
        for (int t = 0; t < nb; t++) idxOf[basis[t]] = t;

        var zz = new List<int>();
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
            if (t2 != t) handled[t2] = true;
            // ket = SE at site i (mask 1<<i); bra = DE at sites j,k (mask (1<<j)|(1<<k)).
            zz.Add(Zz(nBlock, 1 << i) - Zz(nBlock, (1 << j) | (1 << k)));
        }
        return zz.ToArray();
    }

    /// <summary>The FULL (SE,DE) block M(q) (R-even AND R-odd, not the symmetrized sector), dim =
    /// nBlock·C(nBlock,2), in the (SE site i, DE pair (j,k)) basis (i outer, pair inner; same ordering as
    /// the raw L of <see cref="BuildTwoTimesSymBlock"/>): diagonal −2 (overlap, i ∈ {j,k}) / −6 (no overlap),
    /// SE/ket hop −i·2q, DE/bra hop +i·2q. The eigenvalues of this block's R-even orthonormal projection are
    /// the scout's λ (<see cref="BuildTwoTimesSymBlock"/> is the ×2-cleared R-even sector of this). The
    /// substrate for the reflection-parity / odd-N real-q diabolic mechanism: pair it with
    /// <see cref="ReflectionPermutation"/> to split the block into R-even and R-odd sectors.</summary>
    public static Complex[,] BuildFullBlock(int nBlock, Complex q)
    {
        var dePairs = new List<(int J, int K)>();
        for (int a = 0; a < nBlock; a++)
            for (int b = a + 1; b < nBlock; b++)
                dePairs.Add((a, b));
        var pairIndex = new Dictionary<(int, int), int>();
        for (int pj = 0; pj < dePairs.Count; pj++) pairIndex[dePairs[pj]] = pj;

        var basis = new List<(int I, int Pair)>();
        for (int i = 0; i < nBlock; i++)
            for (int pj = 0; pj < dePairs.Count; pj++)
                basis.Add((i, pj));
        int nb = basis.Count;
        var idxOf = new Dictionary<(int, int), int>();
        for (int t = 0; t < nb; t++) idxOf[basis[t]] = t;

        Complex hop = Complex.ImaginaryOne * 2.0 * q;       // +i·2q  (DE/bra hop)
        Complex neg = -hop;                                  // −i·2q  (SE/ket hop)

        var L = new Complex[nb, nb];
        for (int col = 0; col < nb; col++)
        {
            var (i, pj) = basis[col];
            var (j, k) = dePairs[pj];

            foreach (int i2 in new[] { i - 1, i + 1 })          // SE hop (ket)
                if (i2 >= 0 && i2 < nBlock)
                    L[idxOf[(i2, pj)], col] += neg;

            foreach (int nj in new[] { j - 1, j + 1 })          // DE hop on j (bra)
                if (nj >= 0 && nj < nBlock && nj != k)
                {
                    var np = nj < k ? (nj, k) : (k, nj);
                    L[idxOf[(i, pairIndex[np])], col] += hop;
                }
            foreach (int nk in new[] { k - 1, k + 1 })          // DE hop on k (bra)
                if (nk >= 0 && nk < nBlock && nk != j)
                {
                    var np = j < nk ? (j, nk) : (nk, j);
                    L[idxOf[(i, pairIndex[np])], col] += hop;
                }

            L[col, col] += (i == j || i == k)                   // γ = 1 diagonal: −2 (overlap) / −6
                ? new Complex(-2, 0)
                : new Complex(-6, 0);
        }
        return L;
    }

    /// <summary>The site-reflection R: (SE i, DE pair (j,k)) ↦ (nBlock−1−i, sorted(nBlock−1−k, nBlock−1−j)),
    /// as a permutation of the full-block basis index (same ordering as <see cref="BuildFullBlock"/>). An
    /// involution that commutes with the block; its fixed points (perm[t] == t) are the reflection-fixed
    /// singletons, present only at odd nBlock (center SE site × self-mirror DE pair), count (nBlock−1)/2.</summary>
    public static int[] ReflectionPermutation(int nBlock)
    {
        var dePairs = new List<(int J, int K)>();
        for (int a = 0; a < nBlock; a++)
            for (int b = a + 1; b < nBlock; b++)
                dePairs.Add((a, b));
        var pairIndex = new Dictionary<(int, int), int>();
        for (int pj = 0; pj < dePairs.Count; pj++) pairIndex[dePairs[pj]] = pj;

        var basis = new List<(int I, int Pair)>();
        for (int i = 0; i < nBlock; i++)
            for (int pj = 0; pj < dePairs.Count; pj++)
                basis.Add((i, pj));
        int nb = basis.Count;
        var idxOf = new Dictionary<(int, int), int>();
        for (int t = 0; t < nb; t++) idxOf[basis[t]] = t;

        var perm = new int[nb];
        for (int t = 0; t < nb; t++)
        {
            var (i, pj) = basis[t];
            var (j, k) = dePairs[pj];
            int mi = nBlock - 1 - i, mj = nBlock - 1 - j, mk = nBlock - 1 - k;
            var mp = mj < mk ? (mj, mk) : (mk, mj);
            perm[t] = idxOf[(mi, pairIndex[mp])];
        }
        return perm;
    }
}
