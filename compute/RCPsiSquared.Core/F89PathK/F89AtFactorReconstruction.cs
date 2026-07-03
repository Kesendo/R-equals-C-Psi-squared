using System;
using System.Collections.Generic;
using System.Numerics;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>Reconstructs the AT-locked factor (the radical-closed, integrable half of the path-k
/// (SE,DE) spectrum) as an exact Z[i] polynomial, ×2-scaled to match the cleared block from
/// <see cref="F89PathKSeDeBlock"/> (roots 2·λ_AT). Dividing it out of the live full characteristic
/// polynomial isolates the H_B-mixed factor F_d. For path-3 the AT factor is the product of the two
/// explicit F_a/F_b quadratics (<see cref="ForPath3"/>). For path-k≥4 (<see cref="ForPathK"/>) the
/// naive "single-particle F_a + 2-particle DE-Slater F_b" multiset rule FAILS (the F_b modes are not
/// simple Slater sums from path-5 on); the AT factor is instead reconstructed from the rate-confined
/// invariant subspace: with M = D + iK, AT = ∏ over the rate sectors of charpoly(M|W), where W is the
/// largest M-invariant subspace inside that sector (found by the iterative shrink W ← {v : Kv ∈ span W}
/// to avoid the Kᵐ coefficient blow-up), mapped back to Z[i] by the iᵈ·p(−i(λ−r0)) substitution.</summary>
public static class F89AtFactorReconstruction
{
    /// <summary>The ×2-scaled path-3 AT factor AT = F_a·F_b at integer q0 (degree 4, monic, Z[i]).</summary>
    public static GaussianInteger[] ForPath3(int q0)
    {
        // F_a = λ² + (2iq+4)λ + (4q²+4iq+4);  F_b = λ² + (2iq+12)λ + (4q²+12iq+36).
        var fa = new GaussianInteger[]
        {
            new(4L * q0 * q0 + 4, 4L * q0),     // 4q²+4 + 4q·i
            new(4, 2L * q0),                    // 4 + 2q·i
            GaussianInteger.One,
        };
        var fb = new GaussianInteger[]
        {
            new(4L * q0 * q0 + 36, 12L * q0),   // 4q²+36 + 12q·i
            new(12, 2L * q0),                   // 12 + 2q·i
            GaussianInteger.One,
        };
        var at = GaussianPolynomial.Multiply(fa, fb);       // degree 4, roots λ_AT
        var scaled = new GaussianInteger[at.Length];        // roots 2λ_AT: coeff[k]·2^(deg−k)
        for (int k = 0; k < at.Length; k++)
        {
            BigInteger pow = BigInteger.Pow(2, at.Length - 1 - k);
            scaled[k] = new GaussianInteger(at[k].Re * pow, at[k].Im * pow);
        }
        return scaled;
    }

    /// <summary>FULL-D general path-k AT factor (×2-scaled, over Z[i]): the rate-confined
    /// invariant-subspace construction. M = D + iK (D = real diagonal rates ×2 ∈ {−4,−12}, K = the
    /// real-symmetric hopping = Im of the cleared block). For each rate sector, W = the largest
    /// M-invariant subspace inside it (= nullspace of [P_Uc·Kᵐ]); the sector factor is
    /// charpoly(M|W) = iᵈ·p(−i(λ−r0)) where p = charpoly of the rate-restriction K|W. No F_d import.
    /// Proven against the oracle (= C/AT) for path-4/5/6; reference: f89_pathk_galois.py full-d.</summary>
    public static GaussianInteger[] ForPathK(int k)
    {
        var (rate, ksym, n) = RateAndK(k);

        GaussianInteger[] at = { GaussianInteger.One };
        foreach (var r0 in new BigInteger[] { -4, -12 })        // the ×2-cleared sector rates
        {
            var W = SectorBasis(rate, r0, n);
            if (W.GetLength(1) == 0) continue;
            W = LargestInvariantSubspace(ksym, W);              // shrink to the K-invariant part
            int dim = W.GetLength(1);
            if (dim == 0) continue;

            var R = Restrict(ksym, W);                          // K | W  (dim × dim, BigRational)
            var p = BigRationalMatrixCharpoly.Characteristic(R);
            at = GaussianPolynomial.Multiply(at, ISubstitute(p, r0, dim));
        }
        return at;
    }

    /// <summary>The AT-locked invariant subspace ITSELF (the union over rate sectors of the largest
    /// K-invariant subspace), as a real basis of column vectors (SymDim × AtDegree, not orthonormal).
    /// This is the SAME subspace whose characteristic polynomial <see cref="ForPathK"/> returns as the AT
    /// factor (validated against the oracle by <c>F89FullDReconstructionTests</c>); exposing the basis lets
    /// the diabolic scout read the F_d residual roots at ANY q EXACTLY, as the eigenvalues of M(q)
    /// compressed onto the orthogonal complement of this subspace (the block-triangular split: the AT
    /// subspace is M(q)-invariant for every q because the rate diagonal D and the hopping K each preserve
    /// it and K-invariance is scale-invariant in q, so the compression's eigenvalues are exactly the full
    /// spectrum minus the AT part). AT-free by construction, with no nearest-match strand partition (which
    /// floods at the F_53 strand density of path-6/N=7) and no continuity tracking.</summary>
    public static double[,] AtInvariantSubspaceBasis(int k)
    {
        var (rate, ksym, n) = RateAndK(k);
        var cols = new List<BigRational[]>();
        foreach (var r0 in new BigInteger[] { -4, -12 })
        {
            var W = SectorBasis(rate, r0, n);
            if (W.GetLength(1) == 0) continue;
            W = LargestInvariantSubspace(ksym, W);
            for (int c = 0; c < W.GetLength(1); c++)
            {
                var col = new BigRational[n];
                for (int i = 0; i < n; i++) col[i] = W[i, c];
                cols.Add(col);
            }
        }
        var basis = new double[n, cols.Count];
        for (int c = 0; c < cols.Count; c++)
            for (int i = 0; i < n; i++)
                basis[i, c] = (double)cols[c][i].Numerator / (double)cols[c][i].Denominator;
        return basis;
    }

    /// <summary>The R-ODD analog of <see cref="AtInvariantSubspaceBasis"/>: the largest rate-confined
    /// K-invariant subspace of the R-odd sector of the FULL (SE,DE) block (physical scale, rates ∈ {−2,−6};
    /// the R-even story lives in the ×2-cleared S₂-symmetric block, this one in its reflection-antisymmetric
    /// complement). Returned as a real basis of column vectors in the R-ODD 2-CYCLE COORDINATES of
    /// <see cref="F89PathKSeDeBlock.ROddBasis"/> (column c there = coordinate c here; the shared convention),
    /// oddDim × atOddDegree, not orthonormal; atOddDegree may be 0 (then the whole R-odd sector is residual).
    /// Same exact construction as the R-even case (<see cref="AtInvariantSubspaceBasis"/>): per rate sector
    /// the iterative shrink W ← {v ∈ W : Kv ∈ span W}, all in BigRational, so the compression onto the
    /// orthogonal complement reads the R-odd RESIDUAL roots at any q exactly (the sectorbraid arc's R-odd
    /// deep-loci probe; the subspace is L(q)-invariant for every q because D preserves each rate sector and
    /// K-invariance is q-scale-invariant).</summary>
    public static double[,] ROddAtInvariantSubspaceBasis(int k)
    {
        var (rate, ksym, n) = ROddRateAndK(k);
        var cols = new List<BigRational[]>();
        foreach (var r0 in new BigInteger[] { -2, -6 })         // the full-block (physical) sector rates
        {
            var W = SectorBasis(rate, r0, n);
            if (W.GetLength(1) == 0) continue;
            W = LargestInvariantSubspace(ksym, W);
            for (int c = 0; c < W.GetLength(1); c++)
            {
                var col = new BigRational[n];
                for (int i = 0; i < n; i++) col[i] = W[i, c];
                cols.Add(col);
            }
        }
        var basis = new double[n, cols.Count];
        for (int c = 0; c < cols.Count; c++)
            for (int i = 0; i < n; i++)
                basis[i, c] = (double)cols[c][i].Numerator / (double)cols[c][i].Denominator;
        return basis;
    }

    /// <summary>The exact per-sector data of the AT factor of the ×2-cleared (SE,DE) sector block, in the
    /// block's OWN pencil units 2M(q) = A + q·C with C = i·K: for each cleared rate sector r0 ∈ {−4, −12}
    /// (the diagonal of A), the INTEGER characteristic polynomial of K|W on the largest K-invariant subspace
    /// W inside the sector. Because W is q-independent (K-invariance is q-scale-invariant) and
    /// 2M(q)|W = r0·I + i·q·(K|W), the AT strands are exactly q-linear and the bivariate AT factor is
    /// ∏_sectors Σ_m p_m·i^{d−m}·(Λ−r0)^m·q^{d−m} ∈ Z[i][Λ, q]. Integrality of p = charpoly(K|W) is
    /// guaranteed by Gauss's lemma (the AT factor is a monic-in-Λ factor of the block's monic integer
    /// charpoly over the integrally closed Z[i][q]); a guard throws if it fails numerically. rOdd selects
    /// the reflection-antisymmetric block (<see cref="F89PathKSeDeBlock.BuildTwoTimesROddBlock"/>) instead
    /// of the S₂-symmetric one (<see cref="F89PathKSeDeBlock.BuildTwoTimesSymBlock"/>); at q0=2 the sector
    /// product reproduces <see cref="ForPathK"/> exactly (pinned by the FOLDRESULTANT gate). Consumed
    /// per-sample mod p by the fold-resultant certificate (sectorbraid arc, remainder R1).</summary>
    public static IReadOnlyList<AtSector> ClearedAtSectors(int k, bool rOdd)
    {
        int nBlock = k + 1;
        var b0 = rOdd
            ? F89PathKSeDeBlock.BuildTwoTimesROddBlock(0, nBlock)
            : F89PathKSeDeBlock.BuildTwoTimesSymBlock(0, nBlock);
        var b1 = rOdd
            ? F89PathKSeDeBlock.BuildTwoTimesROddBlock(1, nBlock)
            : F89PathKSeDeBlock.BuildTwoTimesSymBlock(1, nBlock);
        int n = b0.GetLength(0);

        var rate = new BigInteger[n];
        var ksym = new BigRational[n, n];                       // K = −i·C = Im(C), real
        for (int i = 0; i < n; i++)
        {
            if (!b0[i, i].Im.IsZero)
                throw new InvalidOperationException("the cleared block's q=0 diagonal must be real.");
            rate[i] = b0[i, i].Re;
            for (int j = 0; j < n; j++)
            {
                if (i != j && !b0[i, j].Equals(GaussianInteger.Zero))
                    throw new InvalidOperationException("the cleared block at q=0 must be diagonal.");
                var c = b1[i, j] - b0[i, j];                    // C = blk(1) − blk(0), pure imaginary
                if (!c.Re.IsZero)
                    throw new InvalidOperationException("the hopping direction C must be pure imaginary.");
                ksym[i, j] = new BigRational(c.Im);
            }
        }

        var sectors = new List<AtSector>();
        foreach (var r0 in new BigInteger[] { -4, -12 })        // the ×2-cleared sector rates
        {
            var w = SectorBasis(rate, r0, n);
            if (w.GetLength(1) == 0) continue;
            w = LargestInvariantSubspace(ksym, w);
            int dim = w.GetLength(1);
            if (dim == 0) continue;

            var p = BigRationalMatrixCharpoly.Characteristic(Restrict(ksym, w));
            var coeffs = new BigInteger[p.Length];
            for (int m = 0; m < p.Length; m++)
            {
                if (!p[m].IsInteger)
                    throw new InvalidOperationException(
                        $"charpoly(K|W) coefficient {m} of sector r0={r0} is not an integer (Gauss's lemma violated?).");
                coeffs[m] = p[m].Numerator;
            }
            sectors.Add(new AtSector(r0, coeffs));
        }
        return sectors;
    }

    /// <summary>The FULL (SE,DE) block at q0=2 projected onto the R-odd 2-cycle basis (columns e_t − e_perm[t],
    /// increasing t, the un-normalized integer form of <see cref="F89PathKSeDeBlock.ROddBasis"/>): the rate
    /// diagonal (each 2-cycle carries one rate ∈ {−2,−6}, R commutes with the dephasing) and the projected
    /// real-symmetric hopping K_odd = (WᵀW)⁻¹Wᵀ·Im(L)·W = ½Wᵀ·Im(L)·W (WᵀW = 2I), which equals UᵀIm(L)U for
    /// the orthonormal U = W/√2, so the coordinates match <see cref="F89PathKSeDeBlock.ROddBasis"/> exactly.
    /// All entries are exact integers at q0=2; a guard throws if a non-integer sneaks in.</summary>
    private static (BigInteger[] rate, BigRational[,] ksym, int n) ROddRateAndK(int k)
    {
        int nBlock = k + 1;
        var L = F89PathKSeDeBlock.BuildFullBlock(nBlock, new Complex(2, 0));
        var perm = F89PathKSeDeBlock.ReflectionPermutation(nBlock);
        var pairs = new List<(int T, int T2)>();
        for (int t = 0; t < perm.Length; t++)
            if (perm[t] > t) pairs.Add((t, perm[t]));

        int n = pairs.Count;
        var rate = new BigInteger[n];
        var ksym = new BigRational[n, n];
        for (int r = 0; r < n; r++)
        {
            rate[r] = ExactInteger(L[pairs[r].T, pairs[r].T].Real);
            for (int c = 0; c < n; c++)
            {
                double v = (L[pairs[r].T, pairs[c].T].Imaginary - L[pairs[r].T, pairs[c].T2].Imaginary
                          - L[pairs[r].T2, pairs[c].T].Imaginary + L[pairs[r].T2, pairs[c].T2].Imaginary) / 2.0;
                ksym[r, c] = new BigRational(ExactInteger(v));
            }
        }
        return (rate, ksym, n);
    }

    /// <summary>Checked double → BigInteger for the exactly-integer entries of the q0=2 full block.</summary>
    private static BigInteger ExactInteger(double v)
    {
        double r = Math.Round(v);
        if (Math.Abs(v - r) > 1e-9)
            throw new InvalidOperationException($"expected an exact integer entry, got {v}");
        return new BigInteger(r);
    }

    /// <summary>The ×2-cleared (SE,DE) block at q0=2 reduced to its integer rate diagonal D (the AT rates
    /// ∈ {−4,−12}) and its real-symmetric hopping K = Im(2M). Shared by <see cref="ForPathK"/> and
    /// <see cref="AtInvariantSubspaceBasis"/> so both read the IDENTICAL invariant subspace.</summary>
    private static (BigInteger[] rate, BigRational[,] ksym, int n) RateAndK(int k)
    {
        var block = F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0: 2, nBlock: k + 1);
        int n = block.GetLength(0);
        var rate = new BigInteger[n];
        var ksym = new BigRational[n, n];                       // K = Im(2M), real symmetric
        for (int i = 0; i < n; i++)
        {
            rate[i] = block[i, i].Re;
            for (int j = 0; j < n; j++) ksym[i, j] = new BigRational(block[i, j].Im);
        }
        return (rate, ksym, n);
    }

    /// <summary>The standard basis of the rate sector { i : rate[i] = r0 } (n × |sector|).</summary>
    private static BigRational[,] SectorBasis(BigInteger[] rate, BigInteger r0, int n)
    {
        var idx = new List<int>();
        for (int i = 0; i < n; i++) if (rate[i] == r0) idx.Add(i);
        var w = new BigRational[n, idx.Count];
        for (int i = 0; i < n; i++) for (int c = 0; c < idx.Count; c++) w[i, c] = BigRational.Zero;
        for (int c = 0; c < idx.Count; c++) w[idx[c], c] = BigRational.One;
        return w;
    }

    /// <summary>The largest K-invariant subspace contained in span(W), by iteratively keeping only
    /// the v ∈ W with Kv ∈ span(W) until stable. Uses K once per step (no Kᵐ blow-up).</summary>
    private static BigRational[,] LargestInvariantSubspace(BigRational[,] k, BigRational[,] w)
    {
        while (true)
        {
            int d = w.GetLength(1);
            var kw = BigRationalLinearAlgebra.Multiply(k, w);                 // K·W
            var wt = BigRationalLinearAlgebra.Transpose(w);
            var coords = BigRationalLinearAlgebra.Solve(
                BigRationalLinearAlgebra.Multiply(wt, w),
                BigRationalLinearAlgebra.Multiply(wt, kw));                   // (WᵀW)⁻¹ Wᵀ K W
            var resid = SubtractMatrices(kw, BigRationalLinearAlgebra.Multiply(w, coords)); // (I−P)K W
            var c = BigRationalLinearAlgebra.Nullspace(resid);               // { c : K(Wc) ∈ span W }
            if (c.Count == d) return w;                                      // stable
            if (c.Count == 0) return new BigRational[w.GetLength(0), 0];
            var cMat = new BigRational[d, c.Count];
            for (int j = 0; j < c.Count; j++) for (int i = 0; i < d; i++) cMat[i, j] = c[j][i];
            w = BigRationalLinearAlgebra.Multiply(w, cMat);                  // W ← W·C
        }
    }

    /// <summary>K | W = (WᵀW)⁻¹ Wᵀ K W (W must be K-invariant).</summary>
    private static BigRational[,] Restrict(BigRational[,] k, BigRational[,] w)
    {
        var wt = BigRationalLinearAlgebra.Transpose(w);
        return BigRationalLinearAlgebra.Solve(
            BigRationalLinearAlgebra.Multiply(wt, w),
            BigRationalLinearAlgebra.Multiply(wt, BigRationalLinearAlgebra.Multiply(k, w)));
    }

    private static BigRational[,] SubtractMatrices(BigRational[,] a, BigRational[,] b)
    {
        int n = a.GetLength(0), m = a.GetLength(1);
        var r = new BigRational[n, m];
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) r[i, j] = a[i, j] - b[i, j];
        return r;
    }

    /// <summary>AT_sector(λ) = det(λI − M|W) = iᵈ·p(−i(λ−r0)) over Z[i], where p = charpoly(K|W)
    /// (real) and M|W = r0·I + i·(K|W). Returns the sector factor as Gaussian-integer coefficients
    /// (lowest-first); throws if any coefficient is not a Gaussian integer.</summary>
    private static GaussianInteger[] ISubstitute(BigRational[] p, BigInteger r0, int d)
    {
        var re = new BigRational[d + 1];
        var im = new BigRational[d + 1];
        for (int l = 0; l <= d; l++) { re[l] = BigRational.Zero; im[l] = BigRational.Zero; }
        BigInteger negR0 = -r0;
        for (int m = 0; m <= d; m++)
        {
            if (p[m].IsZero) continue;
            int e = (d + m) % 4;                                   // i^(d+m): unit (uRe, uIm)
            int uRe = e == 0 ? 1 : e == 2 ? -1 : 0;
            int uIm = e == 1 ? 1 : e == 3 ? -1 : 0;
            int sign = (m % 2 == 0) ? 1 : -1;
            for (int l = 0; l <= m; l++)
            {
                BigInteger scale = Binom(m, l) * BigInteger.Pow(negR0, m - l);   // C(m,l)(−r0)^{m−l}
                BigRational term = p[m] * new BigRational(scale) * sign;
                re[l] += term * uRe;
                im[l] += term * uIm;
            }
        }
        var result = new GaussianInteger[d + 1];
        for (int l = 0; l <= d; l++)
        {
            if (!re[l].IsInteger || !im[l].IsInteger)
                throw new InvalidOperationException($"AT sector coefficient {l} is not a Gaussian integer.");
            result[l] = new GaussianInteger(re[l].Numerator, im[l].Numerator);
        }
        return result;
    }

    private static BigInteger Binom(int m, int l)
    {
        BigInteger r = BigInteger.One;
        for (int i = 0; i < l; i++) r = r * (m - i) / (i + 1);
        return r;
    }
}

/// <summary>One rate sector of the ×2-cleared (SE,DE) block's AT factor: the cleared rate r0 (∈ {−4, −12})
/// and the monic INTEGER characteristic polynomial (lowest-first) of K|W on the sector's largest K-invariant
/// subspace W. The sector's exact bivariate AT contribution is Σ_m KCharpoly[m]·i^{d−m}·(Λ−r0)^m·q^{d−m}
/// (d = degree); its strands are the q-linear Λ = r0 + i·q·κ_j, κ_j the eigenvalues of K|W. Produced by
/// <see cref="F89AtFactorReconstruction.ClearedAtSectors"/>, consumed by the fold-resultant certificate.</summary>
public sealed record AtSector(BigInteger Rate, BigInteger[] KCharpoly);
