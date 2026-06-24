using System;
using System.Collections.Generic;
using System.Numerics;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>Reconstructs the AT-locked factor (the radical-closed, integrable half of the path-k
/// (SE,DE) spectrum) as an exact Z[i] polynomial, ×2-scaled to match the cleared block from
/// <see cref="F89PathKSeDeBlock"/> (roots 2·λ_AT). Dividing it out of the live full characteristic
/// polynomial isolates the H_B-mixed factor F_d. For path-3 the AT factor is the product of the two
/// explicit F_a/F_b quadratics; the general path-k≥4 reconstruction (single-particle F_a +
/// 2-particle DE-Slater F_b multiset) is built separately.</summary>
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
        var block = F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0: 2, nBlock: k + 1);
        int n = block.GetLength(0);
        var rate = new BigInteger[n];
        var ksym = new BigRational[n, n];                       // K = Im(2M), real symmetric
        for (int i = 0; i < n; i++)
        {
            rate[i] = block[i, i].Re;
            for (int j = 0; j < n; j++) ksym[i, j] = new BigRational(block[i, j].Im);
        }

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
