using System;
using System.Collections.Generic;
using System.Numerics;

namespace RCPsiSquared.Core.Numerics;

/// <summary>The result of an eigenvector-frame holonomy loop: the 2×2 frame-monodromy after each loop
/// (row-major [M00,M01,M10,M11]) and the normalized 2D-span residual per loop (how well the transported
/// pair stays in the initial coalescing subspace — small on every loop for a clean isolated EP).</summary>
public sealed class HolonomyResult
{
    /// <summary>LoopMonodromy[n−1] = the 2×2 frame-monodromy after n loops, row-major [M00,M01,M10,M11].</summary>
    public Complex[][] LoopMonodromy { get; }
    /// <summary>SpanResidual[n−1] = ‖V0·Mₙ − transported‖_F / ‖transported‖_F after n loops.</summary>
    public double[] SpanResidual { get; }

    public HolonomyResult(Complex[][] loopMonodromy, double[] spanResidual)
    {
        LoopMonodromy = loopMonodromy;
        SpanResidual = spanResidual;
    }

    /// <summary>The two eigenvalues of the one-loop frame-monodromy M₁ (expect ±i for a defective EP2).</summary>
    public Complex[] M1Eigenvalues => Eig2x2(LoopMonodromy[0]);

    /// <summary>Eigenvalues of a row-major 2×2 [a,b,c,d].</summary>
    public static Complex[] Eig2x2(Complex[] m)
    {
        Complex a = m[0], b = m[1], c = m[2], d = m[3];
        Complex tr = a + d, det = a * d - b * c;
        Complex disc = Complex.Sqrt(tr * tr - 4.0 * det);
        return new[] { (tr + disc) / 2.0, (tr - disc) / 2.0 };
    }
}

/// <summary>EigenVECTOR-frame holonomy around a defective exceptional point — the companion of
/// <see cref="Monodromy"/> (which tracks the eigenVALUE swap). Encircling a defective EP2 in the complex
/// parameter plane, the two coalescing right eigenvectors, transported in the biorthogonal vᵀv gauge
/// (the natural inner product of a complex-symmetric matrix, self-orthogonal vᵀv → 0 at the EP), give a
/// frame-monodromy generator with eigenvalues ±i: M₁ ~ a 90° rotation, M₂ = −I, M₄ = I, so the frame is
/// single-valued only after four loops (i⁴ = 1).
///
/// <para>Gauge note (calculated-or-marked): the gauge-INVARIANT content is the swap — M is anti-diagonal,
/// the two eigenvectors exchanged (the eigenVALUE monodromy is order 2). The specific ±i / mod-4 is the
/// biorthogonal-gauge READING: on the real block a cleanly phase-pinned Hermitian vᴴv gauge gives only
/// mod-2 (the ±1 swap, M₂ = +I) — note a bare 2×2 toy is NOT discriminating (LAPACK's phase convention
/// leaks ±i into the Hermitian gauge there too), while the vᵀv self-orthogonality (vᵀv ~ √ε near the EP, so √(vᵀv) ~ ε^{1/4} winds a
/// quarter-turn per loop) is exactly what upgrades mod-2 to mod-4. So vᵀv is the load-bearing gauge — not
/// because the Hermitian gauge "sees nothing" (it sees the swap), but because only vᵀv promotes the swap
/// to the ±i rotation. (The static |cos| merge of <see cref="EpCharacter"/> is a different, non-holonomy
/// object; do not conflate.)</para>
///
/// <para>Numerics: right eigenvectors via LAPACK zgeev (<see cref="MklDirect"/>.EigenvaluesAndVectorsDirectRaw);
/// the coalescing pair is followed by a 2-strand biorthogonal-overlap bijection (distinct strands), ±1-fixed
/// by Euclidean closeness (essential — it sets the step-count floor, per-step phase < π/2).</para></summary>
public static class EigenvectorHolonomy
{
    /// <summary>Track the frame of the two eigenvectors nearest <paramref name="lambda0"/> around the loop
    /// q(θ) = <paramref name="center"/> + <paramref name="radius"/>·e^{iθ}, for the complex-symmetric matrix
    /// family <paramref name="matrixAt"/> (dim×dim). Returns the frame-monodromy after 1..nLoops loops.</summary>
    public static HolonomyResult FrameMonodromy(
        Func<Complex, Complex[,]> matrixAt, int dim, Complex lambda0,
        Complex center, double radius, int nLoops = 4, int stepsPerLoop = 400)
    {
        if (dim < 2) throw new ArgumentOutOfRangeException(nameof(dim), "dim must be ≥ 2");
        if (radius <= 0) throw new ArgumentOutOfRangeException(nameof(radius), "radius must be positive");

        // reference frame at θ=0, biorthogonally normalized
        var (w0, v0all) = Eig(matrixAt(center + radius), dim);
        int[] idx0 = NearestTwo(w0, lambda0, dim);
        var prev = new Complex[dim, 2];
        for (int t = 0; t < 2; t++)
        {
            var v = BNorm(Col(v0all, idx0[t], dim));
            for (int i = 0; i < dim; i++) prev[i, t] = v[i];
        }
        var V0 = (Complex[,])prev.Clone();

        var loopM = new List<Complex[]>();
        var loopR = new List<double>();
        int total = nLoops * stepsPerLoop;
        for (int n = 1; n <= total; n++)
        {
            double theta = 2.0 * Math.PI * n / stepsPerLoop;
            var q = center + radius * new Complex(Math.Cos(theta), Math.Sin(theta));
            var (_, vall) = Eig(matrixAt(q), dim);
            var next = new Complex[dim, 2];
            int chosen0 = -1;                                   // 2-strand bijection: the two strands must
            for (int t = 0; t < 2; t++)                         // pick DISTINCT eigenvectors (else a degenerate
            {                                                   // frame, det → 0, on a poorly-chosen lambda0)
                // best-matching new eigenvector by |biorthogonal overlap| to prev column t
                int best = -1; double bestMag = -1.0;
                for (int j = 0; j < dim; j++)
                {
                    if (t == 1 && j == chosen0) continue;       // strand 1 may not reuse strand 0's pick
                    double mag = BDot(prev, t, vall, j, dim).Magnitude;
                    if (mag > bestMag) { bestMag = mag; best = j; }
                }
                if (t == 0) chosen0 = best;
                var vj = BNorm(Col(vall, best, dim));
                // Load-bearing, not cosmetic: BNorm's principal Complex.Sqrt(vᵀv) flips sign as vᵀv crosses
                // its branch cut once per loop; fixing the ±1 by Euclidean (Hermitian) closeness to the
                // previous step restores the continuous quarter-turn. The mod-2 → mod-4 promotion itself
                // comes from the vᵀv normalizer (√(vᵀv) ~ ε^{1/4}, see the class docstring); this sign-fix
                // only preserves step-to-step continuity so the numerics reveal that branch instead of
                // hiding it behind the cut. It corrects a ±1 only, so it sets the step-count requirement:
                // the per-step phase must stay < π/2 (fails only very close to the EP).
                Complex herm = Complex.Zero;
                for (int i = 0; i < dim; i++) herm += Complex.Conjugate(prev[i, t]) * vj[i];
                if (herm.Real < 0) for (int i = 0; i < dim; i++) vj[i] = -vj[i];
                for (int i = 0; i < dim; i++) next[i, t] = vj[i];
            }
            prev = next;
            if (n % stepsPerLoop == 0)
            {
                var (m, resid) = FitFrame(V0, prev, dim);
                loopM.Add(m);
                loopR.Add(resid);
            }
        }
        return new HolonomyResult(loopM.ToArray(), loopR.ToArray());
    }

    static (Complex[] w, Complex[] vFlat) Eig(Complex[,] mat, int dim)
    {
        var a = new Complex[(long)dim * dim];
        for (int j = 0; j < dim; j++)
            for (int i = 0; i < dim; i++)
                a[i + (long)j * dim] = mat[i, j];              // column-major (zgeev), destroyed in place
        var (w, v) = MklDirect.EigenvaluesAndVectorsDirectRaw(a, dim);
        return (w, v);                                          // v column-major: v[j*dim + i]
    }

    static int[] NearestTwo(Complex[] w, Complex lambda0, int dim)
    {
        int a = 0, b = 1; double da = double.MaxValue, db = double.MaxValue;
        for (int i = 0; i < dim; i++)
        {
            double d = (w[i] - lambda0).Magnitude;
            if (d < da) { db = da; b = a; da = d; a = i; }
            else if (d < db) { db = d; b = i; }
        }
        return new[] { a, b };
    }

    static Complex[] Col(Complex[] vFlat, int j, int dim)
    {
        var v = new Complex[dim];
        for (int i = 0; i < dim; i++) v[i] = vFlat[(long)j * dim + i];
        return v;
    }

    static Complex[] BNorm(Complex[] v)
    {
        Complex n2 = Complex.Zero;
        foreach (var z in v) n2 += z * z;                      // vᵀv (no conjugate)
        var s = Complex.Sqrt(n2);
        var r = new Complex[v.Length];
        for (int i = 0; i < v.Length; i++) r[i] = v[i] / s;
        return r;
    }

    static Complex BDot(Complex[,] prev, int t, Complex[] vFlat, int j, int dim)
    {
        Complex s = Complex.Zero;
        for (int i = 0; i < dim; i++) s += prev[i, t] * vFlat[(long)j * dim + i];   // prevᵀ · V_j
        return s;
    }

    // Least-squares 2×2 M with V0·M ≈ prev (Hermitian pseudo-inverse M = (V0ᴴV0)⁻¹ V0ᴴ prev),
    // plus the normalized Frobenius residual of the fit.
    static (Complex[] m, double resid) FitFrame(Complex[,] V0, Complex[,] prev, int dim)
    {
        Complex g00 = 0, g01 = 0, g10 = 0, g11 = 0, r00 = 0, r01 = 0, r10 = 0, r11 = 0;
        for (int i = 0; i < dim; i++)
        {
            Complex a0 = Complex.Conjugate(V0[i, 0]), a1 = Complex.Conjugate(V0[i, 1]);
            g00 += a0 * V0[i, 0]; g01 += a0 * V0[i, 1];
            g10 += a1 * V0[i, 0]; g11 += a1 * V0[i, 1];
            r00 += a0 * prev[i, 0]; r01 += a0 * prev[i, 1];
            r10 += a1 * prev[i, 0]; r11 += a1 * prev[i, 1];
        }
        Complex detG = g00 * g11 - g01 * g10;
        Complex i00 = g11 / detG, i01 = -g01 / detG, i10 = -g10 / detG, i11 = g00 / detG;
        Complex m00 = i00 * r00 + i01 * r10;
        Complex m01 = i00 * r01 + i01 * r11;
        Complex m10 = i10 * r00 + i11 * r10;
        Complex m11 = i10 * r01 + i11 * r11;
        double num = 0, den = 0;
        for (int i = 0; i < dim; i++)
        {
            Complex p0 = V0[i, 0] * m00 + V0[i, 1] * m10;
            Complex p1 = V0[i, 0] * m01 + V0[i, 1] * m11;
            num += (p0 - prev[i, 0]).Magnitude * (p0 - prev[i, 0]).Magnitude
                 + (p1 - prev[i, 1]).Magnitude * (p1 - prev[i, 1]).Magnitude;
            den += prev[i, 0].Magnitude * prev[i, 0].Magnitude + prev[i, 1].Magnitude * prev[i, 1].Magnitude;
        }
        return (new[] { m00, m01, m10, m11 }, Math.Sqrt(num / den));
    }
}
