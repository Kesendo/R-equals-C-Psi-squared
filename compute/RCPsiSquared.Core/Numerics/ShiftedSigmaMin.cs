using System;
using System.Numerics;

namespace RCPsiSquared.Core.Numerics;

/// <summary>Smallest singular value of (M − s·I) from ONE dense LU factorization (zgetrf) plus inverse
/// power iteration on [(M−s)ᴴ(M−s)]⁻¹ = (M−s)⁻¹(M−s)⁻ᴴ (two zgetrs solves per step, TRANS='C' then 'N').
/// No full spectrum, no SVD: cost = one O(d³) LU + O(d²) per iteration (~10-30 iterations). The step-3
/// shell-census probe (docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md §6-§7).
///
/// <para>Direction of rigor (the census's honesty note): σ_min(M−s) ≤ |λ−s| for EVERY eigenvalue λ of M
/// (‖(M−s)v‖ = |λ−s|·‖v‖ ≥ σ_min·‖v‖), so a LARGE σ_min excludes eigenvalues near s. The norm-growth
/// estimate 1/√‖w_k‖ is an upper bound on σ_min at EVERY step and decreases monotonically (Cauchy-Schwarz
/// on the moments of the PSD iteration matrix), so the returned value is numerical census evidence
/// converging FROM ABOVE, not a certified lower bound. A small σ_min does NOT by itself certify an
/// eigenvalue near s (non-normal pseudospectrum); membership verdicts lean on the containment corollary.
/// At a member shift the coalescing near-defective PAIR makes the two smallest singular values comparable;
/// the value estimate still lands between them (both ≈ 0), which is all the census reads.</para></summary>
public static class ShiftedSigmaMin
{
    public readonly record struct Result(double SigmaMin, int Iterations, bool Converged, bool ExactlySingular);

    /// <summary>Estimate σ_min(m − shift·I) for a dense square block (small/test path; copies into
    /// column-major flat storage first).</summary>
    public static Result Estimate(Complex[,] m, Complex shift, int maxIter = 200, double relTol = 1e-4, int seed = 12345)
    {
        int d = m.GetLength(0);
        if (m.GetLength(1) != d) throw new ArgumentException("square matrix required", nameof(m));
        var a = new Complex[(long)d * d];
        for (int col = 0; col < d; col++)
            for (int row = 0; row < d; row++)
                a[(long)col * d + row] = m[row, col] - (row == col ? shift : Complex.Zero);
        return EstimateColumnMajor(a, d, maxIter, relTol, seed);
    }

    /// <summary>Estimate σ_min for a column-major flat matrix that ALREADY has the shift subtracted.
    /// 'a' is destroyed (LU in place). The big-block entry point: no copy, so the caller's sector matrix
    /// is the only O(d²) allocation (the LP64 wall d ≤ 46340 is the caller's responsibility).</summary>
    public static Result EstimateColumnMajor(Complex[] a, int d, int maxIter = 200, double relTol = 1e-4, int seed = 12345)
    {
        var ipiv = new int[d];
        try { MklDirect.LuFactorizeRaw(a, d, ipiv); }
        catch (InvalidOperationException) { return new Result(0.0, 0, true, true); }   // U exactly singular ⟹ σ_min = 0

        var rng = new Random(seed);
        var v = new Complex[d];
        for (int i = 0; i < d; i++) v[i] = new Complex(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
        Normalize(v);

        double sigma = double.PositiveInfinity;
        int stallCount = 0;
        for (int it = 1; it <= maxIter; it++)
        {
            // v ← (M−s)⁻¹(M−s)⁻ᴴ v, in place: growth ‖v‖ → 1/σ_min²
            MklDirect.LuSolveRaw(a, d, ipiv, v, conjugateTranspose: true);
            MklDirect.LuSolveRaw(a, d, ipiv, v, conjugateTranspose: false);
            double growth = Norm(v);
            if (growth == 0.0) return new Result(0.0, it, true, false);
            double next = 1.0 / Math.Sqrt(growth);
            for (int i = 0; i < d; i++) v[i] /= growth;
            // require the relative-change criterion TWICE in a row before declaring convergence (guards
            // the near-degenerate member-pair stall); the estimate is an upper bound at every step, so a
            // residual stall is conservative (never turns a true non-member into a false member).
            bool small = Math.Abs(next - sigma) <= relTol * next;
            stallCount = small ? stallCount + 1 : 0;
            sigma = next;
            if (stallCount >= 2 || sigma < 1e-14) return new Result(sigma, it, true, false);
        }
        return new Result(sigma, maxIter, false, false);
    }

    static double Norm(Complex[] v)
    {
        double s = 0;
        foreach (var z in v) s += z.Real * z.Real + z.Imaginary * z.Imaginary;
        return Math.Sqrt(s);
    }

    static void Normalize(Complex[] v)
    {
        double n = Norm(v);
        for (int i = 0; i < v.Length; i++) v[i] /= n;
    }
}
