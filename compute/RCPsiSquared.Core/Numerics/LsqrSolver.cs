using System;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;

namespace RCPsiSquared.Core.Numerics;

/// <summary>
/// LSQR (Paige-Saunders 1982) solve of the shifted system (M − shift·I)·x = rhs on a
/// <see cref="WeightCoherenceSectorCsr.Csr"/> weight-coherence sector. The [N7] fallback inner engine
/// for <see cref="SparseShiftedSigmaMin"/>: Golub-Kahan bidiagonalization works implicitly on the
/// Hermitian PSD operator (M−s)ᴴ(M−s), so convergence is governed by the SINGULAR value distribution
/// only; the non-normal eigenvalue clustering that makes identity-preconditioned BiCGStab diverge at
/// near-spectrum shifts (measured at the N=9 (3,3) × λ_A cell: BiCGStab residual 6.6 after 100k
/// iterations, while the same operator at a generic shift solves in 10) is structurally invisible to
/// LSQR. The residual norm is monotonically nonincreasing, so LSQR has no BiCGStab-style breakdown:
/// budget exhaustion is a genuine conditioning verdict, not an unlucky recurrence.
///
/// <para>Complex Golub-Kahan recurrences: with A = M − s·I,
///   β₁u₁ = b,   α₁v₁ = Aᴴu₁,
///   β_{k+1}u_{k+1} = Av_k − α_k u_k,   α_{k+1}v_{k+1} = Aᴴu_{k+1} − β_{k+1}v_k.
/// The scalars α, β are REAL (they are norms; each step normalizes u and v), so the plane rotations
/// eliminating the lower bidiagonal are the standard real LSQR rotations; only the vectors are complex,
/// and the conjugate transpose Aᴴ appears exactly where the real algorithm has Aᵀ.
/// The adjoint apply uses the same identity <see cref="SparseShiftedSigmaMin"/> documents:
/// (M − s·I)ᴴ = Mᴴ − conj(s)·I, so Aᴴu = CsrOps.MultiplyShifted(mh, conj(shift), u, ·) with
/// mh = CsrOps.HermitianTranspose(m).</para>
///
/// <para>Same honesty contract as <see cref="BiCgStabSolver"/> (which remains a landed, tested kernel):
/// at exit the TRUE relative residual is recomputed with one extra matvec, and
/// <see cref="Outcome.Converged"/> is that true residual being ≤ relTol.</para>
/// </summary>
public static class LsqrSolver
{
    /// <summary>Outcome of a <see cref="Solve(WeightCoherenceSectorCsr.Csr, WeightCoherenceSectorCsr.Csr, Complex, Complex[], Complex[], double, int)"/>
    /// call: whether the recomputed TRUE relative residual met the tolerance, the number of LSQR
    /// iterations taken, and that true relative residual.</summary>
    public readonly record struct Outcome(bool Converged, int Iterations, double RelResidual);

    /// <summary>Convenience overload computing the Hermitian transpose internally (one O(nnz) pass).
    /// Callers that already hold mh (the estimator) use the primary overload to avoid recomputing it.</summary>
    public static Outcome Solve(WeightCoherenceSectorCsr.Csr m, Complex shift, Complex[] rhs, Complex[] x,
        double relTol, int maxIter)
        => Solve(m, CsrOps.HermitianTranspose(m), shift, rhs, x, relTol, maxIter);

    /// <summary>LSQR for (M − shift·I)·x = rhs. mh MUST be CsrOps.HermitianTranspose(m); it is used for
    /// the adjoint apply Aᴴu = (Mᴴ − conj(shift)·I)u. To solve the HERMITIAN system (M − s·I)ᴴ w = z,
    /// call Solve(mh, m, Complex.Conjugate(s), z, w, ...): the forward operator Mᴴ − conj(s)·I then has
    /// adjoint (Mᴴ − conj(s)·I)ᴴ = M − s·I, which is exactly the swapped pair's adjoint apply. Initial
    /// guess x ≡ 0 (x is cleared on entry, no warm start, matching the BiCGStab kernel).</summary>
    public static Outcome Solve(WeightCoherenceSectorCsr.Csr m, WeightCoherenceSectorCsr.Csr mh, Complex shift,
        Complex[] rhs, Complex[] x, double relTol, int maxIter)
    {
        if (m is null) throw new ArgumentNullException(nameof(m));
        if (mh is null) throw new ArgumentNullException(nameof(mh));
        if (rhs is null) throw new ArgumentNullException(nameof(rhs));
        if (x is null) throw new ArgumentNullException(nameof(x));
        if (rhs.Length != m.Dim) throw new ArgumentException($"rhs length {rhs.Length} != dim {m.Dim}", nameof(rhs));
        if (x.Length != m.Dim) throw new ArgumentException($"x length {x.Length} != dim {m.Dim}", nameof(x));
        if (mh.Dim != m.Dim) throw new ArgumentException($"mh dim {mh.Dim} != m dim {m.Dim}", nameof(mh));
        if (relTol <= 0) throw new ArgumentOutOfRangeException(nameof(relTol));
        if (maxIter < 1) throw new ArgumentOutOfRangeException(nameof(maxIter));

        int n = m.Dim;
        Array.Clear(x, 0, n);

        double bNorm = Norm(rhs);
        if (bNorm < 1e-300) return new Outcome(true, 0, 0.0);   // rhs = 0 => x = 0 is exact

        Complex conjShift = Complex.Conjugate(shift);
        var u = new Complex[n];
        var v = new Complex[n];
        var tmp = new Complex[n];

        // β₁u₁ = b
        Array.Copy(rhs, u, n);
        double beta = bNorm;
        Scale(u, 1.0 / beta);

        // α₁v₁ = Aᴴu₁
        CsrOps.MultiplyShifted(mh, conjShift, u, v);
        double alpha = Norm(v);
        if (alpha > 0) Scale(v, 1.0 / alpha);

        var w = (Complex[])v.Clone();
        double phiBar = beta;      // ‖r_k‖ estimate (monotone nonincreasing)
        double rhoBar = alpha;

        int iterations = maxIter;
        for (int k = 1; k <= maxIter; k++)
        {
            // β_{k+1}u_{k+1} = Av_k − α_k u_k
            CsrOps.MultiplyShifted(m, shift, v, tmp);
            for (int i = 0; i < n; i++) u[i] = tmp[i] - alpha * u[i];
            beta = Norm(u);
            if (beta > 0) Scale(u, 1.0 / beta);

            // α_{k+1}v_{k+1} = Aᴴu_{k+1} − β_{k+1}v_k
            CsrOps.MultiplyShifted(mh, conjShift, u, tmp);
            for (int i = 0; i < n; i++) v[i] = tmp[i] - beta * v[i];
            alpha = Norm(v);
            if (alpha > 0) Scale(v, 1.0 / alpha);

            // plane rotation eliminating β from the lower bidiagonal (all scalars real)
            double rho = Math.Sqrt(rhoBar * rhoBar + beta * beta);
            double c = rhoBar / rho;
            double s = beta / rho;
            double theta = s * alpha;
            rhoBar = -c * alpha;
            double phi = c * phiBar;
            phiBar = s * phiBar;

            double t1 = phi / rho;
            double t2 = -theta / rho;
            for (int i = 0; i < n; i++)
            {
                x[i] += t1 * w[i];
                w[i] = v[i] + t2 * w[i];
            }

            if (phiBar <= relTol * bNorm) { iterations = k; break; }
            // exact-arithmetic termination: the Krylov space is exhausted, the LS solution is reached
            if (alpha == 0.0 || beta == 0.0) { iterations = k; break; }
        }

        // TRUE relative residual: one extra matvec, independent of the recurrence estimate above.
        CsrOps.MultiplyShifted(m, shift, x, tmp);
        double num = 0.0;
        for (int i = 0; i < n; i++)
        {
            Complex d = tmp[i] - rhs[i];
            num += d.Real * d.Real + d.Imaginary * d.Imaginary;
        }
        double trueRel = Math.Sqrt(num) / bNorm;
        return new Outcome(trueRel <= relTol, iterations, trueRel);
    }

    private static double Norm(Complex[] v)
    {
        double s = 0;
        foreach (var z in v) s += z.Real * z.Real + z.Imaginary * z.Imaginary;
        return Math.Sqrt(s);
    }

    private static void Scale(Complex[] v, double f)
    {
        for (int i = 0; i < v.Length; i++) v[i] *= f;
    }
}
