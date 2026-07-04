using System;
using System.Numerics;
using RCPsiSquared.Core.BlockSpectrum.JordanWigner;
using RCPsiSquared.Core.F89PathK;

namespace RCPsiSquared.Core.Numerics;

/// <summary>
/// BiCGStab solve of the shifted system (M − shift·I)·x = rhs on a <see cref="WeightCoherenceSectorCsr.Csr"/>
/// weight-coherence sector, with the IDENTITY preconditioner hard-wired. The algorithm body is copied verbatim
/// from <see cref="BlockSpectrum.SparseShiftInvertArnoldi"/>'s private SolveShiftedSystem (a single point of
/// algorithmic truth for BiCGStab in this repo), with the preconditioner fixed to M⁻¹ = I: these coherence
/// blocks are structurally NOT diagonally dominant (the diagonal is only the dissipator −2·n_diff, much smaller
/// than the ±2iq hopping, and many rows carry an exact-zero diagonal), so Jacobi (diagonal) preconditioning
/// amplifies near-zero shifted diagonals to ~1/shift and destabilises BiCGStab. This is the same lesson pinned
/// by the Klein-4-group self-paired sub-block path (see SparseShiftInvertArnoldi's PreconditionerKind docstring
/// and the KleinFourGroupSelfPaired* tests): identity preconditioning preserves the natural conditioning of
/// (M − shift·I).
///
/// <para>At exit the TRUE relative residual ‖(M − shift·I)·x − rhs‖ / ‖rhs‖ is recomputed with one extra matvec
/// (not the recurrence residual, which can drift from the real one on ill-conditioned shifts), and
/// <see cref="Outcome.Converged"/> is that true residual being ≤ relTol. The near-zero rho / omega / rt·v / t·t
/// breakdown guards from the source are preserved as early loop exits; a breakdown that stalls leaves a large
/// true residual, so it reports Converged = false, while a lucky breakdown (r already ≈ 0) reports the achieved
/// small residual honestly.</para>
/// </summary>
public static class BiCgStabSolver
{
    /// <summary>Outcome of a <see cref="Solve"/> call: whether the recomputed TRUE relative residual met the
    /// tolerance, the number of BiCGStab iterations taken (or where the loop exited), and that true relative
    /// residual.</summary>
    public readonly record struct Outcome(bool Converged, int Iterations, double RelResidual);

    /// <summary>Right-preconditioned BiCGStab for (M − shift·I)·x = rhs with M⁻¹ = I (identity preconditioning;
    /// see the class docstring for why Jacobi is unusable on these non-diagonally-dominant blocks). Initial guess
    /// x ≡ 0. rhs and x must be length m.Dim; x is overwritten with the solution. Returns the true-residual-based
    /// <see cref="Outcome"/>.</summary>
    public static Outcome Solve(WeightCoherenceSectorCsr.Csr m, Complex shift, Complex[] rhs, Complex[] x,
        double relTol, int maxIter)
    {
        if (m is null) throw new ArgumentNullException(nameof(m));
        if (rhs is null) throw new ArgumentNullException(nameof(rhs));
        if (x is null) throw new ArgumentNullException(nameof(x));
        if (rhs.Length != m.Dim) throw new ArgumentException($"rhs length {rhs.Length} != dim {m.Dim}", nameof(rhs));
        if (x.Length != m.Dim) throw new ArgumentException($"x length {x.Length} != dim {m.Dim}", nameof(x));
        if (relTol <= 0) throw new ArgumentOutOfRangeException(nameof(relTol));
        if (maxIter < 1) throw new ArgumentOutOfRangeException(nameof(maxIter));

        int n = m.Dim;
        var r = new Complex[n];
        var rt = new Complex[n];
        var p = new Complex[n];
        var v = new Complex[n];
        var s = new Complex[n];
        var t = new Complex[n];
        var z = new Complex[n];   // z = M⁻¹·p ; identity preconditioner => z = p
        var y = new Complex[n];   // y = M⁻¹·s ; identity preconditioner => y = s

        Array.Clear(x, 0, n);
        Array.Copy(rhs, r, n);
        Array.Copy(rhs, rt, n);

        double bNorm = Math.Sqrt(KrylovOps.NormSquared(rhs));
        if (bNorm < 1e-300) return new Outcome(true, 0, 0.0);   // rhs = 0 => x = 0 is exact

        Complex rho = Complex.One, alpha = Complex.One, omega = Complex.One;
        Complex rhoPrev;
        int iterations = maxIter;

        for (int k = 1; k <= maxIter; k++)
        {
            rhoPrev = rho;
            rho = KrylovOps.ConjugateDot(rt, r);
            if (rho.Magnitude < 1e-300) { iterations = k; break; }          // breakdown

            if (k == 1)
            {
                Array.Copy(r, p, n);
            }
            else
            {
                if (omega.Magnitude < 1e-300) { iterations = k; break; }    // breakdown
                Complex beta = (rho / rhoPrev) * (alpha / omega);
                for (int i = 0; i < n; i++) p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }

            for (int i = 0; i < n; i++) z[i] = p[i];                        // identity preconditioner
            CsrOps.MultiplyShifted(m, shift, z, v);
            Complex rtv = KrylovOps.ConjugateDot(rt, v);
            if (rtv.Magnitude < 1e-300) { iterations = k; break; }          // breakdown
            alpha = rho / rtv;

            for (int i = 0; i < n; i++) s[i] = r[i] - alpha * v[i];
            double sNorm = Math.Sqrt(KrylovOps.NormSquared(s));
            if (sNorm / bNorm < relTol)
            {
                for (int i = 0; i < n; i++) x[i] += alpha * z[i];
                iterations = k;
                break;
            }

            for (int i = 0; i < n; i++) y[i] = s[i];                        // identity preconditioner
            CsrOps.MultiplyShifted(m, shift, y, t);
            Complex tt = KrylovOps.ConjugateDot(t, t);
            if (tt.Magnitude < 1e-300) { iterations = k; break; }           // breakdown
            omega = KrylovOps.ConjugateDot(t, s) / tt;

            for (int i = 0; i < n; i++)
            {
                x[i] += alpha * z[i] + omega * y[i];
                r[i] = s[i] - omega * t[i];
            }

            double rNorm = Math.Sqrt(KrylovOps.NormSquared(r));
            if (rNorm / bNorm < relTol) { iterations = k; break; }
            if (omega.Magnitude < 1e-300) { iterations = k; break; }        // breakdown
        }

        // TRUE relative residual: one extra matvec, independent of the recurrence residual above.
        var res = new Complex[n];
        CsrOps.MultiplyShifted(m, shift, x, res);
        double num = 0.0;
        for (int i = 0; i < n; i++)
        {
            Complex d = res[i] - rhs[i];
            num += d.Real * d.Real + d.Imaginary * d.Imaginary;
        }
        double trueRel = Math.Sqrt(num) / bNorm;
        return new Outcome(trueRel <= relTol, iterations, trueRel);
    }
}
