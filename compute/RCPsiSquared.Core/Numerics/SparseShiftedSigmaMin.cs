using System;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;

namespace RCPsiSquared.Core.Numerics;

/// <summary>
/// σ_min(M − s) for a CSR block via inverse power iteration on (M−s)ᴴ(M−s), the two linear
/// solves per step done ITERATIVELY (LSQR). Mirrors the dense ShiftedSigmaMin recipe, but under
/// INEXACT inner solves the result is an ESTIMATE, not a rigorous from-above bound [review N5].
/// Intended for cells where σ_min is LARGE (non-member/exclusion and Bendixson-control cells). At
/// near-defective member cells (σ_min ~ (gap/2)² ~ 1e-13) the inner solves cannot reach meaningful
/// residuals in double precision; the estimator must then report Converged=false (honesty over
/// availability); the member reading is SectorWitnessTransport's job (Task 4).
///
/// <para>Inner engine history: the first engine was identity-preconditioned BiCGStab
/// (<see cref="BiCgStabSolver"/>, which remains a landed, tested kernel) and was FALSIFIED by the
/// N=9 (3,3) × λ_A cell (2026-07-04: every solve diverged, rr 6.6+ after 100k iterations, while the
/// same operator at a generic far shift solved in 10; a non-normal eigenvalue-cluster effect that
/// no small imaginary offset of the shift escapes). <see cref="LsqrSolver"/> is the review-sanctioned
/// replacement [N7]: Golub-Kahan bidiagonalization works implicitly on the Hermitian PSD operator
/// (M−s)ᴴ(M−s), so convergence is governed by the singular value distribution only and the cluster
/// divergence is structurally invisible to it (the same cell: 3329 iterations, 1.5 s).</para>
///
/// <para>Bookkeeping mirror of <see cref="ShiftedSigmaMin.EstimateColumnMajor"/>: v enters each
/// outer step unit-normalized, the two solves apply [(M−s)(M−s)ᴴ]⁻¹, growth ‖w‖ → 1/σ_min², the
/// estimate 1/√growth is an upper bound on σ_min at every step, and convergence is declared only
/// after the dense source's TWO-in-a-row relative-change stall guard [review N8]. LSQR's residual is
/// monotone (no BiCGStab-style breakdown), so an inner solve that exhausts its budget is a genuine
/// hard-cell signal: the estimate and stall counter freeze on such steps (only CLEAN steps, both
/// inner solves converged, confirm convergence), v is still advanced (the exhausted LSQR iterate is
/// the monotone least-squares iterate, dominated by the small-singular directions, so the retry runs
/// a genuinely different RHS), and the run bails to Converged=false after
/// <see cref="MaxConsecutiveInnerFailures"/> consecutive failed steps.</para>
/// </summary>
public static class SparseShiftedSigmaMin
{
    public sealed record Options(int MaxOuter, double RelTol, double InnerRelTol, int InnerMaxIter, int Seed)
    {
        // InnerRelTol 1e-8, NOT 1e-10 [N5]: at the worst trusted (non-member) cells kappa(M−s) ~ 1e4-1e5,
        // so the attainable iterative floor is ~kappa*eps_mach ~ 1e-12-1e-11; 1e-8 clears it with margin
        // while the outer RelTol 1e-3 only needs the applied operator good to ~1e-3. InnerMaxIter 50_000:
        // the measured worst trusted cell (N=9 (3,3) even at lambda_A, d=3536) takes 3329 LSQR iterations
        // to 1e-6 from a random RHS; 50k covers the 1e-8 target with an order of magnitude of headroom.
        public static Options Default { get; } = new(MaxOuter: 200, RelTol: 1e-3, InnerRelTol: 1e-8,
                                                     InnerMaxIter: 50_000, Seed: 12345);
    }

    public readonly record struct Result(double SigmaMin, int Outer, bool Converged, long InnerIterations);

    /// <summary>Consecutive failed outer steps (an inner LSQR solve exhausted its budget) tolerated
    /// before the run bails to Converged=false. LSQR has no breakdown and its residual is monotone, so
    /// a budget exhaustion is a systematic conditioning verdict, not an unlucky recurrence; measured
    /// trusted cells (N=5 strip, N=9 (3,3)) show ZERO inner failures under LSQR, so 3 is pure escape
    /// margin (the advanced RHS gets two retries) while keeping a true member-cell probe cheap.</summary>
    private const int MaxConsecutiveInnerFailures = 3;

    public static Result Estimate(WeightCoherenceSectorCsr.Csr m, Complex shift, Options opts)
    {
        if (m is null) throw new ArgumentNullException(nameof(m));
        if (opts is null) throw new ArgumentNullException(nameof(opts));
        int d = m.Dim;
        if (d == 0) return new Result(double.NaN, 0, false, 0L);

        // (M − s)ᴴ = Mᴴ − conj(s)·I, since (M − sI)ᴴ = Mᴴ − s̄·Iᴴ = Mᴴ − s̄·I. The forward solve uses
        // the pair (m, mh) with shift; the Hermitian solve uses the SWAPPED pair (mh, m) with conj(s),
        // because the operator Mᴴ − conj(s)·I has adjoint (Mᴴ − conj(s)·I)ᴴ = M − s·I. Both directions
        // reuse the one Hermitian transpose computed here.
        var mh = CsrOps.HermitianTranspose(m);

        var rng = new Random(opts.Seed);
        var v = new Complex[d];
        for (int i = 0; i < d; i++) v[i] = new Complex(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
        Normalize(v);

        var z = new Complex[d];
        var w = new Complex[d];

        double sigma = double.PositiveInfinity;
        int stallCount = 0;
        int consecutiveFail = 0;
        long innerTotal = 0;
        for (int outer = 1; outer <= opts.MaxOuter; outer++)
        {
            // v ← (M−s)⁻ᴴ (M−s)⁻¹ v, in place via two iterative solves: growth ‖v‖ → 1/σ_min²
            // (v enters unit-normalized). Iteration matrix (M−s)⁻ᴴ(M−s)⁻¹ = [(M−s)(M−s)ᴴ]⁻¹, whose
            // dominant eigenvalue is 1/σ_min(M−s)², the same value the dense recipe converges to.
            var o1 = LsqrSolver.Solve(m, mh, shift, v, z, opts.InnerRelTol, opts.InnerMaxIter);
            innerTotal += o1.Iterations;
            var o2 = LsqrSolver.Solve(mh, m, Complex.Conjugate(shift), z, w, opts.InnerRelTol, opts.InnerMaxIter);
            innerTotal += o2.Iterations;

            double growth = Norm(w);
            bool advanceable = growth > 0.0 && double.IsFinite(growth);
            // Advance v regardless of inner accuracy: even a budget-exhausted LSQR iterate is dominated
            // by the smallest-singular-vector direction, which is what the power iteration must amplify;
            // this also changes the RHS so the retry after a failed step is not the identical solve.
            if (advanceable)
                for (int i = 0; i < d; i++) v[i] = w[i] / growth;

            if (o1.Converged && o2.Converged && advanceable)
            {
                consecutiveFail = 0;
                double next = 1.0 / Math.Sqrt(growth);
                // Dense-source stall guard: require the relative-change criterion TWICE in a row before
                // declaring convergence (guards the near-degenerate member-pair stall). Only CLEAN steps
                // count toward the stall, so a failed step between two clean small-change steps is
                // bridged, never turning a locked non-member estimate into a false negative.
                bool small = Math.Abs(next - sigma) <= opts.RelTol * next;
                stallCount = small ? stallCount + 1 : 0;
                sigma = next;
                if (stallCount >= 2 || sigma < 1e-14) return new Result(sigma, outer, true, innerTotal);
            }
            else if (++consecutiveFail >= MaxConsecutiveInnerFailures)
            {
                // Persistent inner failure without the estimate ever locking: the near-defective member
                // signature. Report the honest non-answer (Converged=false); see Report for the value.
                return new Result(Report(sigma), outer, false, innerTotal);
            }
        }
        return new Result(Report(sigma), opts.MaxOuter, false, innerTotal);
    }

    // On a non-converged exit the outer estimate is not trustworthy as a verdict: report the last CLEAN
    // estimate if one exists (an honest upper bound so far), else NaN when no clean step ever completed.
    // Either way Converged=false carries the honesty; callers must not read the value as a membership call.
    private static double Report(double sigma) => double.IsInfinity(sigma) ? double.NaN : sigma;

    private static double Norm(Complex[] v)
    {
        double s = 0;
        foreach (var z in v) s += z.Real * z.Real + z.Imaginary * z.Imaginary;
        return Math.Sqrt(s);
    }

    private static void Normalize(Complex[] v)
    {
        double n = Norm(v);
        for (int i = 0; i < v.Length; i++) v[i] /= n;
    }
}
