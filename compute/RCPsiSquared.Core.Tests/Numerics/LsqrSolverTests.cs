using System;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>The LSQR kernel (the [N7] fallback inner engine for the sparse σ_min estimator): the
/// moderate-conditioning solve, the honest-failure pattern on budget exhaustion, and the (m, mh)
/// adjoint-swap convention the estimator's second solve relies on. The decisive test (LSQR converges
/// on the exact N=9 cell that falsified BiCGStab) lives in SparseShiftedSigmaMinTests, which can
/// reach SectorShellCensus.RefineSeed.</summary>
public class LsqrSolverTests
{
    private static WeightCoherenceSectorCsr.Csr SmallBlock()
        => WeightCoherenceSectorCsr.BuildFull(5, 2, 3, new Complex(0.7, 0.0));

    [Fact]
    public void Lsqr_SolvesShiftedSystem_AtModerateConditioning()
    {
        // The same moderate shifted system the BiCGStab kernel test solves.
        var m = SmallBlock();
        var shift = new Complex(-3.0, 0.4);       // generic: far from the spectrum's tight spots
        var rng = new Random(11);
        var rhs = new Complex[m.Dim];
        for (int i = 0; i < m.Dim; i++) rhs[i] = new Complex(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
        var x = new Complex[m.Dim];
        var outcome = LsqrSolver.Solve(m, shift, rhs, x, relTol: 1e-10, maxIter: 20000);
        Assert.True(outcome.Converged, $"stalled at relres {outcome.RelResidual}");
        var r = new Complex[m.Dim];
        CsrOps.MultiplyShifted(m, shift, x, r);
        double num = 0, den = 0;
        for (int i = 0; i < m.Dim; i++) { num += MagSq(r[i] - rhs[i]); den += MagSq(rhs[i]); }
        Assert.True(Math.Sqrt(num / den) <= 1e-8, "true residual must match the reported one");
    }

    [Fact]
    public void Lsqr_ReportsHonestFailure_WhenIterationBudgetExhausted()
    {
        // Same honesty contract as the BiCGStab kernel: a starved budget at a near-singular shift
        // must report Converged=false keyed on the recomputed TRUE residual, never a false success.
        var m = WeightCoherenceSectorCsr.BuildReflectionSector(5, 1, 2, new Complex(0.620878, 0.0), odd: false);
        var shift = new Complex(-3.95, 0.0);
        var rng = new Random(23);
        var rhs = new Complex[m.Dim];
        for (int i = 0; i < m.Dim; i++) rhs[i] = new Complex(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
        var x = new Complex[m.Dim];
        var outcome = LsqrSolver.Solve(m, shift, rhs, x, relTol: 1e-10, maxIter: 5);
        Assert.False(outcome.Converged);
        Assert.True(outcome.RelResidual > 1e-10, $"expected an honest large residual, got {outcome.RelResidual}");
    }

    [Fact]
    public void Lsqr_AdjointSwapConvention_SolvesTheHermitianSystem()
    {
        // The estimator's second solve: (M−s)ᴴ w = z is LSQR on the swapped pair, Solve(mh, m,
        // conj(shift), ...), because the forward operator Mᴴ − conj(s)·I has adjoint M − s·I. This
        // pins the swap convention so a future edit cannot silently break the estimator's step 2.
        var m = SmallBlock();
        var mh = CsrOps.HermitianTranspose(m);
        var shift = new Complex(-3.0, 0.4);
        var rng = new Random(17);
        var z = new Complex[m.Dim];
        for (int i = 0; i < m.Dim; i++) z[i] = new Complex(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
        var w = new Complex[m.Dim];
        var outcome = LsqrSolver.Solve(mh, m, Complex.Conjugate(shift), z, w, relTol: 1e-10, maxIter: 20000);
        Assert.True(outcome.Converged, $"stalled at relres {outcome.RelResidual}");
        // verify against the explicit forward apply of (M−s)ᴴ = Mᴴ − conj(s)·I
        var r = new Complex[m.Dim];
        CsrOps.MultiplyShifted(mh, Complex.Conjugate(shift), w, r);
        double num = 0, den = 0;
        for (int i = 0; i < m.Dim; i++) { num += MagSq(r[i] - z[i]); den += MagSq(z[i]); }
        Assert.True(Math.Sqrt(num / den) <= 1e-8, "true residual of the Hermitian system must be small");
    }

    private static double MagSq(Complex z) => z.Real * z.Real + z.Imaginary * z.Imaginary;
}
