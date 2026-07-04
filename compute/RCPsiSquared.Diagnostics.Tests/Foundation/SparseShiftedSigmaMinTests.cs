using System;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>
/// The sparse σ_min must agree with the dense ShiftedSigmaMin exactly where the sparse
/// path will be TRUSTED: non-member and control cells (σ_min large, conditioning moderate).
/// At member-like shifts (near-defective, σ_min ~ (gap/2)², conditioning ~1e14) the sparse
/// estimator is NOT required to succeed; it is required to be HONEST (Converged=false or an
/// upper-bound-consistent value). The member reading belongs to the witness (Task 4).
/// </summary>
public class SparseShiftedSigmaMinTests
{
    // The recorded 6-digit lambda_A of the N=5 locus-1 seed (RealDefectiveSeeds), refined near this value.
    private const double LambdaAN5Locus1 = -4.618886;

    private readonly Xunit.Abstractions.ITestOutputHelper _out;
    public SparseShiftedSigmaMinTests(Xunit.Abstractions.ITestOutputHelper output) => _out = output;

    [Fact]
    public void Lsqr_Converges_OnTheCellThatFalsifiedBiCgStab()
    {
        // THE decisive inner-engine gate [review N7]: identity-preconditioned BiCGStab diverges on the
        // N=9 (3,3) even sector at the refined lambda_A (measured 2026-07-04: rr=6.6 after 100k iters,
        // while the same operator at a generic far shift solves in 10 iterations, a non-normal
        // eigenvalue-cluster effect). LSQR's convergence is governed by the singular value distribution
        // only, so it must solve this exact system. If this fails, the sparse instrument has no engine.
        var seed = RealDefectiveSeeds.ForN(9).Single(s => Math.Abs(s.QStar - 0.511958) < 1e-6);
        var (qRefined, lambda, _) = SectorShellCensus.RefineSeed(seed);
        var m = WeightCoherenceSectorCsr.BuildReflectionSector(9, 3, 3, new Complex(qRefined, 0), odd: false);
        var mh = CsrOps.HermitianTranspose(m);
        var rng = new Random(99);
        var rhs = new Complex[m.Dim];
        for (int i = 0; i < m.Dim; i++) rhs[i] = new Complex(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
        var x = new Complex[m.Dim];
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var outcome = LsqrSolver.Solve(m, mh, lambda, rhs, x, relTol: 1e-6, maxIter: 200_000);
        sw.Stop();
        _out.WriteLine($"LSQR N=9 (3,3) even d={m.Dim} lambda={lambda.Real:F6}: conv={outcome.Converged} " +
                       $"it={outcome.Iterations} rr={outcome.RelResidual:E3} wall={sw.Elapsed.TotalSeconds:F1}s");
        Assert.True(outcome.Converged, $"LSQR must converge where BiCGStab diverged; rr={outcome.RelResidual:E3} " +
                                       $"after {outcome.Iterations} iterations");
    }

    [Theory]
    [InlineData(2, 2)]   // in-shell non-member under the lambda_A shift (NOTE: (2,2) IS a mu-member; this test must stay on lambda_A)
    [InlineData(1, 4)]   // interior core (window-excluded; trivially conditioned; kept as an easy anchor)
    public void NonMemberCells_MatchDense(int p, int w)
    {
        var seed = RealDefectiveSeeds.ForN(5).Single(s => Math.Abs(s.QStar - 0.620878) < 1e-6);
        var (qRefined, lambda, _) = SectorShellCensus.RefineSeed(seed);

        // Guard: the probe shift MUST be lambda_A (near -4.6189), never mu = -lambda - 2N (~ -5.38).
        // (2,2) is a mu-member; a silent switch to mu would turn this into a member cell and defeat
        // the whole point of the non-member agreement test.
        Assert.True(Math.Abs(lambda.Real - LambdaAN5Locus1) < 5e-3,
            $"shift must be lambda_A ~ {LambdaAN5Locus1}, got {lambda.Real}");
        var q = new Complex(qRefined, 0);

        foreach (bool odd in new[] { false, true })
        {
            // Dense reference: build the column-major sector, subtract the shift on the diagonal in
            // place, run the LU-based estimator (exactly the census's Run recipe).
            var (a, d) = WeightCoherenceBlock.BuildReflectionSectorColumnMajor(5, p, w, q, odd);
            if (d == 0) continue;                                  // no odd sector on all-palindromic blocks
            for (int i = 0; i < d; i++) a[(long)i * d + i] -= lambda;
            var dense = ShiftedSigmaMin.EstimateColumnMajor(a, d);

            var sparse = SparseShiftedSigmaMin.Estimate(
                WeightCoherenceSectorCsr.BuildReflectionSector(5, p, w, q, odd), lambda,
                SparseShiftedSigmaMin.Options.Default);
            _out.WriteLine($"N=5 ({p},{w}) odd={odd} d={d}: dense={dense.SigmaMin:E6} sparse={sparse.SigmaMin:E6} " +
                           $"outer={sparse.Outer} innerIters={sparse.InnerIterations}");

            Assert.True(sparse.Converged,
                $"({p},{w}) odd={odd}: sparse must converge at a non-member cell (relres too high)");
            Assert.True(Math.Abs(sparse.SigmaMin - dense.SigmaMin) <= 0.2 * dense.SigmaMin,
                $"({p},{w}) odd={odd}: sparse {sparse.SigmaMin:E6} vs dense {dense.SigmaMin:E6}");
        }
    }

    [Fact]
    public void NonMemberCell_RepresentativeRegime_N9_MatchesDense()
    {
        // The representative sigma_min ~1e-3 regime cell [review N11]: (3,3) at the N=9 R-odd seed,
        // both parity sectors (dims ~2.3k, dense LU fine).
        var seed = RealDefectiveSeeds.ForN(9).Single(s => Math.Abs(s.QStar - 0.511958) < 1e-6);
        var (qRefined, lambda, _) = SectorShellCensus.RefineSeed(seed);
        var q = new Complex(qRefined, 0);
        const int p = 3, w = 3;

        foreach (bool odd in new[] { false, true })
        {
            var (a, d) = WeightCoherenceBlock.BuildReflectionSectorColumnMajor(9, p, w, q, odd);
            if (d == 0) continue;
            for (int i = 0; i < d; i++) a[(long)i * d + i] -= lambda;
            var dense = ShiftedSigmaMin.EstimateColumnMajor(a, d);

            var sw = System.Diagnostics.Stopwatch.StartNew();
            var sparse = SparseShiftedSigmaMin.Estimate(
                WeightCoherenceSectorCsr.BuildReflectionSector(9, p, w, q, odd), lambda,
                SparseShiftedSigmaMin.Options.Default);
            sw.Stop();
            // the Task-7 Step-0 feasibility datum: LSQR inner iterations + wall time at the N=9 cell
            _out.WriteLine($"N=9 (3,3) odd={odd} d={d}: dense={dense.SigmaMin:E6} sparse={sparse.SigmaMin:E6} " +
                           $"outer={sparse.Outer} innerIters={sparse.InnerIterations} wall={sw.Elapsed.TotalSeconds:F1}s");

            Assert.True(sparse.Converged,
                $"(3,3) odd={odd}: sparse must converge in the sigma_min ~1e-3 regime");
            Assert.True(Math.Abs(sparse.SigmaMin - dense.SigmaMin) <= 0.2 * dense.SigmaMin,
                $"(3,3) odd={odd}: sparse {sparse.SigmaMin:E6} vs dense {dense.SigmaMin:E6}");
        }
    }

    [Fact]
    public void MemberCell_IsHonest_NeverSilentlyWrong()
    {
        // (2,3) x lambda_A, even sector: a member block probed at the near-defective shift. The sparse
        // estimator is not required to converge here; it must be HONEST. On convergence the value must
        // be upper-bound-consistent with the dense reference; on non-convergence it must be an honest
        // non-answer (NaN or a positive last estimate). The test fails ONLY on a silently-wrong value.
        var seed = RealDefectiveSeeds.ForN(5).Single(s => Math.Abs(s.QStar - 0.620878) < 1e-6);
        var (qRefined, lambda, _) = SectorShellCensus.RefineSeed(seed);
        var q = new Complex(qRefined, 0);
        const int p = 2, w = 3;
        const bool odd = false;   // even sector

        var (a, d) = WeightCoherenceBlock.BuildReflectionSectorColumnMajor(5, p, w, q, odd);
        for (int i = 0; i < d; i++) a[(long)i * d + i] -= lambda;
        var dense = ShiftedSigmaMin.EstimateColumnMajor(a, d);

        var r = SparseShiftedSigmaMin.Estimate(
            WeightCoherenceSectorCsr.BuildReflectionSector(5, p, w, q, odd), lambda,
            SparseShiftedSigmaMin.Options.Default);
        _out.WriteLine($"N=5 (2,3) even (member) d={d}: dense={dense.SigmaMin:E6} sparse={r.SigmaMin:E6} " +
                       $"conv={r.Converged} outer={r.Outer} innerIters={r.InnerIterations}");

        if (r.Converged)
            Assert.True(r.SigmaMin >= dense.SigmaMin * 0.5 && r.SigmaMin <= 10 * Math.Sqrt(dense.SigmaMin),
                $"converged value must be upper-bound-consistent: sparse {r.SigmaMin:E6} vs dense {dense.SigmaMin:E6}");
        else
            Assert.True(double.IsNaN(r.SigmaMin) || r.SigmaMin > 0,
                $"an honest non-answer expected: got {r.SigmaMin}");
    }
}
