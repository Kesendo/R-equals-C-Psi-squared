using System;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>
/// TASK 7 STEP 0 — THE MANDATORY at-N=11-scale EXCLUSION CONTROLS (sectorbraid large-N program,
/// docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md §6-§7). The N=9 calibration (SparseCensusCalibrationTests)
/// proved the sparse instruments reproduce the dense census cell for cell — but N=9's largest FD sector is
/// ~7938, so it CANNOT see the 54k-107k invit regime the six deferred N=11 cores live in, and those cores
/// have NO dense ground truth. Before any N=11 money run (Task 7 Step 1) these two controls must pass:
///
/// <para><b>(a) AT-SIZE AGREEMENT.</b> On the LARGEST PROBEABLE N=11 sector ((3,6) even, dim ~38140, just
/// under the managed-LP64 wall) at seed 0.502123, the dense σ_min and the sparse-invit σ_min must agree to
/// |sparse − dense| ≤ 0.2·dense. This is the closest at-size proxy for the deferred cores and the
/// iterations + wall-time extrapolation datum for the 54k-107k cells.</para>
///
/// <para><b>(b) PLANTED-SMALL DETECTION.</b> On a deferred-size MEMBER cell ((4,5)×λ_A, sector ~76230,
/// where the W-transport witness independently proves σ_min is genuinely small), the invit must EITHER
/// converge to a SMALL value (it can detect a hidden small σ_min at this scale and conditioning) OR
/// honestly report Converged=false. It must NOT converge to a LARGE value: that would mean large invit
/// readings cannot be trusted as exclusions, and Step 1 must NOT run (STOP, re-plan the inner engine [N7]:
/// LSQR on (M−s) at effective κ not κ², or shift-invert Lanczos on the Hermitian augmented system).</para>
///
/// <para><b>SLOW.</b> (a)'s dense LU is a 38140² complex factorization (~23 GB, tens of minutes; the test
/// host enables System.GC.AllowVeryLargeObjects for it); (b)'s invit runs to the full iteration budget
/// when it honestly bails (~100 min ceiling). Category SLOW_SHELLCENSUS_N11 — never in a fast sweep; run
/// each fact on its own (<c>--filter FullyQualifiedName~AtSizeAgreement</c> then <c>~PlantedSmall</c>) so a
/// gate failure on (a) is caught before spending (b)'s budget.</para>
/// </summary>
public sealed class SectorShellCensusN11ControlTests
{
    private readonly ITestOutputHelper _out;
    public SectorShellCensusN11ControlTests(ITestOutputHelper o) => _out = o;

    private static RealSeed N11Seed() =>
        RealDefectiveSeeds.ForN(11).Single(s => Math.Abs(s.QStar - 0.502123) < 1e-6);

    // ---- (a) AT-SIZE AGREEMENT: dense vs sparse-invit on the largest probeable N=11 sector ------------
    [Fact]
    [Trait("Category", "SLOW_SHELLCENSUS_N11")]
    public void AtSizeAgreement_LargestProbeableSector_DenseVsSparseInvit()
    {
        const int n = 11, p = 3, w = 6;   // (3,6): full block 76230, EVEN sector ~38140, just under the wall
        var seed = N11Seed();
        var (qR, lambda, pairGap) = SectorShellCensus.RefineSeed(seed);
        var q = new Complex(qR, 0);
        Complex s = lambda;               // λ_A shift; (3,6) is a NON-member here → σ_min large, well-conditioned

        // the even sector is the larger of the two R-parity sectors (dimEven = fixed + pairs ≥ pairs = dimOdd)
        var perm = WeightCoherenceBlock.ReflectionPermutation(n, p, w);
        int fixedCount = 0;
        for (int i = 0; i < perm.Length; i++) if (perm[i] == i) fixedCount++;
        int pairs = (perm.Length - fixedCount) / 2;
        int dimEven = fixedCount + pairs, dimOdd = pairs;
        _out.WriteLine($"(3,6) full={perm.Length} dimEven={dimEven} dimOdd={dimOdd} (probing EVEN); " +
                       $"qR={qR:F9} lambdaA={lambda.Real:F9} pairGap={pairGap:E3}");

        // dense ground truth: build the even sector column-major, shift the diagonal, inverse-power σ_min
        var denseSw = Stopwatch.StartNew();
        var (a, d) = WeightCoherenceBlock.BuildReflectionSectorColumnMajor(n, p, w, q, odd: false);
        for (int i = 0; i < d; i++) a[(long)i * d + i] -= s;
        var dense = ShiftedSigmaMin.EstimateColumnMajor(a, d);
        denseSw.Stop();
        a = Array.Empty<Complex>();   // release the ~23 GB before the sparse leg
        GC.Collect();

        // sparse-invit on the SAME sector at the SAME shift
        var csr = WeightCoherenceSectorCsr.BuildReflectionSector(n, p, w, q, odd: false);
        var sparseSw = Stopwatch.StartNew();
        var sparse = SparseShiftedSigmaMin.Estimate(csr, s, SparseShiftedSigmaMin.Options.Default);
        sparseSw.Stop();

        double rel = Math.Abs(sparse.SigmaMin - dense.SigmaMin) / dense.SigmaMin;
        _out.WriteLine($"  dense  σ_min={dense.SigmaMin:E6} conv={dense.Converged} iters={dense.Iterations} wall={denseSw.Elapsed.TotalMinutes:F1}min");
        _out.WriteLine($"  sparse σ_min={sparse.SigmaMin:E6} conv={sparse.Converged} outer={sparse.Outer} innerIters={sparse.InnerIterations} wall={sparseSw.Elapsed.TotalMinutes:F1}min");
        _out.WriteLine($"  relDiff={rel:F4}  (gate: sparse converged AND relDiff ≤ 0.20)  [extrapolation datum for the 54k-107k cores]");

        Assert.True(dense.Converged, $"dense σ_min did not converge on the at-size sector d={d}");
        Assert.True(sparse.Converged,
            $"[N7] STOP: sparse-invit did NOT converge on the at-size sector d={d} (dense={dense.SigmaMin:E6}); " +
            "the estimator cannot exclude at the 38140 scale, so it cannot exclude the larger deferred cores — re-plan the inner engine before Step 1.");
        Assert.True(rel <= 0.2,
            $"[N7] STOP: at-size DISAGREEMENT — sparse {sparse.SigmaMin:E6} vs dense {dense.SigmaMin:E6}, " +
            $"relDiff {rel:F4} > 0.20; large invit readings are not trustworthy at N=11 scale, do NOT run Step 1.");
    }

    // ---- (b) PLANTED-SMALL DETECTION: invit on a deferred-size member the witness proves small ---------
    [Fact]
    [Trait("Category", "SLOW_SHELLCENSUS_N11")]
    public void PlantedSmall_DeferredMemberCell_InvitDetectsSmallOrReportsFalse()
    {
        const int n = 11, p = 4, w = 5;   // (4,5)×λ_A: a band MEMBER, sector ~76230 (deferred size)
        var seed = N11Seed();
        var (qR, lambda, pairGap) = SectorShellCensus.RefineSeed(seed);
        var q = new Complex(qR, 0);
        Complex s = lambda;
        double memberTol = Math.Max(1e-9, 10.0 * pairGap);   // the census member cut at this seed

        // the witness independently proves σ_min genuinely small (a from-above bound) and names the carried parity
        var bounds = SectorWitnessTransport.MemberUpperBounds(seed, m => _out.WriteLine("[witness] " + m));
        Assert.True(bounds.TryGetValue((p, w, "lambdaA"), out var wb),
            "the witness produced no (4,5)×lambdaA bound — cannot set up the planted-small control");
        _out.WriteLine($"(4,5)×λ_A witness bound={wb.Bound:E4} carriedOdd={wb.CarriedOdd} " +
                       $"qR={qR:F9} lambdaA={lambda.Real:F9} pairGap={pairGap:E3} memberTol={memberTol:E3}");
        Assert.True(wb.Bound < memberTol,
            $"witness bound {wb.Bound:E4} is not below memberTol {memberTol:E3}; this is not a genuine planted-small cell");

        // the invit on the CARRIED-parity sector (where σ_min is genuinely small): the honest reading
        var csr = WeightCoherenceSectorCsr.BuildReflectionSector(n, p, w, q, odd: wb.CarriedOdd);
        var sw = Stopwatch.StartNew();
        var sparse = SparseShiftedSigmaMin.Estimate(csr, s, SparseShiftedSigmaMin.Options.Default);
        sw.Stop();
        _out.WriteLine($"(4,5) carried sector d={csr.Dim}: invit σ_min={sparse.SigmaMin:E6} conv={sparse.Converged} " +
                       $"outer={sparse.Outer} innerIters={sparse.InnerIterations} wall={sw.Elapsed.TotalMinutes:F1}min");

        // THE GATE: never a LARGE converged reading (that would falsely EXCLUDE a real member). Acceptable:
        // an honest Converged=false (the documented near-defective behavior), or a small converged value.
        bool falseLargeExclusion = sparse.Converged && sparse.SigmaMin > memberTol;
        Assert.False(falseLargeExclusion,
            $"[N7] STOP: the invit CONVERGED to a LARGE σ_min {sparse.SigmaMin:E6} > memberTol {memberTol:E3} on a " +
            $"witness-proven-small member (bound {wb.Bound:E4}) — large invit readings cannot be trusted as exclusions; " +
            "do NOT run Step 1, re-plan the inner engine [N7].");

        _out.WriteLine(sparse.Converged
            ? $"  PASS: invit converged SMALL ({sparse.SigmaMin:E6} ≤ memberTol {memberTol:E3}) — it detects the hidden small σ_min"
            : $"  PASS: invit honestly reported Converged=false — the near-defective member is not falsely excluded");
    }
}
