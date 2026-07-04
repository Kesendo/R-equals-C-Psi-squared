using System;
using System.Collections.Generic;
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
/// TASK 5 — THE N=9 SPARSE-VS-DENSE CALIBRATION GATE (sectorbraid large-N program,
/// docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md §6-§7). The pre-registered decision gate the sparse
/// instruments must clear BEFORE any N=11 run: at both defective-seed parities the sparse member
/// instrument (<see cref="SectorWitnessTransport.MemberUpperBounds"/>) and the sparse non-member
/// instrument (<see cref="SparseShiftedSigmaMin"/> over the census CSR) must reproduce the DENSE
/// census (<see cref="SectorShellCensus.Run"/>, ground truth) cell for cell.
///
/// <para><b>DECISION RULE (verbatim).</b> If ANY threshold fails, the sparse path is NOT fit for
/// N=11; stop, report the failing cells and numbers, and re-plan. This is a decision gate, not a
/// tuning knob; a silently loosened threshold is a falsified instrument.</para>
///
/// <para><b>KNOWN BLIND SPOT (verbatim).</b> N=9's largest FD sector is ~7938, so this gate cannot
/// see the 54k-107k invit regime; Task 7 Step 0 exists for that and is mandatory.</para>
///
/// <para><b>Pre-registered thresholds (immutable).</b> Member cells: bound &lt; memberTol/100 with
/// memberTol = 10·pairGap (the census's own adaptive cut); bound &gt;= 0.5·denseCarried with
/// denseCarried the dense entry's σ column in the CARRIED parity; and the carried-parity column is
/// the SMALL one (≤ the other parity, pinning the F7 alternation law against the dense ground truth
/// at both seed parities). Non-member and control cells: sparse.Converged; and |sparse − dense| ≤
/// 0.2·dense per parity, at the entry's own shift value. Per-cell wall time + InnerIterations are
/// recorded (the Task-7 N=11 feasibility extrapolation data); ratios are recorded, never gating.</para>
/// </summary>
public sealed class SparseCensusCalibrationTests
{
    private readonly ITestOutputHelper _out;
    public SparseCensusCalibrationTests(ITestOutputHelper o) => _out = o;

    // The two N=9 calibration seeds (RealDefectiveSeeds.ForN(9)): the R-odd and R-even loci. One
    // Theory row per seed => per-seed failure accumulation (a gate failure reports EVERY failing
    // cell at once, never just the first).
    [Theory]
    [Trait("Category", "SHELLCENSUS_SPARSE_CAL")]
    [InlineData(0.511958, -1)]   // R-odd  calibration seed (lambda_A = -4.1581)
    [InlineData(0.591760, +1)]   // R-even calibration seed (lambda_A = -5.1206)
    public void SparseCensus_ReproducesDense_AtN9(double qStar, int rParity)
    {
        const int n = 9;
        var seed = RealDefectiveSeeds.ForN(n).Single(s => Math.Abs(s.QStar - qStar) < 1e-6 && s.RParity == rParity);

        var gate = Stopwatch.StartNew();
        var failures = new List<string>();

        // ---- (1) DENSE CENSUS = ground truth (default Options; nothing deferred at N=9) ----------
        var opts = new SectorShellCensus.Options { Log = m => _out.WriteLine("[census] " + m) };
        var denseSw = Stopwatch.StartNew();
        var census = SectorShellCensus.Run(seed, opts);
        denseSw.Stop();
        Assert.True(census.SeedUsable,
            $"seed q*={qStar} R={rParity} must refine to a usable pair gap; got {census.PairGap:E3}");

        double qRefined = census.QRefined;                 // Run exposes the refined locus: no second RefineSeed
        Complex lambda = census.RefinedLambdaA;
        double pairGap = census.PairGap;
        double memberTol = 10.0 * pairGap;                 // pre-registered; = the census adaptive cut at N=9
        var q = new Complex(qRefined, 0);
        Complex mu = -lambda - 2.0 * n;

        var summary = census.Summarize();
        _out.WriteLine($"=== SEED q*={qStar} R={rParity} (N={n}) ===");
        _out.WriteLine($"qRefined={qRefined:F9} lambdaA={lambda.Real:F9} mu={mu.Real:F9} pairGap={pairGap:E3} " +
                       $"memberTol={memberTol:E3} censusVerdict={summary.Verdict} denseWall={denseSw.Elapsed.TotalMinutes:F1}min");

        // ---- (2) MEMBER CELLS: the from-above witness bound vs the dense census placement ---------
        var bounds = SectorWitnessTransport.MemberUpperBounds(seed, m => _out.WriteLine("[witness] " + m));
        var expectedMembers = SectorShellCensus.ExpectedMembers(n);
        _out.WriteLine("--- member cells (bound / denseCarried) ---");
        foreach (var kv in bounds.OrderBy(k => k.Key.P).ThenBy(k => k.Key.W).ThenBy(k => k.Key.Shift))
        {
            var key = kv.Key;
            var val = kv.Value;
            var e = census.Entries.FirstOrDefault(x => x.P == key.P && x.W == key.W && x.Shift == key.Shift);
            if (e is null || !e.Probed)
            {
                failures.Add($"member ({key.P},{key.W},{key.Shift}): no probed dense entry found");
                continue;
            }

            double carriedCol = val.CarriedOdd ? e.SigmaMinOdd : e.SigmaMinEven;   // the CARRIED parity column
            double otherCol = val.CarriedOdd ? e.SigmaMinEven : e.SigmaMinOdd;
            double ratio = carriedCol > 0 && double.IsFinite(carriedCol) ? val.Bound / carriedCol : double.NaN;
            _out.WriteLine($"({key.P},{key.W}) {key.Shift} carriedOdd={val.CarriedOdd}: bound={val.Bound:E4} " +
                           $"denseCarried={carriedCol:E4} denseOther={otherCol:E4} ratio={ratio:F3} " +
                           $"memberTol/100={memberTol / 100:E3}");

            // (a) deep from-above reading, below the pre-registered census cut / 100
            if (!(val.Bound < memberTol / 100))
                failures.Add($"member ({key.P},{key.W},{key.Shift}): bound {val.Bound:E4} " +
                             $"NOT < memberTol/100 {memberTol / 100:E4}");
            // (b) from-above, not collapsed to a fabricated 0: within a factor 2 of the dense carried value
            if (!(val.Bound >= 0.5 * carriedCol))
                failures.Add($"member ({key.P},{key.W},{key.Shift}): bound {val.Bound:E4} " +
                             $"NOT >= 0.5*denseCarried {0.5 * carriedCol:E4}");
            // (c) the F7 alternation law: the carried parity is the SMALL dense column
            if (!(carriedCol <= otherCol))
                failures.Add($"member ({key.P},{key.W},{key.Shift}): carried parity (odd={val.CarriedOdd}) " +
                             $"col {carriedCol:E4} NOT <= other {otherCol:E4}");
        }
        // no silently missing member cell: every ExpectedMembers key must have been read by the witness
        foreach (var m in expectedMembers)
            if (!bounds.ContainsKey(m))
                failures.Add($"member {m}: ExpectedMembers key missing from MemberUpperBounds");

        // ---- (3) NON-MEMBER + CONTROL CELLS: sparse LSQR vs the dense census, per parity ----------
        _out.WriteLine("--- non-member / control cells (sparse vs dense, per parity) ---");
        var nonMembers = census.Entries
            .Where(e => e.Probed && (e.Method == "lu-rparity" || e.Method == "control")
                        && !expectedMembers.Contains((e.P, e.W, e.Shift)))
            .OrderBy(e => e.P).ThenBy(e => e.W).ThenBy(e => e.Shift)
            .ToList();

        foreach (var e in nonMembers)
        {
            foreach (bool odd in new[] { false, true })
            {
                double dense = odd ? e.SigmaMinOdd : e.SigmaMinEven;
                if (!double.IsFinite(dense)) continue;                 // empty parity sector: nothing to compare

                var csr = WeightCoherenceSectorCsr.BuildReflectionSector(n, e.P, e.W, q, odd);
                if (csr.Dim == 0) continue;

                var cellSw = Stopwatch.StartNew();
                var sparse = SparseShiftedSigmaMin.Estimate(csr, e.ShiftValue, SparseShiftedSigmaMin.Options.Default);
                cellSw.Stop();

                double rel = dense > 0 ? Math.Abs(sparse.SigmaMin - dense) / dense : double.NaN;
                _out.WriteLine($"({e.P},{e.W}) {e.Shift} odd={odd} d={csr.Dim}: dense={dense:E6} sparse={sparse.SigmaMin:E6} " +
                               $"relDiff={rel:F4} conv={sparse.Converged} outer={sparse.Outer} " +
                               $"innerIters={sparse.InnerIterations} wall={cellSw.Elapsed.TotalSeconds:F1}s");

                if (!sparse.Converged)
                    failures.Add($"nonmember ({e.P},{e.W},{e.Shift}) odd={odd} d={csr.Dim}: sparse did NOT converge " +
                                 $"(dense={dense:E6} sparse={sparse.SigmaMin:E6})");
                else if (!(Math.Abs(sparse.SigmaMin - dense) <= 0.2 * dense))
                    failures.Add($"nonmember ({e.P},{e.W},{e.Shift}) odd={odd} d={csr.Dim}: sparse {sparse.SigmaMin:E6} " +
                                 $"vs dense {dense:E6} relDiff {rel:F4} > 0.20");
            }
        }

        gate.Stop();
        _out.WriteLine($"=== SEED q*={qStar} R={rParity}: gate wall {gate.Elapsed.TotalMinutes:F1}min, " +
                       $"{nonMembers.Count} non-member cells, {failures.Count} failing checks ===");

        // Accumulated assertion: report EVERY failing cell at once (a gate failure is a decision-gate
        // STOP, not a tuning knob; the caller re-plans from the full table, not from the first failure).
        Assert.True(failures.Count == 0,
            $"N=9 SPARSE CALIBRATION GATE FAILED at seed q*={qStar} R={rParity} ({failures.Count} failing checks). " +
            $"DECISION-GATE STOP: the sparse path is NOT fit for N=11 until these are understood.\n  " +
            string.Join("\n  ", failures));
    }
}
