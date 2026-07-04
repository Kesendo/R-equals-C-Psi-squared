using System;
using System.Linq;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>
/// TASK 6 — THE SPARSE PROBE DOCKED INTO THE CENSUS DEFERRED BRANCH (sectorbraid large-N program,
/// docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md §6-§7). The census's deferred branch (a sector past the
/// managed-LP64 wall) no longer emits a bare "deferred" row: with <see cref="SectorShellCensus.Options.SparseForDeferred"/>
/// on (the default) it routes MEMBER cells to the from-above W-transport witness
/// (<see cref="SectorWitnessTransport.MemberUpperBounds"/>, method "sparse-witness") and NON-MEMBER /
/// control cells to the sparse inverse-power estimator (<see cref="SparseShiftedSigmaMin"/>, method
/// "sparse-invit"), honest <c>Converged</c> flags throughout.
///
/// <para>THE LOWERED-WALL TRICK: at N=5 the full answer is a fast dense gate, so we run the SAME seed
/// twice — once dense (MaxSectorDim=46000, nothing defers) and once with MaxSectorDim=40 so the bigger
/// blocks are forced onto the sparse path — and demand the sparse run reproduce the dense one cell for
/// cell within the pre-registered Task-5 thresholds, with the twin verdict PASS (witness-assisted).</para>
///
/// <para>Two guards beyond the twin: (i) with SparseForDeferred=false the old bare "deferred" rows are
/// preserved byte-compatibly (no regression of the landed behavior); (ii) [F1] a probed non-member whose
/// sparse estimator does NOT converge forces PARTIAL with the cell named, NEVER a silent PASS — the
/// exclusion half is the claim's promotion blocker, so converged=0 can never contribute to a PASS.</para>
/// </summary>
public sealed class SectorShellCensusSparseTests
{
    private readonly ITestOutputHelper _out;
    public SectorShellCensusSparseTests(ITestOutputHelper o) => _out = o;

    private static RealSeed N5Seed() =>
        RealDefectiveSeeds.ForN(5).Single(s => Math.Abs(s.QStar - 0.620878) < 1e-6);

    private void DumpMethods(string tag, ShellCensusResult r)
    {
        _out.WriteLine($"--- {tag}: method mix ---");
        foreach (var g in r.Entries.GroupBy(e => e.Method).OrderBy(g => g.Key))
            _out.WriteLine($"  {g.Key}: {g.Count()}  [{string.Join(" ", g.Select(e => $"({e.P},{e.W})x{e.Shift}"))}]");
        var sum = r.Summarize();
        _out.WriteLine($"  verdict={sum.Verdict} found={sum.FoundMembers.Count} unresolved={sum.Unresolved.Count}");
    }

    // ---- (1) THE TWIN: sparse path reproduces the dense census cell for cell -----------------------
    [Fact]
    [Trait("Category", "SHELLCENSUS_SPARSE")]
    public void LoweredWall_SparsePath_ReproducesDense_AtN5()
    {
        var seed = N5Seed();

        var dense = SectorShellCensus.Run(seed, new SectorShellCensus.Options());          // wall 46000: nothing defers
        var sparse = SectorShellCensus.Run(seed,
            new SectorShellCensus.Options { MaxSectorDim = 40, SparseForDeferred = true }); // wall 40: bigger blocks go sparse
        DumpMethods("dense", dense);
        DumpMethods("sparse", sparse);

        // (a) identical classification: both find EXACTLY the containment diamond
        Assert.Equal(SectorShellCensus.ExpectedMembers(5), dense.Summarize().FoundMembers);
        Assert.Equal(SectorShellCensus.ExpectedMembers(5), sparse.Summarize().FoundMembers);

        // (b) the wall actually forced the sparse path AND removed the bare deferred rows
        Assert.Contains(sparse.Entries, e => e.Method == "sparse-witness");
        Assert.Contains(sparse.Entries, e => e.Method == "sparse-invit");
        Assert.DoesNotContain(sparse.Entries, e => e.Method == "deferred");

        // (c) per-cell agreement vs the dense twin, within the pre-registered Task-5 thresholds
        foreach (var se in sparse.Entries.Where(e => e.Method is "sparse-witness" or "sparse-invit"))
        {
            var de = dense.Entries.Single(e => e.P == se.P && e.W == se.W && e.Shift == se.Shift);
            Assert.True(de.Probed, $"dense twin of ({se.P},{se.W})x{se.Shift} not probed");
            if (se.Method == "sparse-witness")
            {
                // classification keys on the carried-parity witness bound: a deep member reading below
                // memberTol/100, not collapsed to a fabricated 0 (within a factor 2 of the dense value)
                Assert.True(se.SigmaMin < sparse.MemberTol / 100,
                    $"witness ({se.P},{se.W})x{se.Shift}: {se.SigmaMin:E3} NOT < memberTol/100 {sparse.MemberTol / 100:E3}");
                Assert.True(se.SigmaMin >= 0.5 * de.SigmaMin,
                    $"witness ({se.P},{se.W})x{se.Shift}: {se.SigmaMin:E3} NOT >= 0.5*dense {0.5 * de.SigmaMin:E3}");
            }
            else
            {
                Assert.True(se.Converged, $"sparse-invit ({se.P},{se.W})x{se.Shift} did NOT converge");
                Assert.True(Math.Abs(se.SigmaMin - de.SigmaMin) <= 0.2 * de.SigmaMin,
                    $"sparse-invit ({se.P},{se.W})x{se.Shift}: {se.SigmaMin:E3} vs dense {de.SigmaMin:E3} > 0.20");
            }
        }

        // (d) the twin verdicts: the sparse run is the witness-assisted twin of the blind dense PASS
        Assert.Equal("PASS", dense.Summarize().Verdict);
        Assert.Equal("PASS (witness-assisted)", sparse.Summarize().Verdict);
    }

    // ---- (2) NO REGRESSION: flag off keeps the bare deferred rows ----------------------------------
    [Fact]
    [Trait("Category", "SHELLCENSUS_SPARSE")]
    public void SparseForDeferredFalse_PreservesTheDeferredRows_AtN5()
    {
        var seed = N5Seed();
        var r = SectorShellCensus.Run(seed,
            new SectorShellCensus.Options { MaxSectorDim = 40, SparseForDeferred = false });
        DumpMethods("no-sparse", r);

        // the old behavior: sectors past the wall are named "deferred", probed=false, sigma NaN
        var deferred = r.Entries.Where(e => e.Method == "deferred").ToList();
        Assert.NotEmpty(deferred);
        Assert.All(deferred, e =>
        {
            Assert.False(e.Probed);
            Assert.True(double.IsNaN(e.SigmaMin));
        });
        Assert.DoesNotContain(r.Entries, e => e.Method is "sparse-witness" or "sparse-invit");
        // deferred members exist => PARTIAL, never PASS
        Assert.Equal("PARTIAL", r.Summarize().Verdict);
    }

    // ---- (3) [F1] a non-converged non-member forces PARTIAL, never a silent PASS -------------------
    [Fact]
    [Trait("Category", "SHELLCENSUS_SPARSE")]
    public void NonConvergedNonMember_ForcesPartial_Named_NeverPass_AtN5()
    {
        var seed = N5Seed();
        // InnerMaxIter=1 starves every sparse-invit inner solve: the non-member cells cannot converge,
        // while the witness members (no inner solve on the carried parity) stay soundly classified.
        var starved = SparseShiftedSigmaMin.Options.Default with { InnerMaxIter = 1 };
        var r = SectorShellCensus.Run(seed,
            new SectorShellCensus.Options { MaxSectorDim = 40, SparseForDeferred = true, Sparse = starved });
        DumpMethods("starved", r);

        var sum = r.Summarize();
        // a probed non-member that did not converge is UNRESOLVED and named, forcing PARTIAL
        Assert.NotEmpty(sum.Unresolved);
        Assert.Contains(r.Entries, e => e.Method == "sparse-invit" && e.Probed && !e.Converged);
        Assert.Equal("PARTIAL", sum.Verdict);
        Assert.NotEqual("PASS", sum.Verdict);
        Assert.NotEqual("PASS (witness-assisted)", sum.Verdict);
    }
}
