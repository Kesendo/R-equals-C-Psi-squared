using System;
using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The REAL SEED CENSUS: does the containment corollary's one per-N input — the existence of a
/// real defective EP on the (1,2) block at odd N — extend past the proven N = 5, 7? The instrument is the
/// PT-breaking COUNT CHANGE (<see cref="PathKMonodromyScout.FindRealDefectiveByCountChange"/>): at odd N the
/// reflection sectors are self-conjugate, so a real defective EP is where two real residual strands merge
/// and leave the axis, and the REAL-root count jumps by 2. Counting is immune to the closest-pair masking
/// that hides √-EPs from the gap-field scans at F_53/F_116 density; the validation gates below include the
/// masking victim itself (the N=7 seed q* = 1.5148, absent from today's R-even gap scan).
///
/// <para>PRE-COMMITTED CRITERIA: a SEED at N = a count-change point at real q in [0.2, 3] whose AT-aware
/// defective-anywhere reading on the raw (1,2) block is Defective. Diabolic tangencies (count changes that
/// read semisimple) are reported, not counted. Coverage is the window, as always. Run:
/// <c>dotnet test "compute/RCPsiSquared.Diagnostics.Tests" --filter "Category=SEEDCENSUS"
/// --logger "console;verbosity=detailed"</c> (N=11 runs under Category=SLOW_SEEDCENSUS).</para></summary>
public sealed class RealSeedCensusTests
{
    private readonly ITestOutputHelper _out;
    public RealSeedCensusTests(ITestOutputHelper o) => _out = o;

    private const double QLo = 0.2, QHi = 3.0, Step = 0.01;

    private List<PathKMonodromyScout.PtBreakPoint> ScanBoth(int k, out List<PathKMonodromyScout.PtBreakPoint> even,
        out List<PathKMonodromyScout.PtBreakPoint> odd)
    {
        even = PathKMonodromyScout.FindRealDefectiveByCountChange(k, rOdd: false, QLo, QHi, Step);
        odd = PathKMonodromyScout.FindRealDefectiveByCountChange(k, rOdd: true, QLo, QHi, Step);
        return even.Concat(odd).OrderBy(p => p.QStar).ToList();
    }

    private void Print(string label, List<PathKMonodromyScout.PtBreakPoint> pts)
    {
        _out.WriteLine($"{label}: {pts.Count} count-change event(s)");
        foreach (var p in pts)
            _out.WriteLine($"  q*={p.QStar:F6}  realCount {p.RealCountBelow}→{p.RealCountAbove}  " +
                           (p.IsDefective
                               ? $"DEFECTIVE  λ={p.Lambda.Real:F4}{(p.Lambda.Imaginary >= 0 ? "+" : "-")}{Math.Abs(p.Lambda.Imaginary):F4}i  gap={p.MinDefectiveGap:E1}"
                               : "semisimple (tangency/graze)"));
    }

    [Fact]
    [Trait("Category", "SEEDCENSUS")]
    public void N5_Validation_KnownSeeds_Recovered()
    {
        var all = ScanBoth(4, out var even, out var odd);
        Print("N=5 R-even", even);
        Print("N=5 R-odd", odd);

        // The two known R-even real defective loci and the known R-odd one (PROOF_CODIM1 / the census).
        Assert.Contains(even, p => p.IsDefective && Math.Abs(p.QStar - 0.620878) < 2e-3);
        Assert.Contains(even, p => p.IsDefective && Math.Abs(p.QStar - 1.077615) < 2e-3);
        Assert.Contains(odd, p => p.IsDefective && Math.Abs(p.QStar - 2.804888) < 2e-3);
    }

    [Fact]
    [Trait("Category", "SEEDCENSUS")]
    public void N7_Validation_MaskingVictim_Found()
    {
        var all = ScanBoth(6, out var even, out var odd);
        Print("N=7 R-even", even);
        Print("N=7 R-odd", odd);

        // The known N=7 seed q* = 1.5148 — the locus the R-even gap-field scan MISSES at F_53 density.
        Assert.Contains(all, p => p.IsDefective && Math.Abs(p.QStar - 1.5148) < 3e-3);
    }

    [Fact]
    [Trait("Category", "SEEDCENSUS")]
    public void N9_TheProbe_SeedExistence()
    {
        var all = ScanBoth(8, out var even, out var odd);
        Print("N=9 R-even", even);
        Print("N=9 R-odd", odd);

        var seeds = all.Where(p => p.IsDefective).ToList();
        _out.WriteLine("");
        _out.WriteLine(seeds.Count > 0
            ? $"VERDICT N=9: {seeds.Count} real defective EP(s) on (1,2) in the window ⟹ the containment " +
              "corollary's per-N seed input EXTENDS to N=9 (the 12-set diamond membership then follows with " +
              "zero new compute, PROOF_CODIM1 §7)."
            : "VERDICT N=9: NO real defective EP found in the window ⟹ the seed input does NOT extend here; " +
              "the corollary stays conditional at N=9 (window-scoped statement).");

        // Instrument sanity: at odd N the real residual population must drain as q grows (the strands leave
        // the axis pairwise), so the scan cannot come back empty. (The documented 0.4755 tangency is a GRAZE:
        // no net count change across the interval, invisible to a grid count by design — see the docstring.)
        Assert.NotEmpty(even);
    }

    [Fact]
    [Trait("Category", "SLOW_SEEDCENSUS")]
    public void N11_TheProbe_SeedExistence()
    {
        var all = ScanBoth(10, out var even, out var odd);
        Print("N=11 R-even", even);
        Print("N=11 R-odd", odd);

        var seeds = all.Where(p => p.IsDefective).ToList();
        _out.WriteLine("");
        _out.WriteLine(seeds.Count > 0
            ? $"VERDICT N=11: {seeds.Count} real defective EP(s) on (1,2) in the window ⟹ the seed input " +
              "EXTENDS to N=11."
            : "VERDICT N=11: NO real defective EP found in the window (window-scoped).");
        Assert.NotEmpty(all);   // at minimum the real-population structure must produce SOME count change
    }
}
