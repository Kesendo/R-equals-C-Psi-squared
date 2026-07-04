using System;
using System.Globalization;
using System.IO;
using System.Linq;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Cli.Commands;

/// <summary>Drive the step-3 shell census (<see cref="SectorShellCensus"/>): at each real defective seed
/// locus of the (1,2) block, probe σ_min(L(p,w) − s) for both shifts s ∈ {λ_A, μ = −λ_A − 2N} on the
/// fundamental-domain strip, R-parity split, window-gated (the window-shell lemma excludes analytically
/// where the Bendixson window misses Re s). Verdict per seed: PASS / PARTIAL (LP64-deferred members
/// listed) / DISAGREE. CSV per seed to --out.
///
/// usage: rcpsi shellcensus --n 9 [--seed 2.137549 | --all-seeds] [--max-sector-dim 46000] [--out dir]</summary>
public static class ShellCensusCommand
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    public static int Run(string[] args)
    {
        var p = new ArgParser(args);
        int n = (int)(p.OptionalDouble("n") ?? 9);
        double? seedQ = p.OptionalDouble("seed");
        bool allSeeds = p.HasFlag("all-seeds");
        int maxSectorDim = (int)(p.OptionalDouble("max-sector-dim") ?? 46000);
        string outDir = p.OptionalString("out") ?? Path.Combine("simulations", "results", "sector_shell_census");

        var seeds = RealDefectiveSeeds.ForN(n).ToList();
        if (seeds.Count == 0)
        {
            Console.Error.WriteLine($"no recorded seeds for N={n} (registry covers N=5,7,9,11)");
            return 2;
        }
        var toRun = allSeeds ? seeds
            : seedQ is double sq ? seeds.Where(s => Math.Abs(s.QStar - sq) < 1e-6).ToList()
            : new() { seeds.First(s => s.RParity == +1) };
        if (toRun.Count == 0)
        {
            Console.Error.WriteLine($"--seed {seedQ?.ToString(Inv)} matches no recorded N={n} seed; recorded: " +
                string.Join(", ", seeds.Select(s => s.QStar.ToString("F6", Inv))));
            return 2;
        }

        Directory.CreateDirectory(outDir);
        Console.WriteLine($"# shellcensus: N={n}, {toRun.Count} seed(s), maxSectorDim={maxSectorDim}, out={outDir}");
        Console.WriteLine($"# expected members (containment corollary, FD strip): " +
            string.Join(" ", SectorShellCensus.ExpectedMembers(n).OrderBy(m => (m.Shift, m.P))
                .Select(m => $"({m.P},{m.W})x{m.Shift}")));

        bool anyDisagree = false;
        foreach (var seed in toRun)   // strictly sequential — never stack LUs
        {
            Console.WriteLine($"\n=== seed q*={seed.QStar.ToString("F6", Inv)} (R-{(seed.RParity > 0 ? "even" : "odd")}, {seed.Origin}) ===");
            var opts = new SectorShellCensus.Options { MaxSectorDim = maxSectorDim, Log = s => Console.WriteLine("  " + s) };
            var result = SectorShellCensus.Run(seed, opts);
            string csv = Path.Combine(outDir, $"shell_census_N{n}_q{seed.QStar.ToString("F6", Inv)}.csv");
            result.WriteCsv(csv);

            var sum = result.Summarize();
            Console.WriteLine($"  refined q*={result.QRefined.ToString("F9", Inv)} lambda_A={result.RefinedLambdaA.Real.ToString("F9", Inv)} " +
                              $"pairGap={result.PairGap.ToString("E2", Inv)} memberTol={sum.MemberTol.ToString("E2", Inv)}");
            Console.WriteLine($"  found members:    {Fmt(sum.FoundMembers)}");
            Console.WriteLine($"  expected(probe):  {Fmt(sum.ExpectedProbeable)}");
            if (sum.DeferredMembers.Count > 0)
                Console.WriteLine($"  DEFERRED members: {Fmt(sum.DeferredMembers)} (LP64 wall; the matrix-free follow-up's work list)");
            if (sum.Ambiguous.Count > 0)
                Console.WriteLine($"  AMBIGUOUS:        {Fmt(sum.Ambiguous)} (sigma_min in [memberTol, excludeTol])");
            Console.WriteLine($"  worst non-member sigma_min = {sum.WorstNonMemberSigma.ToString("E3", Inv)} " +
                              $"(exclude headroom x{(sum.WorstNonMemberSigma / opts.ExcludeTol).ToString("F1", Inv)})");
            Console.WriteLine($"  {sum.Verdict}  ({result.Elapsed.TotalMinutes.ToString("F1", Inv)} min, csv: {csv})");
            anyDisagree |= sum.Verdict == "DISAGREE";
        }
        return anyDisagree ? 2 : 0;
    }

    private static string Fmt(System.Collections.Generic.IEnumerable<(int P, int W, string Shift)> ms) =>
        string.Join(" ", ms.OrderBy(m => (m.Shift, m.P)).Select(m => $"({m.P},{m.W})x{m.Shift}")) is { Length: > 0 } s ? s : "(none)";
}
