using System.Globalization;
using System.Numerics;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Cli.Commands;

/// <summary>The G3 explorer: sweep the complex-q plane for the F89 path-3 octic's branch points (EPs),
/// lasso each from a common base, and assemble their transpositions into the monodromy = Galois group.
/// Parameterised so regions/cell/base can be swept WITHOUT an edit-rebuild cycle (same binary, new args).
///
/// usage: rcpsi gmscan [--re lo,hi] [--im lo,hi] [--cell d] [--q0 re,im]
/// example: rcpsi gmscan --re -1.8,1.8 --im -0.22,0.22 --cell 0.05 --q0 2,0</summary>
public static class GaloisMonodromyScanCommand
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    public static int Run(string[] args)
    {
        var p = new ArgParser(args);
        var (reLo, reHi) = Pair(p.OptionalString("re") ?? "-1.8,1.8");
        var (imLo, imHi) = Pair(p.OptionalString("im") ?? "-0.22,0.22");
        double cell = p.OptionalDouble("cell") ?? 0.05;
        var (q0re, q0im) = Pair(p.OptionalString("q0") ?? "2,0");

        Console.WriteLine($"# gmscan: q0={q0re.ToString("0.##", Inv)}{Sign(q0im)}i, " +
                          $"region re[{reLo.ToString("0.##", Inv)},{reHi.ToString("0.##", Inv)}] " +
                          $"im[{imLo.ToString("0.##", Inv)},{imHi.ToString("0.##", Inv)}], cell={cell.ToString("0.###", Inv)}");

        var r = GaloisMonodromyWitness.Assemble(reLo, reHi, imLo, imHi, cell, new Complex(q0re, q0im));

        Console.WriteLine($"\noctic strands at q0 (index: lambda):");
        for (int i = 0; i < r.OcticRoots.Length; i++)
            Console.WriteLine($"  [{i}] {r.OcticRoots[i].Real.ToString("0.000", Inv)}{Sign(r.OcticRoots[i].Imaginary)}i");

        Console.WriteLine($"\n{r.Eps.Count} branch points (q -> octic transposition):");
        foreach (var (q, a, b, movedAll, movedOctic) in r.Eps.OrderBy(e => e.Q.Real).ThenBy(e => e.Q.Imaginary))
        {
            string verdict = a < 0
                ? $"rejected (moved {movedAll}; octic strands moved: {{{string.Join(",", movedOctic)}}})"
                : movedAll == 2
                    ? $"swap ({string.Join(" ", movedOctic)})"
                    : $"{movedAll}-cycle ({string.Join(" ", movedOctic)})  [connects these strands]";
            Console.WriteLine($"  q={q.Real.ToString("0.000", Inv)}{Sign(q.Imaginary)}i  ->  {verdict}");
        }

        // group strands by component
        var byComp = Enumerable.Range(0, 8).GroupBy(s => r.StrandComponent[s])
            .Select(g => g.OrderBy(x => x).ToList()).OrderByDescending(g => g.Count).ToList();
        Console.WriteLine($"\ntransposition graph: {r.Edges.Count} edges {string.Join(" ", r.Edges.Select(e => $"({e.A} {e.B})"))}");
        Console.WriteLine($"components ({r.Components}), largest = {r.Largest}/8:");
        foreach (var g in byComp)
            Console.WriteLine($"  {{{string.Join(",", g)}}}{(g.Count == 1 ? "   <-- isolated strand (its EP is unfound)" : "")}");

        Console.WriteLine($"\nVERDICT: {(r.Components == 1 ? "CONNECTED -> the transpositions generate S_8 = Gal(F_8), monodromy = Galois from below"
            : $"{r.Components} components, {r.Largest}/8 strands connected. Not yet S_8: widen/retarget the EP search.")}");
        return 0;
    }

    private static (double, double) Pair(string s)
    {
        var parts = s.Split(',');
        if (parts.Length != 2) throw new ArgumentException($"expected 'a,b', got '{s}'");
        return (double.Parse(parts[0], Inv), double.Parse(parts[1], Inv));
    }

    private static string Sign(double x) => x >= 0 ? $"+{x.ToString("0.000", Inv)}" : x.ToString("0.000", Inv);
}
