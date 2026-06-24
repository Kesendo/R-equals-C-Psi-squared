using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using System.Text;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Visualization.Plotters;

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

        if (p.HasFlag("mirror")) PrintMirror(reLo, reHi, imLo, imHi, cell);
        if (p.HasFlag("map")) PrintHeatmap(reLo, reHi, imLo, imHi);

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

        if (p.OptionalString("png") is { } png) SavePng(reLo, reHi, imLo, imHi, r, png);
        if (p.OptionalString("lambda-png") is { } lpng) SaveLambdaPng(p.OptionalDouble("lq") ?? 1.5, lpng);
        return 0;
    }

    // the dual image: the octic spectrum in the λ-plane at a fixed q (time-killed, each mode one λ), in the
    // symphony reel palette so it lays beside the molecule spectra of THE_SHARED_SKELETON. The AT-locked roots
    // sit on the absorption rungs Re λ = −2γ, −6γ (cyan); the H_B-mixed octic (pink) spreads between them.
    private static void SaveLambdaPng(double lq, string lpng)
    {
        var all = GaloisMonodromyWitness.AllRootsAt(new Complex(lq, 0));   // 12 λ at γ = 1, J = lq
        bool IsAt(Complex z) => Math.Abs(z.Real + 2) < 1e-3 || Math.Abs(z.Real + 6) < 1e-3;
        var at = all.Where(IsAt).Select(z => (z.Real, z.Imaginary)).ToArray();
        var oct = all.Where(z => !IsAt(z)).Select(z => (z.Real, z.Imaginary)).ToArray();
        OcticSpectrumPlot.Save(oct, at, -4.0, lq,
            $"F89 path-3 octic spectrum in the lambda-plane at q = J/gamma = {lq.ToString("0.##", Inv)}   (time-killed: each mode = one lambda)",
            lpng);
        Console.WriteLine($"saved {lpng}");
    }

    // render the gap field as a cyberpunk + Matrix PNG (EPs magenta, diabolic cyan) for embedding in docs.
    private static void SavePng(double reLo, double reHi, double imLo, double imHi,
        GaloisMonodromyWitness.ScanResult r, string png)
    {
        int cols = 360;
        double step = (reHi - reLo) / cols;
        int rows = Math.Max(2, (int)((imHi - imLo) / step));
        // sample the gap at CELL CENTRES (reLo + (col+½)·step) so the heatmap aligns pixel-exactly with the
        // EP markers (which sit at the true q); the Extent below is set to exactly the sampled grid, so no
        // half-cell shift and no vertical stretch.
        var gap = GaloisMonodromyWitness.OcticGapField(reLo + step / 2, imLo + step / 2, step, cols, rows);

        var intensity = new double[rows, cols];                 // ScottPlot row 0 = top = max imag
        for (int row = 0; row < rows; row++)
            for (int col = 0; col < cols; col++)
            {
                double g = gap[col, rows - 1 - row];
                intensity[row, col] = double.IsNaN(g) ? 1.0 : Math.Exp(-g / 0.5);   // q=0 super-branch -> bright core
            }
        double reHiGrid = reLo + cols * step, imHiGrid = imLo + rows * step;

        // mark EVERY found branch point (the locus is symmetric); A>=0 vs A<0 is a Galois-assembly
        // distinction, not a property of the EP, so it must not gate the visual (it caused a missing mirror).
        var eps = r.Eps.Select(e => (e.Q.Real, e.Q.Imaginary)).ToArray();
        var diab = new List<(double, double)>();                // the roots of 3q^4 + q^2 - 1 = 0
        double q2a = (-1 + Math.Sqrt(13)) / 6, q2b = (-1 - Math.Sqrt(13)) / 6;
        if (q2a > 0) { double q = Math.Sqrt(q2a); diab.Add((q, 0)); diab.Add((-q, 0)); }
        if (q2b < 0) { double q = Math.Sqrt(-q2b); diab.Add((0, q)); diab.Add((0, -q)); }
        var diabIn = diab.Where(d => d.Item1 >= reLo && d.Item1 <= reHi && d.Item2 >= imLo && d.Item2 <= imHi).ToArray();

        GapFieldPlot.Save(intensity, reLo, reHiGrid, imLo, imHiGrid, eps, diabIn,
            "F89 path-3 octic: branch locus in q = J/gamma   (EP = magenta, diabolic = cyan)", png);
        Console.WriteLine($"saved {png}");
    }

    // the flashlight: an ASCII heatmap of the octic min-gap field. Dark (@) = near a branch point (gap -> 0),
    // space = far. One cheap parallel sweep illuminates the whole EP landscape at a glance.
    private static void PrintHeatmap(double reLo, double reHi, double imLo, double imHi)
    {
        const int cols = 100;
        double step = (reHi - reLo) / cols;
        int rows = Math.Max(2, (int)((imHi - imLo) / step));
        var gap = GaloisMonodromyWitness.OcticGapField(reLo, imLo, step, cols, rows);

        double max = 0;
        for (int ir = 0; ir < cols; ir++)
            for (int ii = 0; ii < rows; ii++)
                if (!double.IsNaN(gap[ir, ii])) max = Math.Max(max, gap[ir, ii]);

        const string ramp = "@%#*+=:-. ";                     // index 0 = smallest gap (a branch point)
        Console.WriteLine($"\noctic min-gap heatmap (@ = at a branch point, space = far); re in [{reLo.ToString("0.##", Inv)},{reHi.ToString("0.##", Inv)}], im top->bottom [{imHi.ToString("0.##", Inv)} .. {imLo.ToString("0.##", Inv)}]:");
        for (int ii = rows - 1; ii >= 0; ii--)
        {
            var sb = new StringBuilder("  ");
            for (int ir = 0; ir < cols; ir++)
            {
                double g = gap[ir, ii];
                if (double.IsNaN(g)) { sb.Append(' '); continue; }
                int idx = max > 0 ? (int)(Math.Sqrt(g / max) * (ramp.Length - 1)) : ramp.Length - 1;
                sb.Append(ramp[Math.Clamp(idx, 0, ramp.Length - 1)]);
            }
            Console.WriteLine(sb.ToString());
        }
    }

    // The gate for "EPs are a mirror pair swapping": map each branch point's collision lambda_EP and test
    // the palindrome. The (SE,DE) spectrum is mirror-symmetric about Re lambda = -4 (the centre between the
    // -2/-6 rungs). If the diabolic sits ON the axis (its own mirror, silent) and the EPs come in mirror
    // pairs about it, the branch locus IS a palindrome: the map of where observer and observed swap.
    private static void PrintMirror(double reLo, double reHi, double imLo, double imHi, double cell)
    {
        (Complex lam, double gap) Collision(Complex q)
        {
            var r = GaloisMonodromyWitness.OcticRootsAt(q);
            int bi = 0, bj = 1; double bg = double.PositiveInfinity;
            for (int i = 0; i < r.Length; i++)
                for (int j = i + 1; j < r.Length; j++)
                {
                    double g = (r[i] - r[j]).Magnitude;
                    if (g < bg) { bg = g; bi = i; bj = j; }
                }
            return ((r[bi] + r[bj]) / 2, bg);
        }
        string F(Complex z) => $"{z.Real.ToString("0.00", Inv)}{Sign(z.Imaginary)}i";

        Console.WriteLine("\n# the palindrome gate: collision lambda_EP per branch point (centre = Re lambda -4):");
        double q2a = (-1 + Math.Sqrt(13)) / 6, q2b = (-1 - Math.Sqrt(13)) / 6;
        var diab = new List<Complex>();
        if (q2a > 0) { double q = Math.Sqrt(q2a); diab.Add(new Complex(q, 0)); diab.Add(new Complex(-q, 0)); }
        if (q2b < 0) { double q = Math.Sqrt(-q2b); diab.Add(new Complex(0, q)); diab.Add(new Complex(0, -q)); }
        foreach (var d in diab)
        {
            var (lam, _) = Collision(d);
            Console.WriteLine($"  diabolic  q={F(d),-14} lambda_EP={F(lam),-14} |Re+4|={Math.Abs(lam.Real + 4).ToString("0.000", Inv)}");
        }
        var eps = GaloisMonodromyWitness.FindBranchPoints(reLo, reHi, imLo, imHi, cell);
        foreach (var ep in eps.OrderBy(e => e.Q.Real).ThenBy(e => e.Q.Imaginary))
        {
            var (lam, _) = Collision(ep.Q);
            Console.WriteLine($"  EP        q={F(ep.Q),-14} lambda_EP={F(lam),-14} |Re+4|={Math.Abs(lam.Real + 4).ToString("0.000", Inv)}");
        }
    }

    private static (double, double) Pair(string s)
    {
        var parts = s.Split(',');
        if (parts.Length != 2) throw new ArgumentException($"expected 'a,b', got '{s}'");
        return (double.Parse(parts[0], Inv), double.Parse(parts[1], Inv));
    }

    private static string Sign(double x) => x >= 0 ? $"+{x.ToString("0.000", Inv)}" : x.ToString("0.000", Inv);
}
