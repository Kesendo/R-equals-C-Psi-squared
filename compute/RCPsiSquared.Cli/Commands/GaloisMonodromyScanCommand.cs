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

        if (p.HasFlag("zeros")) PrintZeros(r);

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
        if (p.OptionalString("graph-png") is { } gpng) SaveGraphPng(gpng);
        return 0;
    }

    // mine the "zeros": the σ_T-fixed strands (Re λ = −4, on the fold), with full-precision frequencies and
    // their braid-graph connectivity. Tom's principle: if the zeros differ, a structure must connect them.
    private static void PrintZeros(GaloisMonodromyWitness.ScanResult r)
    {
        var (sigmaT, _, _, _, _) = GaloisMonodromyWitness.PalindromeStrandPairing();
        var nbr = new Dictionary<int, List<int>>();
        for (int i = 0; i < 8; i++) nbr[i] = new List<int>();
        foreach (var (a, b) in r.Edges)
        {
            if (!nbr[a].Contains(b)) nbr[a].Add(b);
            if (!nbr[b].Contains(a)) nbr[b].Add(a);
        }

        Console.WriteLine("\n# the zeros (σ_T-fixed, on the fold Re λ=−4) and the structure between them:");
        for (int i = 0; i < 8; i++)
        {
            var lam = r.OcticRoots[i];
            string cls = sigmaT[i] == i ? "ZERO (fold)" : $"twin<->{sigmaT[i]}";
            Console.WriteLine($"  [{i}] lambda={lam.Real.ToString("F6", Inv)}{(lam.Imaginary >= 0 ? "+" : "")}{lam.Imaginary.ToString("F6", Inv)}i  " +
                              $"{cls,-12} braids {{{string.Join(",", nbr[i].OrderBy(x => x))}}} (deg {nbr[i].Count})");
        }

        var zeros = Enumerable.Range(0, 8).Where(i => sigmaT[i] == i).ToList();
        var ims = zeros.Select(i => r.OcticRoots[i].Imaginary).OrderBy(x => x).ToList();
        var twins = Enumerable.Range(0, 8).Where(i => sigmaT[i] != i).Select(i => r.OcticRoots[i].Imaginary).OrderBy(x => x).ToList();
        Console.WriteLine($"\n  zero frequencies (Im lambda, sorted): {string.Join(", ", ims.Select(x => x.ToString("F4", Inv)))}");
        Console.WriteLine($"  sum(zero Im) = {ims.Sum().ToString("F4", Inv)};  consecutive gaps: {string.Join(", ", ims.Zip(ims.Skip(1), (a, b) => (b - a).ToString("F4", Inv)))}");
        Console.WriteLine($"  twin frequencies (Im lambda): {string.Join(", ", twins.Select(x => x.ToString("F4", Inv)))}");
        Console.WriteLine($"  sum(all 8 Im) = {(ims.Sum() + twins.Sum()).ToString("F4", Inv)}");

        // THE PATH from one zero to the next: the shortest route through the braid graph. [n]=a zero
        // (self-mirror), (n)=a twin (a +/- mode) crossed on the way. The connection between the zeros, not
        // the zeros themselves: the zeros are linked THROUGH the twins, over the EP transitions.
        Console.WriteLine("\n  the path from one zero to the next (shortest braid route; [z]=zero, (t)=twin crossed):");
        for (int x = 0; x < zeros.Count; x++)
            for (int y = x + 1; y < zeros.Count; y++)
            {
                var path = Bfs(nbr, zeros[x], zeros[y]);
                if (path.Count == 0) { Console.WriteLine($"    {zeros[x]} -> {zeros[y]}: (disconnected)"); continue; }
                string render = string.Join(" - ", path.Select(s => sigmaT[s] == s ? $"[{s}]" : $"({s})"));
                int twinsCrossed = path.Count(s => sigmaT[s] != s);
                Console.WriteLine($"    {zeros[x]} -> {zeros[y]}: {render}   ({twinsCrossed} twin(s) crossed, {path.Count - 1} hops)");
            }
    }

    // shortest path src->dst through the undirected braid graph (empty list if disconnected).
    private static System.Collections.Generic.List<int> Bfs(
        System.Collections.Generic.Dictionary<int, System.Collections.Generic.List<int>> adj, int src, int dst)
    {
        var prev = new System.Collections.Generic.Dictionary<int, int> { { src, -1 } };
        var queue = new System.Collections.Generic.Queue<int>();
        queue.Enqueue(src);
        while (queue.Count > 0)
        {
            int u = queue.Dequeue();
            if (u == dst) break;
            foreach (var v in adj[u])
                if (!prev.ContainsKey(v)) { prev[v] = u; queue.Enqueue(v); }
        }
        var path = new System.Collections.Generic.List<int>();
        if (!prev.ContainsKey(dst)) return path;
        for (int c = dst; c != -1; c = prev[c]) path.Add(c);
        path.Reverse();
        return path;
    }

    // the transposition graph with the σ_T fold-mirror: lay the 8 strands so σ_T is the left-right reflection
    // (fold-fixed on the axis x=0, mirror-twins at ±x, height = Im λ, which twins share), draw the braid edges,
    // and witness the non-invariance with one braid edge + its σ_T-image ghost (NOT a braid edge). The C-T
    // obstruction (σ_T non-central) made visible. Canonical base q0=2 (matches PalindromeStrandPairing's σ_T).
    private static void SaveGraphPng(string gpng)
    {
        var r = GaloisMonodromyWitness.Assemble(-1.80, 1.80, -0.22, 0.22, 0.05, new Complex(2, 0));
        var (sigmaT, _, _, _, _) = GaloisMonodromyWitness.PalindromeStrandPairing();
        int Mn(int x, int y) => System.Math.Min(x, y);
        int Mx(int x, int y) => System.Math.Max(x, y);

        // the witness: the first braid edge whose σ_T-image is NOT a braid edge (σ_T breaks the braiding).
        var edgeSet = new HashSet<(int, int)>(r.Edges.Select(e => (Mn(e.A, e.B), Mx(e.A, e.B))));
        (int a, int b)? witness = null;
        foreach (var (a, b) in r.Edges)
            if (!edgeSet.Contains((Mn(sigmaT[a], sigmaT[b]), Mx(sigmaT[a], sigmaT[b])))) { witness = (a, b); break; }
        if (witness is not { } w) { Console.WriteLine("no σ_T-broken braid edge found (unexpected)"); return; }

        // only the strands of this one witness: the braid (a,c) and its σ_T-image (σa, σc).
        int a0 = w.a, c0 = w.b, a1 = sigmaT[a0], c1 = sigmaT[c0];
        var involved = new List<int> { a0, c0, a1, c1 }.Distinct().ToList();

        // y rows: σ_T-fixed strands on top (apex), twin pairs below; a twin pair shares one row (same Im λ).
        var groups = involved.GroupBy(s => Mn(s, sigmaT[s]))
            .OrderBy(g => g.Any(s => sigmaT[s] == s) ? 0 : 1).ToList();
        var posY = new Dictionary<int, double>();
        double yy = groups.Count;
        foreach (var g in groups) { foreach (var s in g) posY[s] = yy; yy -= 1.0; }

        int n = involved.Count;
        var local = new Dictionary<int, int>();
        var xs = new double[n]; var ys = new double[n]; var labels = new int[n]; var onFold = new bool[n];
        for (int i = 0; i < n; i++)
        {
            int s = involved[i];
            local[s] = i;
            onFold[i] = sigmaT[s] == s;
            xs[i] = onFold[i] ? 0.0 : (s < sigmaT[s] ? -1.0 : 1.0);   // twin goes left if lower-indexed
            ys[i] = posY[s];
            labels[i] = s;
        }

        TranspositionGraphPlot.SaveWitness(
            xs, ys, labels, onFold,
            (local[a0], local[c0]), (local[a1], local[c1]),
            "F89 path-3: the Re λ = −4 fold (σ_T) is not a symmetry of the braiding\n" +
            "a real braid edge reflects across the fold to a non-braid  ⟹  σ_T is non-central",
            gpng);
        Console.WriteLine($"saved {gpng}");
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
            "F89 path-3 octic: branch locus in q = J/gamma   (EP = magenta, diabolic = gold)", png);
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
