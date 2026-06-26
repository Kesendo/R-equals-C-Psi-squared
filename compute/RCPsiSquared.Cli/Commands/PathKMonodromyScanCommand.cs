using System;
using System.Globalization;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Cli.Commands;

/// <summary>Drive the path-k monodromy scout (<see cref="PathKMonodromyScout"/>): generalise the path-3
/// octic Galois = monodromy (S_8) off path-3. Build the path-k (SE,DE) block, remove the AT-locked factor,
/// lasso the residual F_d's branch points from a common base, and assemble the transposition graph (connected
/// ⟺ S_d). Validates against the known S_8 (k=3) / S_18 (k=4); maps the uncomputed path-4 branch geometry.
///
/// usage: rcpsi pkmono [--k 3] [--re lo,hi] [--im lo,hi] [--cell d] [--q0 2,0]
/// note: AT identification is built at q0=2, so keep --q0 2,0 (the canonical base).</summary>
public static class PathKMonodromyScanCommand
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    public static int Run(string[] args)
    {
        var p = new ArgParser(args);
        int k = (int)(p.OptionalDouble("k") ?? 3);
        var (reLo, reHi) = Pair(p.OptionalString("re") ?? "-2,2");
        var (imLo, imHi) = Pair(p.OptionalString("im") ?? "-0.3,0.3");
        double cell = p.OptionalDouble("cell") ?? 0.05;
        var (q0re, q0im) = Pair(p.OptionalString("q0") ?? "2,0");

        Console.WriteLine($"# pkmono: path-{k} monodromy = Galois (generalises the path-3 octic S_8 off path-3)");
        Console.WriteLine($"# q0={q0re.ToString("0.##", Inv)}{Sign(q0im)}i, region re[{reLo.ToString("0.##", Inv)},{reHi.ToString("0.##", Inv)}] " +
                          $"im[{imLo.ToString("0.##", Inv)},{imHi.ToString("0.##", Inv)}], cell={cell.ToString("0.###", Inv)}");

        var r = PathKMonodromyScout.Scan(k, reLo, reHi, imLo, imHi, cell, new Complex(q0re, q0im));

        Console.WriteLine($"# path-{r.K} (N_block={r.NBlock}): block dim {r.BlockDim} = AT {r.AtDim} + residual F_{r.ResidualDim}");
        Console.WriteLine($"\nresidual roots at q0 (index: lambda):");
        for (int i = 0; i < r.ResidualRoots.Length; i++)
            Console.WriteLine($"  [{i,2}] {r.ResidualRoots[i].Real.ToString("0.000", Inv)}{Sign(r.ResidualRoots[i].Imaginary)}i");

        Console.WriteLine($"\n{r.Eps.Count} branch points (q -> residual transposition):");
        foreach (var (q, a, b, movedRes) in r.Eps.OrderBy(e => e.Q.Real).ThenBy(e => e.Q.Imaginary))
        {
            string verdict = a < 0
                ? $"rejected (residual strands moved: {{{string.Join(",", movedRes)}}})"
                : movedRes.Length == 2
                    ? $"swap ({string.Join(" ", movedRes)})"
                    : $"{movedRes.Length}-cycle ({string.Join(" ", movedRes)})  [connects these]";
            Console.WriteLine($"  q={q.Real.ToString("0.000", Inv)}{Sign(q.Imaginary)}i  ->  {verdict}");
        }

        Console.WriteLine($"\ntransposition graph: {r.Edges.Count} edges {string.Join(" ", r.Edges.Select(e => $"({e.A} {e.B})"))}");
        var byComp = Enumerable.Range(0, r.ResidualDim).GroupBy(s => r.StrandComponent[s])
            .Select(g => g.OrderBy(x => x).ToList()).OrderByDescending(g => g.Count).ToList();
        Console.WriteLine($"components ({r.Components}), largest = {r.Largest}/{r.ResidualDim}:");
        foreach (var g in byComp)
            Console.WriteLine($"  {{{string.Join(",", g)}}}{(g.Count == 1 ? "   <-- isolated (its EP is unfound)" : "")}");

        Console.WriteLine($"\nVERDICT: {(r.Components == 1
            ? $"CONNECTED -> the transpositions generate S_{r.ResidualDim} = Gal(F_{r.ResidualDim}), monodromy = Galois from below"
            : $"{r.Components} components, {r.Largest}/{r.ResidualDim} connected. Not yet S_{r.ResidualDim}: widen/retarget the EP search.")}");

        PrintFoldStructure(r);
        return 0;
    }

    // the σ_T fold structure on the residual strands (zeros / within-block twins / cross-block partners) and,
    // when within-block zeros exist (path-3 only), the braid-graph road between them through the twins.
    private static void PrintFoldStructure(PathKMonodromyScout.ScanResult r)
    {
        var fp = r.FoldPartner;
        var zeros = Enumerable.Range(0, r.ResidualDim).Where(i => fp[i] == i).ToList();
        var twins = Enumerable.Range(0, r.ResidualDim).Where(i => fp[i] >= 0 && fp[i] != i).ToList();
        int cross = Enumerable.Range(0, r.ResidualDim).Count(i => fp[i] < 0);
        var twinPairs = twins.Where(i => i < fp[i]).Select(i => $"({i}<->{fp[i]})");

        Console.WriteLine($"\n# the global fold lambda -> -conj(lambda) - 2N (sigma = N = {r.NBlock}): the sigma_T structure on the residual:");
        Console.WriteLine($"  zeros (on-fold, self-mirror Re lambda = -{r.NBlock}): {zeros.Count}  {{{string.Join(",", zeros)}}}");
        Console.WriteLine($"  within-block twins (+/- modes): {twins.Count / 2} pairs  {string.Join(" ", twinPairs)}");
        Console.WriteLine($"  CROSS-block (mirror partner lives in (SE,w_N-2)): {cross}");

        if (zeros.Count == 0)
        {
            Console.WriteLine("  => 0 within-block zeros: the N=4 'zeros + road' structure is N=4-only; at N>=5 the");
            Console.WriteLine("     self-mirror partners are CROSS-block, so the connection between the zeros IS the");
            Console.WriteLine("     cross-block fold (SE,DE)<->(SE,w_N-2) (foldcross), not an intra-block braid road.");
            return;
        }

        // road: shortest braid route between each zero pair; [z]=zero, (t)=twin crossed (the path-3 picture).
        var nbr = new Dictionary<int, List<int>>();
        for (int i = 0; i < r.ResidualDim; i++) nbr[i] = new List<int>();
        foreach (var (a, b) in r.Edges)
        {
            if (!nbr[a].Contains(b)) nbr[a].Add(b);
            if (!nbr[b].Contains(a)) nbr[b].Add(a);
        }
        Console.WriteLine("  the path from one zero to the next (shortest braid route; [z]=zero, (t)=twin crossed):");
        for (int x = 0; x < zeros.Count; x++)
            for (int y = x + 1; y < zeros.Count; y++)
            {
                var path = Bfs(nbr, zeros[x], zeros[y]);
                if (path.Count == 0) { Console.WriteLine($"    {zeros[x]} -> {zeros[y]}: (disconnected)"); continue; }
                string render = string.Join(" - ", path.Select(s => fp[s] == s ? $"[{s}]" : $"({s})"));
                int twinsCrossed = path.Count(s => fp[s] != s);
                Console.WriteLine($"    {zeros[x]} -> {zeros[y]}: {render}   ({twinsCrossed} twin(s) crossed, {path.Count - 1} hops)");
            }
    }

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

    private static (double, double) Pair(string s)
    {
        var parts = s.Split(',');
        if (parts.Length != 2) throw new ArgumentException($"expected 'a,b', got '{s}'");
        return (double.Parse(parts[0], Inv), double.Parse(parts[1], Inv));
    }

    private static string Sign(double x) => x >= 0 ? $"+{x.ToString("0.000", Inv)}" : x.ToString("0.000", Inv);
}
