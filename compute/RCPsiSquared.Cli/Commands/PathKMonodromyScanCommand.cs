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

        if (p.HasFlag("diabolic"))   // the N=4->N=5 forward edge: hunt the residual's diabolic points
        {
            // --residual (REQUIRED from path-5/N=6): scan the residual strands only, tracked from q0=2, so the
            // AT-locked exact degeneracies (dense same-⟨n_XY⟩ crossings) do not flood the gap field. Without it
            // the full-block scan at N>=6 returns hundreds of residual=False gap=0 AT crossings and isolates none.
            PrintDiabolics(k, reLo, reHi, imLo, imHi, cell, p.HasFlag("residual"));
            return 0;
        }

        if (p.HasFlag("delta-flip"))   // the integrability test: does XXZ Delta kill the diabolic?
        {
            PrintDeltaFlip(k, p.OptionalString("q"), p.OptionalString("lam"), p.OptionalString("deltas"));
            return 0;
        }

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

    // the integrability test: track a path-k diabolic as XXZ anisotropy Δ turns on. H(Δ)=J(XX+YY)+JΔ·ZZ;
    // the ZZ is Hermitian so the AT rate is untouched, but it breaks the additivity that makes the crossing
    // semisimple. A true (integrable) diabolic flips DEFECTIVE (geo 2->1) or LIFTS at Δ>0; a defective EP
    // (control) just drifts. Reproduces the committed N=4 table with --k 3 --q 0.658983,0 --lam -4,1.318.
    private static void PrintDeltaFlip(int k, string? qStr, string? lamStr, string? deltasStr)
    {
        int n = k + 1;
        var (qre, qim) = Pair(qStr ?? "0.6407,0.180");
        var (lre, lim) = Pair(lamStr ?? "-4.077,-1.115");
        var deltas = (deltasStr ?? "0,0.02,0.05,0.1,0.2,0.5").Split(',').Select(s => double.Parse(s, Inv)).ToArray();
        var q0 = new Complex(qre, qim); var lam0 = new Complex(lre, lim);

        Console.WriteLine($"\n# DELTA-FLIP path-{k} (N={n}): track the coalescence at q={qre.ToString("0.####", Inv)}{Sign(qim)}i, " +
                          $"lambda={lre.ToString("0.###", Inv)}{Sign(lim)}i under XXZ anisotropy Delta");
        Console.WriteLine("# H(D) = J(XX+YY) + J*D*ZZ; ZZ Hermitian => AT rate untouched; additivity broken => an integrable diabolic dies");
        Console.WriteLine("  Delta   verdict     alg geo    dep       gap        q*");
        bool diabolicAt0 = false, survivesAtPositive = false;
        foreach (var d in deltas)
        {
            var t = XxzCoherenceBlock.TrackDiabolicUnderDelta(n, q0, lam0, d);
            if (d == 0 && t.Verdict == XxzCoherenceBlock.DeltaFlipVerdict.Diabolic) diabolicAt0 = true;
            if (d > 0 && t.Survived) survivesAtPositive = true;
            Console.WriteLine($"  {d.ToString("0.###", Inv),5}  {t.Verdict,-9}  {t.Algebraic}   {t.Geometric}   " +
                              $"{t.Departure.ToString("0.0000", Inv),8}  {t.Gap.ToString("E2", Inv)}  " +
                              $"{t.QStar.Real.ToString("0.0000", Inv)}{Sign(t.QStar.Imaginary)}i");
        }
        Console.WriteLine($"\n# GATE: diabolic at Delta=0? {(diabolicAt0 ? "YES" : "NO")};  survives at Delta>0? " +
                          $"{(survivesAtPositive ? "YES -> refutes integrability-protection" : "NO -> defects/lifts => integrability-protected (DIABOLIC_BY_INTEGRABILITY's gate, off N=4)")}");
    }

    // the diabolic hunt (Q1-Q3 of the forward-edge plan): find the residual's coalescences, classify each
    // (diabolic vs defective), and read the merge Re against the AT rung-2 line (-4) vs the palindrome centre (-N).
    private static void PrintDiabolics(int k, double reLo, double reHi, double imLo, double imHi, double cell, bool residualOnly)
    {
        int nBlock = k + 1;
        Console.WriteLine($"\n# DIABOLIC HUNT path-{k} (N_block={nBlock}){(residualOnly ? " [residual-only: AT-flood excluded, tracked from q0=2]" : "")}: scan q in re[{reLo.ToString("0.##", Inv)},{reHi.ToString("0.##", Inv)}] " +
                          $"im[{imLo.ToString("0.##", Inv)},{imHi.ToString("0.##", Inv)}], cell={cell.ToString("0.###", Inv)}, q=0 mask 0.20");
        // q<->lambda coverage sanity (R-6): the full-block Re(lambda) span at the region-centre q.
        var (a, c) = PathKMonodromyScout.BuildLinear(nBlock);
        var qc = new Complex((reLo + reHi) / 2, (imLo + imHi) / 2);
        var spec = PathKMonodromyScout.AllRootsAt(a, c, qc);
        Console.WriteLine($"# coverage: at region-centre q={qc.Real.ToString("0.##", Inv)}{Sign(qc.Imaginary)}i the block Re(lambda) spans " +
                          $"[{spec.Min(z => z.Real).ToString("0.00", Inv)}, {spec.Max(z => z.Real).ToString("0.00", Inv)}] " +
                          $"(AT rung-2 line = -4; palindrome centre = -N = -{nBlock})");

        var found = PathKMonodromyScout.FindDiabolics(k, reLo, reHi, imLo, imHi, cell, residualOnly: residualOnly);
        Console.WriteLine($"\n{found.Count} coalescence(s):");
        foreach (var d in found.OrderBy(x => x.QValue.Real).ThenBy(x => x.QValue.Imaginary))
            Console.WriteLine($"  q={d.QValue.Real.ToString("0.0000", Inv)}{Sign(d.QValue.Imaginary)}i  lambda={d.MergeLambda.Real.ToString("0.000", Inv)}{Sign(d.MergeLambda.Imaginary)}i  " +
                              $"{(d.IsSemisimple ? "DIABOLIC (semisimple)" : "defective EP")}  gap={d.Gap.ToString("E2", Inv)} " +
                              $"loopId={d.LoopIsIdentity} residual={d.PairIsResidual} gap-exponent={d.GapScalingExponent.ToString("0.00", Inv)}");

        var diab = found.Where(d => d.IsSemisimple && d.PairIsResidual).ToList();
        Console.WriteLine("\n# GATES (zeros_connecting_structure forward edge):");
        if (diab.Count == 0)
        {
            Console.WriteLine("  Q1 EXISTENCE: FAIL-in-region - no semisimple residual diabolic in the scanned region.");
            Console.WriteLine("     This is the R-1/R-2 generically-EXPECTED outcome at N>=5: a diabolic is codim-3-complex");
            Console.WriteLine("     (overdetermined by 4 real conditions in a 2-DOF q-scan); the N=4 self-fold that auto-");
            Console.WriteLine("     satisfied them is gone. NOT a non-existence proof beyond this region+resolution - re-scan");
            Console.WriteLine($"     wider/finer or escalate to Route B. Defective EPs in region: {found.Count(d => !d.IsSemisimple)}.");
            return;
        }
        Console.WriteLine($"  Q1 EXISTENCE: PASS - {diab.Count} semisimple residual diabolic(s).");
        foreach (var d in diab)
        {
            double dist4 = Math.Abs(d.MergeLambda.Real + 4), distN = Math.Abs(d.MergeLambda.Real + nBlock);
            string loc = (dist4 < 1e-2 && distN < 1e-2)
                ? $"Re=-4 = -N: the N=4 TRIPLE COINCIDENCE (AT rung-2 = palindrome centre = self-fold). Expected at path-3."
                : dist4 < 1e-2
                ? "ON the AT rung-2 line Re=-4 but NOT the palindrome centre -N. SURPRISE at N>=5 (the self-fold is gone): a residual symmetry must force it - the AT alone does NOT pin a degeneracy to its midline. Next: find that symmetry."
                : distN < 1e-2 ? $"ON the palindrome centre Re=-N=-{nBlock} (rides the fold, not the AT midline)"
                : $"NEITHER -4 nor -{nBlock}: a new mechanism sets the location.";
            Console.WriteLine($"  Q2 LOCATION: q={d.QValue.Real.ToString("0.0000", Inv)}, Re(lambda)={d.MergeLambda.Real.ToString("0.000", Inv)} -> {loc}");
            Console.WriteLine($"  Q3 CHARACTER: semisimple confirmed (gap-exponent {d.GapScalingExponent.ToString("0.00", Inv)} ~1 linear; identity loop={d.LoopIsIdentity}).");
        }
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
