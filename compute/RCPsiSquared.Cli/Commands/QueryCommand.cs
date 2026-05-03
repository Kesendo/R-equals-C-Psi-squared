using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;

using RCPsiSquared.Core.Knowledge;
namespace RCPsiSquared.Cli.Commands;

/// <summary>The Object Manager query layer: ask typed questions of an
/// <see cref="F86KnowledgeBase"/> tree without writing a tree-walking script.
///
/// <para>Supported queries (via <c>--q</c>):</para>
/// <list type="bullet">
///   <item><c>tier-inventory</c> — counts of claims per tier</item>
///   <item><c>tier=&lt;name&gt;</c> — list claims at a specific tier</item>
///   <item><c>anchors</c> — unique anchor pointers across all claims</item>
///   <item><c>witnesses-at --c &lt;c&gt; --wN &lt;N&gt;</c> — Interior + Endpoint HWHM/Q witnesses at (c, N)</item>
///   <item><c>per-block-qpeak --c &lt;c&gt;</c> — Q_SCALE per-block Q_peak for chromaticity c</item>
///   <item><c>per-bond-qpeak --c &lt;c&gt; --N &lt;N&gt; --bond &lt;Endpoint|Interior&gt;</c> — fine-grid Q_peak</item>
///   <item><c>extract-gamma --j &lt;double&gt; --c &lt;c&gt;</c> — γ₀-extraction protocol from measured J*</item>
///   <item><c>open</c> — open theoretical items</item>
///   <item><c>retracted</c> — retracted claims</item>
///   <item><c>compare</c> — measured (full ResonanceScan) vs predictions (Interior + Endpoint)</item>
/// </list>
///
/// <para>Examples:</para>
/// <code>
///   rcpsi query --N 5 --n 1 --gamma 0.05 --q tier-inventory
///   rcpsi query --N 5 --n 1 --gamma 0.05 --q tier=Tier1Derived
///   rcpsi query --N 7 --n 2 --gamma 0.05 --q witnesses-at --c 3 --wN 7
///   rcpsi query --N 8 --n 1 --gamma 0.05 --q extract-gamma --j 0.09 --c 4
///   rcpsi query --N 5 --n 1 --gamma 0.05 --q compare --q-grid-points 30
/// </code>
/// </summary>
public static class QueryCommand
{
    public static int Run(string[] args)
    {
        var p = new ArgParser(args);
        p.RequireNoPositional();
        int N = p.RequireInt("N");
        int n = p.RequireInt("n");
        double gamma = p.RequireDouble("gamma");
        string query = p.OptionalString("q") ?? throw new ArgumentException("missing --q <query>");
        int? qGridPoints = p.OptionalDouble("q-grid-points") is { } v ? (int)v : null;

        var block = new CoherenceBlock(N, n, gamma);
        var kb = new F86KnowledgeBase(block);

        return ExecuteQuery(kb, query, p, qGridPoints);
    }

    private static int ExecuteQuery(F86KnowledgeBase kb, string query, ArgParser p, int? qGridPoints)
    {
        if (query == "tier-inventory")
        {
            Console.WriteLine($"# {kb.DisplayName}");
            Console.WriteLine($"# {kb.TierInventoryLine()}");
            foreach (var (tier, count) in kb.CountByTier().OrderBy(kv => kv.Key))
                Console.WriteLine($"  {tier.Label(),-30} {count}");
            return 0;
        }

        if (query.StartsWith("tier="))
        {
            string tierName = query["tier=".Length..];
            if (!Enum.TryParse<Tier>(tierName, ignoreCase: true, out var tier))
            {
                Console.Error.WriteLine($"unknown tier: {tierName}; valid: {string.Join(", ", Enum.GetNames<Tier>())}");
                return 2;
            }
            Console.WriteLine($"# {kb.DisplayName} — claims at {tier.Label()}");
            foreach (var c in kb.ClaimsAtTier(tier))
                Console.WriteLine($"  {c.DisplayName}  —  {c.Summary}");
            return 0;
        }

        if (query == "anchors")
        {
            Console.WriteLine($"# {kb.DisplayName} — anchors referenced");
            foreach (var anchor in kb.AnchorsReferenced().OrderBy(a => a))
                Console.WriteLine($"  {anchor}");
            return 0;
        }

        if (query == "witnesses-at")
        {
            int c = p.RequireInt("c");
            int wN = p.RequireInt("wN");
            var (interior, endpoint) = kb.WitnessesAt(c, wN);
            Console.WriteLine($"# witnesses at c={c} N={wN}");
            Console.WriteLine($"  Interior: {(interior is null ? "(no witness)" : interior.Summary)}");
            Console.WriteLine($"  Endpoint: {(endpoint is null ? "(no witness)" : endpoint.Summary)}");
            return 0;
        }

        if (query == "per-block-qpeak")
        {
            int c = p.RequireInt("c");
            var claim = kb.PerBlockQPeak(c);
            if (claim is null) { Console.Error.WriteLine($"no per-block Q_peak claim for c={c}"); return 2; }
            Console.WriteLine($"# {claim.DisplayName}");
            Console.WriteLine($"  {claim.Summary}");
            return 0;
        }

        if (query == "per-bond-qpeak")
        {
            int c = p.RequireInt("c");
            int qN = p.RequireInt("N");
            string bondStr = p.OptionalString("bond") ?? throw new ArgumentException("--bond required");
            if (!Enum.TryParse<BondClass>(bondStr, ignoreCase: true, out var bondClass))
            { Console.Error.WriteLine($"unknown bond: {bondStr}"); return 2; }
            var w = kb.PerBondQPeak(c, qN, bondClass);
            Console.WriteLine($"# per-bond Q_peak at c={c} N={qN} {bondClass}");
            Console.WriteLine($"  {(w is null ? "(no witness)" : w.Summary)}");
            return 0;
        }

        if (query == "extract-gamma")
        {
            double jStar = p.RequireDouble("j");
            int c = p.RequireInt("c");
            var perBlock = kb.PerBlockQPeak(c);
            if (perBlock is null) { Console.Error.WriteLine($"no per-block Q_peak for c={c}"); return 2; }
            double extracted = perBlock.ExtractGammaZero(jStar);
            Console.WriteLine($"# γ₀ extraction: J* = {jStar:G6}, c = {c}");
            Console.WriteLine($"  γ₀ = J* / Q_peak(c={c}) = {jStar:G6} / {perBlock.QPeakValue:F2} = {extracted:G6}");
            return 0;
        }

        if (query == "open")
        {
            Console.WriteLine($"# {kb.DisplayName} — open theoretical items");
            foreach (var oq in kb.OpenItems())
            {
                Console.WriteLine($"  {oq.Name}");
                Console.WriteLine($"    {oq.Description}");
                Console.WriteLine($"    approach: {oq.Approach}");
            }
            return 0;
        }

        if (query == "retracted")
        {
            Console.WriteLine($"# {kb.DisplayName} — retracted claims");
            foreach (var r in kb.Retractions())
            {
                Console.WriteLine($"  [RETRACTED] {r.Name}");
                Console.WriteLine($"    was: {r.PreviousFormula}");
                Console.WriteLine($"    refutation: {r.Refutation}");
            }
            return 0;
        }

        if (query == "compare")
        {
            double[] qGrid = qGridPoints is { } np ? ResonanceScan.LinearQGrid(0.20, 4.00, np) : ResonanceScan.DefaultQGrid();
            var curve = new ResonanceScan(kb.Block).ComputeKCurve(qGrid);
            Console.WriteLine($"# {kb.DisplayName} — measured ({qGrid.Length} Q points) vs predicted");
            foreach (var match in kb.ComparePredictions(curve))
                Console.WriteLine($"  {match.Prediction.BondClass}: {match.Description} → {(match.Within ? "within" : "OUTSIDE")} tolerance");
            return 0;
        }

        Console.Error.WriteLine($"unknown query: {query}");
        Console.Error.WriteLine("known queries: tier-inventory, tier=<name>, anchors, witnesses-at, " +
                                "per-block-qpeak, per-bond-qpeak, extract-gamma, open, retracted, compare");
        return 2;
    }

}
