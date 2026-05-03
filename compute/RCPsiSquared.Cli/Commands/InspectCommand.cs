using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Decomposition.Views;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;
using RCPsiSquared.Visualization.Inspection;

namespace RCPsiSquared.Cli.Commands;

/// <summary>The terminal-side Object Manager: walks an <see cref="IInspectable"/> root and
/// emits an indented tree, JSON file, or both. Roots are selected via <c>--root</c>; depth
/// of expansion via <c>--max-depth</c>; an optional <c>--q-sweep</c> attaches the L_eff
/// Q-sweep view (3D / 4D structure with EVD per Q) under the root for free.
/// </summary>
public static class InspectCommand
{
    public static int Run(string[] args)
    {
        var p = new ArgParser(args);
        p.RequireNoPositional();
        int N = p.RequireInt("N");
        int n = p.RequireInt("n");
        double gamma = p.RequireDouble("gamma");
        string rootKind = p.OptionalString("root") ?? "fourmode";
        int maxDepth = (int)(p.OptionalDouble("max-depth") ?? 4);
        bool withQSweep = p.HasFlag("q-sweep");
        string? exportJson = p.OptionalString("export-json");

        bool withMeasured = p.HasFlag("with-measured");
        int? qGridPoints = p.OptionalDouble("q-grid-points") is { } v ? (int)v : null;

        var block = new CoherenceBlock(N, n, gamma);
        IInspectable root = rootKind switch
        {
            "fourmode" => BuildFourModeRoot(block, withQSweep, qGridPoints),
            "f86" => BuildF86Root(block, withMeasured, qGridPoints),
            _ => throw new ArgumentException($"unknown root: {rootKind}; known: fourmode, f86"),
        };

        bool wroteSomething = false;

        if (exportJson is not null)
        {
            InspectionJsonExporter.WriteToFile(root, exportJson);
            Console.Error.WriteLine($"# wrote JSON to {exportJson}");
            wroteSomething = true;
        }

        // Print tree to stdout if no export, OR if both export and tree are wanted (export
        // doesn't suppress the tree by default; --json-only suppresses).
        bool jsonOnly = p.HasFlag("json-only");
        if (!jsonOnly)
        {
            Console.WriteLine(ConsoleTreeRenderer.Render(root, maxDepth));
            wroteSomething = true;
        }

        return wroteSomething ? 0 : 2;
    }

    private static IInspectable BuildF86Root(CoherenceBlock block, bool withMeasured, int? qGridPoints)
    {
        var kb = new F86KnowledgeBase(block);
        if (!withMeasured) return kb;

        double[] qGrid = qGridPoints is { } np ? ResonanceScan.LinearQGrid(0.20, 4.00, np) : ResonanceScan.DefaultQGrid();
        var measuredCurve = new ResonanceScan(block).ComputeKCurve(qGrid);
        var matches = kb.CompareTo(measuredCurve);
        return new InspectableNode(
            displayName: kb.DisplayName + " + measured comparison",
            summary: $"{kb.Summary}; scan over {qGrid.Length} Q points",
            children: new IInspectable[]
            {
                kb,
                InspectableNode.Group("measured vs predicted",
                    matches.Cast<IInspectable>().ToArray()),
            });
    }

    private static IInspectable BuildFourModeRoot(CoherenceBlock block, bool withQSweep, int? qGridPoints)
    {
        if (block.C < 2)
            throw new InvalidOperationException(
                $"--root fourmode requires chromaticity ≥ 2; got c={block.C} for N={block.N}, n={block.LowerPopcount}.");
        var eff = FourModeEffective.Build(block);
        if (!withQSweep) return eff;

        double[] qGrid = qGridPoints is { } np ? ResonanceScan.LinearQGrid(0.20, 4.00, np) : ResonanceScan.DefaultQGrid();
        var sweep = new LEffSweepView(eff, qGrid);
        return new InspectableNode(
            displayName: $"FourModeEffective+QSweep (c={block.C}, N={block.N}, n={block.LowerPopcount}, γ₀={block.GammaZero:G3})",
            summary: $"4-mode block + L_eff sweep over {qGrid.Length} Q points",
            children: new IInspectable[] { eff, sweep });
    }

}
