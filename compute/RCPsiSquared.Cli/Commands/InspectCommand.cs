using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Decomposition.Views;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;
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
        string rootKind = p.OptionalString("root") ?? "fourmode";
        int maxDepth = (int)(p.OptionalDouble("max-depth") ?? 4);
        bool withQSweep = p.HasFlag("q-sweep");
        string? exportJson = p.OptionalString("export-json");
        bool withMeasured = p.HasFlag("with-measured");
        int? qGridPoints = p.OptionalDouble("q-grid-points") is { } v ? (int)v : null;

        IInspectable root = rootKind switch
        {
            "f71" => new F71KnowledgeBase(N),
            "f1" => BuildF1Root(p, N),
            "f87" => BuildF87Root(p, N),
            "pi2" => BuildPi2Root(p, N),
            "fourmode" => BuildFourModeRoot(BuildCoherenceBlock(p, N), withQSweep, qGridPoints),
            "f86" => BuildF86Root(BuildCoherenceBlock(p, N), withMeasured, qGridPoints),
            "c2hwhm" => C2HwhmRatio.Build(BuildCoherenceBlock(p, N), BuildOptionalQGrid(p, qGridPoints)),
            _ => throw new ArgumentException($"unknown root: {rootKind}; known: fourmode, f86, c2hwhm, f71, f1, f87, pi2"),
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

    private static F1KnowledgeBase BuildF1Root(ArgParser p, int N)
    {
        int? bondCount = p.OptionalDouble("bond-count") is { } bv ? (int)bv : null;
        int? d2 = p.OptionalDouble("degree-squared-sum") is { } dv ? (int)dv : null;
        return new F1KnowledgeBase(N, bondCount, d2);
    }

    private static F87KnowledgeBase BuildF87Root(ArgParser p, int N)
    {
        double j = p.OptionalDouble("J") ?? 1.0;
        double gamma = p.OptionalDouble("gamma") ?? 0.05;
        var chain = new ChainSystem(N, J: j, GammaZero: gamma);
        return new F87KnowledgeBase(chain);
    }

    private static Pi2KnowledgeBase BuildPi2Root(ArgParser p, int N)
    {
        double j = p.OptionalDouble("J") ?? 1.0;
        double gamma = p.OptionalDouble("gamma") ?? 0.05;
        var chain = new ChainSystem(N, J: j, GammaZero: gamma);
        return new Pi2KnowledgeBase(chain);
    }

    private static CoherenceBlock BuildCoherenceBlock(ArgParser p, int N)
    {
        int n = p.RequireInt("n");
        double gamma = p.RequireDouble("gamma");
        return new CoherenceBlock(N, n, gamma);
    }

    /// <summary>Returns an explicit Q-grid only when the caller passed at least one of
    /// <c>--q-lo</c>, <c>--q-hi</c>, <c>--q-grid-points</c>; otherwise null lets the underlying
    /// primitive use its own canonical default. At N≥9 the default upper bound 4.0 clips
    /// flanking-Interior peaks (Q_peak sticks at the edge) — pass <c>--q-hi 6.0</c> or higher.
    /// </summary>
    private static IReadOnlyList<double>? BuildOptionalQGrid(ArgParser p, int? qGridPoints)
    {
        double? lo = p.OptionalDouble("q-lo");
        double? hi = p.OptionalDouble("q-hi");
        if (lo is null && hi is null && qGridPoints is null) return null;
        double gridLo = lo ?? 0.20;
        double gridHi = hi ?? 4.00;
        int gridPoints = qGridPoints ?? 153;
        return ResonanceScan.LinearQGrid(gridLo, gridHi, gridPoints);
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
