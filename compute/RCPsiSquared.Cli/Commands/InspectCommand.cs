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
/// of expansion via <c>--max-depth</c> (default 4, or 1 for <c>--root f86</c> to avoid
/// firing heavy witness computes at full breadth); an optional <c>--q-sweep</c> attaches
/// the L_eff Q-sweep view (3D / 4D structure with EVD per Q) under the root for free.
/// Pass <c>--claim &lt;ClassName&gt;</c> instead of <c>--root</c> to render a single
/// registered Claim from the typed-knowledge registry without instantiating any root
/// knowledge base.
/// </summary>
public static class InspectCommand
{
    public static int Run(string[] args)
    {
        var p = new ArgParser(args);
        p.RequireNoPositional();

        // --claim short-circuit: render a single registered Claim, no root knowledge base.
        string? claimName = p.OptionalString("claim");
        if (claimName is not null) return RunClaim(p, claimName);

        int N = p.RequireInt("N");
        string rootKind = p.OptionalString("root") ?? "fourmode";
        int maxDepth = (int)(p.OptionalDouble("max-depth") ?? DefaultDepthForRoot(rootKind));
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
            "c2cpsi" => BuildC2CpsiRoot(BuildCoherenceBlock(p, N), p),
            "c2cpsi-scan" => BuildC2CpsiScanRoot(BuildCoherenceBlock(p, N), p),
            _ => throw new ArgumentException($"unknown root: {rootKind}; known: fourmode, f86, c2hwhm, c2cpsi, c2cpsi-scan, f71, f1, f87, pi2"),
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

    /// <summary>Builds <see cref="BlockCpsiTrajectory"/> at the given <c>--Q</c> across a
    /// time grid spanning [0, <c>--t-max</c>] with <c>--t-points</c> samples (defaults
    /// t_max = 4·t_peak = 1/γ₀, t-points = 41). The trajectory tests Question B from
    /// the 2026-05-07 open-questions memory: does CΨ_block cross 1/4 under L_block(Q)
    /// evolution from the maximally-coherent initial state?</summary>
    private static IInspectable BuildC2CpsiRoot(CoherenceBlock block, ArgParser p)
    {
        double q = p.OptionalDouble("Q") ?? 1.0;
        double tMax = p.OptionalDouble("t-max") ?? (1.0 / block.GammaZero);
        int tPoints = p.OptionalDouble("t-points") is { } np ? (int)np : 41;
        var timeGrid = new double[tPoints];
        for (int i = 0; i < tPoints; i++) timeGrid[i] = tMax * i / (tPoints - 1);
        return BlockCpsiTrajectory.Build(block, q, timeGrid);
    }

    /// <summary>Builds <see cref="C2BlockCpsiQScan"/> with single-bond perturbation across
    /// a Q range. Args: <c>--bond N</c> (perturbed bond index, 0..NumBonds-1),
    /// <c>--delta D</c> (perturbation magnitude in units of γ₀, default 0.5),
    /// <c>--snap-t T</c> (snapshot time, default t_peak = 1/(4γ₀)),
    /// <c>--q-lo</c>, <c>--q-hi</c>, <c>--q-grid-points</c> (defaults 0.2, 3.0, 29).</summary>
    private static IInspectable BuildC2CpsiScanRoot(CoherenceBlock block, ArgParser p)
    {
        int bond = p.OptionalDouble("bond") is { } b ? (int)b : 0;
        double deltaUnits = p.OptionalDouble("delta") ?? 0.5;
        double delta = deltaUnits * block.GammaZero;
        double snapT = p.OptionalDouble("snap-t") ?? (1.0 / (4.0 * block.GammaZero));
        double qLo = p.OptionalDouble("q-lo") ?? 0.2;
        double qHi = p.OptionalDouble("q-hi") ?? 3.0;
        int qPoints = p.OptionalDouble("q-grid-points") is { } np ? (int)np : 29;
        var qGrid = ResonanceScan.LinearQGrid(qLo, qHi, qPoints);
        return C2BlockCpsiQScan.Build(block, bond, delta, snapT, qGrid);
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

    private static int RunClaim(ArgParser p, string claimName)
    {
        int maxDepth = (int)(p.OptionalDouble("max-depth") ?? 4);
        string? exportJson = p.OptionalString("export-json");
        var registry = KnowledgeCommand.BuildRegistry();
        var claim = registry.All().FirstOrDefault(c => c.GetType().Name == claimName);
        if (claim is null)
        {
            Console.Error.WriteLine($"error: no registered claim with type name '{claimName}'");
            return 1;
        }

        bool wroteSomething = false;
        if (exportJson is not null)
        {
            InspectionJsonExporter.WriteToFile(claim, exportJson);
            Console.Error.WriteLine($"# wrote JSON to {exportJson}");
            wroteSomething = true;
        }
        bool jsonOnly = p.HasFlag("json-only");
        if (!jsonOnly)
        {
            Console.WriteLine(ConsoleTreeRenderer.Render(claim, maxDepth));
            wroteSomething = true;
        }
        return wroteSomething ? 0 : 2;
    }

    /// <summary>Default <c>--max-depth</c> when the user did not pass one. F86's full
    /// knowledge base, expanded at depth ≥ 2, accesses Summary/Payload on
    /// <see cref="F86.UniversalShapeWitness"/> instances and triggers their
    /// <see cref="Resonance.WitnessCache"/>-backed full-block scans across the entire
    /// (c, N) anchor grid — OOM-prone on workstation memory. Default depth 1 keeps the
    /// render to the F86KB root + Tier groups only; pass <c>--max-depth N</c> for deeper.</summary>
    private static int DefaultDepthForRoot(string rootKind) =>
        rootKind == "f86" ? 1 : 4;
}
