using System.Globalization;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Confirmations;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Decomposition.Views;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.OpenArcs;
using RCPsiSquared.Core.Resonance;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Visualization.Inspection;

namespace RCPsiSquared.Cli.Commands;

/// <summary>The terminal-side Object Manager: walks an <see cref="IInspectable"/> root and
/// emits an indented tree, JSON file, or both. Roots are selected via <c>--root</c>; depth
/// of expansion via <c>--max-depth</c> (default 4, or 1 for <c>--root f86</c> to avoid
/// firing heavy witness computes at full breadth); an optional <c>--q-sweep</c> attaches
/// the L_eff Q-sweep view (3D / 4D structure with EVD per Q) under the root for free.
/// Pass <c>--claim &lt;ClassName&gt;</c> instead of <c>--root</c> to render a single
/// registered Claim from the typed-knowledge registry without instantiating any root
/// knowledge base. Pass <c>--draw</c> to draw leaf payloads as ASCII art (vector bars,
/// matrix heatmap, curve plot) under each node instead of only their shape tag.
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

        string rootKind = p.OptionalString("root") ?? "fourmode";
        var entry = Catalog.FirstOrDefault(e => e.Name == rootKind)
            ?? throw new ArgumentException(
                $"unknown root: {rootKind}; known: {string.Join(", ", Catalog.Select(e => e.Name))}");

        // N is required only for roots that consume the global --N (fourmode, f71, mirror, …);
        // the qudit and world roots carry their own dimensions, so they run without --N.
        // If such a root is handed an explicit --N it is silently ignored; warn once on stderr
        // so a caller does not believe the flag did anything (ArgParser.HasFlag tells us the
        // flag was actually passed, as opposed to its synthesized default). A handful of
        // RequiresN:false roots DO honour an optionally-passed --N (e.g. decoder, whose witness
        // size changes with N); those are flagged HonorsOptionalN and stay silent — the value
        // was used, so claiming it was ignored would be a lie.
        if (!entry.RequiresN && !entry.HonorsOptionalN && p.HasFlag("N"))
            Console.Error.WriteLine(
                $"note: --N is not used by root '{entry.Name}' (it carries its own dimensions; ignored).");
        int N = entry.RequiresN ? p.RequireInt("N") : (p.OptionalDouble("N") is { } nv ? (int)nv : 1);
        bool withQSweep = p.HasFlag("q-sweep");
        string? exportJson = p.OptionalString("export-json");
        bool withMeasured = p.HasFlag("with-measured");
        int? qGridPoints = p.OptionalDouble("q-grid-points") is { } v ? (int)v : null;
        var ctx = new InspectRootContext(p, N, withQSweep, withMeasured, qGridPoints);

        int maxDepth = (int)(p.OptionalDouble("max-depth") ?? entry.DefaultDepth);
        IInspectable root = entry.Factory(ctx);

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
            Console.WriteLine(ConsoleTreeRenderer.Render(root, maxDepth, drawPayloads: p.HasFlag("draw")));
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

    /// <summary>The live-data root: one <see cref="MirrorSystem"/> (the conductor's stand) built
    /// from <c>--N --J --gamma --htype XY|Heisenberg --topology chain|star|ring</c>, surfacing its
    /// voices (slow modes with Δ-portfolios, the F1 palindrome, the clock's two hands) as a tree.
    /// Not a Claim from the registry; built on demand. N is capped at 7 (joint-popcount blocked
    /// spectrum). Pair with <c>--draw</c> to draw the mode portfolios as bars.</summary>
    private static IInspectable BuildMirrorRoot(ArgParser p, int N)
    {
        if (N < 1 || N > 7)
            throw new ArgumentException($"--root mirror needs N in 1..7 (joint-popcount blocked spectrum); got {N}");
        double j = p.OptionalDouble("J") ?? 1.0;
        double gamma = p.OptionalDouble("gamma") ?? 0.1;
        var htype = (p.OptionalString("htype") ?? "XY").ToLowerInvariant() switch
        {
            "heisenberg" or "heis" => HamiltonianType.Heisenberg,
            _ => HamiltonianType.XY,
        };
        var topo = (p.OptionalString("topology") ?? "chain").ToLowerInvariant() switch
        {
            "star" => TopologyKind.Star,
            "ring" => TopologyKind.Ring,
            _ => TopologyKind.Chain,
        };
        var hamiltonian = new ChainSystem(N, j, gamma, htype, topo).BuildHamiltonian();
        var channels = Enumerable.Range(0, N).Select(l => new ChannelRate($"q{l}", gamma)).ToList();
        return new MirrorSystem(N, hamiltonian, channels);
    }

    /// <summary>The post-EP flow GameObject: a single excitation evolved across a Q-grid, with
    /// per-site occupation ⟨n_site⟩(τ) relaxing to 1/N. Args: <c>--N 1..6</c>,
    /// <c>--q-list 0.5,1.0,1.5,2.5</c>, <c>--t-max 6.0</c>, <c>--t-points 60</c> (Python defaults).
    /// Pair with <c>--draw</c> to plot the trajectory curves.</summary>
    private static IInspectable BuildFlowRoot(ArgParser p, int N)
    {
        string qListStr = p.OptionalString("q-list") ?? "0.5,1.0,1.5,2.5";
        var qGrid = qListStr.Split(',')
            .Select(s => double.Parse(s.Trim(), CultureInfo.InvariantCulture)).ToArray();
        double tMax = p.OptionalDouble("t-max") ?? 6.0;
        int tPoints = p.OptionalDouble("t-points") is { } np ? (int)np : 60;
        var tauGrid = new double[tPoints];
        for (int i = 0; i < tPoints; i++) tauGrid[i] = tMax * i / (tPoints - 1);

        double[]? profile = null;
        string? profileStr = p.OptionalString("gamma-profile");
        if (profileStr is not null)
        {
            profile = profileStr.Split(',')
                .Select(s => double.Parse(s.Trim(), CultureInfo.InvariantCulture)).ToArray();
            if (p.HasFlag("fix-total"))
                profile = PostEpFlowField.NormalizeToTotal(profile, N);
        }
        return new PostEpFlowField(N, qGrid, tauGrid, profile);
    }

    /// <summary>The in-between navigator (the Object Manager telescope): sweeps a dimension's
    /// parameter and reads the marks against the in-between. Six axes:
    /// <c>--axis crossover</c> (default; sweeps the bond angle θ, frozen marks + pure rotation),
    /// <c>--axis jdefect</c> (sweeps a single bond's δJ, palindrome held but spectrum moving +
    /// eigenvector mixing; extra args <c>--defect-bond</c> default 0, <c>--delta-j-max</c> default 0.1),
    /// <c>--axis interior</c> (the ¼-to-½ interior read as a horizon: the heading θ → 0 from the
    /// interior, the Mandelbrot recursion crawling at the cusp ¼ where time stops, the slowing-is-ours
    /// seam, and the γ-invariant dwell carrying Bell+ through the fold; extra args <c>--eps-lo</c>
    /// default 1e-4, <c>--eps-hi</c> default 0.25, <c>--eps-points</c> default 13, <c>--tol</c> default
    /// 1e-12, <c>--rel-k</c> default 1e-3; N-free, uses <c>--gamma</c> only),
    /// <c>--axis spiral</c> (the interior axis in 2D: the cusp ¼ as a circle |CΨ|=¼ that every spiral
    /// crosses, the crossing angle the only free thing, the one Kingston steered; extra args
    /// <c>--omega</c> default 0.4, <c>--phi0</c> default 0, <c>--omega-points</c> default 9,
    /// <c>--tmax-factor</c> default 4; N-free, uses <c>--gamma</c> only),
    /// <c>--axis approach</c> (the family of approach shapes: the partial-entanglement start
    /// |ψ(α)⟩=cosα|00⟩+sinα|11⟩ swept over s=sin2α; CΨ(0)=s/3, crosses ¼ iff s>3/4, harmonic fraction
    /// s²/2, every member shares the carrier 4γ; extra args <c>--s-lo</c> default 0.3, <c>--s-hi</c>
    /// default 1.0, <c>--s-points</c> default 8, <c>--tmax-factor</c> default 6; N-free, uses
    /// <c>--gamma</c> only),
    /// and <c>--axis ep</c> (the exceptional point: sweeps Q across Q_EP=2/g_eff and reads the birth of the
    /// rotation, the two real decay channels coalescing defectively at −4γ₀ [the Takt pins], the Rotation
    /// angle lifting off [the F95 angle], the eigenvector overlap min(x,1/x)→1 [the defective pinch], and the
    /// IBM Kingston onset; extra args <c>--g-eff</c> default 4/3, <c>--q-lo</c> default 0.3, <c>--q-hi</c>
    /// default 4, <c>--q-points</c> default 41; N-free, uses <c>--gamma</c> only). Shared args (crossover /
    /// jdefect): <c>--gamma</c>, <c>--theta-points</c>, <c>--slow-count</c>; those two are N capped 1..6
    /// (dense Liouvillian + eigenvectors). Pair with <c>--draw</c> to plot the curves and heatmaps.</summary>
    private static IInspectable BuildBetweenRoot(ArgParser p, int N)
    {
        if (N < 1 || N > 6)
            throw new ArgumentException($"--root between needs N in 1..6 (dense Liouvillian + eigenvectors); got {N}");
        double gamma = p.OptionalDouble("gamma") ?? 0.5;
        int points = p.OptionalDouble("theta-points") is { } tp ? (int)tp : 25;
        int slowCount = p.OptionalDouble("slow-count") is { } sc ? (int)sc : 16;
        string axisName = (p.OptionalString("axis") ?? "crossover").ToLowerInvariant();

        if (axisName == "jdefect")
        {
            int defectBond = p.OptionalDouble("defect-bond") is { } db ? (int)db : 0;
            double deltaJMax = p.OptionalDouble("delta-j-max") ?? 0.1;
            return new JDefectField(N, gamma, defectBond, deltaJMax, points, slowCount);
        }

        if (axisName == "interior")
        {
            double epsLo = p.OptionalDouble("eps-lo") ?? 1e-4;
            double epsHi = p.OptionalDouble("eps-hi") ?? 0.25;
            int epsPoints = p.OptionalDouble("eps-points") is { } ep ? (int)ep : 13;
            double tol = p.OptionalDouble("tol") ?? 1e-12;
            double relK = p.OptionalDouble("rel-k") ?? 1e-3;
            return new InteriorHorizonField(epsLo, epsHi, epsPoints, tol, relK, gamma);
        }

        if (axisName == "spiral")
        {
            double omega = p.OptionalDouble("omega") ?? 0.4;
            double phi0 = p.OptionalDouble("phi0") ?? 0.0;
            int omegaPoints = p.OptionalDouble("omega-points") is { } op ? (int)op : 9;
            double tMaxFactor = p.OptionalDouble("tmax-factor") ?? 4.0;
            return new ComplexCuspSpiralField(gamma, omega, phi0, omegaPoints, tMaxFactor);
        }

        if (axisName == "approach")
        {
            double sLo = p.OptionalDouble("s-lo") ?? 0.3;
            double sHi = p.OptionalDouble("s-hi") ?? 1.0;
            int sPoints = p.OptionalDouble("s-points") is { } sp ? (int)sp : 8;
            double tMaxFactor = p.OptionalDouble("tmax-factor") ?? 6.0;
            return new ApproachFamilyField(gamma, sLo, sHi, sPoints, tMaxFactor);
        }

        if (axisName == "ep")
        {
            double gEff = p.OptionalDouble("g-eff") ?? 4.0 / 3.0;
            double qLo = p.OptionalDouble("q-lo") ?? 0.3;
            double qHi = p.OptionalDouble("q-hi") ?? 4.0;
            int qPoints = p.OptionalDouble("q-points") is { } qp ? (int)qp : 41;
            return new EpField(gamma, gEff, qLo, qHi, qPoints);
        }

        DimensionAxis axis = axisName switch
        {
            "crossover" => DimensionAxis.Crossover(N, gamma, points),
            _ => throw new ArgumentException($"unknown axis: {axisName}; known: crossover, jdefect, interior, spiral, approach, ep"),
        };
        return new DimensionField(axis, slowCount);
    }

    /// <summary>The F121 live lab: builds a <see cref="QuditPartialPalindromeWitness"/> that
    /// materialises the qudit full-Cartan dephasing dissipator at inspect time and recomputes
    /// the partial-palindrome ceiling, the product cap, and the non-product remainder from the
    /// live spectrum. Args: <c>--qudit-d</c> (local dimension, default 3), <c>--qudit-n</c>
    /// (sites, default 2), <c>--gamma</c> (default 0.05). Guarded at d^(2N) ≤ 1024 so
    /// (3,2)=81, (4,2)=256, (3,3)=729 are admitted and (4,3)=4096 is not. Pair with
    /// <c>--draw</c> to plot the live spectrum.</summary>
    private static IInspectable BuildQuditRoot(ArgParser p)
    {
        int d = p.OptionalDouble("qudit-d") is { } dv ? (int)dv : 3;
        int n = p.OptionalDouble("qudit-n") is { } nv ? (int)nv : 2;
        double gamma = p.OptionalDouble("gamma") ?? 0.05;
        return new QuditPartialPalindromeWitness(d, n, gamma);
    }

    /// <summary>The zoom-out, the Symphony: one open quantum system, one time evolution, every lens on a
    /// shared timeline, plus the cross-lens events axis. Args (honored like <c>mirror</c>):
    /// <c>--N 2..5</c> (default 3), <c>--J</c> (default 1), <c>--gamma</c> (default 0.1),
    /// <c>--htype XY|Heisenberg</c>, <c>--topology chain|star|ring</c>, <c>--initial bell|excitation|bonding</c>
    /// (default bell), <c>--t-max</c> (default 1/γ), <c>--t-points</c> (default 60). Pair with
    /// <c>--draw</c> to plot the CΨ, K, and light curves.
    ///
    /// <para>The painters' movement (PTF): pass <c>--defect-bond &lt;int&gt;</c> (and optionally
    /// <c>--delta-j</c>, default 0.02) to grow the "movement: painters" child — the piece played a
    /// second (defected) and a guard third (δJ/2) time, the per-site α_i field, the closure Σ ln α,
    /// and the live K₁ chiral mirror. Canonical to the XY chain with a real initial state; declines
    /// honestly otherwise.</para></summary>
    private static IInspectable BuildSymphonyRoot(ArgParser p, int N)
    {
        if (N < 2 || N > Symphony.MaxN)
            throw new ArgumentException($"--root symphony needs N in 2..{Symphony.MaxN} (dense d²×d² Liouvillian); got {N}");
        double j = p.OptionalDouble("J") ?? 1.0;
        double gamma = p.OptionalDouble("gamma") ?? 0.1;
        var htype = (p.OptionalString("htype") ?? "XY").ToLowerInvariant() switch
        {
            "heisenberg" or "heis" => HamiltonianType.Heisenberg,
            _ => HamiltonianType.XY,
        };
        var topo = (p.OptionalString("topology") ?? "chain").ToLowerInvariant() switch
        {
            "star" => TopologyKind.Star,
            "ring" => TopologyKind.Ring,
            _ => TopologyKind.Chain,
        };
        var initial = (p.OptionalString("initial") ?? "bell").ToLowerInvariant() switch
        {
            "excitation" or "single" or "se" => InitialStateKind.SingleExcitation,
            "bonding" or "sine" => InitialStateKind.BondingMode,
            _ => InitialStateKind.BellPair,
        };
        double tMax = p.OptionalDouble("t-max") ?? double.NaN;
        int tPoints = p.OptionalDouble("t-points") is { } np ? (int)np : 60;
        int? defectBond = p.OptionalDouble("defect-bond") is { } db ? (int)db : null;
        double deltaJ = p.OptionalDouble("delta-j") ?? 0.02;
        double? tempoRatio = p.OptionalDouble("tempo-ratio");
        return new Symphony(N, j, gamma, htype, topo, initial, tMax, tPoints, defectBond, deltaJ, carrierPair: null, tempoRatio: tempoRatio);
    }

    /// <summary>The F116 live lab: builds a <see cref="GoldenRouterWitness"/> that re-runs the soft-certifier
    /// router machinery at inspect time, watches the §7.12 ceiling close on the Z-middle cases (the 2 → 0
    /// step, the window-summed golden router), and root-finds the metallic frame ratio r(c) from the router
    /// itself, comparing it against the closed-form metallic mean. Args: <c>--router-c</c> (comma-separated
    /// weights, default 0,1,2,3). Pair with <c>--draw</c> to plot the residual ‖{W, S}‖_F over r at c = 1
    /// (dipping to zero at φ).</summary>
    private static IInspectable BuildRouterRoot(ArgParser p)
    {
        IReadOnlyList<double>? weights = null;
        if (p.OptionalString("router-c") is { } raw)
        {
            var parsed = raw.Split(',', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)
                .Select(s => double.Parse(s, System.Globalization.CultureInfo.InvariantCulture))
                .ToList();
            if (parsed.Count > 0) weights = parsed;
        }
        return new GoldenRouterWitness(weights);
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
        var registry = KnowledgeRegistryFactory.BuildDefault();
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
            Console.WriteLine(ConsoleTreeRenderer.Render(claim, maxDepth, drawPayloads: p.HasFlag("draw")));
            wroteSomething = true;
        }
        return wroteSomething ? 0 : 2;
    }

    /// <summary>The single source of truth for inspect roots: name, one-line description, the
    /// factory that builds the <see cref="IInspectable"/> from a parsed context, and the default
    /// <c>--max-depth</c> for a bare render. The <c>Run</c> switch is a lookup into this list;
    /// the world root walks it. F86 keeps default depth 1 because its full knowledge base,
    /// expanded at depth ≥ 2, triggers <see cref="Resonance.WitnessCache"/>-backed full-block
    /// scans across the entire (c, N) anchor grid — OOM-prone on workstation memory. The world
    /// root keeps default depth 2 so a bare <c>--root world</c> never enumerates a catalog node's
    /// children (and so never fires any heavy factory). Pass <c>--max-depth N</c> for deeper.
    /// </summary>
    public static readonly IReadOnlyList<RootCatalogEntry> Catalog = new RootCatalogEntry[]
    {
        new("fourmode", "the live 4-mode effective block (default root)",
            c => BuildFourModeRoot(BuildCoherenceBlock(c.Parser, c.N), c.WithQSweep, c.QGridPoints)),
        new("f86", "F86 universal-shape knowledge base (Q-peak, t-peak, EP)",
            c => BuildF86Root(BuildCoherenceBlock(c.Parser, c.N), c.WithMeasured, c.QGridPoints),
            DefaultDepth: 1),
        new("c2hwhm", "C₂ HWHM-ratio across a Q-grid",
            c => C2HwhmRatio.Build(BuildCoherenceBlock(c.Parser, c.N), BuildOptionalQGrid(c.Parser, c.QGridPoints))),
        new("c2cpsi", "block CΨ trajectory through the ¼ cusp under L_block(Q)",
            c => BuildC2CpsiRoot(BuildCoherenceBlock(c.Parser, c.N), c.Parser)),
        new("c2cpsi-scan", "block CΨ Q-scan with a single-bond perturbation",
            c => BuildC2CpsiScanRoot(BuildCoherenceBlock(c.Parser, c.N), c.Parser)),
        new("f71", "F71 chain-mirror knowledge base",
            c => new F71KnowledgeBase(c.N)),
        new("f1", "F1 palindrome knowledge base",
            c => BuildF1Root(c.Parser, c.N)),
        new("f87", "F87 trichotomy knowledge base",
            c => BuildF87Root(c.Parser, c.N)),
        new("pi2", "Π² polarity knowledge base",
            c => BuildPi2Root(c.Parser, c.N)),
        new("mirror", "the live MirrorSystem (slow modes, palindrome, the clock's two hands)",
            c => BuildMirrorRoot(c.Parser, c.N)),
        new("symphony", "the zoom-out: one system, one evolution, every lens on a shared timeline",
            c => BuildSymphonyRoot(c.Parser, c.N)),
        new("envelope", "the Envelope Theorem checked live: global CΨ peaks non-increasing (theorem) vs the carrier pair's beating rise (freedom)",
            c => new EnvelopeTheoremWitness(c.N)),
        new("clock", "the two clocks: the γ-protected band-edge ladder (N≥3) and the γ-pulled exceptional point (N=2)",
            c => new ClockHandLadderWitness(
                c.Parser.OptionalDouble("J") ?? ClockHandLadderWitness.DefaultJ,
                c.Parser.OptionalDouble("gamma") ?? ClockHandLadderWitness.DefaultGamma),
            RequiresN: false),
        new("horizon", "the coherence horizon Q*(N): where the slowest mode stops oscillating = the carbon Frost-Hückel coherent↔incoherent threshold",
            _ => new CoherenceHorizonWitness(), RequiresN: false),
        new("decoder", "reading power measured live: Fisher information vs Q per readout basis - resolution grows with the Q-factor; the exceptional point reads worst",
            c => new ReadingPowerWitness(c.Parser.HasFlag("N") ? c.N : 4),
            RequiresN: false, HonorsOptionalN: true),
        new("flow", "the post-EP single-excitation flow to 1/N",
            c => BuildFlowRoot(c.Parser, c.N)),
        new("surface", "the sterile<->birth-canal boundary computed as a live surface in gamma-profile " +
            "space (N=5): the light-freeze mechanism read through every lens (the whole surface, not the " +
            "s*=0.709 line)",
            c => new BirthCanalSurfaceWitness(
                    grid: c.Parser.OptionalDouble("grid") is { } g ? (int)g : 9),
            RequiresN: false),
        new("reduction", "the birth-canal boundary as a Liouville sector reduction: the |1-exc><vac| (0,1) " +
            "block (N-dim, validated vs the full witness at N=5, runs past it) and the {0,2} junction where it " +
            "crosses to the coherence-horizon mode at N>=6",
            c => new SectorReductionWitness(
                    c.Parser.HasFlag("N") ? c.N : 5,
                    c.Parser.OptionalDouble("gamma") ?? 0.5,
                    TopologyKind.Chain),
            RequiresN: false, HonorsOptionalN: true),
        new("survivor", "the dynamic survival probe: WHERE the longest-lived dissipative mode lives across the " +
            "three physically-grounded topologies - the interior incompleteness (C=0.5) coherence on DISPERSIVE " +
            "extended matter (chain: polyenes/spin-chains/proton-wire; ring: aromatics/light-harvesting), the " +
            "boundary hub coherence on the hub-localized central-spin STAR (NV/quantum-dot/mediator, the " +
            "counterexample); lifetime <n_XY> ~ Q^2/N^2, ring/chain -> 4 (model-independent)",
            c => new IncompletenessSurvivorWitness(
                    c.Parser.HasFlag("N") ? c.N : 6,
                    c.Parser.OptionalDouble("q") ?? 1.5),
            RequiresN: false, HonorsOptionalN: true),
        new("between", "the in-between navigator (six axes: crossover/jdefect/interior/spiral/approach/ep)",
            c => BuildBetweenRoot(c.Parser, c.N)),
        new("qudit", "F121 qudit partial palindrome, recomputed live",
            c => BuildQuditRoot(c.Parser), RequiresN: false),
        new("router", "F116 golden/metallic router, ceiling closed live",
            c => BuildRouterRoot(c.Parser), RequiresN: false),
        new("arcs", "the open-arcs ledger: started, not finished, not forgotten",
            _ => OpenArcsInspectableNode.Build(), RequiresN: false),
        new("glossary", "the house language: load-bearing terms in plain words for a stranger (start here)",
            _ => GlossaryInspectableNode.Build(), RequiresN: false),
        new("diagonal", "the one diagonal, recomputed live: the rungs k=popcount(i^j) (rate -2gk), the three " +
            "readings (D-fix/R-anti/judge = the mirror group D4 within a diagonal), the basis-S3 orbit " +
            "{Q_X,Q_Y,Q_Z}, and the L_H even-step dynamics (k-parity conserved) - the whole functioning, " +
            "symmetry + dynamics, S3 |x| D4 (typed: ThreeDephasingDiagonalsOrbitClaim)",
            c => new DiagonalWitness(c.Parser.HasFlag("N") ? c.N : 3,
                                     c.Parser.OptionalDouble("gamma") ?? 0.05),
            RequiresN: false, HonorsOptionalN: true),
        new("mirrorgroup", "the mirror group D₄'s operator algebra, live (the factorization side of " +
            "--root mirror's spectral palindrome): Π_Z = R·D, the dihedral inversion D·Π_Z·D = Π_Y " +
            "(reflections/D_PI_Z_EQUALS_PI_Y), the order-8 closure, and the §3 palindrome split " +
            "generator-by-generator WITH the −2Σγ shift the diagonal witness recenters away (typed: " +
            "MirrorGroupD4Claim)",
            c => new MirrorGroupWitness(c.Parser.HasFlag("N") ? c.N : 3,
                                        c.Parser.OptionalDouble("gamma") ?? 0.05),
            RequiresN: false, HonorsOptionalN: true),
        new("ladders", "the three-ladder hinge, live: girth(ℓ)/rung(k)/moment(j) are not three orthogonal " +
            "axes but the two factors of one F87-hardness coefficient on M = A + γQ, hinged by Q (its spectrum " +
            "is the rung k = N−2k; its action Σ Z_l⊗Z_l projects A's closed walks onto the girth moments). " +
            "P_{m,1} = m·Tr(Q·A^{m−1}) = the girth moments, at every rung (typed: ThreeLadderHingeClaim)",
            c => new LadderHingeWitness(c.Parser.HasFlag("N") ? c.N : 3),
            RequiresN: false, HonorsOptionalN: true),
        new("pascalgram", "F117 Pascal-Gram positivity, live: for a non-bipartite windowed pair every #Q class " +
            "at m* is the equal-leg-total sum of squares P_{m*,d} = (m*/d)·Σ_l⃗ Σ_k⃗ |U|² ≥ 0, so p_{m*}(γ) > 0 " +
            "for every γ>0 (hard at one γ ⟹ hard at all γ). Recomputes the five canonical branch cases " +
            "(d=1/3/5) live from H and reproduces the exact CRT coefficients (typed: WindowedConverseAllGammaClaim)",
            c => new PascalGramPositivityWitness(),
            RequiresN: false),
        new("zeroimmune", "Zero-Sector Immunity, live: a random parity-violating 2-body H gives M ≈ 0 on the " +
            "w=0 ({I,Z}^⊗N) and w=N ({X,Y}^⊗N) Pauli blocks while ‖M‖ > 0 (the non-trivial gate); the classical " +
            "extreme is immune to every 2-body coupling, the palindrome-breaking lives in 0<w<N " +
            "(typed: ZeroSectorImmunityClaim)",
            c => new ZeroSectorImmunityWitness(c.Parser.HasFlag("N") ? c.N : 3,
                                               c.Parser.OptionalDouble("gamma") ?? 0.05),
            RequiresN: false, HonorsOptionalN: true),
        new("world", "the whole Object Manager: every root, the typed claims, the hardware confirmations, the open-arcs ledger, and the glossary (try --root glossary first if the language is new)",
            BuildWorldRoot, DefaultDepth: 2, RequiresN: false),
    };

    /// <summary>The world root: one tree over the entire Object Manager. The "roots" group lists
    /// every catalog entry (lazily, so a shallow render fires no heavy factory), then the typed
    /// claim registry grouped by tier, then the hardware confirmations. Default render depth 2
    /// (see <see cref="Catalog"/>) keeps a bare <c>--root world</c> to section headers only.
    /// </summary>
    private static IInspectable BuildWorldRoot(InspectRootContext ctx)
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var rootNodes = Catalog
            .Where(e => e.Name != "world")
            .Select(e => (IInspectable)new LazyInspectableNode(
                e.Name, e.Description, () => e.Factory(ctx)))
            .ToArray();
        var rootsGroup = new InspectableNode(
            displayName: "roots",
            summary: $"{rootNodes.Length} inspect root(s)",
            children: rootNodes);

        var claimsNode = ClaimRegistryInspectableNode.Build(registry);
        var confirmationsNode = ConfirmationsInspectableNode.Build();
        var arcsNode = OpenArcsInspectableNode.Build();

        return new InspectableNode(
            displayName: "world",
            summary: $"{rootNodes.Length} roots, {registry.Count} claims, " +
                     $"{ConfirmationsRegistry.All.Count} confirmations, " +
                     $"{OpenArcsRegistry.OpenCount} open arcs",
            children: new[] { (IInspectable)rootsGroup, claimsNode, confirmationsNode, arcsNode });
    }
}
