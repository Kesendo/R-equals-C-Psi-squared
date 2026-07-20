using System.Globalization;
using System.IO;
using System.Text;
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
using RCPsiSquared.Core.Knowledge;
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

    /// <summary>The self-mirror fixed point AS AN OBJECT: build the system (the x/y/z frame) exactly as
    /// --root mirror does, then return the <see cref="SelfMirrorObject"/> that INHERITS it. The jump:
    /// System (bears x/y/z) → Object (inherits the frame, is the F1 fixed point Re λ = −σ).</summary>
    private static IInspectable BuildSelfMirrorRoot(ArgParser p, int N)
    {
        if (N < 1 || N > 7)
            throw new ArgumentException($"--root selfmirror needs N in 1..7 (joint-popcount blocked spectrum); got {N}");
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
        return new SelfMirrorObject(new MirrorSystem(N, hamiltonian, channels));   // the object inherits the system
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
    /// and <c>--axis ep</c> (the toy 2×2 rate-channel exceptional point: sweeps Q across Q_EP=2/g_eff and
    /// reads the birth of the rotation in the 2-level reduction, the toy's two real decay channels coalescing
    /// defectively at −4γ₀ [the Takt pins], the Rotation angle lifting off [the F95 angle], the toy eigenvector
    /// overlap min(x,1/x)→1 [the toy EP pinch], and the IBM Kingston single-excitation-walk overdamped→revival
    /// handover at Q≈1.5. The defectiveness is the toy reduction's; the physical (n,n+1) chain block is
    /// non-normal near Q_peak; its own real-axis defective EPs are F89's scattered seeds, not this clean
    /// Q_EP pinch [F86a-retraction corrected 2026-07-07, see LocalGlobalEpLink / PROOF_F86A section The real-axis EP]. Extra args <c>--g-eff</c> default 4/3, <c>--q-lo</c> default
    /// 0.3, <c>--q-hi</c> default 4, <c>--q-points</c> default 41; N-free, uses <c>--gamma</c> only). Shared args (crossover /
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

    /// <summary>The eigenVECTOR holonomy around the (1,2)-block defective seed (i⁴=1 memory loop): encircle
    /// the EP q* in the complex plane and recompute the frame-monodromy of the two coalescing eigenvectors
    /// in the biorthogonal vᵀv gauge. Arg: <c>--N</c> (default 5; the full 50-dim block tracks cleanly).
    /// <c>--N 9</c> = the living seed studied 2026-07-07 (the full block leaks on odd loops — the per-loop
    /// span residual makes it visible). Pair with <c>--draw</c> to plot the span residual per loop.</summary>
    private static IInspectable BuildHolonomyRoot(ArgParser p)
    {
        int n = p.OptionalDouble("N") is { } nv ? (int)nv : 5;
        return new SeedHolonomyWitness(n);
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
    /// honestly otherwise.</para>
    /// <para>The seam movement (movement 4, the calibration topology): pass <c>--calibrate</c> (and
    /// optionally <c>--lab-unit &lt;name&gt;</c>, default "model-unit") to grow "movement: the seam" —
    /// the γ-anchor (γ₀ = gap/2) and J-anchor (J from the coherence hand) recover (γ₀, J), and the
    /// over-determination gate certifies they reconcile through Q (firing outside the XY, Q ≥ Q*(N)
    /// domain). Spectrum-only; the converse of the clock movement (--tempo-ratio).</para></summary>
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
        bool calibrate = p.HasFlag("calibrate");
        string? labUnit = p.OptionalString("lab-unit");
        var sym = new Symphony(N, j, gamma, htype, topo, initial, tMax, tPoints, defectBond, deltaJ,
            carrierPair: null, tempoRatio: tempoRatio, calibrate: calibrate, labUnit: labUnit);
        if (p.HasFlag("export"))
            ExportSymphonyCsv(sym, p.OptionalString("export-name"));
        return sym;
    }

    /// <summary>Writes the Symphony's already-computed reel to CSV for an external plotter (--export):
    /// the Liouvillian spectrum (every mode's whole life as one point λ) and the film curves (global CΨ,
    /// local CΨ, light over the shared timeline). Lets simulations/reel_and_projector.py draw the
    /// "with t-axis" film and the "without t-axis" time-erased spectrum from the LIVE LAB, not the
    /// deprecated Python framework. I/O lives here in the Cli layer; the Symphony witness stays pure.</summary>
    private static void ExportSymphonyCsv(Symphony s, string? exportName = null)
    {
        string baseDir = Path.Combine(RepoRootLocator.Require(), "simulations", "results", "symphony_reel");
        string dir = exportName is { Length: > 0 } ? Path.Combine(baseDir, exportName) : baseDir;
        Directory.CreateDirectory(dir);
        var inv = CultureInfo.InvariantCulture;
        double sigma = s.N * s.Gamma;

        var ev = new StringBuilder();
        ev.AppendLine($"# N={s.N} J={s.J.ToString("R", inv)} gamma={s.Gamma.ToString("R", inv)} "
            + $"Q={(s.J / s.Gamma).ToString("R", inv)} sigma={sigma.ToString("R", inv)} "
            + $"center={(-sigma).ToString("R", inv)}");
        ev.AppendLine("Re,Im");
        foreach (var lam in s.LiouvillianEigenvalues)
            ev.AppendLine($"{lam.Real.ToString("R", inv)},{lam.Imaginary.ToString("R", inv)}");
        File.WriteAllText(Path.Combine(dir, "symphony_eigenvalues.csv"), ev.ToString());

        var cu = new StringBuilder();
        cu.AppendLine("t,K,global_CPsi,local_CPsi,light");
        var t = s.TimeGrid;
        var st = s.States;
        for (int i = 0; i < t.Count; i++)
            cu.AppendLine(
                $"{t[i].ToString("R", inv)},{(s.Gamma * t[i]).ToString("R", inv)},"
                + $"{Symphony.Cpsi(st[i]).ToString("R", inv)},{s.LocalCpsi(st[i]).ToString("R", inv)},"
                + $"{Symphony.LightContent(st[i]).ToString("R", inv)}");
        File.WriteAllText(Path.Combine(dir, "symphony_curves.csv"), cu.ToString());

        System.Console.WriteLine($"[export] symphony reel written to {dir}");
        System.Console.WriteLine($"[export]   symphony_eigenvalues.csv  ({s.LiouvillianEigenvalues.Count} modes)");
        System.Console.WriteLine($"[export]   symphony_curves.csv       ({t.Count} time points)");
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
        new("time-exclusion", "TIME_IRREVERSIBILITY_EXCLUSION, the typed argument behind F49's cross-term value: the anticommutator {L_H, L_Dc} of the Hamiltonian and F1-centered dephasing Liouvillians vanishes ONLY at N=2 (Frobenius-orthogonality, the Pythagorean split L_c²=L_H²+L_Dc²), growing for N>2 as R(N)=√((N−2)/(N·4^(N−1))); but the COMMUTATOR [L_H, L_Dc] ≠ 0 already at N=2, so the vanishing is orthogonality, NOT a separability/reversibility criterion — the naive arrow-of-time reading is EXCLUDED. Builds both norms across N via LindbladianBuilder and cross-checks the anticommutator against the F49 closed form fed a live per-bond norm (two computations meeting). --N sweeps up to that N (≤5)",
            c => new TimeIrreversibilityExclusionWitness(
                c.Parser.HasFlag("N") ? c.N : 4,
                c.Parser.OptionalDouble("gamma") ?? 0.05,
                c.Parser.OptionalDouble("J") ?? 1.0),
            RequiresN: false, HonorsOptionalN: true),
        new("quarter-uniqueness", "UNIQUENESS_PROOF, the typed argument behind the ¼ value: α=2 is the UNIQUE Rényi order whose fold threshold CΨ*_α=(α−1)^(α−1)/(α^α·Ψ^(α−2)) is state-independent (the Ψ^(α−2) factor vanishes only there), where it equals ¼; and the α=2 fixed-point discriminant D=1−4CΨ has its single zero at ¼, the unique bifurcation boundary. Sweeps α at two probe states (spread=0 only at α=2) + the discriminant across CΨ; elementary arithmetic, exact. --psi-low / --psi-high set the probe states",
            c => new QuarterBoundaryUniquenessWitness(
                c.Parser.OptionalDouble("psi-low") ?? 0.3,
                c.Parser.OptionalDouble("psi-high") ?? 0.7),
            RequiresN: false),
        new("noise-origin", "INCOMPLETENESS_PROOF, the noise-origin 5-candidate elimination: dephasing noise cannot originate WITHIN the d(d−2)=0 ontology, so it comes from OUTSIDE (the incompleteness, the V-Effect one dimension up). Candidate 5 (the dimension algebra d²−2d=0 ⟹ d∈{0,2}, d=1 and d≥3 excluded) is recomputed live; candidate 1 is the typed F63 [Π²,L]=0 constraint; candidates 2-4 (single-qubit decay γ_eff=0, the bath's infinite regress, d=0's property-lessness) are surfaced from the proof, the heavier process-tomography compute deferred. --max-d sets the dimension sweep",
            c => new NoiseOriginExclusionWitness(
                c.Parser.OptionalDouble("max-d") is { } m ? (int)m : NoiseOriginExclusionWitness.DefaultMaxD),
            RequiresN: false),
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
        new("bandedge", "the topology-general band edge = J × adjacency spectral radius; the law + the gap-dominance map (chain/star/ring)",
            c => new TopologyBandEdgeWitness(
                c.Parser.OptionalDouble("J") ?? TopologyBandEdgeWitness.DefaultJ,
                c.Parser.OptionalDouble("gamma") ?? TopologyBandEdgeWitness.DefaultGamma),
            RequiresN: false),
        new("ceiling", "the structural ceiling closed forms: g2(K_N)=4/N, g2(star_N)=4/(N−1), the N=4 (2,2) outlier 2−2/√3, from the commutant rep structure",
            _ => new StructuralCeilingWitness(), RequiresN: false),
        new("protection", "the qudit mirror-protection scaling law: the per-site product mirror protects (2d)^N / d^{2N} = (2/d)^N of the coherence space, decaying exponentially in d; = 1 ⟺ d=2 (the qubit is the unique full-mirror dimension). Rates stay d-independent.",
            _ => new QuditMirrorProtectionWitness(), RequiresN: false),
        new("horizon", "the coherence horizon Q*(N): where the slowest mode stops oscillating = the carbon Frost-Hückel coherent↔incoherent threshold",
            _ => new CoherenceHorizonWitness(), RequiresN: false),
        new("epcharacter", "the artifact-free EP-character diagnostic (Riesz ‖P‖ / departure-from-normality / geo-vs-alg): the non-eig sibling of PhaseRigidity that confirms the coherence-horizon √-EP is genuinely DEFECTIVE (a Jordan block, dep≈4, geo 1<alg 2), not a diabolic degeneracy. Gate-first (toy Jordan→DEFECTIVE, diag→DIABOLIC); the family that misfired in the F86a retraction is here corroborating-but-not-load-bearing",
            _ => new EpCharacterWitness(), RequiresN: false),
        new("f89octic", "the F89 path-3 octic EP-character (live EpCharacter): DIABOLIC (semisimple — eigenvalues coalesce, eigenvectors independent; geo=alg=2, dep≈0), NOT a defective EP. Grid-free root cause: the (3q⁴+q²−1) discriminant factor has even multiplicity 2. The diabolic sibling of the coherence-horizon defective √-EP (--root epcharacter)",
            _ => new F89OcticCharacterWitness(), RequiresN: false),
        new("crosstriple", "the F127/F128/F129/F130/F133/F134 family, live: the F127 GF(p) variety slice + core identity + closed form + sheet lattice, the F128 flip lemma (exact ℤ) + factorization + sharper locus, the F129 exact ℤ[ζ_2n] level-collision census (n ≤ 60) + family-inventory sum tie, the F130 exact collision decoupling (every proof cell + controls), the F133 symplectic closed form of W (sin-s lemma + read-off spot + GF(p) certificate), and the F134 two-row reflection law (table + read-off + the l = 2 fence). Each check carries its own discriminating control; the proof objects stay the 527/527 grid+CRT wall and the committed f12*/f13* gates. Seconds to ~a minute (the F133 halves build), fully deterministic",
            _ => new CrossTripleOrthogonalityWitness(), RequiresN: false),
        new("holonomy", "the eigenVECTOR holonomy around the (1,2)-block defective seed = the mod-4 memory loop i⁴=1: encircling the EP, the coalescing eigenvector frame (transported in the biorthogonal vᵀv gauge) rotates 90°/loop (M₁ eig ±i), M₂=−I, M₄=+I — the operator NinetyDegreeMirrorMemory at the LIVING seed, the eigenVECTOR-phase companion of the eigenVALUE swap (Numerics/Monodromy). Default N=5 (full block clean, ~1s); --N 9 = the 2026-07-07 living seed (~2 min; the witness auto-selects the R=+1 sector, where the full block leaks). Live vᵀv-gauge frame transport",
            c => BuildHolonomyRoot(c.Parser), RequiresN: false, HonorsOptionalN: true),
        new("f89galois", "the F89 path-k H_B-mixed Galois groups (k=3..6): Gal(F_d/Q(i)(q)) = S_d for d=8/18/32/53, NON-SOLVABLE ⟹ the decay rates λ_k(q) have no radical closure in q=J/γ. all four paths fully live (block → Berkowitz over Z[i] → isolate F_d by dividing out the reconstructed AT factor, validation triple → Frobenius ⟹ S_d); the AT factor is rebuilt from the rate-confined invariant subspace, never imported (the committed F_d literals are only a test cross-check). The Galois sibling of the EP-character witness (--root f89octic)",
            _ => new F89PathKGaloisWitness(), RequiresN: false),
        new("topowritability", "topology controls the relaxation's Galois writability: the wiring's automorphism group caps the (SE,DE) factor degrees N-independently when it is large (complete K_N → cap 4, radically writable for all N; star K_{1,N−1} → cap 9, a fixed S_9 scramble) and lets them grow when it is small (ring D_N, chain S_2 → S_8/18/32/53). Live: builds the (SE,DE) block, verifies [L,ρ(g)]=0 for Aut(G), and recomputes the cap = the standard-rep multiplicity via the character sum. The classification sibling of --root f89galois",
            _ => new TopologyGaloisWritabilityWitness(), RequiresN: false),
        new("galoischaos", "Galois vs spectral chaos: does the H_B-mixed half (chain Galois S_8/18/32/53, no radical closure) read as dissipative quantum chaos (GinUE) at fixed q? A clean NULL — it reads Poisson-like/sub-Poisson, still on the integrable frequency lattice. Algebraic chaos (Galois over the q=J/γ field) ≠ spectral chaos (RMT at fixed q); the former does not imply the latter. Live: the shared (SE,DE) block → MathNet EVD → complex spacing ratio (Sá-Ribeiro-Prosen), whose own GinUE reference confirms it CAN see chaos. The sector-resolved test the global RMT scan left open; the q-parametric monodromy (--root f89octic discriminant/EP loci) is where the Galois structure does live. The spectral-statistics sibling of --root topowritability",
            _ => new GaloisSpectralChaosWitness(), RequiresN: false),
        new("fillcsr", "the F89 Door-C decisive follow-up: is dissipative quantum chaos (GinUE) a FILLING threshold, not an integrability one? Door-C found the DILUTE (SE,DE)=(1,2) coherence block stays Poisson / non-GinUE under every integrability-breaking knob (Δ, random Z-field) because a 2-excitation sector cannot thermalize. This witness builds the GENERAL (wKet,wBra) block at EXTENSIVE filling (near N/2) and re-runs the same disordered CSR: the DENSE block IS chaotic — ⟨|z|⟩ at the GinUE value and ⟨cosθ⟩ negative, climbing toward GinUE with N (≈−0.09/−0.13/−0.16 at N=6/7/8 = 43/56/67% of the size-matched GinUE angle), while the dilute (1,2) stays flat at ⟨cosθ⟩≈0. So the null is a FILLING threshold: chaos switches on with extensive excitation content, not with breaking the Galois/Hamiltonian integrability. Class A licensed by the unequal weight (p,p+1) (Π maps it to the conjugate (p+1,p) block, not itself; conjugation-match ≈0 under disorder). Live: general WeightCoherenceBlock + random field → MathNet EVD → complex spacing ratio, pooled per-spectrum with finite-size-matched refs. The decisive sequel to --root galoischaos",
            _ => new FillingThresholdWitness(), RequiresN: false),
        new("galoismonodromy", "where the F89 octic's Galois structure lives spectrally: the q-parametric MONODROMY (the sequel to the --root galoischaos fixed-q null). As q=J/γ loops the complex plane the 8 octic roots braid, and that braiding IS the Galois group (monodromy = Galois over the function field). G2 live: a loop around the diabolic point q_EP≈0.659 returns the identity (the coalescing pair does not swap, a double discriminant zero = transversal crossing), confirming --root f89octic semisimple-not-defective by an independent route. The gateway to the S_8-generation gate (the simple zeros of P_10 carry transpositions that generate Gal(F_8)=S_8). Built on the trusted 12×12 block + the validated Monodromy tracker",
            _ => new GaloisMonodromyWitness(), RequiresN: false),
        new("branchpalindrome", "the F89 octic branch locus is a palindrome (the typed home F89BranchLocusPalindromeClaim): the EP/diabolic collisions are mirror-symmetric about Re λ = −σ = −4, FORCED by the F1 palindrome carried antiunitarily on the (SE,DE) block (T L(q) T⁻¹ = −L(q̄)−2σ, a same-q fold at real q; all-q holomorphic identity F8(λ,q)=F8(−λ−8,−q); spectral action λ↦−λ̄−2σ). Two-sided gate live: the octic closes under the antiunitary mirror, not the linear one; every EP on the line or in a mirror pair, no orphan. The line is the palindrome's gift, the diabolic's silence is integrability's (separate, XXZ Δ≠0 stays on-line yet defects). The q-direction sibling of --root galoismonodromy; reading reflections/ON_WHO_WATCHES_WHOM.md",
            _ => new BranchLocusPalindromeWitness(), RequiresN: false),
        new("diabolicparity", "the odd-N real-q diabolic onset, grounded from below: the dimension-mismatch / sector-swap. The site reflection R splits the (SE,DE) block into R-even/R-odd; the realness antiunitarity Σ L Σ = L† maps even↔odd at EVEN N (σ_even = conj σ_odd EXACTLY, neither sector self-conjugate ⟹ a real-axis collision in R-even is the generic pseudo-Hermitian DEFECTIVE EP), but the odd-N reflection-fixed central site makes dim(even)−dim(odd) = (N−1)/2 ≠ 0, forbidding the cross-pairing and forcing the antiunitarity within each sector (self-conjugate ⟹ R-even carries real eigenvalues, two of which cross semisimply = the real-q diabolic). Live: builds the full block, splits by R, measures self/cross-conjugacy across N=5..9, and reproduces the C# scout's N=7 (λ=−4.942) and N=9 (λ=−5.424) diabolics in R-even. The from-below grounding of the diabolic_over_higher_n arc; reading experiments/F89_PATH_K_DIABOLIC.md and SLOW_MODE_R_PARITY.md",
            _ => new DiabolicReflectionParityWitness(), RequiresN: false),
        new("crossfold", "Move 4, answered: the (SE,DE) diabolics PAIR across the (SE,DE)↔(SE,w_{N−2}) cross-block fold, because that fold is an EXACT antiunitary similarity. The branch-locus palindrome's bra bit-flip ρ[a,b]→ρ[a,b̄] (F89c lemma, n_diff(a,b̄)=N−n_diff) maps the (SE,DE)=(w1,w2) block to the (SE,w_{N−2})=(w1,N−2) block; the witness builds both blocks + the bra-complement permutation P and checks L(1,N−2)(q̄) = −P·conj(L(1,2)(q))·Pᵀ − 2N·I to machine zero (N=4..9, every q). STRONGER than foldcross's spectrum match: an antiunitary similarity preserves Jordan structure, so a diabolic at (q,λ) maps to a diabolic at (q̄,−λ̄−2N) with identical gap and character. N=4 is the degenerate partner=self within-block self-fold; for N≥5 the partner is a different block and the N=4 on-line 'zeros' become cross-block mirror partners. Reproduces the N=7 real-q diabolic pairing (λ=−4.942 ↔ −9.058). The Move-4 grounding of the diabolic_over_higher_n arc; reading experiments/F89_PATH_K_DIABOLIC.md and F89_BRANCH_LOCUS_PALINDROME.md",
            _ => new CrossFoldSimilarityWitness(), RequiresN: false),
        new("sectorbraid", "the Multi-Sector Monodromy verdict, live: is the F89 (1,2) octic's S₈ braid CONFINED to the (1,2) coherence orbit or SHARED across the joint-popcount sectors of L(q)? N-DEPENDENT (the finding): CONFINED to the D₄ orbit {(1,2),(2,1),(2,3),(3,2)} at N=4 (the dense half-filled core (2,2) is braid-free, a Door-C echo), but SPREADS at N=5 to a symmetric 12-sector diamond INCLUDING the dense core (2,2), splitting into two cross-fold-conjugate families of 6 that carry a BYTE-IDENTICAL shared eigenvalue λ (same branch point). Mechanism = free-fermion/AT additivity: the |bra−ket|=1 SE-DE rung's EP is invariant under a diagonal mode-spectator (so (1,2)≡(2,3)≡(3,4) byte-identically) and Family B is the F89d cross-fold image λ↦−λ̄−2N; membership = {|bra−ket|=1, popcounts∈[1,N−1]} ∪ cross-fold. The λ value is a γ-driven EP OFF the AT rate lines (Re λ not −2·integer at real loci); its single-particle-level construction is the open piece. Live: MultiSectorMonodromyCensus.Run(N) (~15s at N=4, ~2min at N=5) + the additivity-embedding gate at a real reference locus. Default N=4; --N 5 for the generic spread. Reading experiments/F89_MULTI_SECTOR_MONODROMY.md; sibling of --root crossfold (F89d) and --root galoismonodromy",
            c => new SectorBraidWitness(c.Parser.HasFlag("N") ? c.N : 4), RequiresN: false, HonorsOptionalN: true),
        new("monodromymirror", "how the F89 palindrome reaches the octic monodromy: the mirror SPLITS at the Galois boundary. C-K: the q↦−q̄ reflection (from L(q)*=L(−q̄)) intertwines the monodromy (σ_K = identity in the aligned labelling, each cluster EP carries the same braid as its mirror; forced by L(q)*=L(−q̄)), the branch-locus palindrome lifted from the seams' positions to their braids. C-T: the Re=−4 spectral fold induces a non-central strand involution σ_T (four fixed on the fold + two mirror-twin 2-cycles) that does NOT commute with the monodromy (commuting with the full S_8 would force it central, but Z(S_8)=1), so it is not a loop-independent symmetry of the braiding (conjugation by it is still an inner automorphism of S_8). The braid-level sibling of --root galoismonodromy and --root branchpalindrome; the from-below ground of reflections/ON_WHO_WATCHES_WHOM.md",
            _ => new MonodromyMirrorWitness(), RequiresN: false),
        new("secondclock", "the stitch: the {0,2}/half-filling coherence (the second clock) is ONE mode whose regime = map(band degeneracy, dispersion) — joining the EP horizon, the structural ceiling, and the gradual star into one node (live N=4 full-Liouvillian gate)",
            _ => new SecondClockRegimeWitness(), RequiresN: false),
        new("starseam", "the star's frozen seam: the longest-lived coherence never un-freezes (N≥5) — its survivor is the [H,A]=0 commutant (1,1) coherence, frozen by construction, the survivor iff g2=4/(N−1)≤1; N=4 (4/3>1) un-freezes (the (2,2)/K₄ outlier). The third member of chain(SE-EP)/ring(frozen crossing)/star(frozen commutant): the structural ceiling read dynamically",
            _ => new StarFrozenSeamWitness(), RequiresN: false),
        new("niven", "the Niven root: Niven's theorem on the SE cyclotomic angle π/(N+1) is the number-theoretic ceiling on the spectrum's closed forms (three faces; N=4 = first golden, band edge = φ); the arithmetic root of the small-N specials",
            _ => new NivenRationalityRootWitness(), RequiresN: false),
        new("transition", "F124 the band-edge transition invariant: the full bond-transition matrix M[b,k]=⟨ψ_k|V_b|ψ_1⟩ (all N modes) has ‖M‖_F² + λ_min(MMᵀ) = z = 2 exactly (‖M‖_F²=2−E, λ_min=E=(4/(N+1))sin²(π/(N+1))). The real content λ_min=E is the Dirichlet-edge coupling (an SSH/Peierls edge effect); frame reading λ_min=σ_min²=the lower frame bound, kernel = the K-partner ψ_N. Only the band-edge carrier makes staggered the genuine minimum (interior carrier → sum<2); the location dictionary k=2..N gives λ_min=0",
            _ => new BandEdgeTransitionInvariantWitness(), RequiresN: false),
        new("resolution", "F124's conditioning read as a defect-localization RESOLUTION LIMIT: σ_min=√E the reconstruction floor (the lower frame bound), κ=λ_max/λ_min ~ N² the noise amplification, contrast σ_max/σ_min=√κ ~ N (a staggered q=π zone-boundary defect √κ ~ N times harder to localize, matched-filter SNR); the worst direction is F124's staggered λ_min eigenvector (the q=π diffraction limit); the floor σ_min ~ (N+1)^(−3/2), E·(N+1)³ → 4π². One object in three trades (inverse problem / observability Gramian / optics MTF). NOT the decoder's 1.5 ambiguity",
            _ => new BandEdgeResolutionLimitWitness(), RequiresN: false),
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
        new("blockspectrum", "the joint-popcount block spectrum, live: the (N+1)² sector decomposition " +
            "(halved by X⊗N, quartered by the F1 Π orbit -- the 100->50->25 story at N=9), the F1 palindrome " +
            "{λ}={−2σ−λ} reconstructed sector-by-sector (full at N≤7), the (0,1) band-edge Absorption floor " +
            "Re=−2γ, and the N=9 banked headline read live from chain_N9.json -- the browsable face of the " +
            "SLOW_N9 result (arc block_spectrum_n9)",
            c => new BlockSpectrumWitness(
                    c.Parser.HasFlag("N") ? c.N : 6,
                    c.Parser.OptionalDouble("gamma") ?? 0.5,
                    c.Parser.OptionalDouble("J") ?? 1.0),
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
        new("stone", "THE STONE (felt_time arc B): the PTF painter closure Sum_i ln(alpha_i), via the CANONICAL " +
            "Symphony FitAlpha on the mode-isolating probe rho_0 = I/d + eps*Herm(mode), reads the mode's first-order " +
            "RATE shift - OUT + sign-coherent (rate-shift) for the soft survivor interior (2,2), IN (frozen) for the " +
            "rigid (0,1) band edge: the TRAJECTORY-level dual of the eigenvalue value/vector split. Probe-state-specific " +
            "(review-pinned, not a universal law), the rate shift certified by sign-coherence; N in 4..5",
            c => new StoneSurvivorClosureWitness(
                    c.Parser.HasFlag("N") ? c.N : 4,
                    c.Parser.OptionalDouble("q") ?? 1.5,
                    c.Parser.OptionalDouble("delta-j") ?? 0.02),
            RequiresN: false, HonorsOptionalN: true),
        new("gradient", "(D) THE CLOSURE FUNCTIONAL (felt_time arc D): the survivor's first-order bond rate shift " +
            "dRe(b) ~ (density-mode gradient at bond b)^2 - the diffusion Rayleigh quotient (amplitude^2). The slow " +
            "survivor is a DENSITY/diffusion mode; a delta-J defect perturbs the local diffusion coefficient, so dRe ~ " +
            "(n(j)-n(j+1))^2, ~0 at the no-flux chain ends, mirror-symmetric, Q-invariant. The eigenvalue-level dual of " +
            "the PTF closure (inspect --root stone). dRe/grad^2 bond-independent, log-log slope ~2; N in 4..5",
            c => new SurvivorDiffusionGradientWitness(
                    c.Parser.HasFlag("N") ? c.N : 4,
                    c.Parser.OptionalDouble("q") ?? 1.5),
            RequiresN: false, HonorsOptionalN: true),
        new("renewal", "THE DEPHASING-FRONT RENEWAL REPRESENTATION (proof PROOF_DEPHASING_FRONT_RENEWAL.md): the " +
            "exact solution of the WATCHED walk. The watched single excitation is the unwatched wave repeatedly " +
            "caught and released, P_n(t) = e^{−Γt}·S_n(t) with the Volterra refill ladder (★), Γ = 4γ. The j=0 " +
            "term is the coherent front, the j≥1 halo the incoherent refill; the ladder closes to " +
            "Ŝ(p,z) = 1/(√(z²+a²)−Γ), a = 4J·sin(p/2), with the diffusive pole D = 2J²/Γ (the F123 sibling). Six " +
            "from-below checks: renewal-vs-RK4, probability conservation, the coherent-front Bessel identity, the " +
            "Γ=0 clean-wave limit, the Haken-Strobl plateau, the I₁ Airy constant. Typed: DephasingFrontRenewalClaim",
            c => new DephasingFrontRenewalWitness(
                    c.Parser.HasFlag("N") ? c.N : 27,
                    seed: null,                                       // interior seed = N/2
                    j: c.Parser.OptionalDouble("J") ?? 1.0,
                    gamma: c.Parser.OptionalDouble("gamma") ?? 0.15),
            RequiresN: false, HonorsOptionalN: true),
        new("record", "THE RECORD LAWS (F135 + F136, proofs PROOF_RECORD_PARITY_LAW.md + PROOF_RECORD_LETTER_LAW.md): " +
            "who records, and WHAT. Under H = Σ Δ·ZZ from |+⟩^⊗N every pair page is closed-form for all t (the " +
            "Absorption substrate); at t* = π/4 the watcher-ratio parities decide (even → perfect record, odd → " +
            "blind, non-integer → 1−h₂((1+β)/2)), and the shared dressers' parity picks the family: all even + " +
            "write bond → pointer record of Z_S, all odd → Bell record (letter = dresser parity, m odd YY / m even " +
            "XX, zero pointer content), mixed → dark. Bell pays both sites, pointer only the witness, watching the " +
            "writers is free; a pendant S role-swaps; anti-pointer redundancy is not deg-bounded (K_{R+1,2}); " +
            "fully-witnessed worlds = stars + K_N, girth ≥ 5 + leafless = dark. The battery recomputes every case " +
            "closed-form-vs-full-state at inspect time. Typed: RecordParityLawClaim + RecordLetterLawClaim",
            c => new RecordLawWitness(c.Parser.OptionalDouble("gamma") ?? 0.05),
            RequiresN: false),
        new("trichotomy",
            "the chain/ring/star survivor trichotomy as one sweep: carbon un-freeze read (RouteSweep) + absolute Δn-seam read",
            c => new TrichotomyWitness(
                    c.Parser.HasFlag("N") ? c.N : 6,
                    c.Parser.OptionalDouble("q") ?? 1.5),
            DefaultDepth: 2, RequiresN: false, HonorsOptionalN: true),
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
        new("directsum", "the direct-sum decomposition live: L = L_even ⊕ L_odd by n_XY parity (off-parity block " +
            "exactly 0, sectors 2^(2N−1) each), Π sector map checked column-complete (odd N exchanges, even N " +
            "preserves), sector-restricted palindrome machine-zero, plus the selective-breaking cross as controls " +
            "(T1 breaks the mirror not the wall; a transverse field the wall not the mirror) " +
            "(typed: DirectSumDecompositionClaim)",
            c => new DirectSumDecompositionWitness(c.Parser.HasFlag("N") ? c.N : 3,
                                                   c.Parser.OptionalDouble("gamma") ?? 0.05),
            RequiresN: false, HonorsOptionalN: true),
        new("zeroimmune", "Zero-Sector Immunity, live: a random parity-violating 2-body H gives M ≈ 0 on the " +
            "w=0 ({I,Z}^⊗N) and w=N ({X,Y}^⊗N) Pauli blocks while ‖M‖ > 0 (the non-trivial gate); the classical " +
            "extreme is immune to every 2-body coupling, the palindrome-breaking lives in 0<w<N " +
            "(typed: ZeroSectorImmunityClaim)",
            c => new ZeroSectorImmunityWitness(c.Parser.HasFlag("N") ? c.N : 3,
                                               c.Parser.OptionalDouble("gamma") ?? 0.05),
            RequiresN: false, HonorsOptionalN: true),
        new("selfmirror", "the self-mirror fixed point AS AN OBJECT (the jump: the System bears x/y/z, the Object inherits it). It does NOT own x/y/z (z=watched Re λ, x-y=motion Im λ, σ=center, all the system's); its only delta is BEING the F1 fixed point Re λ = −σ, the k=N/2 self-mirror rung (even-N populated, odd-N empty)",
            c => BuildSelfMirrorRoot(c.Parser, c.Parser.HasFlag("N") ? c.N : 4), RequiresN: false, HonorsOptionalN: true),
        new("seedcount", "the seed-existence counting theorem live: r(0⁺) − r(∞) = N − 1 (odd N) on the (1,2) " +
            "block via the three counting lemmas: (N2) the −2 rung = N−1 paths of N vertices, (FF) nullity(C) " +
            "= the fusion-resonance count, (N1′) the ordering-sector theorem (K₆ = three no-passing sectors " +
            "gauged to −H₃, spec(K₆) = 3×{−(λa+λb+λc)}). SVD nullities + combinatorial counts + the exact-zero " +
            "cross-sector/gauge gates + two-sided nonzero controls, recomputed at inspect time; odd N ≤ 9 " +
            "(typed: SeedExistenceCountingClaim)",
            c => new SeedExistenceCountingWitness(c.Parser.HasFlag("N") ? c.N : 5),
            RequiresN: false, HonorsOptionalN: true),
        new("betaexotic", "the β-exotic excluded exactly at N = 5 and N = 7, both R-parities: the certified " +
            "squarefree-layer reading of disc_Λ(F_res) has maximum root multiplicity 2 off q = 0, and a " +
            "Puiseux-3/2 point (the β-exotic) would need 3. Re-runs the D-only certificate at inspect time " +
            "(~0.4 s per parity at N = 5, ~2.5 min at N = 7, so N = 5 is the default; pass --N 7 to pay for it) " +
            "and reads MaxDiscMultiplicity + DiscLayersCertified; the one-way lift needs a prime good at BOTH " +
            "ends of the q-axis. Per-N, not a law: N = 9 is out of reach by this route, the all-N scalar s₆ ≠ 0 " +
            "stays open (typed: BetaExoticPerNExclusionClaim)",
            c => new BetaExoticExclusionWitness(c.Parser.HasFlag("N") ? c.N : 5),
            RequiresN: false, HonorsOptionalN: true),
        new("label", "the label layer, typed: the watcher is its letter - the 4^N Pauli strings are ONE shared " +
            "eigenbasis of all three letter dissipators with three price lists (rate -2g*n_anti(S,P), the " +
            "disagreement with the held letter alone); the letter swap (Klein V4 / basis-S3) relocates which " +
            "cells pay, entry-exactly; only the identity rides free under every watcher. The exact core of " +
            "docs/quantum LABELS_TRANSLATED s2 and DEPHASING_TRANSLATED s4: even the environment routes by a " +
            "label (typed: WatchedLetterRoutingClaim)",
            c => new WatchedLetterRoutingWitness(c.Parser.HasFlag("N") ? c.N : 3,
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
