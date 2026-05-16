using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>Bond Hamiltonian scenario behind the F87-trichotomy IBM tomography
/// snapshots, plus the Heisenberg control. Used as the typed key into the
/// <see cref="IbmBlockCpsiHardwareTable"/> witness rows.</summary>
public enum BlockCpsiScenario
{
    /// <summary>Heisenberg control: H = (J/4) Σ_b (X_b X_{b+1} + Y_b Y_{b+1} + Z_b Z_{b+1}); Π²-even, full SU(2)-symmetric.</summary>
    Heisenberg,
    /// <summary>F87 truly: H = (J/4) Σ_b (X_b X_{b+1} + Y_b Y_{b+1}); Π²-even bilinears.</summary>
    Truly,
    /// <summary>F87 soft: H = (J/4) Σ_b (X_b Y_{b+1} + Y_b X_{b+1}); Π²-odd bilinears.</summary>
    Soft,
    /// <summary>F87 hard: H = (J/4) Σ_b (X_b X_{b+1} + X_b Y_{b+1}); mixed Π²-even + Π²-odd.</summary>
    Hard,
}

/// <summary>F86 Tier-2-Verified table: state-level C_block readouts on the
/// 2026-04-26 IBM framework_snapshots (Marrakesh / Kingston / Fez hardware plus the
/// Aer Trotter-noiseless simulator baseline) through Theorem 2's universal 1/4 ceiling.
///
/// <para>Source pipeline: <c>simulations/_block_cpsi_lens_ibm_snapshots.py</c>; pinned
/// values from <c>simulations/results/_block_cpsi_lens_ibm_snapshots.txt</c>. The 2-qubit
/// reduced state ρ on (q0, q2) of the N=3 |+−+⟩ chain at t=0.8, J=1.0 is reconstructed from
/// the 16-Pauli expectations in the snapshot JSON, and <see cref="BlockCoherenceContent.Compute"/>
/// is applied to the (popcount-0, popcount-1) and (popcount-1, popcount-2) coherence
/// blocks for each of the four scenarios (Heisenberg / Truly / Soft / Hard).</para>
///
/// <para><b>Theorem 2 satisfaction:</b> every witness has C_block ≤ 1/4 + 1e-12 (the
/// universal Mandelbrot-cardioid ceiling per <see cref="Symmetry.QuarterAsBilinearMaxvalClaim"/>);
/// asserted in the matching test. The lens has order-of-magnitude resolution at the
/// hardware level: the actual min/max range across the 24 ibm_* witnesses is exposed
/// live in the inspect tree (the "hardware resolution" child computes from the witness
/// list). The state in question (|+−+⟩ at N=3 reduced to (q0, q2)) is structurally far
/// from the canonical Dicke maximiser <c>(|D_n⟩+|D_{n+1}⟩)/√2</c>, so the lens cannot
/// approach the 1/4 ceiling on this data; a future Dicke-prepared run is needed for that.</para>
///
/// <para><b>Π²-odd asymmetry signature (ideal continuous evolution):</b></para>
/// <list type="bullet">
///   <item>Heisenberg block(0,1) ≡ block(1,2) → ratio 1.000 (Π²-symmetric)</item>
///   <item>Truly      block(0,1) ≡ block(1,2) → ratio 1.000 (Π²-symmetric)</item>
///   <item>Soft       block(0,1) &gt; block(1,2) → ratio &gt; 1 (Π²-odd-broken)</item>
///   <item>Hard       block(0,1) &gt; block(1,2) → ratio &gt; 1 (Π²-odd-broken, weaker than Soft)</item>
/// </list>
/// <para>The Π²-odd-driven Hamiltonians (Soft, Hard) break the popcount-block symmetry
/// that Truly and Heisenberg preserve. Whether this ratio is derivable from F88
/// popcount-coherence remains an OpenQuestion (see <see cref="OpenQuestionPi2OddAsymmetryDerivation"/>);
/// the live numerical values come from <see cref="IdealAsymmetryRatios"/>.</para>
///
/// <para><b>Trotter dominates hardware noise on Soft block(0,1):</b> the continuous-time
/// ideal does not survive the 3-step Trotter approximation; Aer noiseless already loses
/// most of the Soft regime before any hardware noise enters. Visible directly in the
/// witness rows (compare CBlockMeasured for Aer-Soft-block(0,1) against CBlockIdealContinuous).</para>
/// </summary>
public sealed class IbmBlockCpsiHardwareTable : Claim
{
    /// <summary>1/4 = (1/2)² bilinear-apex maxval — the typed Pi2-Foundation
    /// parent that grounds the universal Theorem 2 ceiling against which every
    /// row in <see cref="Witnesses"/> is asserted. Each <c>FractionOfQuarter</c>
    /// reading is literally <c>CBlockMeasured / Quarter.Value</c>. Added
    /// 2026-05-16 as a typed ctor parent — the safest F86-Sammelbecken edge,
    /// promoting the Theorem-2 reading from registration-discard / docstring
    /// mention into the actual ancestor graph (per Wave 5 of the 2026-05-16
    /// inheritance-map sweep).</summary>
    public Symmetry.QuarterAsBilinearMaxvalClaim Quarter { get; }

    public IbmBlockCpsiHardwareTable(Symmetry.QuarterAsBilinearMaxvalClaim quarter)
        : base(
            name: "IBM 2026-04-26 framework_snapshots through Theorem 2 C_block lens",
            tier: Tier.Tier2Verified,
            anchor: "simulations/_block_cpsi_lens_ibm_snapshots.py + simulations/results/_block_cpsi_lens_ibm_snapshots.txt + docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md (Theorem 2) + compute/RCPsiSquared.Core/Symmetry/QuarterAsBilinearMaxvalClaim (typed parent — Theorem 2 ceiling)")
    {
        Quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
    }

    public override string DisplayName =>
        "IBM 2026-04-26 C_block snapshots vs 1/4 ceiling (4 backends × 4 scenarios × 2 blocks)";

    public override string Summary =>
        "Tier 2 Verified: state-level C_block on (q0, q2) reduced ρ from N=3 |+−+⟩ at t=0.8, J=1.0 across Aer / Marrakesh / Kingston / Fez. All 32 witnesses satisfy Theorem 2's C_block ≤ 1/4. The (|+−+⟩, q0/q2) state is structurally far from the canonical Dicke maximiser, so the lens does not approach the ceiling on this data. Π²-odd asymmetry block(0,1)/block(1,2) > 1 for Soft and Hard scenarios at the continuous-time ideal level, =1 for Heisenberg + Truly.";

    /// <summary>Run date of the IBM tomography snapshots that produced these witnesses.</summary>
    public const string RunDate = "2026-04-26";

    /// <summary>Initial state of the N=3 chain whose (q0, q2) reduced ρ is the lens target.</summary>
    public const string InitialState = "|+−+⟩ at N=3";

    /// <summary>Evaluation time of the underlying tomography circuits.</summary>
    public const double EvalTime = 0.8;

    /// <summary>Coupling J of the bond Hamiltonian for the four scenarios.</summary>
    public const double CouplingJ = 1.0;

    /// <summary>Trotter steps used in the circuit assembly (matches the JSON's parameters).</summary>
    public const int TrotterSteps = 3;

    /// <summary>The pinned 32-row table: (backend, scenario, block n) → measured C_block
    /// + the noiseless continuous-time-evolution ideal at the same parameters.</summary>
    public IReadOnlyList<BlockCpsiSnapshotWitness> Witnesses => _witnesses;

    /// <summary>Π²-odd asymmetry ratio: ideal continuous-time block(0,1) divided by
    /// block(1,2) per scenario, computed from <see cref="Witnesses"/>. Symmetric scenarios
    /// (Heisenberg, Truly) sit at 1.0; Π²-odd scenarios (Soft, Hard) lift above 1.0
    /// because Π²-odd dynamics break the popcount-block symmetry.</summary>
    public IReadOnlyDictionary<BlockCpsiScenario, double> IdealAsymmetryRatios => _idealAsymmetryRatios.Value;

    private readonly Lazy<IReadOnlyDictionary<BlockCpsiScenario, double>> _idealAsymmetryRatios = new(() =>
    {
        var dict = new Dictionary<BlockCpsiScenario, double>();
        foreach (var scen in Scenarios)
        {
            double c01 = _witnesses.First(w => w.Scenario == scen && w.BlockN == 0).CBlockIdealContinuous;
            double c12 = _witnesses.First(w => w.Scenario == scen && w.BlockN == 1).CBlockIdealContinuous;
            dict[scen] = c12 == 0.0 ? double.PositiveInfinity : c01 / c12;
        }
        return dict;
    });

    /// <summary>Open question: is the Π²-odd block(0,1) / block(1,2) asymmetry ratio
    /// (Soft and Hard each > 1, Heisenberg and Truly each = 1) derivable in closed form
    /// from F88 popcount-coherence at the operator level? Closure path: match against
    /// <see cref="Symmetry.PopcountCoherencePi2Odd"/> at the matching (n_p, n_q, HD)
    /// shape, plus Π²-odd integral on |+−+⟩.</summary>
    public string OpenQuestionPi2OddAsymmetryDerivation =>
        "Is the Π²-odd block(0,1)/block(1,2) asymmetry ratio (Soft and Hard each > 1, " +
        "Heisenberg and Truly each = 1) derivable in closed form from F88 popcount-coherence " +
        "+ Π²-odd integral on |+−+⟩? Tier-1-promotion path for the empirical pattern documented here.";

    public static IReadOnlyList<BlockCpsiScenario> Scenarios { get; } = new[]
    {
        BlockCpsiScenario.Heisenberg,
        BlockCpsiScenario.Truly,
        BlockCpsiScenario.Soft,
        BlockCpsiScenario.Hard,
    };

    public static IReadOnlyList<string> Backends { get; } = new[]
    {
        "aer (Trotter noiseless)", "ibm_marrakesh", "ibm_kingston", "ibm_fez",
    };

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("source pipeline",
                summary: "simulations/_block_cpsi_lens_ibm_snapshots.py (one-shot exploration; reuses reconstruct_2qubit_rho from _f88_lens_ibm_framework_snapshots)");
            yield return new InspectableNode("run date",
                summary: $"{RunDate}; initial = {InitialState}, t = {EvalTime}, J = {CouplingJ}, n_trotter = {TrotterSteps}");
            yield return new InspectableNode("Theorem 2 ceiling",
                summary: $"every witness satisfies C_block ≤ 1/4 = {BlockCoherenceContent.Quarter} (PROOF_BLOCK_CPSI_QUARTER)");
            yield return BuildHardwareResolutionNode();
            yield return InspectableNode.Group("ideal asymmetry block(0,1)/block(1,2) per scenario",
                BuildAsymmetryNodes().Cast<IInspectable>().ToArray());
            yield return new InspectableNode("OpenQuestion: Π²-odd asymmetry derivation",
                summary: OpenQuestionPi2OddAsymmetryDerivation);
            yield return InspectableNode.Group("witnesses (32 rows: 4 backends × 4 scenarios × 2 blocks)",
                _witnesses.Cast<IInspectable>().ToArray());
        }
    }

    private InspectableNode BuildHardwareResolutionNode()
    {
        var ibm = _witnesses.Where(w => w.Backend.StartsWith("ibm_")).ToList();
        var lo = ibm.OrderBy(w => w.CBlockMeasured).First();
        var hi = ibm.OrderByDescending(w => w.CBlockMeasured).First();
        return new InspectableNode("hardware resolution",
            summary:
                $"C_block range across the {ibm.Count} hardware witnesses: " +
                $"{lo.CBlockMeasured:F4} ({lo.Backend.Replace("ibm_", "")} {lo.Scenario} block({lo.BlockN},{lo.BlockN + 1})) " +
                $"to {hi.CBlockMeasured:F4} ({hi.Backend.Replace("ibm_", "")} {hi.Scenario} block({hi.BlockN},{hi.BlockN + 1})), " +
                $"i.e. {lo.FractionOfQuarter * 100:F1}% to {hi.FractionOfQuarter * 100:F1}% of the 1/4 ceiling");
    }

    private IEnumerable<InspectableNode> BuildAsymmetryNodes()
    {
        foreach (var scen in Scenarios)
        {
            double r = IdealAsymmetryRatios[scen];
            yield return new InspectableNode($"[{scen}] block(0,1)/(1,2)",
                summary: $"ratio = {r:F3} ({(Math.Abs(r - 1.0) < 0.01 ? "Π²-symmetric" : "Π²-odd-broken")})");
        }
    }

    private static readonly BlockCpsiSnapshotWitness[] _witnesses = new[]
    {
        new BlockCpsiSnapshotWitness("aer (Trotter noiseless)", BlockCpsiScenario.Heisenberg, 0, 0.0378, 0.0642),
        new BlockCpsiSnapshotWitness("aer (Trotter noiseless)", BlockCpsiScenario.Heisenberg, 1, 0.0280, 0.0642),
        new BlockCpsiSnapshotWitness("aer (Trotter noiseless)", BlockCpsiScenario.Truly,      0, 0.0192, 0.1064),
        new BlockCpsiSnapshotWitness("aer (Trotter noiseless)", BlockCpsiScenario.Truly,      1, 0.0219, 0.1064),
        new BlockCpsiSnapshotWitness("aer (Trotter noiseless)", BlockCpsiScenario.Soft,       0, 0.0279, 0.1838),
        new BlockCpsiSnapshotWitness("aer (Trotter noiseless)", BlockCpsiScenario.Soft,       1, 0.0596, 0.0499),
        new BlockCpsiSnapshotWitness("aer (Trotter noiseless)", BlockCpsiScenario.Hard,       0, 0.0879, 0.1707),
        new BlockCpsiSnapshotWitness("aer (Trotter noiseless)", BlockCpsiScenario.Hard,       1, 0.0422, 0.0785),
        new BlockCpsiSnapshotWitness("ibm_marrakesh",           BlockCpsiScenario.Heisenberg, 0, 0.0412, 0.0642),
        new BlockCpsiSnapshotWitness("ibm_marrakesh",           BlockCpsiScenario.Heisenberg, 1, 0.0388, 0.0642),
        new BlockCpsiSnapshotWitness("ibm_marrakesh",           BlockCpsiScenario.Truly,      0, 0.0386, 0.1064),
        new BlockCpsiSnapshotWitness("ibm_marrakesh",           BlockCpsiScenario.Truly,      1, 0.0217, 0.1064),
        new BlockCpsiSnapshotWitness("ibm_marrakesh",           BlockCpsiScenario.Soft,       0, 0.0310, 0.1838),
        new BlockCpsiSnapshotWitness("ibm_marrakesh",           BlockCpsiScenario.Soft,       1, 0.0753, 0.0499),
        new BlockCpsiSnapshotWitness("ibm_marrakesh",           BlockCpsiScenario.Hard,       0, 0.1099, 0.1707),
        new BlockCpsiSnapshotWitness("ibm_marrakesh",           BlockCpsiScenario.Hard,       1, 0.0543, 0.0785),
        new BlockCpsiSnapshotWitness("ibm_kingston",            BlockCpsiScenario.Heisenberg, 0, 0.0453, 0.0642),
        new BlockCpsiSnapshotWitness("ibm_kingston",            BlockCpsiScenario.Heisenberg, 1, 0.0502, 0.0642),
        new BlockCpsiSnapshotWitness("ibm_kingston",            BlockCpsiScenario.Truly,      0, 0.0328, 0.1064),
        new BlockCpsiSnapshotWitness("ibm_kingston",            BlockCpsiScenario.Truly,      1, 0.0329, 0.1064),
        new BlockCpsiSnapshotWitness("ibm_kingston",            BlockCpsiScenario.Soft,       0, 0.0409, 0.1838),
        new BlockCpsiSnapshotWitness("ibm_kingston",            BlockCpsiScenario.Soft,       1, 0.0920, 0.0499),
        new BlockCpsiSnapshotWitness("ibm_kingston",            BlockCpsiScenario.Hard,       0, 0.1134, 0.1707),
        new BlockCpsiSnapshotWitness("ibm_kingston",            BlockCpsiScenario.Hard,       1, 0.0460, 0.0785),
        new BlockCpsiSnapshotWitness("ibm_fez",                 BlockCpsiScenario.Heisenberg, 0, 0.0140, 0.0642),
        new BlockCpsiSnapshotWitness("ibm_fez",                 BlockCpsiScenario.Heisenberg, 1, 0.0084, 0.0642),
        new BlockCpsiSnapshotWitness("ibm_fez",                 BlockCpsiScenario.Truly,      0, 0.0060, 0.1064),
        new BlockCpsiSnapshotWitness("ibm_fez",                 BlockCpsiScenario.Truly,      1, 0.0166, 0.1064),
        new BlockCpsiSnapshotWitness("ibm_fez",                 BlockCpsiScenario.Soft,       0, 0.0132, 0.1838),
        new BlockCpsiSnapshotWitness("ibm_fez",                 BlockCpsiScenario.Soft,       1, 0.0303, 0.0499),
        new BlockCpsiSnapshotWitness("ibm_fez",                 BlockCpsiScenario.Hard,       0, 0.0845, 0.1707),
        new BlockCpsiSnapshotWitness("ibm_fez",                 BlockCpsiScenario.Hard,       1, 0.0347, 0.0785),
    };
}

/// <summary>One row of the IBM C_block snapshot table: a frozen empirical readout for
/// (backend, scenario, block n) at the 2026-04-26 |+−+⟩ run.</summary>
/// <param name="Backend">Source backend label (e.g. "ibm_marrakesh").</param>
/// <param name="Scenario">Hamiltonian scenario from <see cref="BlockCpsiScenario"/>.</param>
/// <param name="BlockN">Block index n; the block is (popcount-n, popcount-(n+1)). Only n ∈ {0, 1} for the 2-qubit reduced ρ.</param>
/// <param name="CBlockMeasured">Measured C_block content from the snapshot JSON's 16-Pauli reconstruction.</param>
/// <param name="CBlockIdealContinuous">Reference: C_block from the noiseless continuous-time ideal at the same parameters.</param>
public sealed record BlockCpsiSnapshotWitness(
    string Backend,
    BlockCpsiScenario Scenario,
    int BlockN,
    double CBlockMeasured,
    double CBlockIdealContinuous
) : IInspectable
{
    /// <summary>Fraction of the universal 1/4 ceiling: <c>CBlockMeasured / 0.25</c>.</summary>
    public double FractionOfQuarter => CBlockMeasured / BlockCoherenceContent.Quarter;

    /// <summary>Deviation from the noiseless continuous-time ideal:
    /// <c>CBlockMeasured − CBlockIdealContinuous</c>. Mixes Trotter discretisation error
    /// and (for ibm_* backends) hardware noise.</summary>
    public double DeltaFromIdeal => CBlockMeasured - CBlockIdealContinuous;

    public string DisplayName =>
        $"[{Backend}] {Scenario} block({BlockN},{BlockN + 1})";

    public string Summary =>
        $"C_block = {CBlockMeasured:F4} ({FractionOfQuarter * 100:F1}% of 1/4); " +
        $"ideal = {CBlockIdealContinuous:F4}; Δ = {DeltaFromIdeal:+0.0000;-0.0000;+0.0000}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("backend", summary: Backend);
            yield return new InspectableNode("scenario", summary: Scenario.ToString());
            yield return new InspectableNode("block (n, n+1)", summary: $"({BlockN}, {BlockN + 1})");
            yield return InspectableNode.RealScalar("CBlockMeasured", CBlockMeasured, "F4");
            yield return InspectableNode.RealScalar("CBlockIdealContinuous", CBlockIdealContinuous, "F4");
            yield return InspectableNode.RealScalar("FractionOfQuarter", FractionOfQuarter, "F3");
            yield return InspectableNode.RealScalar("DeltaFromIdeal", DeltaFromIdeal, "F4");
        }
    }

    public InspectablePayload Payload =>
        new InspectablePayload.Real($"C_block ({Backend}, {Scenario}, ({BlockN},{BlockN + 1}))", CBlockMeasured, "F4");
}
