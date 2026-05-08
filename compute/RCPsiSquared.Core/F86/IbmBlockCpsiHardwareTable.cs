using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 Tier-2-Verified table: state-level C_block readouts on the
/// 2026-04-26 IBM framework_snapshots (Marrakesh / Kingston / Fez hardware plus the
/// Aer Trotter-noiseless simulator baseline) through Theorem 2's universal 1/4 ceiling.
///
/// <para>Source pipeline: <c>simulations/_block_cpsi_lens_ibm_snapshots.py</c>; pinned
/// values from <c>simulations/results/_block_cpsi_lens_ibm_snapshots.txt</c>. The 2-qubit
/// reduced state ρ on (q0, q2) of the N=3 |+−+⟩ chain at t=0.8, J=1.0 is reconstructed from
/// the 16-Pauli expectations in the snapshot JSON, and <see cref="BlockCoherenceContent.Compute"/>
/// is applied to the (popcount-0, popcount-1) and (popcount-1, popcount-2) coherence
/// blocks for each of the four scenarios (heisenberg / truly / soft / hard).</para>
///
/// <para><b>Theorem 2 satisfaction:</b> every witness has C_block ≤ 1/4 + 1e-12 (the
/// universal Mandelbrot-cardioid ceiling per <see cref="Symmetry.QuarterAsBilinearMaxvalClaim"/>);
/// asserted in the matching test. Hardware values range 2.4% to 45.4% of the ceiling: the
/// lens has order-of-magnitude resolution, refuting the "all noisy hardware sits near 1/4"
/// concern. Decoherence pulls states AWAY from the ceiling (toward 0), not toward it.</para>
///
/// <para><b>Π²-odd asymmetry signature (ideal continuous evolution):</b></para>
/// <list type="bullet">
///   <item>heisenberg block(0,1) ≡ block(1,2) = 0.0642 (ratio 1.000, Π²-symmetric)</item>
///   <item>truly      block(0,1) ≡ block(1,2) = 0.1064 (ratio 1.000, Π²-symmetric)</item>
///   <item>soft       block(0,1) = 0.1838, block(1,2) = 0.0499 (ratio 3.683, Π²-odd)</item>
///   <item>hard       block(0,1) = 0.1707, block(1,2) = 0.0785 (ratio 2.175, mixed)</item>
/// </list>
/// <para>The Π²-odd-driven Hamiltonians (soft, hard) break the popcount-block symmetry
/// that truly + heisenberg preserve. Whether this ratio is derivable from F88
/// popcount-coherence remains an OpenQuestion (see <see cref="OpenQuestionPi2OddAsymmetryDerivation"/>).</para>
///
/// <para><b>Trotter dominates hardware noise on soft block(0,1):</b> continuous → Aer
/// drops from 0.184 to 0.028 (−0.156); Aer → hw mean adds only +0.001. The Trotterised
/// approximation washes most of the soft-regime ideal before hardware noise enters.</para>
/// </summary>
public sealed class IbmBlockCpsiHardwareTable : Claim
{
    public IbmBlockCpsiHardwareTable()
        : base(
            name: "IBM 2026-04-26 framework_snapshots through Theorem 2 C_block lens",
            tier: Tier.Tier2Verified,
            anchor: "simulations/_block_cpsi_lens_ibm_snapshots.py + simulations/results/_block_cpsi_lens_ibm_snapshots.txt + docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md (Theorem 2)")
    { }

    public override string DisplayName =>
        "IBM 2026-04-26 C_block snapshots vs 1/4 ceiling (4 backends × 4 scenarios × 2 blocks)";

    public override string Summary =>
        "Tier 2 Verified: state-level C_block on (q0, q2) reduced ρ from N=3 |+−+⟩ at t=0.8, J=1.0 across Aer / Marrakesh / Kingston / Fez. All 32 witnesses satisfy Theorem 2's C_block ≤ 1/4. Hardware range: 2.4% – 45.4% of the 1/4 ceiling. Π²-odd asymmetry block(0,1):block(1,2) ≈ 3.7 (soft) / 2.2 (hard), 1.0 (truly + heisenberg) at the continuous-time ideal level.";

    /// <summary>Run date of the IBM tomography snapshots that produced these witnesses.
    /// All 32 measurements come from the same 2026-04-26 jobs (file timestamp 105948 for
    /// the three IBM hardware backends; Aer baseline file 012516).</summary>
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
    /// (heisenberg, truly) sit at 1.0; Π²-odd scenarios (soft, hard) lift above 1.0
    /// because Π²-odd dynamics break the popcount-block symmetry.</summary>
    public IReadOnlyDictionary<string, double> IdealAsymmetryRatios => _idealAsymmetryRatios.Value;

    private readonly Lazy<IReadOnlyDictionary<string, double>> _idealAsymmetryRatios = new(() =>
    {
        var dict = new Dictionary<string, double>();
        foreach (string scen in Scenarios)
        {
            double c01 = _witnesses.First(w => w.Scenario == scen && w.BlockN == 0).CBlockIdealContinuous;
            double c12 = _witnesses.First(w => w.Scenario == scen && w.BlockN == 1).CBlockIdealContinuous;
            dict[scen] = c12 == 0.0 ? double.PositiveInfinity : c01 / c12;
        }
        return dict;
    });

    /// <summary>Open question: is the Π²-odd block(0,1) / block(1,2) asymmetry ratio
    /// (≈ 3.683 for soft, ≈ 2.175 for hard at the ideal continuous level) derivable
    /// in closed form from F88 popcount-coherence at the operator level? Closure path:
    /// match against <see cref="Symmetry.PopcountCoherencePi2Odd"/> at the matching
    /// (n_p, n_q, HD) shape, plus Π²-odd integral on |+−+⟩.</summary>
    public string OpenQuestionPi2OddAsymmetryDerivation =>
        "Is the Π²-odd block(0,1)/block(1,2) asymmetry ratio (3.683 soft, 2.175 hard) " +
        "derivable in closed form from F88 popcount-coherence + Π²-odd integral on |+−+⟩? " +
        "Tier-1-promotion path for the empirical pattern documented here.";

    public static IReadOnlyList<string> Scenarios { get; } = new[]
    {
        "heisenberg", "truly", "soft", "hard",
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
            yield return new InspectableNode("hardware resolution",
                summary: "C_block range across the 24 hardware witnesses: 0.0060 (Fez truly block(0,1)) to 0.1134 (Kingston hard block(0,1)), i.e. 2.4% to 45.4% of the 1/4 ceiling");
            foreach (string scen in Scenarios)
            {
                double r = IdealAsymmetryRatios[scen];
                yield return new InspectableNode($"ideal asymmetry [{scen}] block(0,1)/(1,2)",
                    summary: $"ratio = {r:F3} ({(Math.Abs(r - 1.0) < 0.01 ? "Π²-symmetric" : "Π²-odd-broken")})");
            }
            yield return new InspectableNode("OpenQuestion: Π²-odd asymmetry derivation",
                summary: OpenQuestionPi2OddAsymmetryDerivation);
            yield return InspectableNode.Group("witnesses (32 rows: 4 backends × 4 scenarios × 2 blocks)",
                _witnesses.Cast<IInspectable>().ToArray());
        }
    }

    private static readonly BlockCpsiSnapshotWitness[] _witnesses = new[]
    {
        // Aer (Trotter noiseless)
        new BlockCpsiSnapshotWitness("aer (Trotter noiseless)", "heisenberg", 0, 0.0378, 0.0642),
        new BlockCpsiSnapshotWitness("aer (Trotter noiseless)", "heisenberg", 1, 0.0280, 0.0642),
        new BlockCpsiSnapshotWitness("aer (Trotter noiseless)", "truly",      0, 0.0192, 0.1064),
        new BlockCpsiSnapshotWitness("aer (Trotter noiseless)", "truly",      1, 0.0219, 0.1064),
        new BlockCpsiSnapshotWitness("aer (Trotter noiseless)", "soft",       0, 0.0279, 0.1838),
        new BlockCpsiSnapshotWitness("aer (Trotter noiseless)", "soft",       1, 0.0596, 0.0499),
        new BlockCpsiSnapshotWitness("aer (Trotter noiseless)", "hard",       0, 0.0879, 0.1707),
        new BlockCpsiSnapshotWitness("aer (Trotter noiseless)", "hard",       1, 0.0422, 0.0785),
        // ibm_marrakesh
        new BlockCpsiSnapshotWitness("ibm_marrakesh", "heisenberg", 0, 0.0412, 0.0642),
        new BlockCpsiSnapshotWitness("ibm_marrakesh", "heisenberg", 1, 0.0388, 0.0642),
        new BlockCpsiSnapshotWitness("ibm_marrakesh", "truly",      0, 0.0386, 0.1064),
        new BlockCpsiSnapshotWitness("ibm_marrakesh", "truly",      1, 0.0217, 0.1064),
        new BlockCpsiSnapshotWitness("ibm_marrakesh", "soft",       0, 0.0310, 0.1838),
        new BlockCpsiSnapshotWitness("ibm_marrakesh", "soft",       1, 0.0753, 0.0499),
        new BlockCpsiSnapshotWitness("ibm_marrakesh", "hard",       0, 0.1099, 0.1707),
        new BlockCpsiSnapshotWitness("ibm_marrakesh", "hard",       1, 0.0543, 0.0785),
        // ibm_kingston
        new BlockCpsiSnapshotWitness("ibm_kingston", "heisenberg", 0, 0.0453, 0.0642),
        new BlockCpsiSnapshotWitness("ibm_kingston", "heisenberg", 1, 0.0502, 0.0642),
        new BlockCpsiSnapshotWitness("ibm_kingston", "truly",      0, 0.0328, 0.1064),
        new BlockCpsiSnapshotWitness("ibm_kingston", "truly",      1, 0.0329, 0.1064),
        new BlockCpsiSnapshotWitness("ibm_kingston", "soft",       0, 0.0409, 0.1838),
        new BlockCpsiSnapshotWitness("ibm_kingston", "soft",       1, 0.0920, 0.0499),
        new BlockCpsiSnapshotWitness("ibm_kingston", "hard",       0, 0.1134, 0.1707),
        new BlockCpsiSnapshotWitness("ibm_kingston", "hard",       1, 0.0460, 0.0785),
        // ibm_fez
        new BlockCpsiSnapshotWitness("ibm_fez", "heisenberg", 0, 0.0140, 0.0642),
        new BlockCpsiSnapshotWitness("ibm_fez", "heisenberg", 1, 0.0084, 0.0642),
        new BlockCpsiSnapshotWitness("ibm_fez", "truly",      0, 0.0060, 0.1064),
        new BlockCpsiSnapshotWitness("ibm_fez", "truly",      1, 0.0166, 0.1064),
        new BlockCpsiSnapshotWitness("ibm_fez", "soft",       0, 0.0132, 0.1838),
        new BlockCpsiSnapshotWitness("ibm_fez", "soft",       1, 0.0303, 0.0499),
        new BlockCpsiSnapshotWitness("ibm_fez", "hard",       0, 0.0845, 0.1707),
        new BlockCpsiSnapshotWitness("ibm_fez", "hard",       1, 0.0347, 0.0785),
    };
}

/// <summary>One row of the IBM C_block snapshot table: a frozen empirical readout for
/// (backend, scenario, block n) at the 2026-04-26 |+−+⟩ run.</summary>
/// <param name="Backend">Source backend label (e.g. "ibm_marrakesh").</param>
/// <param name="Scenario">Hamiltonian scenario: "heisenberg", "truly", "soft", or "hard".</param>
/// <param name="BlockN">Block index n; the block is (popcount-n, popcount-(n+1)). Only n ∈ {0, 1} for the 2-qubit reduced ρ.</param>
/// <param name="CBlockMeasured">Measured C_block content from the snapshot JSON's 16-Pauli reconstruction.</param>
/// <param name="CBlockIdealContinuous">Reference: C_block from the noiseless continuous-time ideal at the same parameters.</param>
public sealed record BlockCpsiSnapshotWitness(
    string Backend,
    string Scenario,
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
            yield return new InspectableNode("scenario", summary: Scenario);
            yield return new InspectableNode("block (n, n+1)", summary: $"({BlockN}, {BlockN + 1})");
            yield return InspectableNode.RealScalar("CBlockMeasured", CBlockMeasured, "F4");
            yield return InspectableNode.RealScalar("CBlockIdealContinuous", CBlockIdealContinuous, "F4");
            yield return InspectableNode.RealScalar("FractionOfQuarter", FractionOfQuarter, "F3");
            yield return InspectableNode.RealScalar("DeltaFromIdeal", DeltaFromIdeal, "F4");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
