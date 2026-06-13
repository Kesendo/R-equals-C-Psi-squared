using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The birth-canal boundary recomputed as a Liouville SECTOR reduction, live
/// (<c>inspect --root reduction</c>). The slowest-rate mode of the boundary is the ODD,
/// number-changing |1-excitation><vacuum| coherence, which lives in the (PCol=0, PRow=1) sector -
/// an N-dim block. This witness builds that block via <see cref="PerBlockLiouvillianBuilder"/>,
/// validates it bit-for-bit against <see cref="PostEpFlowField"/> at N=5, and goes PAST N=5 (the
/// dense witness caps at N=6) to surface the junction: at N>=6 the global slowest can cross to the
/// EVEN {0,2}-coherence in the 2-excitation density block (the coherence-horizon mode; see the
/// <c>birth_canal_horizon_junction</c> arc). Reuses the per-sector machinery, no full 4^N.
///
/// <para>Convention: H_unit*Q (coeff 1, single-particle off-diag 2) = ChainSystem(N, J=2Q).
/// BuildHamiltonian() (XY = (J/2)(XX+YY), off-diag J); the (0,1) block's -2*gamma_l-per-bit
/// dephasing matches PostEpFlowField, pinned by the validation test.</para>
///
/// <para>Honest scope: this banks the N=5 birth-canal result and the analytic flat-gamma blindness.
/// It does NOT assert the V-Effect (w=N/2) self-pair IS the {0,2}-coherence (different
/// decompositions, identity open), and carries no aromaticity / 4n-vs-4n+2 thesis (open, the C8
/// case breaks the naive reading).</para></summary>
public sealed class SectorReductionWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double KernelTol = 1e-7;

    public int N { get; }
    public double Gamma { get; }
    public TopologyKind Topology { get; }

    public SectorReductionWitness(int n = 5, double gamma = 0.5, TopologyKind topology = TopologyKind.Chain)
    {
        if (n < 2 || n > 10) throw new ArgumentOutOfRangeException(nameof(n), n, "N in 2..10 for the sector reduction");
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "gamma must be > 0");
        N = n; Gamma = gamma; Topology = topology;
    }

    /// <summary>H = Q*H_unit = Q*Sum_b(X_bX_{b+1}+Y_bY_{b+1}). ChainSystem XY builds (J/2)(XX+YY),
    /// so J=2Q gives coeff Q. The GammaZero arg is unused for H (dephasing lives in the dissipator);
    /// pass any &gt;0.</summary>
    private static ComplexMatrix QHUnit(int n, double q, TopologyKind topology) =>
        new ChainSystem(n, 2.0 * q, 1.0, HamiltonianType.XY, topology).BuildHamiltonian();

    /// <summary>The flat indices of the (PCol, PRow) sector.</summary>
    private static int[] SectorFlat(int n, int pCol, int pRow)
    {
        var decomp = JointPopcountSectorBuilder.Build(n);
        var s = decomp.SectorRanges.First(r => r.PCol == pCol && r.PRow == pRow);
        var flat = new int[s.Size];
        for (int k = 0; k < s.Size; k++) flat[k] = decomp.Permutation[s.Offset + k];
        return flat;
    }

    /// <summary>The slowest non-kernel rate of the (PCol,PRow) sector block: -max{ Re lambda :
    /// |lambda| &gt; tol }.</summary>
    public static double SectorSlowest(int n, double q, IReadOnlyList<double> gammaProfile,
        int pCol, int pRow, TopologyKind topology)
    {
        var H = QHUnit(n, q, topology);
        var block = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaProfile, SectorFlat(n, pCol, pRow));
        var ev = block.Evd().EigenValues;
        double maxRe = double.NegativeInfinity;
        foreach (var z in ev)
            if (z.Magnitude > KernelTol && z.Real > maxRe) maxRe = z.Real;
        return double.IsNegativeInfinity(maxRe) ? 0.0 : -maxRe;
    }

    /// <summary>The |1-exc><vac| (PCol=0, PRow=1) block slowest rate - the birth-canal boundary's
    /// mode, reproducing PostEpFlowField at N=5 (pinned by test).</summary>
    public static double VacBlockSlowest(int n, double q, IReadOnlyList<double> gammaProfile, TopologyKind topology) =>
        SectorSlowest(n, q, gammaProfile, pCol: 0, pRow: 1, topology);

    private double Deviation(int n, double[] profile, TopologyKind topo) =>
        VacBlockSlowest(n, 1000.0, profile, topo) - VacBlockSlowest(n, 1.5, profile, topo);

    private InspectableNode TheVacReductionNode()
    {
        var kids = new List<IInspectable>();
        // flat-gamma blindness (uniform): rate == 2*gamma, Q-invariant (analytic).
        var uni = Enumerable.Repeat(Gamma, N).ToArray();
        double uLo = VacBlockSlowest(N, 1.5, uni, Topology), uHi = VacBlockSlowest(N, 1000.0, uni, Topology);
        kids.Add(new InspectableNode("flat-gamma blindness",
            summary: $"uniform gamma={Gamma.ToString("0.###", Inv)}: rate {uLo.ToString("0.0000", Inv)} (Q=1.5) -> " +
                     $"{uHi.ToString("0.0000", Inv)} (Q=1000) = 2*gamma, Q-invariant (analytic: -iQh anti-Hermitian)"));
        // N=5 validation breadcrumb (the canal anchor reproduces PostEpFlowField; see test).
        if (N == 5)
        {
            var canal = new[] { 0.25, 1.5, 1.5, 1.5, 0.25 };
            kids.Add(new InspectableNode("N=5 validation (canal anchor)",
                summary: $"(0,1)-block rate {VacBlockSlowest(5, 1.5, canal, Topology).ToString("0.0000", Inv)} -> " +
                         $"{VacBlockSlowest(5, 1000.0, canal, Topology).ToString("0.0000", Inv)} == PostEpFlowField " +
                         "(bit-identical, see SectorReductionWitnessTests); the N-dim block IS the full boundary at N=5"));
        }
        // past N=5: the block builds at this N regardless (the dense witness cannot beyond N=6).
        var deep = DeepEdge(N);
        kids.Add(InspectableNode.RealScalar($"deep-edge deviation (N={N})", Deviation(N, deep, Topology), "0.000000"));
        return new InspectableNode("the |1-exc><vac| reduction",
            summary: $"the (0,1) sector block (N={N}, {Topology}): an N-dim invariant block whose slowest mode is " +
                     "the birth-canal boundary. Reproduces the full 4^N witness at N=5, runs past it here.",
            children: kids);
    }

    /// <summary>The general-N "deep-edge" profile: edges depressed to 0.25, the rest flat to sum N
    /// (the sibling of the N=5 flat-bulk-edge canal anchor).</summary>
    private static double[] DeepEdge(int n)
    {
        double edge = 0.25, rest = (n - 2 * edge) / (n - 2);
        var p = new double[n];
        for (int i = 0; i < n; i++) p[i] = (i == 0 || i == n - 1) ? edge : rest;
        return p;
    }

    /// <summary>The slowest non-kernel mode of the (pc,pr) density block, with its n_diff histogram
    /// (the {0,2} signature) and phase rigidity (the EP detector).</summary>
    private (double Rate, double Rigidity, IReadOnlyDictionary<int,double> Hist) DensityMode(int n, double q, double[] profile, int pc, int pr)
    {
        var H = QHUnit(n, q, Topology);
        var block = PerBlockLiouvillianBuilder.BuildBlockZ(H, profile, SectorFlat(n, pc, pr));
        var modes = PhaseRigidity.Compute(block).Where(m => m.Lambda.Magnitude > KernelTol).ToList();
        if (modes.Count == 0) return (0.0, 1.0, new Dictionary<int,double>());
        var slow = modes.OrderByDescending(m => m.Lambda.Real).First();
        var (_, hist) = LiouvilleOperatorContent.NDiffHistogram(slow.Right, n);
        return (-slow.Lambda.Real, slow.Rigidity, hist);
    }

    private InspectableNode TheJunctionNode()
    {
        if (N < 6 || N > 8)
            return new InspectableNode("the {0,2} junction",
                summary: "the odd<->even crossing is visible at N=6..8 (the (2,2) density block, " +
                         "C(N,2)^2-dim); set N in 6..8 to read it. At N=5 the (0,1) mode always wins.");
        var deep = DeepEdge(N);
        double vacLo = VacBlockSlowest(N, 1.5, deep, Topology);
        var (densLo, rigLo, histLo) = DensityMode(N, 1.5, deep, 2, 2);
        bool crosses = densLo < vacLo;
        string h0 = histLo.GetValueOrDefault(0).ToString("0.00", Inv), h2 = histLo.GetValueOrDefault(2).ToString("0.00", Inv);
        return new InspectableNode("the {0,2} junction",
            summary: $"deep-edge, N={N}, Q=1.5: (0,1) odd rate {vacLo.ToString("0.000", Inv)} vs (2,2) density rate " +
                     $"{densLo.ToString("0.000", Inv)} -> {(crosses ? "the {0,2} density mode WINS (the crossing): " : "(0,1) still wins: ")}" +
                     $"n_diff hist {{0:{h0}, 2:{h2}}}, rigidity {rigLo.ToString("0.000", Inv)}. This is where birth_canal " +
                     "meets coherence_horizon (the {0,2}-coherence = its EP mode). Identity with the V-Effect (w=N/2) " +
                     "self-pair is co-located, OPEN, not claimed.");
    }

    public string DisplayName => $"SectorReductionWitness (N={N}, {Topology}, gamma={Gamma.ToString("0.###", Inv)})";
    public string Summary =>
        "the birth-canal boundary as a Liouville sector reduction: the slowest mode is the |1-exc><vac| " +
        "(0,1) block (N-dim), reproducing the full-4^N witness at N=5 and running past it; at N>=6 the " +
        "global slowest crosses to the {0,2}-coherence (the coherence-horizon mode). Reuses the per-sector " +
        "builder. No V-Effect identity, no aromaticity thesis (both open).";
    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return TheVacReductionNode();
            yield return TheJunctionNode();
            // chain/ring in a later task
        }
    }
    public InspectablePayload Payload => InspectablePayload.Empty;
}
