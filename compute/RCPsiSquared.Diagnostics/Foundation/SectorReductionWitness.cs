using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
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

    public string DisplayName => $"SectorReductionWitness (N={N}, {Topology}, gamma={Gamma.ToString("0.###", Inv)})";
    public string Summary =>
        "the birth-canal boundary as a Liouville sector reduction: the slowest mode is the |1-exc><vac| " +
        "(0,1) block (N-dim), reproducing the full-4^N witness at N=5 and running past it; at N>=6 the " +
        "global slowest crosses to the {0,2}-coherence (the coherence-horizon mode). Reuses the per-sector " +
        "builder. No V-Effect identity, no aromaticity thesis (both open).";
    public IEnumerable<IInspectable> Children => Array.Empty<IInspectable>();   // nodes in later tasks
    public InspectablePayload Payload => InspectablePayload.Empty;
}
