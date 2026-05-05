using RCPsiSquared.Core.ChainSystems;

namespace RCPsiSquared.Core.Calibration;

/// <summary>A calibration-selected physical-qubit chain on an IBM-style backend, ready
/// to be lowered into a logical <see cref="ChainSystem"/> for analysis or circuit
/// construction. Bridges <see cref="IbmCalibration"/> to the existing logical framework.
///
/// <para><see cref="Qubits"/> are the per-site <see cref="QubitData"/> records in chain
/// order; <see cref="QubitIds"/> is the convenience projection
/// (<c>Qubits[i].Qubit</c> for each i) for circuit-builder layout. <see cref="Score"/>
/// is the path score from <see cref="IbmCalibration.BestChain"/>: per-qubit
/// <see cref="IbmCalibration.ScoreQubit"/> sum plus −1000·CZ-error per bond.</para>
///
/// <para>On the 2026-04-25 Marrakesh CSV the best 3-chain is [4, 3, 2] (score ≈ 867)
/// and the best 5-chain is [1, 2, 3, 4, 5] (score ≈ 1247). The same scoring ranks the
/// soft_break path [48, 49, 50] (score ≈ 682) above the framework_snapshots path
/// [0, 1, 2] (score ≈ 597) by ~14%; the corresponding state-level truly-baseline
/// downstream of F88-Lens is 23× cleaner on [48, 49, 50] (see
/// <c>project_f88_lens_ibm_marrakesh.md</c>). The score gap and the F88 ratio are
/// related but not equal: better calibration → cleaner truly-baseline, with the F88
/// reading amplifying the per-qubit difference at the state level.</para>
/// </summary>
public sealed record CalibrationChain(
    double Score,
    IReadOnlyList<QubitData> Qubits)
{
    public IReadOnlyList<int> QubitIds { get; } = Qubits.Select(q => q.Qubit).ToArray();
    public int Length => Qubits.Count;

    /// <summary>Lower this physical-qubit chain to a logical <see cref="ChainSystem"/>
    /// of size <see cref="Length"/>. The logical bonds (0-1, 1-2, …) act as Chain
    /// topology in qubit-index space; the physical mapping back to hardware lives in
    /// <see cref="QubitIds"/>.</summary>
    public ChainSystem ToChainSystem(
        double J,
        double gammaZero,
        HamiltonianType hType = HamiltonianType.XY) =>
        new(N: Length, J: J, GammaZero: gammaZero, HType: hType, Topology: TopologyKind.Chain);

    /// <summary>Lower this chain to a <see cref="ChainSystem"/> with a per-site γ list
    /// derived from each physical qubit's T2: γ_l = 1 / (2·T2_l) (units μs⁻¹, since T2
    /// is read in μs from the IBM CSV; pass <paramref name="J"/> in μs⁻¹ for a
    /// dimensionally consistent Q = J/γ). The returned <see cref="ChainSystem"/>
    /// carries the **mean** γ as <see cref="ChainSystem.GammaZero"/> for convenience
    /// (lossy when γ varies across sites); the per-site array is returned alongside
    /// for callers that want non-uniform γ via
    /// <c>PauliDephasingDissipator.Build(H, gammaPerSite)</c>.</summary>
    public (ChainSystem Chain, double[] GammaPerSite) ToChainSystemWithT2Gamma(
        double J,
        HamiltonianType hType = HamiltonianType.XY)
    {
        double[] gammas = new double[Length];
        double sum = 0.0;
        for (int i = 0; i < Length; i++)
        {
            double t2 = Qubits[i].T2Us;
            double g = t2 > 0 ? 1.0 / (2.0 * t2) : 0.0;
            gammas[i] = g;
            sum += g;
        }
        double meanGamma = Length > 0 ? sum / Length : 0.0;
        var chain = new ChainSystem(N: Length, J: J, GammaZero: meanGamma, HType: hType, Topology: TopologyKind.Chain);
        return (chain, gammas);
    }
}
