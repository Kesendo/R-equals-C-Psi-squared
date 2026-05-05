using RCPsiSquared.Core.ChainSystems;

namespace RCPsiSquared.Core.Calibration;

/// <summary>A calibration-selected physical-qubit chain on an IBM-style backend, ready
/// to be lowered into a logical <see cref="ChainSystem"/> for analysis or circuit
/// construction. Bridges <see cref="IbmCalibration"/> to the existing logical framework.
///
/// <para><see cref="QubitIds"/> are the physical IDs in chain order (the same list a
/// circuit builder would use for layout), with <see cref="Length"/> = QubitIds.Count.
/// <see cref="Score"/> is the path score from <see cref="IbmCalibration.BestChain"/>:
/// per-qubit <see cref="IbmCalibration.ScoreQubit"/> sum plus −1000·CZ-error per bond.</para>
///
/// <para>On the 2026-04-25 Marrakesh CSV the canonical 3-chain is [4, 3, 2] (score ≈ 867)
/// and the canonical 5-chain is [1, 2, 3, 4, 5] (score ≈ 1247); the best 3-chain
/// outscores the documented [48, 49, 50] path used in the soft_break replication
/// (score ≈ 682), which itself outscores the original [0, 1, 2] run (score ≈ 597) by
/// ~14%. See <c>simulations/_f88_lens_ibm_data_extended.py</c> for the corresponding
/// 23× cleaner truly-baseline at the F87-trichotomy state-level reading.</para>
/// </summary>
public sealed record CalibrationChain(
    IReadOnlyList<int> QubitIds,
    double Score,
    IReadOnlyList<QubitData> Qubits)
{
    public int Length => QubitIds.Count;

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
    /// derived from each physical qubit's T2: γ_l = 1 / (2·T2_l) (microsecond inverse).
    /// Returns the chain plus the γ array for callers that want to use
    /// <c>PauliDephasingDissipator.Build(H, gammaPerSite)</c> with non-uniform γ.</summary>
    public (ChainSystem Chain, double[] GammaPerSite) ToChainSystemWithT2Gamma(
        double J,
        HamiltonianType hType = HamiltonianType.XY)
    {
        double[] gammas = new double[Length];
        for (int i = 0; i < Length; i++)
        {
            double t2 = Qubits[i].T2Us;
            gammas[i] = t2 > 0 ? 1.0 / (2.0 * t2) : 0.0;
        }
        double meanGamma = gammas.Average();
        var chain = new ChainSystem(N: Length, J: J, GammaZero: meanGamma, HType: hType, Topology: TopologyKind.Chain);
        return (chain, gammas);
    }
}
