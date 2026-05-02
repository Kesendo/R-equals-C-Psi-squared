using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Decomposition;

/// <summary>SVD of the inter-channel coupling block of M_H_total.
///
/// <para>For HD₁ = 2k−1 and HD₂ = 2k+1 channels of a coherence block, V_inter =
/// P_{HD₁}^† · M_H_total · P_{HD₂} is the H matrix elements between the two HD subspaces.
/// Its top singular vectors |u_0⟩ ∈ HD₁ and |v_0⟩ ∈ HD₂ (lifted to the full block) are the
/// modes that maximally couple under H — the heuristic "EP-partner" pair from the F86
/// 4-mode picture.</para>
///
/// <para>Empirical findings (c=2 chains N=5..8): σ_0 → 2√2 asymptotically; the probe (Dicke
/// state) is orthogonal to |u_0⟩, |v_0⟩; the 2-level same-sign-imaginary off-diagonal form
/// of the proof is realised here as a ±σ_0 real off-diagonal pattern (different phase
/// convention).</para>
/// </summary>
public sealed class InterChannelSvd
{
    public CoherenceBlock Block { get; }
    public int HammingDistance1 { get; }
    public int HammingDistance2 { get; }
    public IReadOnlyList<double> SingularValues { get; }
    public ComplexVector U0InFullBlock { get; }
    public ComplexVector V0InFullBlock { get; }

    public double Sigma0 => SingularValues[0];

    private InterChannelSvd(CoherenceBlock block, int hd1, int hd2,
        double[] singular, ComplexVector u0, ComplexVector v0)
    {
        Block = block;
        HammingDistance1 = hd1;
        HammingDistance2 = hd2;
        SingularValues = singular;
        U0InFullBlock = u0;
        V0InFullBlock = v0;
    }

    public static InterChannelSvd Build(CoherenceBlock block, int hd1, int hd2)
    {
        var pHd1 = HdSubspaceProjector.Build(block, hd1);
        var pHd2 = HdSubspaceProjector.Build(block, hd2);

        var vInter = pHd1.ConjugateTranspose() * block.Decomposition.MhTotal * pHd2; // (n_hd1 × n_hd2)

        var svd = vInter.Svd();
        var singular = svd.S.Select(z => z.Real).ToArray();
        // Top singular vectors
        var uTop = svd.U.Column(0);  // length n_hd1, in HD1 reduced basis
        var vTop = svd.VT.ConjugateTranspose().Column(0);  // length n_hd2, in HD2 reduced basis

        // Lift to full block basis
        var u0Full = pHd1 * uTop;
        var v0Full = pHd2 * vTop;

        return new InterChannelSvd(block, hd1, hd2, singular, u0Full, v0Full);
    }
}
