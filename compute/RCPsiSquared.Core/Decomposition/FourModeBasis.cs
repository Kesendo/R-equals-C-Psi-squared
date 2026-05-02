using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Decomposition;

/// <summary>The minimal orthonormal 4-mode basis in which the F86 K_CC_pr observable
/// dynamics can be modelled: {|c_1⟩, |c_3⟩, |u_0⟩, |v_0⟩}.
///
/// <para>|c_1⟩, |c_3⟩ are the channel-uniform vectors (where the Dicke probe lives).
/// |u_0⟩, |v_0⟩ are the SVD-top inter-channel coupling vectors (the EP-partner modes).
/// At c=2 these four vectors are mutually orthonormal (verified numerically across
/// N=5..8 in `_eq022_b1_step_i_svd_inter_channel.py`).</para>
///
/// <para>The 2-level reduction of the proof's Statement 1 lives in span{|u_0⟩, |v_0⟩};
/// the probe lives in span{|c_1⟩, |c_3⟩}; the K observable couples them via
/// ∂L/∂J = M_H_per_bond[b] which has bond-position-dependent matrix elements between
/// the two 2D subspaces. The bond-class universality of HWHM_left/Q_peak is the
/// fingerprint of these cross-coupling patterns.</para>
/// </summary>
public sealed class FourModeBasis
{
    public CoherenceBlock Block { get; }
    public int Hd1 { get; }
    public int Hd2 { get; }
    public ComplexMatrix BasisMatrix { get; }   // Mtot × 4 columns: c_1, c_3, u_0, v_0

    /// <summary>Maximum |off-diagonal| of the Gram matrix B^† B (should be ~0 for orthonormality).</summary>
    public double OffOrthonormalityResidual { get; }

    private FourModeBasis(CoherenceBlock block, int hd1, int hd2, ComplexMatrix basis, double residual)
    {
        Block = block;
        Hd1 = hd1;
        Hd2 = hd2;
        BasisMatrix = basis;
        OffOrthonormalityResidual = residual;
    }

    public static FourModeBasis Build(CoherenceBlock block, int hd1 = 1, int hd2 = 3)
    {
        var hdChannel = HdChannelBasis.Build(block);
        int k1 = -1, k2 = -1;
        for (int i = 0; i < hdChannel.HammingDistances.Count; i++)
        {
            if (hdChannel.HammingDistances[i] == hd1) k1 = i;
            if (hdChannel.HammingDistances[i] == hd2) k2 = i;
        }
        if (k1 < 0 || k2 < 0)
            throw new ArgumentException($"HD={hd1} or HD={hd2} not present in chromaticity={block.C}.");

        var c1 = hdChannel.P.Column(k1);
        var c3 = hdChannel.P.Column(k2);

        var svd = InterChannelSvd.Build(block, hd1, hd2);
        var u0 = svd.U0InFullBlock;
        var v0 = svd.V0InFullBlock;

        var basis = ComplexMatrix.Build.Dense(block.Basis.MTotal, 4);
        basis.SetColumn(0, c1);
        basis.SetColumn(1, c3);
        basis.SetColumn(2, u0);
        basis.SetColumn(3, v0);

        var gram = basis.ConjugateTranspose() * basis;
        double residual = 0;
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                if (i == j) continue;
                double mag = gram[i, j].Magnitude;
                if (mag > residual) residual = mag;
            }
        }

        return new FourModeBasis(block, hd1, hd2, basis, residual);
    }

    /// <summary>Project a full-block matrix onto the 4-mode basis: B^† · M · B (4×4).</summary>
    public ComplexMatrix Project(ComplexMatrix m) =>
        BasisMatrix.ConjugateTranspose() * m * BasisMatrix;

    /// <summary>Project a full-block vector onto the 4-mode basis: B^† · v (4-vector).</summary>
    public MathNet.Numerics.LinearAlgebra.Vector<Complex> Project(MathNet.Numerics.LinearAlgebra.Vector<Complex> v) =>
        BasisMatrix.ConjugateTranspose() * v;
}
