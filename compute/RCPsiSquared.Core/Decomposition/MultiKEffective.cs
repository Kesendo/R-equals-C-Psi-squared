using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Probes;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Decomposition;

/// <summary>F86 Item 2: full effective model for chromaticity c ≥ 2 — extension of
/// <see cref="FourModeEffective"/> via <see cref="MultiKBasis"/> (3c−2 orthonormal modes).
/// All operators (D, M_h_per_bond[b], M_h_total, probe, S_kernel) projected onto the
/// 3c−2-dim multi-k basis. For c=2 this reduces to the 4-mode effective.
/// </summary>
public sealed class MultiKEffective
{
    public CoherenceBlock Block { get; }
    public MultiKBasis Basis { get; }
    public ComplexMatrix DEff { get; }
    public IReadOnlyList<ComplexMatrix> MhPerBondEff { get; }
    public ComplexMatrix MhTotalEff { get; }
    public ComplexVector ProbeEff { get; }
    public ComplexMatrix SKernelEff { get; }

    public int Dim => Basis.TotalModes;

    private MultiKEffective(CoherenceBlock block, MultiKBasis basis,
        ComplexMatrix dEff, IReadOnlyList<ComplexMatrix> mhPerBond, ComplexMatrix mhTotal,
        ComplexVector probeEff, ComplexMatrix sKernelEff)
    {
        Block = block;
        Basis = basis;
        DEff = dEff;
        MhPerBondEff = mhPerBond;
        MhTotalEff = mhTotal;
        ProbeEff = probeEff;
        SKernelEff = sKernelEff;
    }

    public static MultiKEffective Build(CoherenceBlock block)
    {
        var basis = MultiKBasis.Build(block);
        var decomp = block.Decomposition;
        var probe = DickeBlockProbe.Build(block);
        var sKernel = SpatialSumKernel.Build(block);

        var dEff = basis.Project(decomp.D);
        var mhPerBondEff = decomp.MhPerBond.Select(m => basis.Project(m)).ToList();
        var mhTotalEff = basis.Project(decomp.MhTotal);
        var probeEff = basis.Project(probe);
        var sKernelEff = basis.Project(sKernel);

        return new MultiKEffective(block, basis, dEff, mhPerBondEff, mhTotalEff, probeEff, sKernelEff);
    }

    /// <summary>L_eff(Q) = D_eff + (Q·γ₀)·M_h_total_eff. Uniform-J across all bonds.</summary>
    public ComplexMatrix LEffAtQ(double q)
    {
        double j = q * Block.GammaZero;
        return DEff + (Complex)j * MhTotalEff;
    }
}
