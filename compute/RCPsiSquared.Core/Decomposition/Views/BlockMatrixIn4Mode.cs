using RCPsiSquared.Core.Inspection;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Decomposition.Views;

/// <summary>F86 inspection wrapper for any 4×4 matrix in the 4-mode basis decomposed into
/// the three structural sub-blocks: channel-uniform (rows/cols [0:2]), SVD-top (rows/cols
/// [2:4]), and the cross-block (rows [0:2], cols [2:4] — the bond-position-dependent
/// coupling that splits Endpoint vs Interior).
///
/// <para>Used by <c>BondCouplingIn4Mode</c> for per-bond M_h, and by the wrapper for
/// M_h_total_eff and S_kernel_eff. Computes Frobenius norms per block and the cross-block
/// fraction — the latter is the F86 algebraic fingerprint that decides Endpoint vs Interior
/// shape.</para>
/// </summary>
public sealed class BlockMatrixIn4Mode : IInspectable
{
    public string Label { get; }
    public ComplexMatrix Matrix { get; }
    public Block2x2 ChannelUniformBlock { get; }
    public Block2x2 SvdTopBlock { get; }
    public Block2x2 CrossBlock { get; }
    public Block2x2 CrossBlockReverse { get; }    // rows [2:4] × cols [0:2]

    private readonly Lazy<double> _frobenius;

    public double Frobenius => _frobenius.Value;
    public double CrossBlockFrobenius => CrossBlock.Frobenius;
    public double CrossFrobeniusFraction => Frobenius > 0
        ? (CrossBlock.Frobenius * CrossBlock.Frobenius +
           CrossBlockReverse.Frobenius * CrossBlockReverse.Frobenius)
          / (Frobenius * Frobenius)
        : 0;

    public BlockMatrixIn4Mode(string label, ComplexMatrix matrix)
    {
        if (matrix.RowCount != 4 || matrix.ColumnCount != 4)
            throw new ArgumentException($"BlockMatrixIn4Mode expects 4×4; got {matrix.RowCount}×{matrix.ColumnCount}.");
        Label = label;
        Matrix = matrix;
        ChannelUniformBlock = Block2x2.Extract(matrix, 0, 0, "channel-uniform 2×2",
            FourModeNames.ChannelUniform, FourModeNames.ChannelUniform);
        SvdTopBlock = Block2x2.Extract(matrix, 2, 2, "SVD-top 2×2",
            FourModeNames.SvdTop, FourModeNames.SvdTop);
        CrossBlock = Block2x2.Extract(matrix, 0, 2, "cross (CU → SVD-top) 2×2",
            FourModeNames.ChannelUniform, FourModeNames.SvdTop);
        CrossBlockReverse = Block2x2.Extract(matrix, 2, 0, "cross (SVD-top → CU) 2×2",
            FourModeNames.SvdTop, FourModeNames.ChannelUniform);
        _frobenius = new Lazy<double>(() => matrix.FrobeniusNorm());
    }

    public string DisplayName => Label;
    public string Summary =>
        $"‖·‖_F = {Frobenius:F4}, cross fraction = {CrossFrobeniusFraction:P1}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return ChannelUniformBlock;
            yield return SvdTopBlock;
            yield return CrossBlock;
            yield return CrossBlockReverse;
            yield return InspectableNode.RealScalar("‖·‖_F (full)", Frobenius, "F6");
            yield return InspectableNode.RealScalar("cross-block fraction", CrossFrobeniusFraction, "P3");
        }
    }

    public InspectablePayload Payload =>
        new InspectablePayload.MatrixView(Label, Matrix, FourModeNames.All, FourModeNames.All);
}
