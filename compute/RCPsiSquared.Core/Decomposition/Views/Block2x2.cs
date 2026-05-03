using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Decomposition.Views;

/// <summary>A 2×2 complex sub-block extracted from a larger matrix, with computation:
/// Frobenius norm, max element magnitude, off-diagonal symmetry. Used by the F86 4-mode
/// view types to expose channel-uniform / SVD-top / cross-block decompositions.
/// </summary>
public sealed class Block2x2 : IInspectable
{
    public string Label { get; }
    public ComplexMatrix Matrix { get; }
    public IReadOnlyList<string> RowLabels { get; }
    public IReadOnlyList<string> ColumnLabels { get; }

    private readonly Lazy<double> _frobenius;
    private readonly Lazy<double> _maxMagnitude;
    private readonly Lazy<double> _offDiagonalAsymmetry;

    public double Frobenius => _frobenius.Value;
    public double MaxElementMagnitude => _maxMagnitude.Value;

    /// <summary>|M[0,1]| − |M[1,0]| as a measure of off-diagonal asymmetry. Zero for normal
    /// 2×2 sub-blocks (e.g. Hermitian or anti-Hermitian).</summary>
    public double OffDiagonalAsymmetry => _offDiagonalAsymmetry.Value;

    public Block2x2(string label, ComplexMatrix matrix, IReadOnlyList<string>? rowLabels = null,
        IReadOnlyList<string>? columnLabels = null)
    {
        if (matrix.RowCount != 2 || matrix.ColumnCount != 2)
            throw new ArgumentException($"Block2x2 expects 2×2 matrix; got {matrix.RowCount}×{matrix.ColumnCount}.");
        Label = label;
        Matrix = matrix;
        RowLabels = rowLabels ?? Array.Empty<string>();
        ColumnLabels = columnLabels ?? Array.Empty<string>();
        _frobenius = new Lazy<double>(() => Matrix.FrobeniusNorm());
        _maxMagnitude = new Lazy<double>(() => MatrixUtilities.MaxElementMagnitude(Matrix));
        _offDiagonalAsymmetry = new Lazy<double>(() =>
            Matrix[0, 1].Magnitude - Matrix[1, 0].Magnitude);
    }

    public string DisplayName => Label;
    public string Summary => $"‖·‖_F = {Frobenius:F4}, max |·| = {MaxElementMagnitude:F4}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("Frobenius", Frobenius, "F6");
            yield return InspectableNode.RealScalar("max |element|", MaxElementMagnitude, "F6");
            yield return InspectableNode.RealScalar("|M[0,1]| − |M[1,0]|", OffDiagonalAsymmetry, "F6");
        }
    }

    public InspectablePayload Payload =>
        new InspectablePayload.MatrixView(Label, Matrix,
            RowLabels.Count > 0 ? RowLabels : null,
            ColumnLabels.Count > 0 ? ColumnLabels : null);

    /// <summary>Extract a 2×2 sub-block from a larger matrix at (rowOffset, colOffset).</summary>
    public static Block2x2 Extract(ComplexMatrix source, int rowOffset, int colOffset, string label,
        IReadOnlyList<string>? rowLabels = null, IReadOnlyList<string>? colLabels = null) =>
        new(label, source.SubMatrix(rowOffset, 2, colOffset, 2), rowLabels, colLabels);
}
