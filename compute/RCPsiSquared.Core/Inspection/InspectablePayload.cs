using System.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Inspection;

/// <summary>Tagged union describing the leaf payload of an <see cref="IInspectable"/> node.
/// The Object Manager renders each variant with the appropriate widget (label / stem plot /
/// heatmap / line plot).
///
/// <para>Use <see cref="None"/> for pure containers (an inner node whose value is the union
/// of its children, with no leaf data of its own).</para>
/// </summary>
public abstract record InspectablePayload
{
    /// <summary>No leaf payload — this node is a pure container.</summary>
    public sealed record None : InspectablePayload;

    /// <summary>A real-valued scalar (Frobenius norm, peak position, fraction, …).</summary>
    public sealed record Real(string Label, double Value, string? Format = null) : InspectablePayload;

    /// <summary>A complex-valued scalar.</summary>
    public sealed record Scalar(string Label, Complex Value) : InspectablePayload;

    /// <summary>A complex-valued vector — rendered as a stem plot or magnitude bar chart.</summary>
    public sealed record Vector(string Label, ComplexVector Values, IReadOnlyList<string>? ComponentLabels = null)
        : InspectablePayload;

    /// <summary>A complex-valued matrix — rendered as a magnitude / phase heatmap.</summary>
    public sealed record MatrixView(string Label, ComplexMatrix Values,
        IReadOnlyList<string>? RowLabels = null, IReadOnlyList<string>? ColumnLabels = null)
        : InspectablePayload;

    /// <summary>A 1D curve (Q vs K, t vs ρ_diag, …) — rendered as a line plot.</summary>
    public sealed record Curve(string Label, IReadOnlyList<double> X, IReadOnlyList<double> Y,
        string? XLabel = null, string? YLabel = null) : InspectablePayload;

    public static readonly InspectablePayload Empty = new None();
}
