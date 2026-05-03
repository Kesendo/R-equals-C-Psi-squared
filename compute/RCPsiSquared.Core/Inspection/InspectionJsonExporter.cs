using System.Numerics;
using System.Text.Encodings.Web;
using System.Text.Json;
using System.Text.Json.Serialization;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Inspection;

/// <summary>JSON exporter for any <see cref="IInspectable"/> tree. Handles arbitrary depth
/// (the tree is naturally 3D for per-bond × 4×4 structures, 4D when composed with a Q-sweep
/// or EVD, etc.) — recursive walk + variant-typed payload serialization.
///
/// <para>Every node becomes <c>{ displayName, summary, payload, children }</c>. Payload is
/// a tagged JSON object whose <c>kind</c> field is one of <c>none / real / scalar / vector /
/// matrix / curve</c>, with the variant-specific fields filled in. Complex numbers serialize
/// as <c>{ "re": …, "im": … }</c>.</para>
///
/// <para>The output is consumption-oriented (browse / diff / pipe to other tools). It is not
/// designed to be deserialized back into <see cref="IInspectable"/> objects — the round-trip
/// would lose the computational identity of each wrapper. Re-derive from the source instead.</para>
/// </summary>
public static class InspectionJsonExporter
{
    private static readonly JsonSerializerOptions DefaultOptions = new()
    {
        WriteIndented = true,
        DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull,
        // Keep ⟩, γ, ⊥, ‖·‖, ², · etc. readable in the export instead of Unicode-escaping them.
        // Output is consumption-oriented (browse / diff / pipe), not embedded into a webpage.
        Encoder = JavaScriptEncoder.UnsafeRelaxedJsonEscaping,
    };

    public static string ToJson(IInspectable root, JsonSerializerOptions? options = null) =>
        JsonSerializer.Serialize(ToDto(root), options ?? DefaultOptions);

    public static void WriteToFile(IInspectable root, string path, JsonSerializerOptions? options = null) =>
        File.WriteAllText(path, ToJson(root, options));

    public static InspectionNodeDto ToDto(IInspectable node) => new()
    {
        DisplayName = node.DisplayName,
        Summary = node.Summary,
        Payload = SerializePayload(node.Payload),
        Children = node.Children.Select(ToDto).ToArray(),
    };

    private static PayloadDto SerializePayload(InspectablePayload p) => p switch
    {
        InspectablePayload.None => new PayloadDto { Kind = "none" },
        InspectablePayload.Real r => new PayloadDto
        {
            Kind = "real",
            Label = r.Label,
            RealValue = r.Value,
            Format = r.Format,
        },
        InspectablePayload.Scalar s => new PayloadDto
        {
            Kind = "scalar",
            Label = s.Label,
            ComplexValue = ToCplx(s.Value),
        },
        InspectablePayload.Vector v => new PayloadDto
        {
            Kind = "vector",
            Label = v.Label,
            ComponentLabels = v.ComponentLabels?.ToArray(),
            VectorValues = ToCplxArray(v.Values),
        },
        InspectablePayload.MatrixView m => new PayloadDto
        {
            Kind = "matrix",
            Label = m.Label,
            RowLabels = m.RowLabels?.ToArray(),
            ColumnLabels = m.ColumnLabels?.ToArray(),
            MatrixValues = ToCplxMatrix(m.Values),
        },
        InspectablePayload.Curve c => new PayloadDto
        {
            Kind = "curve",
            Label = c.Label,
            XLabel = c.XLabel,
            YLabel = c.YLabel,
            X = c.X.ToArray(),
            Y = c.Y.ToArray(),
        },
        _ => new PayloadDto { Kind = "unknown" },
    };

    private static CplxDto ToCplx(Complex z) => new() { Re = z.Real, Im = z.Imaginary };

    private static CplxDto[] ToCplxArray(ComplexVector v)
    {
        var a = new CplxDto[v.Count];
        for (int i = 0; i < v.Count; i++) a[i] = ToCplx(v[i]);
        return a;
    }

    private static CplxDto[][] ToCplxMatrix(ComplexMatrix m)
    {
        var a = new CplxDto[m.RowCount][];
        for (int i = 0; i < m.RowCount; i++)
        {
            a[i] = new CplxDto[m.ColumnCount];
            for (int j = 0; j < m.ColumnCount; j++) a[i][j] = ToCplx(m[i, j]);
        }
        return a;
    }
}

public sealed record InspectionNodeDto
{
    [JsonPropertyName("displayName")] public string DisplayName { get; init; } = "";
    [JsonPropertyName("summary")] public string Summary { get; init; } = "";
    [JsonPropertyName("payload")] public PayloadDto Payload { get; init; } = new() { Kind = "none" };
    [JsonPropertyName("children")] public InspectionNodeDto[] Children { get; init; } = Array.Empty<InspectionNodeDto>();
}

public sealed record PayloadDto
{
    [JsonPropertyName("kind")] public string Kind { get; init; } = "none";
    [JsonPropertyName("label")] public string? Label { get; init; }
    [JsonPropertyName("realValue")] public double? RealValue { get; init; }
    [JsonPropertyName("format")] public string? Format { get; init; }
    [JsonPropertyName("complexValue")] public CplxDto? ComplexValue { get; init; }
    [JsonPropertyName("componentLabels")] public string[]? ComponentLabels { get; init; }
    [JsonPropertyName("vectorValues")] public CplxDto[]? VectorValues { get; init; }
    [JsonPropertyName("rowLabels")] public string[]? RowLabels { get; init; }
    [JsonPropertyName("columnLabels")] public string[]? ColumnLabels { get; init; }
    [JsonPropertyName("matrixValues")] public CplxDto[][]? MatrixValues { get; init; }
    [JsonPropertyName("xLabel")] public string? XLabel { get; init; }
    [JsonPropertyName("yLabel")] public string? YLabel { get; init; }
    [JsonPropertyName("x")] public double[]? X { get; init; }
    [JsonPropertyName("y")] public double[]? Y { get; init; }
}

public sealed record CplxDto
{
    [JsonPropertyName("re")] public double Re { get; init; }
    [JsonPropertyName("im")] public double Im { get; init; }
}
