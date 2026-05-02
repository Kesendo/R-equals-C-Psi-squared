using System.Text.Json.Serialization;

namespace RCPsiSquared.Core.Tests.Fixtures;

/// <summary>DTO matching one entry in the Python step_e/step_f summary.json output.
///
/// <para>Read-only contract for validation tests. Mirrors `simulations/results/eq022_*/summary.json`
/// produced by `_eq022_b1_step_e_resonance_shape.py` and `_eq022_b1_step_f_universality_extension.py`.</para>
/// </summary>
public sealed record PythonStepSummaryEntry
{
    [JsonPropertyName("label")] public string? Label { get; init; }
    [JsonPropertyName("c")] public int C { get; init; }
    [JsonPropertyName("N")] public int N { get; init; }
    [JsonPropertyName("n")] public int LowerPopcount { get; init; }
    [JsonPropertyName("block_dim")] public int BlockDim { get; init; }
    [JsonPropertyName("gamma_0")] public double GammaZero { get; init; } = 0.05;
    [JsonPropertyName("Q_peak_interior")] public double? QPeakInterior { get; init; }
    [JsonPropertyName("K_max_interior")] public double? KMaxInterior { get; init; }
    [JsonPropertyName("hwhm_left_interior")] public double? HwhmLeftInterior { get; init; }
    [JsonPropertyName("Q_peak_endpoint")] public double? QPeakEndpoint { get; init; }
    [JsonPropertyName("K_max_endpoint")] public double? KMaxEndpoint { get; init; }
    [JsonPropertyName("hwhm_left_endpoint")] public double? HwhmLeftEndpoint { get; init; }
}
