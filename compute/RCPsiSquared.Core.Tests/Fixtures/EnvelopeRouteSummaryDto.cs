using System.Text.Json.Serialization;

namespace RCPsiSquared.Core.Tests.Fixtures;

/// <summary>DTO for one entry of <c>Fixtures/python_envelope_route_summary.json</c>: the F86 class
/// peaks of one (c, N) block on the envelope max_t K_b(Q, t), computed by the Python second route
/// <c>simulations/f86_envelope_second_route.py</c>.</summary>
public sealed record EnvelopeRouteEntry
{
    [JsonPropertyName("label")] public string? Label { get; init; }
    [JsonPropertyName("c")] public int C { get; init; }
    [JsonPropertyName("N")] public int N { get; init; }
    [JsonPropertyName("n")] public int LowerPopcount { get; init; }
    [JsonPropertyName("gamma_0")] public double GammaZero { get; init; }
    [JsonPropertyName("Q_peak_interior")] public double QPeakInterior { get; init; }
    [JsonPropertyName("K_max_interior")] public double KMaxInterior { get; init; }
    [JsonPropertyName("hwhm_left_over_q_peak_interior")] public double HwhmRatioInterior { get; init; }
    [JsonPropertyName("Q_peak_endpoint")] public double QPeakEndpoint { get; init; }
    [JsonPropertyName("K_max_endpoint")] public double KMaxEndpoint { get; init; }
    [JsonPropertyName("hwhm_left_over_q_peak_endpoint")] public double HwhmRatioEndpoint { get; init; }

    public override string ToString() => Label ?? $"c{C}_N{N}";
}
