namespace RCPsiSquared.Core.Inspection;

/// <summary>Generic factory <see cref="IInspectable"/> for ad-hoc tree nodes (group headers,
/// computed scalars, etc.) where defining a dedicated wrapper type is overkill.
///
/// <para>For domain objects (D_eff, M_h_per_bond, probe, …) prefer a typed wrapper class
/// that implements <see cref="IInspectable"/> directly — that way the wrapper carries
/// computation, not just data.</para>
/// </summary>
public sealed class InspectableNode : IInspectable
{
    public string DisplayName { get; }
    public string Summary { get; }
    public IEnumerable<IInspectable> Children { get; }
    public InspectablePayload Payload { get; }

    public InspectableNode(
        string displayName,
        string summary = "",
        IEnumerable<IInspectable>? children = null,
        InspectablePayload? payload = null)
    {
        DisplayName = displayName;
        Summary = summary;
        Children = children ?? Array.Empty<IInspectable>();
        Payload = payload ?? InspectablePayload.Empty;
    }

    /// <summary>Group header: a container with no payload and a fixed list of children.</summary>
    public static InspectableNode Group(string displayName, params IInspectable[] children) =>
        new(displayName, summary: $"{children.Length} item(s)", children: children);

    /// <summary>Group header with lazy children stream.</summary>
    public static InspectableNode Group(string displayName, IEnumerable<IInspectable> children, int? count = null) =>
        new(displayName, summary: count is null ? "" : $"{count.Value} item(s)", children: children);

    public static InspectableNode RealScalar(string displayName, double value, string? format = null) =>
        new(displayName,
            summary: format is null ? value.ToString("G6") : value.ToString(format),
            payload: new InspectablePayload.Real(displayName, value, format));
}
