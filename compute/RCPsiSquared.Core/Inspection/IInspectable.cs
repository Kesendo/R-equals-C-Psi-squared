namespace RCPsiSquared.Core.Inspection;

/// <summary>An object that can present itself as a node in an inspection tree (Object Manager).
///
/// <para>Every implementer is a computational unit, not a data carrier: the properties
/// surfaced through <see cref="Children"/> and <see cref="Payload"/> are computed lazily
/// from the underlying state, and the per-node <see cref="Summary"/> is a computed
/// fingerprint (norm, support pattern, peak position, …) that lets a tree renderer show
/// meaningful information without unfolding the full subtree.</para>
///
/// <para>Children are <see cref="IEnumerable{T}"/> so implementers can stream lazily;
/// the Object Manager only enumerates when a node is expanded.</para>
/// </summary>
public interface IInspectable
{
    /// <summary>Short label shown in the tree.</summary>
    string DisplayName { get; }

    /// <summary>One-line computed fingerprint shown next to the label.</summary>
    string Summary { get; }

    /// <summary>Sub-nodes — lazy. Each child is itself a computational wrapper.</summary>
    IEnumerable<IInspectable> Children { get; }

    /// <summary>Leaf payload for the renderer (scalar, vector, matrix, curve, or none for pure
    /// containers).</summary>
    InspectablePayload Payload { get; }
}
