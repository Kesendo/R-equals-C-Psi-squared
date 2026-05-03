namespace RCPsiSquared.Core.Inspection;

/// <summary>Query primitives for any <see cref="IInspectable"/> tree. The <see cref="Walk"/>
/// method does a depth-first enumeration; everything else is composable LINQ on the result.
///
/// <para>Use this layer to ask "what is in this object manager"-style questions without
/// hand-walking the tree:</para>
/// <code>
///   kb.Walk().OfType&lt;F86Claim&gt;().Where(c => c.Tier == Tier.Tier1Derived)
///   kb.Walk().Where(n => n.DisplayName.Contains("Q_EP"))
///   kb.Walk().Count(n => n.Payload is InspectablePayload.Real)
/// </code>
/// </summary>
public static class InspectionQuery
{
    /// <summary>Depth-first walk of the IInspectable tree starting at <paramref name="root"/>,
    /// including the root itself. Lazy — does not materialise the children list per node.
    /// </summary>
    public static IEnumerable<IInspectable> Walk(this IInspectable root)
    {
        yield return root;
        foreach (var child in root.Children)
            foreach (var descendant in Walk(child))
                yield return descendant;
    }

    /// <summary>Walk to a maximum depth (root = depth 0, its children = depth 1, …).
    /// Useful when the tree is unbounded or self-referential.</summary>
    public static IEnumerable<IInspectable> Walk(this IInspectable root, int maxDepth)
    {
        if (maxDepth < 0) yield break;
        yield return root;
        if (maxDepth == 0) yield break;
        foreach (var child in root.Children)
            foreach (var descendant in Walk(child, maxDepth - 1))
                yield return descendant;
    }

    /// <summary>All nodes whose payload is a real-valued scalar — convenient for extracting
    /// numerical data from the tree.</summary>
    public static IEnumerable<(IInspectable Node, double Value)> WalkRealScalars(this IInspectable root) =>
        root.Walk()
            .Where(n => n.Payload is InspectablePayload.Real)
            .Select(n => (n, ((InspectablePayload.Real)n.Payload).Value));

    /// <summary>All nodes whose <see cref="IInspectable.DisplayName"/> contains the given
    /// substring (case-sensitive).</summary>
    public static IEnumerable<IInspectable> WhereDisplayNameContains(this IEnumerable<IInspectable> nodes, string needle) =>
        nodes.Where(n => n.DisplayName.Contains(needle));

    /// <summary>Return the first node matching the predicate, or null. Useful for direct
    /// lookups (e.g. "find the t_peak law node").</summary>
    public static IInspectable? FirstOrDefault(this IEnumerable<IInspectable> nodes,
        Func<IInspectable, bool> predicate) => System.Linq.Enumerable.FirstOrDefault(nodes, predicate);

    /// <summary>Group nodes by their <see cref="IInspectable.Payload"/> kind (none / real /
    /// scalar / vector / matrix / curve).</summary>
    public static IDictionary<string, int> CountByPayloadKind(this IEnumerable<IInspectable> nodes) =>
        nodes
            .GroupBy(n => PayloadKindName(n.Payload))
            .ToDictionary(g => g.Key, g => g.Count());

    private static string PayloadKindName(InspectablePayload p) => p switch
    {
        InspectablePayload.None => "none",
        InspectablePayload.Real => "real",
        InspectablePayload.Scalar => "scalar",
        InspectablePayload.Vector => "vector",
        InspectablePayload.MatrixView => "matrix",
        InspectablePayload.Curve => "curve",
        _ => "unknown",
    };
}
