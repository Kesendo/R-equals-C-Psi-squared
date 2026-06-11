using RCPsiSquared.Core.OpenArcs;

namespace RCPsiSquared.Core.Inspection;

/// <summary>Surfaces the <see cref="OpenArcsRegistry"/> as an inspection subtree: one child
/// per arc, each leaf carrying where the arc parked and its next concrete move in the summary
/// (retired arcs carry the retirement reason instead). Lives in Core next to the registry
/// whose data it presents; the world root mounts it as the "arcs" section so the Object
/// Manager displays its own unfinished business.</summary>
public static class OpenArcsInspectableNode
{
    public static IInspectable Build()
    {
        var entries = OpenArcsRegistry.All;
        int open = entries.Count(a => a.Status == OpenArcStatus.Open);
        int retired = entries.Count - open;
        return new InspectableNode(
            displayName: "arcs",
            summary: $"{open} open arc(s), {retired} retired",
            children: entries.Select(Entry).ToArray());
    }

    private static IInspectable Entry(OpenArc a) =>
        new InspectableNode(
            displayName: a.Name,
            summary: a.Status == OpenArcStatus.Retired
                ? $"retired: {a.RetiredReason}"
                : $"{a.Opened} parked: {a.ParkedAt} -> next: {a.NextStep}");
}
