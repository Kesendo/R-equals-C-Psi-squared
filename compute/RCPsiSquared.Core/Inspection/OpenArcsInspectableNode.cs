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
        // Both counts are PREDICATES. The retired count used to be the complement of the open
        // one, which silently absorbs any arc that is neither into "retired".
        int open = entries.Count(a => a.Status == OpenArcStatus.Open);
        int retired = entries.Count(a => a.Status == OpenArcStatus.Retired);
        int parked = entries.Count(a => a.Status == OpenArcStatus.Open && a.ParkedReason is not null);
        return new InspectableNode(
            displayName: "arcs",
            summary: $"{open} open arc(s) ({parked} of them parked), {retired} retired",
            children: entries.Select(Entry).ToArray());
    }

    private static IInspectable Entry(OpenArc a) =>
        new InspectableNode(
            displayName: a.Name,
            summary: a.Status == OpenArcStatus.Retired
                ? $"retired: {a.RetiredReason}"
                : a.ParkedReason is not null
                    ? $"{a.Opened} PARKED: {a.ParkedReason} -> next when it wakes: {a.NextStep}"
                    : $"{a.Opened} parked: {a.ParkedAt} -> next: {a.NextStep}");
}
