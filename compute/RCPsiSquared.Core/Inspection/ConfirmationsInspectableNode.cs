using RCPsiSquared.Core.Confirmations;

namespace RCPsiSquared.Core.Inspection;

/// <summary>Surfaces the hardware <see cref="ConfirmationsRegistry"/> as an inspection subtree:
/// one child per confirmed prediction, each leaf carrying date, machine, observable, and the
/// predicted-vs-measured comparison in its summary. Lives in Core next to the registry whose
/// data it presents; the world root mounts it as the "confirmations" section.</summary>
public static class ConfirmationsInspectableNode
{
    public static IInspectable Build()
    {
        var entries = ConfirmationsRegistry.All;
        return new InspectableNode(
            displayName: "confirmations",
            summary: $"{entries.Count} hardware-confirmed prediction(s)",
            children: entries.Select(Entry).ToArray());
    }

    private static IInspectable Entry(Confirmation c) =>
        new InspectableNode(
            displayName: c.Name,
            summary: $"{c.Date} {c.Machine} {c.Observable}: predicted {c.PredictedValue} vs measured {c.MeasuredValue}");
}
