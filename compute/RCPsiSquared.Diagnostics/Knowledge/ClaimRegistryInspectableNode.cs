using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Knowledge;

/// <summary>Surfaces a <see cref="ClaimRegistry"/> as an inspection subtree: one group per
/// <see cref="Tier"/> (in declaration order), the group summary the claim count, the group
/// children the claims themselves (each already renders as <see cref="IInspectable"/> via the
/// existing --claim path). Lives in Diagnostics, the layer that builds the registry; the world
/// root mounts it as the "claims" section.</summary>
public static class ClaimRegistryInspectableNode
{
    public static IInspectable Build(ClaimRegistry registry)
    {
        var tierGroups = Enum.GetValues<Tier>()
            .Select(t => (Tier: t, Claims: registry.AllOfTier(t)))
            .Where(g => g.Claims.Count > 0)
            .Select(g => InspectableNode.Group(
                $"{g.Tier.Label()} ({g.Claims.Count})",
                g.Claims.Cast<IInspectable>(),
                count: g.Claims.Count))
            .Cast<IInspectable>()
            .ToArray();

        return new InspectableNode(
            displayName: "claims",
            summary: $"{registry.Count} typed claim(s) across {tierGroups.Length} tier(s)",
            children: tierGroups);
    }
}
