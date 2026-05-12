using RCPsiSquared.Core.SymmetryFamily;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="ZGlobalMirrorRefinement"/>; one typed parent
/// (<see cref="SymmetryFamilyInventory"/>).</summary>
public static class ZGlobalMirrorRefinementRegistration
{
    public static ClaimRegistryBuilder RegisterZGlobalMirrorRefinement(
        this ClaimRegistryBuilder builder) =>
        builder.Register<ZGlobalMirrorRefinement>(b =>
        {
            var inv = b.Get<SymmetryFamilyInventory>();
            return new ZGlobalMirrorRefinement(inv);
        });
}
