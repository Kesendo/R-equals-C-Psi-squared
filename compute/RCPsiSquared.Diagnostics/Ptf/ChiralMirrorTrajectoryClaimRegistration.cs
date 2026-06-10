using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Ptf;

/// <summary>Schicht-1 wiring of <see cref="ChiralMirrorTrajectoryClaim"/> (the Tier1Derived
/// PTF K₁ site-wise trajectory identity). Single typed parent <see cref="ChiralKClaim"/>
/// (Tier1Derived, the eigenvalue side of the same sublattice chirality). Must register after
/// <c>RegisterChiralK</c>.</summary>
public static class ChiralMirrorTrajectoryClaimRegistration
{
    public static ClaimRegistryBuilder RegisterChiralMirrorTrajectoryClaim(
        this ClaimRegistryBuilder builder)
    {
        builder.Register<ChiralMirrorTrajectoryClaim>(b =>
        {
            _ = b.Get<ChiralKClaim>(); // typed parent edge (parent before child)
            return new ChiralMirrorTrajectoryClaim();
        });
        return builder;
    }
}
