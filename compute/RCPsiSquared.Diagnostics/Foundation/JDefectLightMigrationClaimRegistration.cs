using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Schicht-1 wiring of <see cref="JDefectLightMigrationClaim"/> (the jdefect axis's
/// light-migration mechanism, Tier1Derived). Two typed parents, both Tier1Derived: <see
/// cref="AbsorptionTheoremClaim"/> (the per-mode absorption identity, exact for any Hermitian
/// H, hence δJ-pointwise along the axis) + <see cref="F1PalindromeIdentity"/> (the palindrome
/// Π·L·Π⁻¹ = −L − 2Σγ·I, Π-invariant along the axis, forcing partner light complementarity).
/// Must register after <c>RegisterAbsorptionTheoremClaim</c> and <c>RegisterF1Family</c>.</summary>
public static class JDefectLightMigrationClaimRegistration
{
    public static ClaimRegistryBuilder RegisterJDefectLightMigrationClaim(
        this ClaimRegistryBuilder builder)
    {
        builder.Register<JDefectLightMigrationClaim>(b =>
        {
            _ = b.Get<AbsorptionTheoremClaim>();   // typed parent edge (the identity)
            _ = b.Get<F1PalindromeIdentity>();     // typed parent edge (the palindrome)
            return new JDefectLightMigrationClaim();
        });
        return builder;
    }
}
