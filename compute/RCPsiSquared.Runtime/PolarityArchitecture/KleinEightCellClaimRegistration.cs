using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="KleinEightCellClaim"/> (Stage 2b
/// of the cubic-unpacking arc; first Cubic3-axis Claim).
///
/// <para>Standalone Claim: no constructor parents. Registered into the typed-knowledge
/// graph so it is visible to <see cref="PolarityCubeMap"/> aggregation and the
/// inspector. PolarityCubeMap.Cubic3Claims grows from 0 to 1.</para></summary>
public static class KleinEightCellClaimRegistration
{
    public static ClaimRegistryBuilder RegisterKleinEightCellClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<KleinEightCellClaim>(_ => new KleinEightCellClaim());
}
