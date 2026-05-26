using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="KleinEightCellClaim"/> (Stage 2b
/// of the cubic-unpacking arc; first Cubic3-axis Claim).
///
/// <para>Typed Klein2 parent: <see cref="KleinFourCellClaim"/>. The Cubic3 8-cell
/// structure lifts the Klein2 4-cell structure by the third Z₂ axis (y_par);
/// wired 2026-05-26 to make the quadratic ↔ cubic inheritance edge explicit.
/// Registered into the typed-knowledge graph so it is visible to
/// <see cref="PolarityCubeMap"/> aggregation and the inspector.</para></summary>
public static class KleinEightCellClaimRegistration
{
    public static ClaimRegistryBuilder RegisterKleinEightCellClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<KleinEightCellClaim>(b =>
            new KleinEightCellClaim(b.Get<KleinFourCellClaim>()));
}
