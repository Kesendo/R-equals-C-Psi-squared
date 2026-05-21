using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="Pi2KleinBilinearTable"/>: the self-computing
/// 9-Pauli-pair-bilinear by 4-Klein-cell assignment table. Single parent edge:
///
/// <list type="bullet">
///   <item><see cref="KleinFourCellClaim"/>: the F88a two-axis Π² decomposition into the
///         four Klein cells (Pp / Pm / Mp / Mm) that this table assigns each bilinear to.</item>
/// </list>
///
/// <para>Tier consistency: <see cref="Pi2KleinBilinearTable"/> is Tier2Empirical
/// (self-computing per bilinear), under the Tier1Derived <see cref="KleinFourCellClaim"/>
/// parent.</para>
///
/// <para>Requires upstream registration: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (registers KleinFourCellClaim).</para></summary>
public static class Pi2KleinBilinearTableRegistration
{
    public static ClaimRegistryBuilder RegisterPi2KleinBilinearTable(
        this ClaimRegistryBuilder builder) =>
        builder.Register<Pi2KleinBilinearTable>(b =>
        {
            _ = b.Get<KleinFourCellClaim>();
            return new Pi2KleinBilinearTable();
        });
}
