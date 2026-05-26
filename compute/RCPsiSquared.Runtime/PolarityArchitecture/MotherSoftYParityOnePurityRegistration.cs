using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="MotherSoftYParityOnePurity"/> (F109).
///
/// <para>Typed Cubic3 parent: <see cref="KleinEightCellClaim"/>. F109 pins
/// the y_par=1 purity of mother Klein (0,0) soft cells across the Z₂³ 8-cell
/// decomposition. Wired 2026-05-26.</para></summary>
public static class MotherSoftYParityOnePurityRegistration
{
    public static ClaimRegistryBuilder RegisterMotherSoftYParityOnePurity(
        this ClaimRegistryBuilder builder) =>
        builder.Register<MotherSoftYParityOnePurity>(b =>
            new MotherSoftYParityOnePurity(b.Get<KleinEightCellClaim>()));
}
