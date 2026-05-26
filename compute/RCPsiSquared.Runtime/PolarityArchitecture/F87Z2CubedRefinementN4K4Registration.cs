using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="F87Z2CubedRefinementN4K4"/> (F106).
///
/// <para>Typed Cubic3 parent: <see cref="KleinEightCellClaim"/>. F106 enumerates
/// 4248 Pauli pairs at N=4 k=4 across the Z₂³ 8-cell decomposition (new
/// enumeration vs F103's 294 at k=3). Wired 2026-05-26.</para></summary>
public static class F87Z2CubedRefinementN4K4Registration
{
    public static ClaimRegistryBuilder RegisterF87Z2CubedRefinementN4K4(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F87Z2CubedRefinementN4K4>(b =>
            new F87Z2CubedRefinementN4K4(b.Get<KleinEightCellClaim>()));
}
