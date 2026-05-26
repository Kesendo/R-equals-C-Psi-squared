using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="F87Z2CubedRefinementN4K3"/> (F103).
///
/// <para>Typed Cubic3 parent: <see cref="KleinEightCellClaim"/>. F103 enumerates
/// 294 Pauli pairs at N=4 k=3 across the Z₂³ 8-cell decomposition; the parent
/// edge makes the (Klein × y_par) enumeration anchor explicit. Wired
/// 2026-05-26.</para></summary>
public static class F87Z2CubedRefinementN4K3Registration
{
    public static ClaimRegistryBuilder RegisterF87Z2CubedRefinementN4K3(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F87Z2CubedRefinementN4K3>(b =>
            new F87Z2CubedRefinementN4K3(b.Get<KleinEightCellClaim>()));
}
