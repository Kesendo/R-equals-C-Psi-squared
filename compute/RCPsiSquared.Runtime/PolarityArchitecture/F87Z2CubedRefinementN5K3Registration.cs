using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="F87Z2CubedRefinementN5K3"/> (F105).
///
/// <para>Typed Cubic3 parent: <see cref="KleinEightCellClaim"/>. F105 enumerates
/// 294 Pauli pairs at N=5 k=3 across the Z₂³ 8-cell decomposition (same
/// enumeration as F103, larger chain). Wired 2026-05-26 to make the Klein2 ↔
/// Cubic3 ↔ YParity chain typed.</para></summary>
public static class F87Z2CubedRefinementN5K3Registration
{
    public static ClaimRegistryBuilder RegisterF87Z2CubedRefinementN5K3(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F87Z2CubedRefinementN5K3>(b =>
            new F87Z2CubedRefinementN5K3(b.Get<KleinEightCellClaim>()));
}
