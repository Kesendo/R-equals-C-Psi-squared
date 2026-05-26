using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="YParityIndependenceAtK3"/> (F102).
///
/// <para>Typed Cubic3 parent: <see cref="KleinEightCellClaim"/>. F102 names the
/// y_par axis where the third Z₂ classifier becomes independent at k_body ≥ 3;
/// KleinEightCellClaim is the structural anchor for the 8-cell Z₂³ decomposition
/// where this independence is enumerated. Wired 2026-05-26.</para></summary>
public static class YParityIndependenceAtK3Registration
{
    public static ClaimRegistryBuilder RegisterYParityIndependenceAtK3(
        this ClaimRegistryBuilder builder) =>
        builder.Register<YParityIndependenceAtK3>(b =>
            new YParityIndependenceAtK3(b.Get<KleinEightCellClaim>()));
}
