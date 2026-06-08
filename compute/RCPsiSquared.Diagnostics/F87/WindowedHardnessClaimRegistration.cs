using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 wiring of <see cref="WindowedHardnessClaim"/> (F115). Single typed parent:
/// <see cref="F87DiagonalCellBipartiteWitnessSet"/>, the F103 §7 diagonal-K (bipartite-chirality)
/// criterion. F115 is the mask-combinatorial reading of that same diagonal-cell hard/soft line, so
/// it sits beneath the bipartite witness set rather than beside it. No chain primitive is needed:
/// the (1 + x)-valuation criterion is pure GF(2) bit arithmetic on the term masks.
///
/// <para>Requires: <see cref="F87DiagonalCellBipartiteWitnessSetRegistration.RegisterF87DiagonalCellBipartiteWitnessSet"/>
/// to have run (the builder errors with <c>MissingParent</c> otherwise).</para>
///
/// <para>Tier consistency: WindowedHardnessClaim is Tier1Derived ← F87DiagonalCellBipartiteWitnessSet
/// (Tier1Derived). The strength-inheritance check is parent ≥ child, i.e. 5 ≥ 5, which passes; both
/// were promoted from candidate in the 2026-06-08 formal promotion pass once the §7.5/§7.6 converse
/// closed (derived modulo standard PT; per §7.10 only the (1 + x)-valuation, not the obstruction
/// size, reaches the spectrum).</para></summary>
public static class WindowedHardnessClaimRegistration
{
    public static ClaimRegistryBuilder RegisterWindowedHardnessClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<WindowedHardnessClaim>(b =>
        {
            _ = b.Get<F87DiagonalCellBipartiteWitnessSet>();
            return new WindowedHardnessClaim();
        });
}
