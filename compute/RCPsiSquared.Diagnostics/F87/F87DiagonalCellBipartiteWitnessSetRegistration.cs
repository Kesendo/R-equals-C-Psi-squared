using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 wiring of <see cref="F87DiagonalCellBipartiteWitnessSet"/> for a given
/// chain. Parameterised by <see cref="ChainSystemPrimitive"/> (so the runtime instance reads N +
/// γ from the registered chain) and depends on two typed parent claims:
/// <see cref="F87TrichotomyClassification"/> (the F87 verdict the §7 criterion is checked
/// against) and <see cref="ChiralKClaim"/> (the chiral sublattice K with KHK=−H that certifies
/// bipartite ⟹ soft).
///
/// <para>Requires: <see cref="Runtime.F1Family.F1FamilyRegistration.RegisterF1Family"/> for
/// <see cref="ChainSystemPrimitive"/>, <see cref="F87FamilyRegistration.RegisterF87Family"/> for
/// <see cref="F87TrichotomyClassification"/>, and the ChiralK registration for
/// <see cref="ChiralKClaim"/>; the builder errors with <c>MissingParent</c> if any are absent.</para>
///
/// <para>Tier consistency: F87DiagonalCellBipartiteWitnessSet is Tier1Candidate (lazy
/// classification per witness) ← F87TrichotomyClassification (Tier1Derived) and ← ChiralKClaim
/// (Tier1Derived). Inheritance check passes (4 ≥ 4). Each inner
/// <see cref="F87DiagonalCellBipartiteWitness"/> self-classifies when accessed; criterion vs F87
/// mismatches surface via <see cref="F87DiagonalCellBipartiteWitness.Matches"/>.</para></summary>
public static class F87DiagonalCellBipartiteWitnessSetRegistration
{
    public static ClaimRegistryBuilder RegisterF87DiagonalCellBipartiteWitnessSet(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F87DiagonalCellBipartiteWitnessSet>(b =>
        {
            _ = b.Get<F87TrichotomyClassification>();
            _ = b.Get<ChiralKClaim>();
            var chain = b.Get<ChainSystemPrimitive>();
            return new F87DiagonalCellBipartiteWitnessSet(chain.System);
        });
}
