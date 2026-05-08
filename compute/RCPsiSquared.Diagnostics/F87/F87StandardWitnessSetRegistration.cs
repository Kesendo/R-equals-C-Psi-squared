using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 wiring of <see cref="F87StandardWitnessSet"/> for a given chain.
/// Parameterised by <see cref="ChainSystemPrimitive"/> (so the runtime instance reads N + γ
/// from the registered chain) and depends on <see cref="F87TrichotomyClassification"/> as
/// its parent claim (the Tier1Derived classification law).
///
/// <para>Requires: <see cref="F1FamilyRegistration.RegisterF1Family"/> for
/// <see cref="ChainSystemPrimitive"/>, plus <see cref="F87FamilyRegistration.RegisterF87Family"/>
/// for <see cref="F87TrichotomyClassification"/>. The latter in turn requires
/// <see cref="Pi2FamilyRegistration"/> for <see cref="PolarityLayerOriginClaim"/>; the
/// builder errors with <c>MissingParent</c> if any of the three families are missing.</para>
///
/// <para>Tier consistency: F87StandardWitnessSet is Tier2Empirical (lazy classification per
/// witness) ← F87TrichotomyClassification (Tier1Derived). Inheritance check passes (5 ≥ 2).
/// Each inner <see cref="F87CanonicalWitness"/> is itself Tier2Empirical and self-classifies
/// when accessed; classification mismatches surface via <see cref="F87CanonicalWitness.Matches"/>.</para></summary>
public static class F87StandardWitnessSetRegistration
{
    public static ClaimRegistryBuilder RegisterF87StandardWitnessSet(this ClaimRegistryBuilder builder) =>
        builder.Register<F87StandardWitnessSet>(b =>
        {
            _ = b.Get<F87TrichotomyClassification>();
            var chain = b.Get<ChainSystemPrimitive>();
            return new F87StandardWitnessSet(chain.System);
        });
}
