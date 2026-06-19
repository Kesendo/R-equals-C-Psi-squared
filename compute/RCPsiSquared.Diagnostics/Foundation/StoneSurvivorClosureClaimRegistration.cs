using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Ptf;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="StoneSurvivorClosureClaim"/> (THE STONE, felt_time arc B; Tier1Candidate).
/// Three typed parents, all registered earlier in the BuildDefault chain:
/// <see cref="AbsorptionTheoremClaim"/> (the -2gamma rate) + <see cref="SurvivalIncompletenessMirrorClaim"/>
/// (the (A) value/vector survivor this confirms) + <see cref="Ptf.ChiralMirrorTrajectoryClaim"/> (the PTF
/// painter closure law). Must register after <c>RegisterSurvivalIncompletenessMirrorClaim</c>.</summary>
public static class StoneSurvivorClosureClaimRegistration
{
    public static ClaimRegistryBuilder RegisterStoneSurvivorClosureClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<StoneSurvivorClosureClaim>(b =>
        {
            var rateLaw = b.Get<AbsorptionTheoremClaim>();                  // the -2gamma rate, a_0
            var valueVector = b.Get<SurvivalIncompletenessMirrorClaim>();    // the (A) survivor this confirms
            var painterClosure = b.Get<ChiralMirrorTrajectoryClaim>();       // the PTF painter closure law (EQ-014)
            return new StoneSurvivorClosureClaim(rateLaw, valueVector, painterClosure);
        });
}
