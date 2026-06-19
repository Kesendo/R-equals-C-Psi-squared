using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="SurvivorDiffusionGradientClaim"/> ((D) the closure functional,
/// felt_time arc D follow-up; Tier1Candidate). Two typed parents, both registered earlier in the
/// BuildDefault chain: <see cref="AbsorptionTheoremClaim"/> (the -2gamma rate the diffusion mode decays at)
/// + <see cref="SurvivalIncompletenessMirrorClaim"/> (the (A) survivor whose density gradient this is).
/// Must register after <c>RegisterSurvivalIncompletenessMirrorClaim</c> (and naturally alongside the sibling
/// <see cref="StoneSurvivorClosureClaim"/>, the trajectory-level dual).</summary>
public static class SurvivorDiffusionGradientClaimRegistration
{
    public static ClaimRegistryBuilder RegisterSurvivorDiffusionGradientClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<SurvivorDiffusionGradientClaim>(b =>
        {
            var rateLaw = b.Get<AbsorptionTheoremClaim>();                 // the -2gamma rate, a_0
            var survivor = b.Get<SurvivalIncompletenessMirrorClaim>();      // the (A) survivor (density mode)
            return new SurvivorDiffusionGradientClaim(rateLaw, survivor);
        });
}
