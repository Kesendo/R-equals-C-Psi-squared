using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="SurvivalIncompletenessMirrorClaim"/> (the survival law and the
/// V-Effect/incompleteness are Pi2-ladder inversion-mirror partners a_0&lt;-&gt;a_2; Tier1Candidate).
/// Typed parents: <see cref="AbsorptionTheoremClaim"/> (a_0, registered above) +
/// <see cref="HalfAsStructuralFixedPointClaim"/> (a_2, a Pi2-foundation root constructed fresh - it is
/// not itself a top-level registered claim). Must register after <c>RegisterAbsorptionTheoremClaim</c>.</summary>
public static class SurvivalIncompletenessMirrorClaimRegistration
{
    public static ClaimRegistryBuilder RegisterSurvivalIncompletenessMirrorClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<SurvivalIncompletenessMirrorClaim>(b =>
        {
            var survival = b.Get<AbsorptionTheoremClaim>();             // typed parent a_0 (the survival law)
            var incompleteness = new HalfAsStructuralFixedPointClaim();  // a_2 anchor (the V-Effect/incompleteness)
            return new SurvivalIncompletenessMirrorClaim(survival, incompleteness);
        });
}
