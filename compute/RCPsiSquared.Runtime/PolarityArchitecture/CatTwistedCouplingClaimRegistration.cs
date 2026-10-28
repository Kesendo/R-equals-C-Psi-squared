using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="CatTwistedCouplingClaim"/> (F168 under a twisted coupling). One parent
/// edge, resolved earlier in <c>BuildDefault</c>:
///
/// <list type="bullet">
///   <item><see cref="CatTurnedFieldSelfEncounterClaim"/>: F168 itself, which this widens; where the
///   angles follow the twist it is F168 in the turned frame with Nφ replaced by Σφ_l.</item>
/// </list>
///
/// <para>Requires: RegisterCatTurnedFieldSelfEncounterClaim.</para></summary>
public static class CatTwistedCouplingClaimRegistration
{
    public static ClaimRegistryBuilder RegisterCatTwistedCouplingClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<CatTwistedCouplingClaim>(b =>
            new CatTwistedCouplingClaim(b.Get<CatTurnedFieldSelfEncounterClaim>()));
}
