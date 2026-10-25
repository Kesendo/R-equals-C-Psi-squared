using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="CatTurnedFieldSelfEncounterClaim"/> (F168). Two parent edges, both
/// resolved earlier in <c>BuildDefault</c>:
///
/// <list type="bullet">
///   <item><see cref="F63LCommutesPi2Pi2Inheritance"/>: the parity of the count of Y and Z letters
///   that keeps the cat's two kinds of words apart (step 3 of the proof).</item>
///   <item><see cref="XGlobalEigenstateMirrorPi2Inheritance"/>: the X⊗N-eigenstate endpoint, the two
///   points where the property holds (step 4, F99's α = 0).</item>
/// </list>
///
/// <para>Requires: RegisterF63LCommutesPi2Pi2Inheritance, RegisterXGlobalEigenstateMirrorPi2Inheritance.</para></summary>
public static class CatTurnedFieldSelfEncounterClaimRegistration
{
    public static ClaimRegistryBuilder RegisterCatTurnedFieldSelfEncounterClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<CatTurnedFieldSelfEncounterClaim>(b =>
        {
            var f63 = b.Get<F63LCommutesPi2Pi2Inheritance>();
            var xGlobal = b.Get<XGlobalEigenstateMirrorPi2Inheritance>();
            return new CatTurnedFieldSelfEncounterClaim(f63, xGlobal);
        });
}
