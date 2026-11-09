using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="OneDephasedSeatClaim"/> (Tier1Derived). Three typed parents, chosen among the
/// results the proof builds on: <see cref="PalindromeTwoEndCountClaim"/> (F158, both ends on spans of Pauli
/// strings and an invertible far element as a carrier), <see cref="PalindromeComplementConnectionClaim"/> (Theorem 6's
/// far equations, written before its spanning-tree hypothesis and solved in the section on planar fields outside it) and <see cref="SeatCutBlindnessClaim"/> (F157, the blind
/// seat). F152 and the colouring page have no claim of their own; F63, typed, is cited only for the three turns
/// the sheet paragraph names and stays a reference rather than an edge. Resolution is topological.</summary>
public static class OneDephasedSeatClaimRegistration
{
    public static ClaimRegistryBuilder RegisterOneDephasedSeatClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<OneDephasedSeatClaim>(b =>
            new OneDephasedSeatClaim(
                b.Get<PalindromeTwoEndCountClaim>(),
                b.Get<PalindromeComplementConnectionClaim>(),
                b.Get<SeatCutBlindnessClaim>()));
}
