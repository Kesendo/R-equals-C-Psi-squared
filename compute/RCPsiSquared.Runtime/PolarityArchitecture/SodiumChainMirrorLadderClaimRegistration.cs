using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="SodiumChainMirrorLadderClaim"/> (F167, the chain). Two parent edges,
/// resolved topologically from <c>BuildDefault</c>:
///
/// <list type="bullet">
///   <item><see cref="Spin32QuadrupolarMirrorsClaim"/>: the single nucleus, the n = 1 case of the ladder.</item>
///   <item><see cref="PalindromeTwoEndCountClaim"/> (F158): the route of the full dipolar coupling's global break.</item>
/// </list>
///
/// <para>Requires: RegisterSpin32QuadrupolarMirrorsClaim, RegisterPalindromeTwoEndCountClaim.</para></summary>
public static class SodiumChainMirrorLadderClaimRegistration
{
    public static ClaimRegistryBuilder RegisterSodiumChainMirrorLadderClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<SodiumChainMirrorLadderClaim>(b =>
            new SodiumChainMirrorLadderClaim(b.Get<Spin32QuadrupolarMirrorsClaim>(), b.Get<PalindromeTwoEndCountClaim>()));
}
