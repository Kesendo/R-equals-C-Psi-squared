using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="MirrorLocalityBoundClaim"/>. Two parent edges, both resolved earlier in
/// <c>BuildDefault</c>:
///
/// <list type="bullet">
///   <item><see cref="PalindromeTwoEndCountClaim"/>: F158, whose §(f8) similarity gives the element's existence.</item>
///   <item><see cref="F84ThermalAmplitudeDampingPi2Inheritance"/>: the thermal channel σ⁻/σ⁺ at unequal rates as a typed object.</item>
/// </list>
///
/// <para>Requires: RegisterPalindromeTwoEndCountClaim, RegisterF84ThermalAmplitudeDampingPi2Inheritance.</para></summary>
public static class MirrorLocalityBoundClaimRegistration
{
    public static ClaimRegistryBuilder RegisterMirrorLocalityBoundClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<MirrorLocalityBoundClaim>(b =>
        {
            var f158 = b.Get<PalindromeTwoEndCountClaim>();
            var thermal = b.Get<F84ThermalAmplitudeDampingPi2Inheritance>();
            return new MirrorLocalityBoundClaim(f158, thermal);
        });
}
