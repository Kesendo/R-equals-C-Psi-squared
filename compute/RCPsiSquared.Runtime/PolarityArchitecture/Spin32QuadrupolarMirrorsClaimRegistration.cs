using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="Spin32QuadrupolarMirrorsClaim"/> (F167). One parent edge, resolved
/// topologically from <c>BuildDefault</c>:
///
/// <list type="bullet">
///   <item><see cref="PalindromeTwoEndCountClaim"/> (F158): the second route to the global-centre
///   theorem for non-negative densities, the normalised jumps being Hermitian strings that square to 1.</item>
/// </list>
///
/// <para>Requires: RegisterPalindromeTwoEndCountClaim.</para></summary>
public static class Spin32QuadrupolarMirrorsClaimRegistration
{
    public static ClaimRegistryBuilder RegisterSpin32QuadrupolarMirrorsClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<Spin32QuadrupolarMirrorsClaim>(b =>
            new Spin32QuadrupolarMirrorsClaim(b.Get<PalindromeTwoEndCountClaim>()));
}
