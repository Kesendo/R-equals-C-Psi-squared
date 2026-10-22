using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="Spin32QuadrupolarMirrorsClaim"/> (F167). Two parent edges, resolved
/// topologically from <c>BuildDefault</c>:
///
/// <list type="bullet">
///   <item><see cref="PalindromeTwoEndCountClaim"/> (F158): the second route to the global-centre
///   theorem for non-negative densities, the normalised jumps being Hermitian strings that square to 1.</item>
///   <item><see cref="TwoByTwoExceptionalPointClaim"/>, for one part of its content: the bright pairs' exceptional
///   point is a setting of the one closed 2 × 2.</item>
/// </list>
///
/// <para>Requires: RegisterPalindromeTwoEndCountClaim, RegisterTwoByTwoExceptionalPointClaim.</para></summary>
public static class Spin32QuadrupolarMirrorsClaimRegistration
{
    public static ClaimRegistryBuilder RegisterSpin32QuadrupolarMirrorsClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<Spin32QuadrupolarMirrorsClaim>(b =>
        {
            _ = b.Get<TwoByTwoExceptionalPointClaim>();   // partial scope: only the bright pairs' exceptional point is a setting of the one closed block
            return new Spin32QuadrupolarMirrorsClaim(b.Get<PalindromeTwoEndCountClaim>());
        });
}
