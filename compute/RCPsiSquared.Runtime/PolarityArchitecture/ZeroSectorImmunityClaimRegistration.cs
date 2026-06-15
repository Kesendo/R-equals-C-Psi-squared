using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="ZeroSectorImmunityClaim"/>. Four parent edges, all resolved
/// earlier in <c>BuildDefault</c>:
///
/// <list type="bullet">
///   <item><see cref="F1PalindromeIdentity"/>: the global palindrome M = 0 this refines on the
///         extreme weight blocks (H-independently).</item>
///   <item><see cref="F61BitAParityPi2Inheritance"/> + <see cref="F63LCommutesPi2Pi2Inheritance"/>:
///         the bit_a/bit_b C₂×C₂ sector machinery that blocks the operator space.</item>
///   <item><see cref="AbsorptionTheoremClaim"/>: the −2Σγ on the w=N corner that Π carries back.</item>
/// </list>
///
/// <para>Requires: RegisterF1Family, RegisterF61BitAParityPi2Inheritance,
/// RegisterF63LCommutesPi2Pi2Inheritance, RegisterAbsorptionTheoremClaim.</para></summary>
public static class ZeroSectorImmunityClaimRegistration
{
    public static ClaimRegistryBuilder RegisterZeroSectorImmunityClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<ZeroSectorImmunityClaim>(b =>
        {
            var f1 = b.Get<F1PalindromeIdentity>();
            var f61 = b.Get<F61BitAParityPi2Inheritance>();
            var f63 = b.Get<F63LCommutesPi2Pi2Inheritance>();
            var abs = b.Get<AbsorptionTheoremClaim>();
            return new ZeroSectorImmunityClaim(f1, f61, f63, abs);
        });
}
