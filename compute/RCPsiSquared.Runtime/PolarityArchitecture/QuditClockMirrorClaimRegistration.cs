using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="QuditClockMirrorClaim"/> (F172). One parent edge, resolved topologically from
/// <c>BuildDefault</c>: <see cref="PalindromeTwoEndCountClaim"/> (F158), whose count the page carries from Hermitian
/// involutions to unitary jumps with a self-adjoint dissipator.
///
/// <para>Requires: RegisterPalindromeTwoEndCountClaim.</para></summary>
public static class QuditClockMirrorClaimRegistration
{
    public static ClaimRegistryBuilder RegisterQuditClockMirrorClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<QuditClockMirrorClaim>(b => new QuditClockMirrorClaim(b.Get<PalindromeTwoEndCountClaim>()));
}
