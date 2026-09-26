using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="QuditMirrorProtectionScalingClaim"/> (2026-06-17): the corollary of the
/// F121 product-mirror cap that the palindrome-protected fraction is P(d, N) / d^{2N} = (2/d)^N for
/// d ≤ 5, the N-th power of 2/d, and = 1 only at the qubit d = 2; from d = 6 the
/// parent <see cref="QuditProductMirrorCap"/> is larger than (2d)^N and the fraction stays exponential.
///
/// <para>Requires <see cref="QuditProductMirrorCapRegistration.RegisterQuditProductMirrorCap"/> earlier in
/// the builder pipeline.</para></summary>
public static class QuditMirrorProtectionScalingClaimRegistration
{
    public static ClaimRegistryBuilder RegisterQuditMirrorProtectionScalingClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<QuditMirrorProtectionScalingClaim>(b =>
            new QuditMirrorProtectionScalingClaim(b.Get<QuditProductMirrorCap>()));
}
