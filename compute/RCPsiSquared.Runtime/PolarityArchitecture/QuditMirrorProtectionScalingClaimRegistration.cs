using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="QuditMirrorProtectionScalingClaim"/> (2026-06-17): the corollary of the
/// F121 product-mirror cap that the palindrome-protected fraction is (2d)^N / d^{2N} = (2/d)^N, decaying
/// exponentially in the local dimension d and = 1 only at the qubit d = 2. A single typed parent edge to
/// the Tier1Derived <see cref="QuditProductMirrorCap"/> (the cap (2d)^N); the total d^{2N} is the
/// Liouville-space dimension and the fraction is exact arithmetic, so the claim is genuinely Tier1Derived.
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
