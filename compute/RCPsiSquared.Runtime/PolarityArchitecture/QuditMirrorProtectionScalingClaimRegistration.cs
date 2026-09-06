using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="QuditMirrorProtectionScalingClaim"/> (2026-06-17): the corollary of the
/// F121 explicit Π_d construction: its shift-aligned coverage is
/// (2d)^N / d^{2N} = (2/d)^N. This is not a universal product optimum; the parent
/// <see cref="QuditProductMirrorCap"/> carries the d=6,N=2 cap retraction.
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
