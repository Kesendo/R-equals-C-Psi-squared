using RCPsiSquared.Core.F1;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="TimeIrreversibilityExclusionClaim"/>. One parent edge:
/// <see cref="F49NonUniformCrossTermClaim"/> (the cross-term value whose (N−2) factor makes the
/// anticommutator vanish at N=2), registered earlier in <c>BuildDefault</c> via
/// <c>RegisterF1Family</c>.
///
/// <para>Requires: RegisterF1Family (which registers F49NonUniformCrossTermClaim).</para></summary>
public static class TimeIrreversibilityExclusionClaimRegistration
{
    public static ClaimRegistryBuilder RegisterTimeIrreversibilityExclusionClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<TimeIrreversibilityExclusionClaim>(b =>
            new TimeIrreversibilityExclusionClaim(b.Get<F49NonUniformCrossTermClaim>()));
}
