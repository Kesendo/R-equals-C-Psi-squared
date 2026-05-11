using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89Path3OcticEpClaim"/>: the path-3
/// octic EP at q ≈ 0.659 with merged eigenvalue −4γ + 2iJ. Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="F89TopologyOrbitClosure"/>: orbit-closure framework parent.</item>
///   <item><see cref="F89PathKAtLockMechanismClaim"/>: AT-rate-midpoint
///         interpretation (Re(λ_EP) = −4γ between rates 2γ and 6γ) requires
///         the AT-lock mechanism context.</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 derived (analytical from disc factor +
/// machine-precision numerical verification). Both parents Tier 1 derived.</para>
///
/// <para>Requires: <see cref="F89TopologyOrbitClosureRegistration.RegisterF89TopologyOrbitClosure"/>
/// and <see cref="F89PathKAtLockMechanismClaimRegistration.RegisterF89PathKAtLockMechanismClaim"/>.</para></summary>
public static class F89Path3OcticEpClaimRegistration
{
    public static ClaimRegistryBuilder RegisterF89Path3OcticEpClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89Path3OcticEpClaim>(b =>
        {
            var f89 = b.Get<F89TopologyOrbitClosure>();
            var atLock = b.Get<F89PathKAtLockMechanismClaim>();
            return new F89Path3OcticEpClaim(f89, atLock);
        });
}
