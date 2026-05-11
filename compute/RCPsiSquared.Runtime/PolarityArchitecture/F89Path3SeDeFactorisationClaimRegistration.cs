using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89Path3SeDeFactorisationClaim"/>:
/// the path-3 (SE, DE) S_2-sym deg-2·deg-2·deg-8 char poly factorisation.
/// Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="F89TopologyOrbitClosure"/>: orbit-closure framework parent.</item>
///   <item><see cref="F89PathKAtLockMechanismClaim"/>: AT-lock mechanism is the
///         structural reason F_a, F_b are at exact AT rates 2γ, 6γ.</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 derived (sympy-verified factorisation,
/// numerical eigvec match for both parent claims also Tier 1 derived).</para>
///
/// <para>Requires: <see cref="F89TopologyOrbitClosureRegistration.RegisterF89TopologyOrbitClosure"/>
/// and <see cref="F89PathKAtLockMechanismClaimRegistration.RegisterF89PathKAtLockMechanismClaim"/>.</para></summary>
public static class F89Path3SeDeFactorisationClaimRegistration
{
    public static ClaimRegistryBuilder RegisterF89Path3SeDeFactorisationClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89Path3SeDeFactorisationClaim>(b =>
        {
            var f89 = b.Get<F89TopologyOrbitClosure>();
            var atLock = b.Get<F89PathKAtLockMechanismClaim>();
            return new F89Path3SeDeFactorisationClaim(f89, atLock);
        });
}
