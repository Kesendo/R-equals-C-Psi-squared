using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89UnifiedFaClosedFormClaim"/>: the
/// unified F_a AT-locked amplitude closed-form template across path-3..6.
/// Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="F89TopologyOrbitClosure"/>: orbit-closure framework parent.</item>
///   <item><see cref="F89PathKAtLockMechanismClaim"/>: AT-lock mechanism is the
///         structural foundation (overlap-only / no-overlap-only support).</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 derived (bit-exact across N=5..20 and q=0.5..3
/// for path-3..6, machine precision diffs ~10⁻¹⁷). Both parents Tier 1 derived.</para>
///
/// <para>Requires: <see cref="F89TopologyOrbitClosureRegistration.RegisterF89TopologyOrbitClosure"/>
/// and <see cref="F89PathKAtLockMechanismClaimRegistration.RegisterF89PathKAtLockMechanismClaim"/>.</para></summary>
public static class F89UnifiedFaClosedFormClaimRegistration
{
    public static ClaimRegistryBuilder RegisterF89UnifiedFaClosedFormClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89UnifiedFaClosedFormClaim>(b =>
        {
            var f89 = b.Get<F89TopologyOrbitClosure>();
            var atLock = b.Get<F89PathKAtLockMechanismClaim>();
            return new F89UnifiedFaClosedFormClaim(f89, atLock);
        });
}
