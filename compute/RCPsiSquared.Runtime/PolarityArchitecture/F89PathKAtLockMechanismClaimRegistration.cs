using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89PathKAtLockMechanismClaim"/>:
/// the universal AT-lock mechanism for path-k (SE, DE) sub-blocks. One typed
/// parent edge:
///
/// <list type="bullet">
///   <item><see cref="F89TopologyOrbitClosure"/>: orbit-closure framework
///         parent (AT-lock is a structural property of the path-k orbit class).</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 derived (numerically verified across
/// path-3..6 with 100% overlap/no-overlap support, F_b sigs at machine zero).
/// Parent F89 also Tier 1 derived.</para>
///
/// <para>Requires: <see cref="F89TopologyOrbitClosureRegistration.RegisterF89TopologyOrbitClosure"/>.</para></summary>
public static class F89PathKAtLockMechanismClaimRegistration
{
    public static ClaimRegistryBuilder RegisterF89PathKAtLockMechanismClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89PathKAtLockMechanismClaim>(b =>
        {
            var f89 = b.Get<F89TopologyOrbitClosure>();
            return new F89PathKAtLockMechanismClaim(f89);
        });
}
