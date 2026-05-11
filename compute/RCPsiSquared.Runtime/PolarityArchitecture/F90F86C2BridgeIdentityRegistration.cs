using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F90F86C2BridgeIdentity"/>: the algebraic
/// identification F86 c=2 K_b ↔ F89 path-(N−1) per-bond Hellmann-Feynman. Two typed
/// parent edges:
///
/// <list type="bullet">
///   <item><see cref="F89TopologyOrbitClosure"/>: F89 framework parent (the bridge
///         resolves F86 c=2 into F89 territory).</item>
///   <item><see cref="F89PathKAtLockMechanismClaim"/>: AT-lock mechanism is the
///         structural reason F86's HWHM_left/Q_peak universal constants exist (4-mode
///         floor 0.6715 = AT-locked F_a/F_b contribution per F89's overlap-only support).</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 derived (algebraic identity by construction +
/// numerical bit-exact verification at 27/29 bonds across N=5..8). Both parents
/// Tier 1 derived.</para>
///
/// <para>F86-side anchor (C2HwhmRatio etc.) lives in the F86 family registry; this
/// bridge claim does not declare F86 as a typed parent because F86 is a collection
/// of partial results not yet closed (per <c>feedback_f86_is_collection_basin</c>
/// memory). The bridge stands on the F89 side: a Tier-1 statement about how F89's
/// machinery encompasses F86 c=2's K_b observable.</para>
///
/// <para>Requires: <see cref="F89TopologyOrbitClosureRegistration.RegisterF89TopologyOrbitClosure"/>
/// and <see cref="F89PathKAtLockMechanismClaimRegistration.RegisterF89PathKAtLockMechanismClaim"/>.</para></summary>
public static class F90F86C2BridgeIdentityRegistration
{
    public static ClaimRegistryBuilder RegisterF90F86C2BridgeIdentity(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F90F86C2BridgeIdentity>(b =>
        {
            var f89 = b.Get<F89TopologyOrbitClosure>();
            var atLock = b.Get<F89PathKAtLockMechanismClaim>();
            return new F90F86C2BridgeIdentity(f89, atLock);
        });
}
