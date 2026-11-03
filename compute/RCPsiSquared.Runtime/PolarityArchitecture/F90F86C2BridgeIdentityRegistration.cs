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
///         structural reason F86's per-bond-class HWHM_left/Q_peak values exist (the bare
///         doubled-PTF floor, 0.7265 on the envelope, is the 2×2 model of the AT-locked F_a/F_b
///         pair per F89's overlap-only support; the per-bond ratio differs across the F71
///         sub-classes, two bands rather than one universal value).</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 derived (algebraic identity by construction +
/// numerical verification at all 22 bonds across N=5..8 to 10⁻¹⁰). Both parents
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
