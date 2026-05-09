using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F71MirrorSymmetryPi2Inheritance"/>:
/// F71's c_1 mirror symmetry as the first F-formula direct edge for
/// <see cref="HalfIntegerMirrorClaim"/>. The N-parity classification HalfIntegerMirror
/// makes (odd N → half-integer w_XY axis, no Pauli string on axis; even N →
/// integer w_XY axis, strings on axis exist) IS the same N-parity F71 makes
/// at the bond-index level (odd N → all bonds in disjoint pairs; even N →
/// self-paired center bond). One parent edge (registration discard, since
/// F71's universal logic does not require N injection):
///
/// <list type="bullet">
///   <item><see cref="HalfIntegerMirrorClaim"/>: registration discard. The
///         claim is parameterised by N (registered at chain N=5 by default
///         in the cockpit); F71's universal MirrorPair logic does not need N
///         at construction time. Before F71 wiring, HalfIntegerMirror had
///         0 descendants — the N-parity Pi2-Foundation primitive was
///         registered but unused (Tom 2026-05-09 mirror-map check).</item>
/// </list>
///
/// <para>Tier consistency: F71 is Tier 1 proven kinematic
/// (PROOF_C1_MIRROR_SYMMETRY); Pi2-Foundation anchoring is documentation-only
/// (no number flows through). Both claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> +
/// <see cref="HalfIntegerMirrorRegistration.RegisterHalfIntegerMirror"/>(N)
/// in the builder pipeline.</para></summary>
public static class F71MirrorSymmetryPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF71MirrorSymmetryPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F71MirrorSymmetryPi2Inheritance>(b =>
        {
            _ = b.Get<HalfIntegerMirrorClaim>();   // direct edge: N-parity classification
            return new F71MirrorSymmetryPi2Inheritance();
        });
}
