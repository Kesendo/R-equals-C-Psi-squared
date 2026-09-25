using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="TransitionBridgeF95SiblingClaim"/>: two distinct
/// positive-b quadratic applications reuse F95's coordinate, one at b=½ in z_rec and the
/// genuine F86 toy 2x2 exceptional point at b=4γ₀ in
/// the decay variable z=−λ, with angle zero
/// at each boundary. The separate FRAGILE_BRIDGE threshold (two qubits per chain) is generically its own EP2 on the real γ axis, at zero
/// decay, not this quadratic. One typed parent,
/// <see cref="F95AngleAtQuadraticZeroPi2Inheritance"/>, the shared lens.
///
/// <para>Requires only the standalone, parentless
/// <see cref="F95AngleAtQuadraticZeroPi2InheritanceRegistration.RegisterF95AngleAtQuadraticZeroPi2Inheritance"/>
/// registration.</para></summary>
public static class TransitionBridgeF95SiblingClaimRegistration
{
    public static ClaimRegistryBuilder RegisterTransitionBridgeF95SiblingClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<TransitionBridgeF95SiblingClaim>(b =>
        {
            var f95 = b.Get<F95AngleAtQuadraticZeroPi2Inheritance>();
            return new TransitionBridgeF95SiblingClaim(f95);
        });
}
