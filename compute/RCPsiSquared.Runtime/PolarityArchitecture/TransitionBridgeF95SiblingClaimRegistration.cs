using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="TransitionBridgeF95SiblingClaim"/>: the cusp CΨ = ¼ (our
/// TransitionBridge) and the F86 exceptional point are F95 siblings, two instances of the angle at a
/// quadratic's discriminant zero, at two anchors b (the cusp at ½, the EP at 4γ₀). One typed parent,
/// <see cref="F95AngleAtQuadraticZeroPi2Inheritance"/>, the shared lens.
///
/// <para>Requires <see cref="F95AngleAtQuadraticZeroPi2InheritanceRegistration.RegisterF95AngleAtQuadraticZeroPi2Inheritance"/>
/// (which itself requires the Pi2 family).</para></summary>
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
