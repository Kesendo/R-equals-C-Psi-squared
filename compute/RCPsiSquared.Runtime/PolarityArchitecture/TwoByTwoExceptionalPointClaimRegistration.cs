using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="TwoByTwoExceptionalPointClaim"/>, the one closed 2 × 2 behind many exceptional
/// points (<c>docs/THE_DOUBLE_ROOT.md</c>, section "One 2 × 2 in many settings"). Parentless: the law is the
/// algebra of a traceless 2 × 2. Its children declare the edge in their own registrations: <c>QEpLaw</c>,
/// <c>LEffMirrorAxisClaim</c>, <c>CoherenceHorizonClaim</c>, <c>TransitionBridgeF95SiblingClaim</c>,
/// <c>Spin32QuadrupolarMirrorsClaim</c>, <c>SeedHolonomyClaim</c> and <c>ExceptionalCouplingSetClaim</c>, each owning one
/// setting of it.</summary>
public static class TwoByTwoExceptionalPointClaimRegistration
{
    public static ClaimRegistryBuilder RegisterTwoByTwoExceptionalPointClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<TwoByTwoExceptionalPointClaim>(_ => new TwoByTwoExceptionalPointClaim());
}
