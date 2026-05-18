using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F1DepolResidualClosedFormPi2Inheritance"/>:
/// F1's depolarizing-residual closed-form coefficients (16/9, 16) as Pi2-Foundation
/// inheritance, parallel to <see cref="F5DepolarizingErrorPi2Inheritance"/>'s wiring of
/// the linear scalar 2·(N−2)/3. Three parent edges:
///
/// <list type="bullet">
///   <item><see cref="F1PalindromeIdentity"/>: the F1 palindrome identity whose
///         depolarizing residual ‖M‖² is being decomposed (registered via
///         <see cref="F1Family.F1FamilyRegistration.RegisterF1Family"/>).</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides the ladder index for
///         <c>a_{−1} = d²</c> (the multiplier in the cross-site coefficient).</item>
///   <item><see cref="Pi2OperatorSpaceMirrorClaim"/>: pins <c>d² = 4</c> at N=1
///         via <see cref="Pi2OperatorSpaceMirrorClaim.PairAt"/>(1).OperatorSpace,
///         from which <c>d² − 1 = 3</c> (and hence 16/9 = (4/3)²) is derived.</item>
/// </list>
///
/// <para>Tier consistency: <see cref="F1DepolResidualClosedFormPi2Inheritance"/> is
/// Tier1Derived, matching its Tier1Derived parents (5 ≥ 5). Unlike the T1 analogue
/// <see cref="F1T1AmplitudeDampingPi2Inheritance"/> (Tier1Candidate because T1's
/// c_1 = 3 is not Pi2-anchored), both depol coefficients reduce cleanly to algebra
/// over <c>d²</c> and <c>d² − 1</c>.</para>
///
/// <para>Requires upstream registrations:
/// <see cref="F1Family.F1FamilyRegistration.RegisterF1Family"/> +
/// <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F88bPopcountCoherenceRegistration.RegisterF88bPopcountCoherence"/> +
/// <see cref="F88bStaticDyadicAnchorRegistration.RegisterF88bStaticDyadicAnchor"/> +
/// <see cref="Pi2OperatorSpaceMirrorRegistration.RegisterPi2OperatorSpaceMirror"/>
/// (same prerequisite set as <see cref="F5DepolarizingErrorPi2InheritanceRegistration"/>).</para></summary>
public static class F1DepolResidualClosedFormPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF1DepolResidualClosedFormPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F1DepolResidualClosedFormPi2Inheritance>(b =>
        {
            _ = b.Get<F1PalindromeIdentity>();
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var mirror = b.Get<Pi2OperatorSpaceMirrorClaim>();
            return new F1DepolResidualClosedFormPi2Inheritance(ladder, mirror);
        });
}
