using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F1T1AmplitudeDampingPi2Inheritance"/>: F1's
/// T1-amplitude-damping ‖M‖² closed-form factors as Pi2-Foundation inheritance. Three
/// parent edges:
///
/// <list type="bullet">
///   <item><see cref="F1PalindromeIdentity"/>: the F1 palindrome identity whose
///         residual ‖M‖² is being scaled.</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_{−(N+1)}</c> (H-part)
///         and <c>a_{3−2N}</c> (T1-part) and <c>a_{−1} = 4</c> (multiplier).</item>
///   <item><see cref="Pi2OperatorSpaceMirrorClaim"/>: pinning <c>4^(N−1) = d²</c>
///         für (N−1) qubits cross-check.</item>
/// </list>
///
/// <para>Tier consistency: F1T1AmplitudeDampingPi2Inheritance is Tier1Candidate
/// (analytic derivation open per F1 open-question Item 2). All three parents are
/// Tier1Derived. TierStrength inheritance trivially passes (5 ≥ 4).</para>
///
/// <para>Requires upstream registrations: <see cref="F1Family.F1FamilyRegistration.RegisterF1Family"/>
/// + <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F88PopcountCoherenceRegistration.RegisterF88PopcountCoherence"/> +
/// <see cref="F88StaticDyadicAnchorRegistration.RegisterF88StaticDyadicAnchor"/> +
/// <see cref="Pi2OperatorSpaceMirrorRegistration.RegisterPi2OperatorSpaceMirror"/>.</para></summary>
public static class F1T1AmplitudeDampingPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF1T1AmplitudeDampingPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F1T1AmplitudeDampingPi2Inheritance>(b =>
        {
            _ = b.Get<F1PalindromeIdentity>();    // parent F-formula identity
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var mirror = b.Get<Pi2OperatorSpaceMirrorClaim>();
            return new F1T1AmplitudeDampingPi2Inheritance(ladder, mirror);
        });
}
