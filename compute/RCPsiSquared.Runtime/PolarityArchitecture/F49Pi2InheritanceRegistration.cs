using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F49Pi2Inheritance"/>: F49's "<c>4^(N−2)</c>"
/// scaling factor as Pi2-Foundation operator-space inheritance. Three parent edges,
/// all Tier1Derived:
///
/// <list type="bullet">
///   <item><see cref="PalindromeResidualScalingClaim"/>: the F49 closed form
///         (parameterised by N and Hamiltonian class).</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_{5−2N} = 4^(N−2)</c>.</item>
///   <item><see cref="Pi2OperatorSpaceMirrorClaim"/>: provides the per-qubit-count
///         pinned mirror pairs that match <c>4^(N−2) = d²</c> for (N−2) qubits.</item>
/// </list>
///
/// <para>Tier consistency: all four Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="F1Family.F1FamilyRegistration.RegisterF1Family"/> for
/// PalindromeResidualScalingClaim + <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>
/// + <see cref="Pi2OperatorSpaceMirrorRegistration.RegisterPi2OperatorSpaceMirror"/>
/// (which itself transitively requires F88StaticDyadicAnchor).</para></summary>
public static class F49Pi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF49Pi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F49Pi2Inheritance>(b =>
        {
            _ = b.Get<PalindromeResidualScalingClaim>();   // F49 closed form parent
            var ladder = b.Get<Pi2DyadicLadderClaim>();    // provides a_{5-2N}
            var mirror = b.Get<Pi2OperatorSpaceMirrorClaim>(); // provides per-qubit-count pinning
            return new F49Pi2Inheritance(ladder, mirror);
        });
}
