using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F5DepolarizingErrorPi2Inheritance"/>:
/// F5's depolarizing error closed-form coefficients <c>2·(N−2)/3</c> as
/// Pi2-Foundation inheritance. Three parent edges:
///
/// <list type="bullet">
///   <item><see cref="F1PalindromeIdentity"/>: the F1 palindrome whose
///         "Breaks for: depolarizing" clause is what F5 quantifies.</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>
///         (DCoefficient) and the (N−2) chain-shift via the same derivation
///         overhead as F49.</item>
///   <item><see cref="Pi2OperatorSpaceMirrorClaim"/>: provides
///         <c>a_{−1} = 4 = d²</c> for N=1, from which DSquaredMinusOne = 3
///         is derived.</item>
/// </list>
///
/// <para>Tier consistency: F5 is Tier 1 proven; the inheritance claim is
/// Tier1Derived; all parents are Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires upstream registrations: F1Family, Pi2Family, Pi2DyadicLadder,
/// F88PopcountCoherence + F88StaticDyadicAnchor + Pi2OperatorSpaceMirror.</para></summary>
public static class F5DepolarizingErrorPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF5DepolarizingErrorPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F5DepolarizingErrorPi2Inheritance>(b =>
        {
            _ = b.Get<F1PalindromeIdentity>();
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var mirror = b.Get<Pi2OperatorSpaceMirrorClaim>();
            return new F5DepolarizingErrorPi2Inheritance(ladder, mirror);
        });
}
