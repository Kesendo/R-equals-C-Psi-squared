using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.F1Family;

/// <summary>Builder extension that registers the F1 family: <see cref="ChainSystemPrimitive"/>,
/// then <see cref="F1PalindromeIdentity"/>, then four Tier-1-derived F1 children that the
/// <see cref="F1KnowledgeBase"/> exposes as top-level properties (the KB's fifth Tier-1
/// child <c>SingleBodyScaling</c> is the omitted one; see the SingleBody paragraph below):
///
/// <list type="bullet">
///   <item><see cref="PalindromeResidualScalingClaim"/> (<see cref="HamiltonianClass.Main"/>):
///         residual ‖M‖² scaling on the main class. The <paramref name="hClass"/> parameter
///         on <c>RegisterF1Family</c> selects which class this single registration carries.</item>
///   <item><see cref="F1T1ResidualClosedForm"/>: ‖M(T1)‖² = 4^(N−1)·[3·Σγ² + 4·(Σγ)²].</item>
///   <item><see cref="F1T1ResidualPi2Decomposition"/>: Π²-orthogonal Pythagorean split of
///         the T1 residual into (anti, sym) parts; depends on F1T1ResidualClosedForm.</item>
///   <item><see cref="F1DepolResidualClosedForm"/>: ‖M(depol)‖² = 4^(N−1)·[(16/9)·Σγ² + 16·(Σγ)²].</item>
/// </list>
///
/// <para>Dependency edges (parent → child):</para>
/// <list type="bullet">
///   <item>ChainSystemPrimitive → F1PalindromeIdentity (topological-order hint;
///         F1PalindromeIdentity has no Core-side parameter but the chain parameterises
///         the surrounding KB context).</item>
///   <item>ChainSystemPrimitive → PalindromeResidualScalingClaim (the N from the chain
///         is what the scaling claim is built against).</item>
///   <item>F1PalindromeIdentity → PalindromeResidualScalingClaim.</item>
///   <item>F1PalindromeIdentity → F1T1ResidualClosedForm.</item>
///   <item>F1PalindromeIdentity → F1T1ResidualPi2Decomposition.</item>
///   <item>F1T1ResidualClosedForm → F1T1ResidualPi2Decomposition (the Π²-decomposition
///         closes the parent total Pythagorically: (anti) + (sym) = (3·Σγ² + 4·(Σγ)²)).</item>
///   <item>F1PalindromeIdentity → F1DepolResidualClosedForm.</item>
/// </list>
///
/// <para><b>SingleBody scaling: deliberately omitted.</b> The
/// <see cref="F1KnowledgeBase"/> also exposes a <c>SingleBodyScaling</c> Tier-1 property
/// (the same <see cref="PalindromeResidualScalingClaim"/> type with
/// <see cref="HamiltonianClass.SingleBody"/>). The runtime
/// <see cref="ClaimRegistryBuilder"/> is type-keyed (one factory per
/// <see cref="System.Type"/>); attempting to register a second
/// <see cref="PalindromeResidualScalingClaim"/> with a different class throws
/// <c>InvariantViolationException(rule: "DuplicateRegistration")</c>. No native
/// discriminator parameter or generic wrapper key exists on the builder, so we surface
/// only one scaling claim per registry build. The <paramref name="hClass"/> parameter
/// lets the caller choose which class this is (Main by default, matching the original
/// behaviour). The KB still exposes both classes for inspection.</para></summary>
public static class F1FamilyRegistration
{
    /// <summary>Register the six F1-family Claims (ChainSystemPrimitive +
    /// F1PalindromeIdentity + PalindromeResidualScalingClaim + F1T1ResidualClosedForm +
    /// F1T1ResidualPi2Decomposition + F1DepolResidualClosedForm) for a given
    /// <paramref name="chain"/>. Default Hamiltonian class for the scaling claim is
    /// <see cref="HamiltonianClass.Main"/>; chain bond count and degree-squared sum
    /// default to <c>null</c> (use <c>FactorChain</c>).</summary>
    public static ClaimRegistryBuilder RegisterF1Family(
        this ClaimRegistryBuilder builder,
        ChainSystem chain,
        HamiltonianClass hClass = HamiltonianClass.Main)
    {
        return builder
            .Register<ChainSystemPrimitive>(_ => new ChainSystemPrimitive(chain))
            .Register<F1PalindromeIdentity>(b =>
            {
                _ = b.Get<ChainSystemPrimitive>();
                return new F1PalindromeIdentity();
            })
            .Register<PalindromeResidualScalingClaim>(b =>
            {
                var primitive = b.Get<ChainSystemPrimitive>();
                _ = b.Get<F1PalindromeIdentity>();
                return new PalindromeResidualScalingClaim(
                    N: primitive.System.N,
                    hClass: hClass);
            })
            .Register<F1T1ResidualClosedForm>(b =>
            {
                _ = b.Get<F1PalindromeIdentity>();
                return new F1T1ResidualClosedForm();
            })
            .Register<F1T1ResidualPi2Decomposition>(b =>
            {
                _ = b.Get<F1PalindromeIdentity>();
                _ = b.Get<F1T1ResidualClosedForm>();
                return new F1T1ResidualPi2Decomposition();
            })
            .Register<F1DepolResidualClosedForm>(b =>
            {
                _ = b.Get<F1PalindromeIdentity>();
                return new F1DepolResidualClosedForm();
            });
    }
}
