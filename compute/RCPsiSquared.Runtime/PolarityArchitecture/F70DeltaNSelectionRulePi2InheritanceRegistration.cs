using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F70DeltaNSelectionRulePi2Inheritance"/>:
/// F70's <c>|ΔN| ≤ k</c> selection rule for k-local partial trace as
/// Pi2-Foundation inheritance. F70 is the kinematic foundation cited by
/// F71 (mirror symmetry of c₁) and F72 (block-diagonal DD⊕CC structure).
/// One parent edge:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_1 = 1</c>
///         (single-site threshold, F77 sibling) and <c>a_0 = 2</c>
///         (pair-local threshold, F1/F66/F86 Q_EP sibling).</item>
/// </list>
///
/// <para>Tier consistency: F70 is Tier 1 proven kinematic; verified at N=5
/// with 9 |ΔN| ≥ 2 pairs giving zero to machine precision. Both claims
/// Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>
/// in the builder pipeline.</para></summary>
public static class F70DeltaNSelectionRulePi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF70DeltaNSelectionRulePi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F70DeltaNSelectionRulePi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            return new F70DeltaNSelectionRulePi2Inheritance(ladder);
        });
}
