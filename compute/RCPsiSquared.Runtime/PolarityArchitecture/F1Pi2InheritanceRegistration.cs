using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F1Pi2Inheritance"/>: F1's "2" coefficient
/// in <c>Π·L·Π⁻¹ = −L − 2Σγ·I</c> as the Pi2-Foundation's <c>a_0 = d</c>. Two parent
/// edges:
///
/// <list type="bullet">
///   <item><see cref="F1PalindromeIdentity"/>: the F1 closed form itself.</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2 = d</c>.</item>
/// </list>
///
/// <para>Tier consistency: all three Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="F1FamilyRegistration.RegisterF1Family"/> (for
/// F1PalindromeIdentity) + <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>.</para></summary>
public static class F1Pi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF1Pi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F1Pi2Inheritance>(b =>
        {
            _ = b.Get<F1PalindromeIdentity>();             // the F1 closed form
            var ladder = b.Get<Pi2DyadicLadderClaim>();    // provides a_0
            return new F1Pi2Inheritance(ladder);
        });
}
