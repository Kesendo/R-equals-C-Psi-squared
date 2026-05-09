using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 wiring of <see cref="F87Pi2Inheritance"/>: F87 ←(F1 residual)←
/// Pi2-Foundation. Two parent edges:
///
/// <list type="bullet">
///   <item><see cref="F87TrichotomyClassification"/>: the F87 closed-form trichotomy
///         itself (Tier1Derived in <c>RCPsiSquared.Diagnostics.F87</c>).</item>
///   <item><see cref="F1Pi2Inheritance"/>: provides the transitive "2" = a_0 anchor
///         through F1's residual M.</item>
/// </list>
///
/// <para>F87 ≠ Klein cells: the earlier KleinFour 4-cell parent edge has been removed.
/// F87 (14/19/3 trichotomy) was developed pre-Klein from Pauli-pair combinatorics;
/// the 4-way Π²-eigenspace decomposition is a distinct primitive that does not feed
/// F87's discriminator. See class docstring.</para>
///
/// <para>Tier consistency: all three Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires upstream registrations: <see cref="F87FamilyRegistration.RegisterF87Family"/>
/// (for F87TrichotomyClassification) plus the Pi2-Foundation registrations
/// <c>RegisterPi2Family + RegisterPi2DyadicLadder + RegisterF1Pi2Inheritance</c>.</para></summary>
public static class F87Pi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF87Pi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F87Pi2Inheritance>(b =>
        {
            _ = b.Get<F87TrichotomyClassification>();
            var f1 = b.Get<F1Pi2Inheritance>();
            return new F87Pi2Inheritance(f1);
        });
}
