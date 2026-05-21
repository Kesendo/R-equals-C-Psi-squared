using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 wiring of <see cref="F89F87TrulyInheritance"/>: the F89-to-F87 bridge
/// claim that F89's bond Hamiltonian H_b = J·(XX+YY) classifies as F87-Truly, so the F89
/// AT-lock Re(λ_n) = −2γ₀ on F_a modes is the n_diff=1 instance of F87-Truly's bit-exact
/// Absorption Theorem. Single parent edge:
///
/// <list type="bullet">
///   <item><see cref="F87TrichotomyClassification"/>: the Tier1Derived F87 trichotomy
///         classification law, the discriminator that places XX+YY in the Truly cell.</item>
/// </list>
///
/// <para>Tier consistency: <see cref="F89F87TrulyInheritance"/> and its parent are both
/// Tier1Derived.</para>
///
/// <para>Requires upstream registration:
/// <see cref="F87FamilyRegistration.RegisterF87Family"/> (registers
/// F87TrichotomyClassification).</para></summary>
public static class F89F87TrulyInheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF89F87TrulyInheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89F87TrulyInheritance>(b =>
        {
            var trichotomy = b.Get<F87TrichotomyClassification>();
            return new F89F87TrulyInheritance(trichotomy);
        });
}
