using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F80PiCommutatorAnticommutatorIdentity"/>:
/// the proven F80 Step 5 identity <c>Π·[H,·]·Π⁻¹ = s·{H,·}</c>, s = −ε_P·ε_Q. One
/// parent edge, <see cref="RCPsiSquared.Core.F1.F1PalindromeIdentity"/>: the Π that
/// conjugates the commutator here is the same order-4 Π of F1's palindrome identity.
///
/// <para>Tier consistency: parent and claim are both Tier1Derived.</para>
///
/// <para>Requires <c>RegisterF1PalindromeIdentity</c> in the builder pipeline.</para></summary>
public static class F80PiCommutatorAnticommutatorIdentityRegistration
{
    public static ClaimRegistryBuilder RegisterF80PiCommutatorAnticommutatorIdentity(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F80PiCommutatorAnticommutatorIdentity>(b =>
        {
            var f1 = b.Get<F1PalindromeIdentity>();                 // typed parent: the Π is F1's Π
            return new F80PiCommutatorAnticommutatorIdentity(f1);
        });
}
