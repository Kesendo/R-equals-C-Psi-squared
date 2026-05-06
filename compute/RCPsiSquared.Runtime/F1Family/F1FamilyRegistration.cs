using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.F1Family;

/// <summary>Schicht 1 registration: ChainSystemPrimitive then F1PalindromeIdentity then
/// PalindromeResidualScalingClaim. F1PalindromeIdentity has no Core-side parameter (it is
/// the project's master palindrome); we declare an edge from ChainSystemPrimitive anyway so
/// the topological order surfaces the parameterisation before the abstract identity.
/// PalindromeResidualScalingClaim depends on both: the closed form is the F1 residual's
/// Frobenius scaling for a given (N, HamiltonianClass) pair on a given chain.</summary>
public static class F1FamilyRegistration
{
    /// <summary>Register the three F1 family Claims for a given <paramref name="chain"/>.
    /// Default Hamiltonian class is <see cref="HamiltonianClass.Main"/>; chain bond count
    /// and degree-squared sum default to <c>null</c> (use <c>FactorChain</c>).</summary>
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
            });
    }
}
