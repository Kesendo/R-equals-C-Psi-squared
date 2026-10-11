using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="CorrelatedZDephasingPriceClaim"/>. Two parent edges, both resolved
/// earlier in <c>BuildDefault</c>:
///
/// <list type="bullet">
///   <item><see cref="F1PalindromeIdentity"/>: the identity this generalises (K = σ·1 its case).</item>
///   <item><see cref="AbsorptionTheoremClaim"/>: the cell rates the shared field reads (proof §2).</item>
/// </list>
///
/// <para>Requires: RegisterF1Family, RegisterAbsorptionTheoremClaim.</para></summary>
public static class CorrelatedZDephasingPriceClaimRegistration
{
    public static ClaimRegistryBuilder RegisterCorrelatedZDephasingPriceClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<CorrelatedZDephasingPriceClaim>(b =>
        {
            var f1 = b.Get<F1PalindromeIdentity>();
            var abs = b.Get<AbsorptionTheoremClaim>();
            return new CorrelatedZDephasingPriceClaim(f1, abs);
        });
}
