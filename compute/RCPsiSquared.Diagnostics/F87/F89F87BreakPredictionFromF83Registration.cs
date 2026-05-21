using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 wiring of <see cref="F89F87BreakPredictionFromF83"/>: the bridge claim
/// that the F87 palindrome-residual norm ‖M(H_F89 + δH)‖_F is predicted bit-exactly by F83's
/// closed-form Π-decomposition of the perturbation δH alone. Single parent edge:
///
/// <list type="bullet">
///   <item><see cref="F89F87TrulyInheritance"/>: establishes that H_F89 = J·(XX+YY) is
///         F87-Truly and contributes 0 to both non-truly norms, so the predicted break norm
///         comes entirely from δH.</item>
/// </list>
///
/// <para>Tier consistency: <see cref="F89F87BreakPredictionFromF83"/> and its parent are
/// both Tier1Derived.</para>
///
/// <para>Requires upstream registration:
/// <see cref="F89F87TrulyInheritanceRegistration.RegisterF89F87TrulyInheritance"/>.</para></summary>
public static class F89F87BreakPredictionFromF83Registration
{
    public static ClaimRegistryBuilder RegisterF89F87BreakPredictionFromF83(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89F87BreakPredictionFromF83>(b =>
        {
            var f89F87Truly = b.Get<F89F87TrulyInheritance>();
            return new F89F87BreakPredictionFromF83(f89F87Truly);
        });
}
