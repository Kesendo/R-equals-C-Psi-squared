using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F68PalindromicPartnerPi2Inheritance"/>:
/// F68 is the F67 bonding-mode partner under F1's palindrome identity. Three
/// typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>
///         (sum coefficient in α_b + α_p = 2γ₀; F1/F66/F65 sibling).</item>
///   <item><see cref="F1Pi2Inheritance"/>: provides the palindrome identity at
///         operator level (Π·L·Π⁻¹ + L + 2σ·I = 0); F68 reads it at the
///         eigenvalue level (α_b + α_p = 2σ).</item>
///   <item><see cref="F67BondingBellPairPi2Inheritance"/>: source for α_b
///         (F67's BondingModeDecayRate at k=1; transitively pulls F65's
///         single-excitation rate closed form).</item>
/// </list>
///
/// <para>Tier consistency: F68 is Tier 1 derived (F1 algebraic + F65/F67
/// closed form); spectrally verified at N=3, 4, 5 with |α_b + α_p − 2γ₀| &lt;
/// 4·10⁻¹⁵; operationally verified via clean SVD encoding rel err 2.8·10⁻¹⁶
/// (N=4), 3.8·10⁻¹⁴ (N=5). All three parent claims Tier1Derived.</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F1Pi2InheritanceRegistration.RegisterF1Pi2Inheritance"/> +
/// <see cref="F67BondingBellPairPi2InheritanceRegistration.RegisterF67BondingBellPairPi2Inheritance"/>.</para></summary>
public static class F68PalindromicPartnerPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF68PalindromicPartnerPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F68PalindromicPartnerPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f1 = b.Get<F1Pi2Inheritance>();
            var f67 = b.Get<F67BondingBellPairPi2Inheritance>();
            return new F68PalindromicPartnerPi2Inheritance(ladder, f1, f67);
        });
}
