using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F44CrooksLikeRateIdentityPi2Inheritance"/>:
/// F44's algebraic rate identity ln(d_fast/d_slow) = 2·artanh(Δd/(2·Σγ)) for
/// any palindromic Liouvillian pair. Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c> twice —
///         as the SumCoefficient (d_fast + d_slow = 2·Σγ) and as the
///         ArTanhCoefficient (in 2·artanh(...) from ln((1+x)/(1−x)) =
///         2·artanh(x)).</item>
///   <item><see cref="F1Pi2Inheritance"/>: provides the operator-level
///         palindrome identity Π·L·Π⁻¹ + L + 2σ·I = 0; F44 reads its
///         eigenvalue-level consequence d_fast + d_slow = 2·Σγ as the
///         pre-condition for the artanh closed form.</item>
/// </list>
///
/// <para>F44 generalizes F68: F68 (α_b + α_p = 2γ₀) is F44 applied to the
/// bonding-mode pair under endpoint Z-dephasing where Σγ = γ₀. F44's artanh
/// reading adds the ratio (not just the sum) of the pair rates as a closed form.</para>
///
/// <para><b>Important:</b> F44 is NOT a Crooks fluctuation theorem. The
/// Jarzynski equality fails (⟨exp(−Δd)⟩ ≈ 0.93, not 1). The palindrome has
/// the FORM of detailed balance without BEING detailed balance. F44 is purely
/// algebraic; the typed wiring reflects this by inheriting from F1 (algebraic
/// palindrome) only, not from any thermodynamic primitive.</para>
///
/// <para>Tier consistency: F44 is Tier 1 proven D08 (one-line algebraic
/// identity from artanh definition); valid for any palindromic Liouvillian,
/// all N. Both parent claims Tier1Derived.</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F1Pi2InheritanceRegistration.RegisterF1Pi2Inheritance"/>.</para></summary>
public static class F44CrooksLikeRateIdentityPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF44CrooksLikeRateIdentityPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F44CrooksLikeRateIdentityPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f1 = b.Get<F1Pi2Inheritance>();
            return new F44CrooksLikeRateIdentityPi2Inheritance(ladder, f1);
        });
}
