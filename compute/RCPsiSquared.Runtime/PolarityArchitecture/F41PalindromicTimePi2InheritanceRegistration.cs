using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F41PalindromicTimePi2Inheritance"/>:
/// F41's palindromic-time corollary of D10. Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>
///         (HoppingCoefficient in 2·J denominator and 2π full-period; same
///         anchor as F1's TwoFactor, F66's UpperPoleCoefficient, F50's
///         DegeneracyFactor).</item>
///   <item><see cref="F1Pi2Inheritance"/>: provides the operator-level
///         palindrome identity Π·L·Π⁻¹ + L + 2σ·I = 0; F41 reads its
///         time-domain consequence (the SFF has a periodic palindromic
///         modulation; F41's t_Pi is the period of the slowest such
///         modulation).</item>
/// </list>
///
/// <para>F41 is a clean F1-corollary at the time-domain level. F44 reads F1's
/// rate identity at the eigenvalue level (ratio); F68 reads F1's pair sum at
/// the eigenvalue level (sum); F41 reads F1's pair difference at the time
/// domain (period). All three are F1-spectral-consequences with distinct
/// closed-form readings.</para>
///
/// <para>Tier consistency: F41 is Tier 1 corollary of D10's w=1 dispersion
/// derivation; FFT-confirmed &lt;1% at N=2..4, 6 in
/// <c>experiments/SPECTRAL_FORM_FACTOR.md</c>. Both parent claims Tier1Derived.</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F1Pi2InheritanceRegistration.RegisterF1Pi2Inheritance"/>.</para></summary>
public static class F41PalindromicTimePi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF41PalindromicTimePi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F41PalindromicTimePi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f1 = b.Get<F1Pi2Inheritance>();
            return new F41PalindromicTimePi2Inheritance(ladder, f1);
        });
}
