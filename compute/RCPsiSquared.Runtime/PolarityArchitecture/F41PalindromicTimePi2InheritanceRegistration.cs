using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F41PalindromicTimePi2Inheritance"/>:
/// F41's palindromic-pair trace-amplitude period from D10 + F1. Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>
///         (HoppingCoefficient in 2·J denominator and 2π full-period; same
///         anchor as F1's TwoFactor, F66's UpperPoleCoefficient, F50's
///         DegeneracyFactor).</item>
///   <item><see cref="F1Pi2Inheritance"/>: provides the operator-level
///         palindrome identity Π·L·Π⁻¹ + L + 2σ·I = 0; F1 supplies the
///         frequency-negated partner of D10's ω_min, producing
///         <c>2cos(ω_min t)</c> in the trace amplitude. The complete SFF is
///         the squared modulus of the full amplitude and has additional
///         doubled and cross frequencies.</item>
/// </list>
///
/// <para>D10 supplies the exact frequency and F1 supplies its partner. F41
/// records the resulting pair-amplitude period. It does not define a physical
/// Heisenberg time or a short/long-time boundary.</para>
///
/// <para>Tier consistency: F41 is Tier 1 corollary of D10's (0,1) coherence block dispersion
/// derivation. The separate producer finds a matching FFT candidate at N=2..4
/// and N=6, but not at N=5 or N=7. Both parent claims are Tier1Derived.</para>
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
