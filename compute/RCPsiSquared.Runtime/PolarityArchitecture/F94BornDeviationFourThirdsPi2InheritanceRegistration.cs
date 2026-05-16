using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F94BornDeviationFourThirdsPi2Inheritance"/>:
/// F94's Δ_|00⟩ = (4/3)·Q²·K³ closed form, derived from Dyson sym3 = 8 bit-exact
/// for the specific setup |0+0+⟩ N=4 Heisenberg ring + Z-dephasing, pair (0,2),
/// |00⟩ outcome.
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_{−1} = 4</c> — the
///         "4" in the coefficient 4/3 (= a_{−1}/3 after the Dyson sym3 ÷ Taylor 3!
///         reduction). Same "4" as F86 t_peak = 1/(4γ₀) and F77's correction
///         denominator.</item>
/// </list>
///
/// <para>Tier consistency: F94 is Tier 1 derived (bit-exact symbolic via Dyson
/// sym3 + numerical verification at 0.3% over 16 samples). The dyadic-ladder
/// parent is Tier 1 derived. Tier composition is well-formed.</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>
/// in the builder pipeline.</para></summary>
public static class F94BornDeviationFourThirdsPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF94BornDeviationFourThirdsPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F94BornDeviationFourThirdsPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            return new F94BornDeviationFourThirdsPi2Inheritance(ladder);
        });
}
