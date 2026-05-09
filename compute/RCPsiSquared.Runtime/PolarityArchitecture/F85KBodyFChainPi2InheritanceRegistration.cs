using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F85KBodyFChainPi2Inheritance"/>:
/// F85 k-body F-chain extension. Generalizes F49 (k=2) and F83 (anti-fraction)
/// to k-body. Three parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_{−1} = 4</c>
///         (Π²-odd factor), <c>a_{−2} = 8</c> (Π²-even non-truly factor),
///         <c>a_0 = 2</c> (Π²-odd count denominator).</item>
///   <item><see cref="F49Pi2Inheritance"/>: F49's 2-body formula is the
///         k = 2 base case of F85.</item>
///   <item><see cref="F83AntiFractionPi2Inheritance"/>: F85's coefficients 4, 8
///         match F83's MNormCoefficientForOdd/ForEvenNontruly; F83 generalizes
///         to k-body with Π²-class grouping per F85.</item>
/// </list>
///
/// <para>Tier consistency: F85 is Tier 1 verified bit-exact k=2,3,4 across
/// 108 explicit Pauli tuple cases. All four claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F49Pi2InheritanceRegistration.RegisterF49Pi2Inheritance"/> +
/// <see cref="F83AntiFractionPi2InheritanceRegistration.RegisterF83AntiFractionPi2Inheritance"/>.</para></summary>
public static class F85KBodyFChainPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF85KBodyFChainPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F85KBodyFChainPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f49 = b.Get<F49Pi2Inheritance>();
            var f83 = b.Get<F83AntiFractionPi2Inheritance>();
            return new F85KBodyFChainPi2Inheritance(ladder, f49, f83);
        });
}
