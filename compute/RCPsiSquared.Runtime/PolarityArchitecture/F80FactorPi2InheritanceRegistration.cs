using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F80FactorPi2Inheritance"/>: F80's "±2i"
/// closed-form factor decomposed into Pi2-Foundation pieces. Two parent edges, both
/// Tier1Derived:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides the "2" via <c>Term(0) = a_0 = d</c>.</item>
///   <item><see cref="Pi2I4MemoryLoopClaim"/>: provides the "i" via <c>PowerOfI(1)</c>.</item>
/// </list>
///
/// <para>The instance is built live with both parent claims so the composition is
/// drift-checked end-to-end through the registry.</para>
///
/// <para>Tier consistency: all three Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>,
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>,
/// <see cref="Pi2I4MemoryLoopRegistration.RegisterPi2I4MemoryLoop"/> in the builder
/// pipeline.</para></summary>
public static class F80FactorPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF80FactorPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F80FactorPi2Inheritance>(b =>
        {
            var f1 = b.Get<RCPsiSquared.Core.F1.F1PalindromeIdentity>();   // typed parent: F80's M is F1's residual
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var loop = b.Get<Pi2I4MemoryLoopClaim>();
            return new F80FactorPi2Inheritance(f1, ladder, loop);
        });
}
