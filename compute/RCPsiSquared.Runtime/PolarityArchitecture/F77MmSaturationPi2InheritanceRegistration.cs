using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F77MmSaturationPi2Inheritance"/>:
/// F77's <c>MM(0) → 1 bit</c> saturation as the first F-formula sitting on the
/// Pi2 dyadic ladder's self-mirror pivot a_1. Two parent edges (one explicit,
/// one as registration discard):
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_1 = 1</c> (saturation
///         bit), <c>a_0 = 2</c>, <c>a_2 = 1/2</c>, <c>a_{-1} = 4</c>.</item>
///   <item><see cref="HalfAsStructuralFixedPointClaim"/>: registration discard
///         to document that the "1/2" in F77's mechanism IS the structural
///         fixed point. Transitively reachable through the ladder; the discard
///         makes the dependency explicit at registration time.</item>
/// </list>
///
/// <para>"ZERO IS THE MIRROR" reading (Tom 2026-05-09): the inheritance graph's
/// self-mirror pivot at n=1 was previously empty; the Manager-query showed
/// Pi2I4MemoryLoop had only 4 descendants, and review of the ladder anchor
/// distribution surfaced F77 (originally formal F77 = MM(0) saturation, distinct
/// from the trichotomy that took the F77 name in early code; per memory
/// feedback_F77_F87_rename) as the missed wiring.</para>
///
/// <para>Tier consistency: F77 is Tier 1 asymptotic proven; Pi2-Foundation
/// anchoring is algebraic-trivial composition. All three claims Tier1Derived
/// (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (registers HalfAsStructuralFixedPointClaim) +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>.</para></summary>
public static class F77MmSaturationPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF77MmSaturationPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F77MmSaturationPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            _ = b.Get<HalfAsStructuralFixedPointClaim>();   // documents 1/2 = structural fixed point
            return new F77MmSaturationPi2Inheritance(ladder);
        });
}
