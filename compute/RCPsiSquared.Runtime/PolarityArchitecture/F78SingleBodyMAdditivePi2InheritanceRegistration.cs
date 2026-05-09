using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F78SingleBodyMAdditivePi2Inheritance"/>:
/// F78 single-body M = Σ_l M_l ⊗ I_(others). Three typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>
///         (EigenvalueCoefficient in ±2c_l·γ·i) and <c>a_{−1} = 4</c>
///         (PerSiteDimension).</item>
///   <item><see cref="Pi2I4MemoryLoopClaim"/>: provides <c>i = i^1</c>
///         (ImaginaryUnit in ±2c_l·γ·i; same anchor as F80's IFactor).</item>
///   <item><see cref="F1Pi2Inheritance"/>: M IS F1's residual operator;
///         F78 reads its single-body additive structure.</item>
/// </list>
///
/// <para>Tier consistency: F78 is Tier 1 proven (Master Lemma + per-site
/// additive structure + direct M_l matrix computation).</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="Pi2I4MemoryLoopRegistration.RegisterPi2I4MemoryLoop"/> +
/// <see cref="F1Pi2InheritanceRegistration.RegisterF1Pi2Inheritance"/>.</para></summary>
public static class F78SingleBodyMAdditivePi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF78SingleBodyMAdditivePi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F78SingleBodyMAdditivePi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var loop = b.Get<Pi2I4MemoryLoopClaim>();
            var f1 = b.Get<F1Pi2Inheritance>();
            return new F78SingleBodyMAdditivePi2Inheritance(ladder, loop, f1);
        });
}
