using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F92Pi2Inheritance"/>: F92's anti-palindromic
/// J orbit (J_b + J_{N−2−b} = 2·J_avg) as the parameter-side J-axis instance of the
/// Pi2-Z₄ rotational structure. One parent edge:
///
/// <list type="bullet">
///   <item><see cref="Pi2I4MemoryLoopClaim"/>: provides the i⁴ = 1 closure that
///         generates the Z₄. F92's 90°-rotation J ↦ 2·J_avg − F71(J) closes at
///         order 4, structurally the same Z₄ that
///         <see cref="NinetyDegreeMirrorMemoryClaim"/> types on the
///         operator-quaternion side.</item>
/// </list>
///
/// <para>Tier consistency: both Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2I4MemoryLoopRegistration.RegisterPi2I4MemoryLoop"/>
/// in the builder pipeline.</para></summary>
public static class F92Pi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF92Pi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F92Pi2Inheritance>(b =>
        {
            var memoryLoop = b.Get<Pi2I4MemoryLoopClaim>();   // Z₄ closure i⁴ = 1
            return new F92Pi2Inheritance(memoryLoop);
        });
}
