using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F91Pi2Inheritance"/>: F91's anti-palindromic
/// γ orbit (γ_l + γ_{N−1−l} = 2·γ_avg) as the parameter-side γ-axis instance of the
/// Pi2-Z₄ rotational structure. One parent edge:
///
/// <list type="bullet">
///   <item><see cref="Pi2I4MemoryLoopClaim"/>: provides the i⁴ = 1 closure that
///         generates the Z₄. F91's 90°-rotation γ ↦ 2·γ_avg − F71(γ) closes at
///         order 4, structurally the same Z₄ that
///         <see cref="NinetyDegreeMirrorMemoryClaim"/> types on the
///         operator-quaternion side.</item>
/// </list>
///
/// <para>Tier consistency: both Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2I4MemoryLoopRegistration.RegisterPi2I4MemoryLoop"/>
/// in the builder pipeline.</para></summary>
public static class F91Pi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF91Pi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F91Pi2Inheritance>(b =>
        {
            var memoryLoop = b.Get<Pi2I4MemoryLoopClaim>();   // Z₄ closure i⁴ = 1
            return new F91Pi2Inheritance(memoryLoop);
        });
}
