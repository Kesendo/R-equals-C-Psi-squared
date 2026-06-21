using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F93Pi2Inheritance"/>: F93's anti-palindromic
/// h orbit (h_l + h_{N−1−l} = 2·h_avg) as the parameter-side h-axis instance of the
/// Pi2 rotational structure (a Klein V₄ on parameters, the order-2 shadow of the
/// operator-side Z₄). One parent edge:
///
/// <list type="bullet">
///   <item><see cref="Pi2I4MemoryLoopClaim"/>: provides the operator-side i⁴ = 1
///         closure (the genuine order-4 Z₄). F93's reshuffle R_{90} : h ↦ 2·h_avg − F71(h)
///         is an INVOLUTION (R_{90}² = identity, order 2) that inherits from that
///         operator-side Z₄ which
///         <see cref="NinetyDegreeMirrorMemoryClaim"/> types on the
///         operator-quaternion side.</item>
/// </list>
///
/// <para>Tier consistency: both Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2I4MemoryLoopRegistration.RegisterPi2I4MemoryLoop"/>
/// in the builder pipeline.</para></summary>
public static class F93Pi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF93Pi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F93Pi2Inheritance>(b =>
        {
            var memoryLoop = b.Get<Pi2I4MemoryLoopClaim>();   // Z₄ closure i⁴ = 1
            return new F93Pi2Inheritance(memoryLoop);
        });
}
