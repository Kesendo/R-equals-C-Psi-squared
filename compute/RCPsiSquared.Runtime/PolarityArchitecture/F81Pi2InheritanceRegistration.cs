using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F81Pi2Inheritance"/>: F81's "2" and "1/2"
/// coefficients as Pi2-Foundation mirror-pair (a_0, a_2), the "Mirror Space"
/// connection through Pi2OperatorSpaceMirror, and the Z₄ generator reading of the
/// Π operator. Four parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides a_0 = 2 and a_2 = 1/2.</item>
///   <item><see cref="HalfAsStructuralFixedPointClaim"/>: typed 1/2 anchor (the
///         50/50 split's structural fixed point).</item>
///   <item><see cref="Pi2OperatorSpaceMirrorClaim"/>: the operator-space dimension
///         d² = 4^N where F81's Π·M·Π⁻¹ action lives, the "Mirror Space" connection
///         Tom 2026-05-08 asked about.</item>
///   <item><see cref="Pi2I4MemoryLoopClaim"/>: Z₄ closure Π⁴ = I; F81's Π is the
///         Z₄ generator that order-4 closes, same algebra as F38 / F80. This edge
///         makes the Π-operator-family inheritance consistent with F38 and F80.</item>
/// </list>
///
/// <para>Tier consistency: all five Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F88PopcountCoherenceRegistration.RegisterF88PopcountCoherence"/> +
/// <see cref="F88StaticDyadicAnchorRegistration.RegisterF88StaticDyadicAnchor"/> +
/// <see cref="Pi2OperatorSpaceMirrorRegistration.RegisterPi2OperatorSpaceMirror"/> +
/// <see cref="Pi2I4MemoryLoopRegistration.RegisterPi2I4MemoryLoop"/>
/// in the builder pipeline.</para></summary>
public static class F81Pi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF81Pi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F81Pi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            _ = b.Get<HalfAsStructuralFixedPointClaim>();          // typed 1/2 anchor
            var mirror = b.Get<Pi2OperatorSpaceMirrorClaim>();      // Mirror Space connection
            var memoryLoop = b.Get<Pi2I4MemoryLoopClaim>();         // Z₄ closure of Π
            return new F81Pi2Inheritance(ladder, mirror, memoryLoop);
        });
}
