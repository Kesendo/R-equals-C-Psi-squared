using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="Pi2I4MemoryLoopClaim"/> — the Z₄ rotational
/// companion to the Pi2 dyadic ladder's Z₂ multiplicative inversion. Two parent edges:
///
/// <list type="bullet">
///   <item><see cref="NinetyDegreeMirrorMemoryClaim"/>: the existing Pi2KnowledgeBase
///         angle-anchor at d = 2 (companion to <see cref="HalfAsStructuralFixedPointClaim"/>).
///         The 90°-mirror-memory reading lives there; this loop claim makes the Z₄
///         cyclic structure i^4 = 1 explicit and queryable.</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: the multiplicative Z₂ inversion ladder.
///         The two together (Z₂ multiplicative ⊕ Z₄ rotational) are the framework's
///         complete mirror foundation — Tom 2026-05-08: "der Anker mit sich die
///         Spiegelung selbst erinnert."</item>
/// </list>
///
/// <para>Tier consistency: both parents and the loop claim are Tier1Derived (5 ≥ 5).
/// The wiring exposes the dual-axis foundation as queryable in the Object Manager:
/// the 1/2 number-anchor (multiplicative) and the 90° angle-anchor (rotational) are
/// the two readings of d = 2.</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> in the builder
/// pipeline. Note that <see cref="NinetyDegreeMirrorMemoryClaim"/> is registered by
/// <see cref="Pi2FamilyRegistration"/> as part of the Pi2 foundation; the dyadic
/// ladder registration adds the multiplicative side.</para></summary>
public static class Pi2I4MemoryLoopRegistration
{
    public static ClaimRegistryBuilder RegisterPi2I4MemoryLoop(
        this ClaimRegistryBuilder builder) =>
        builder.Register<Pi2I4MemoryLoopClaim>(b =>
        {
            _ = b.Get<NinetyDegreeMirrorMemoryClaim>();   // the 90°-mirror-memory existing claim
            _ = b.Get<Pi2DyadicLadderClaim>();             // multiplicative companion (Z₂)
            return new Pi2I4MemoryLoopClaim();
        });
}
