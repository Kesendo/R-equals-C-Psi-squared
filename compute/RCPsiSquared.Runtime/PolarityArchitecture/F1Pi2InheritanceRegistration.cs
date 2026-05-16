using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F1Pi2Inheritance"/>: F1's "2" coefficient
/// in <c>Π·L·Π⁻¹ = −L − 2Σγ·I</c> as the Pi2-Foundation's <c>a_0 = d</c>, plus the
/// Z₄ memory-loop reading of the "−1" sign flip (i² from
/// <see cref="Pi2I4MemoryLoopClaim"/>). Three parent edges:
///
/// <list type="bullet">
///   <item><see cref="F1PalindromeIdentity"/>: the F1 closed form itself.</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2 = d</c>.</item>
///   <item><see cref="Pi2I4MemoryLoopClaim"/>: Z₄ closure of Π; the "−1" sign flip
///         in F1's "−L" reading IS i² = two 90° rotations summed. F1 is named
///         Layer 1 in the Pi2I4MemoryLoop docstring; this edge makes the
///         documented inheritance typed.</item>
/// </list>
///
/// <para>Tier consistency: all four Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="F1FamilyRegistration.RegisterF1Family"/> (for
/// F1PalindromeIdentity) +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="Pi2I4MemoryLoopRegistration.RegisterPi2I4MemoryLoop"/>.</para></summary>
public static class F1Pi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF1Pi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F1Pi2Inheritance>(b =>
        {
            var f1 = b.Get<F1PalindromeIdentity>();
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var memoryLoop = b.Get<Pi2I4MemoryLoopClaim>();
            return new F1Pi2Inheritance(f1, ladder, memoryLoop);
        });
}
