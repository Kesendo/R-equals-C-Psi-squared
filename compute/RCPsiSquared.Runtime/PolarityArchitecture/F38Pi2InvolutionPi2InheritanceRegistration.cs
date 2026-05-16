using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F38Pi2InvolutionPi2Inheritance"/>: F38's
/// <c>Π² = (−1)^w_YZ</c> master operator identity as Pi2-Foundation
/// operator-space bisection. Four typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 · a_{3−2N} =
///         2 · 4^(N−1)</c> = each Π² eigenspace dimension.</item>
///   <item><see cref="Pi2OperatorSpaceMirrorClaim"/>: cross-pinning of <c>4^N</c>
///         as d² for N qubits.</item>
///   <item><see cref="Pi2I4MemoryLoopClaim"/>: Z₄ closure Π⁴ = I → Π²
///         involutivity → Π²-eigenvalues {+1, −1}.</item>
///   <item><see cref="HalfAsStructuralFixedPointClaim"/>: typed parent
///         grounding F38's <c>HalfHalfBalance = EigenspaceDim / FullDim = 1/2</c>
///         exactly. Π² partitions 4^N symmetrically because <c>a_2 = 1/2</c>
///         is the self-dual cell on the dyadic ladder. Added 2026-05-16
///         (previously prose-only docstring mention).</item>
/// </list>
///
/// <para>Tier consistency: F38 is Tier 1 proven (PT_SYMMETRY_ANALYSIS,
/// PROOF_BIT_B_PARITY_SYMMETRY); the Pi2-Foundation anchoring is
/// algebraic-trivial composition. All four claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F88PopcountCoherenceRegistration.RegisterF88PopcountCoherence"/> +
/// <see cref="F88StaticDyadicAnchorRegistration.RegisterF88StaticDyadicAnchor"/> +
/// <see cref="Pi2OperatorSpaceMirrorRegistration.RegisterPi2OperatorSpaceMirror"/> +
/// <see cref="Pi2I4MemoryLoopRegistration.RegisterPi2I4MemoryLoop"/>
/// in the builder pipeline.</para></summary>
public static class F38Pi2InvolutionPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF38Pi2InvolutionPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F38Pi2InvolutionPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var mirror = b.Get<Pi2OperatorSpaceMirrorClaim>();
            var memoryLoop = b.Get<Pi2I4MemoryLoopClaim>();
            var half = b.Get<HalfAsStructuralFixedPointClaim>();
            return new F38Pi2InvolutionPi2Inheritance(ladder, mirror, memoryLoop, half);
        });
}
