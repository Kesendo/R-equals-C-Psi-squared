using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F49cShadowCrossingPi2Inheritance"/>: F49c's
/// <c>R(N) = √((N − 1) / (N · 4^(N−1)))</c> closed form's denominator factor as
/// Pi2-Foundation operator-space inheritance. Two parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_{3−2N} = 4^(N−1)</c>
///         (the denominator factor in R(N)²).</item>
///   <item><see cref="Pi2OperatorSpaceMirrorClaim"/>: cross-pinning of <c>4^(N−1)</c>
///         as d² for (N − 1) qubits.</item>
/// </list>
///
/// <para>Tier consistency: F49c is Tier 1 proven; the Pi2-Foundation anchoring
/// is algebraic-trivial composition. All three claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F88PopcountCoherenceRegistration.RegisterF88PopcountCoherence"/> +
/// <see cref="F88StaticDyadicAnchorRegistration.RegisterF88StaticDyadicAnchor"/> +
/// <see cref="Pi2OperatorSpaceMirrorRegistration.RegisterPi2OperatorSpaceMirror"/>
/// in the builder pipeline.</para></summary>
public static class F49cShadowCrossingPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF49cShadowCrossingPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F49cShadowCrossingPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var mirror = b.Get<Pi2OperatorSpaceMirrorClaim>();
            return new F49cShadowCrossingPi2Inheritance(ladder, mirror);
        });
}
