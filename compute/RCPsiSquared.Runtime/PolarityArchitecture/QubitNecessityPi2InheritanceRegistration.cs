using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="QubitNecessityPi2Inheritance"/>: the
/// per-site Pauli basis split <c>4 = 2 + 2</c> (the typed reading of d² − 2d = 0
/// at d=2) inherits from four Pi2-Foundation claims:
///
/// <list type="bullet">
///   <item><see cref="PolynomialFoundationClaim"/>: the d² − 2d = 0 trunk that
///         enforces the bijection condition.</item>
///   <item><see cref="HalfAsStructuralFixedPointClaim"/>: the C = 0.5 balanced
///         fraction at the per-site basis level.</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c> (immune count)
///         and <c>a_{-1} = 4</c> (total Pauli count).</item>
///   <item><see cref="Pi2OperatorSpaceMirrorClaim"/>: provides the d² = 4 for N=1
///         pinning, cross-verifying the ladder reading.</item>
/// </list>
///
/// <para>Tier consistency: all five Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F88PopcountCoherenceRegistration.RegisterF88PopcountCoherence"/> +
/// <see cref="F88StaticDyadicAnchorRegistration.RegisterF88StaticDyadicAnchor"/> +
/// <see cref="Pi2OperatorSpaceMirrorRegistration.RegisterPi2OperatorSpaceMirror"/>
/// in the builder pipeline.</para></summary>
public static class QubitNecessityPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterQubitNecessityPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<QubitNecessityPi2Inheritance>(b =>
        {
            _ = b.Get<PolynomialFoundationClaim>();             // d² − 2d = 0 trunk
            _ = b.Get<HalfAsStructuralFixedPointClaim>();       // C = 0.5 balanced
            var ladder = b.Get<Pi2DyadicLadderClaim>();         // a_0, a_{-1}
            var mirror = b.Get<Pi2OperatorSpaceMirrorClaim>();   // d² for N=1 pinning
            return new QubitNecessityPi2Inheritance(ladder, mirror);
        });
}
