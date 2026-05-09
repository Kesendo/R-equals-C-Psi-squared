using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F49bCenteredDissipatorPi2Inheritance"/>:
/// F49b's <c>‖L_Dc‖² = γ²·4^N·N</c> closed form's "4^N" factor as Pi2-Foundation
/// operator-space inheritance. Two parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_{1−2N} = 4^N</c>.</item>
///   <item><see cref="Pi2OperatorSpaceMirrorClaim"/>: direct cross-pinning of
///         <c>4^N = d²</c> for N qubits (no derivation overhead, the most direct
///         operator-space anchor of any F-formula in the inheritance graph).</item>
/// </list>
///
/// <para>Tier consistency: F49b is Tier 1 proven (Lemma 1 of PROOF_CROSS_TERM_FORMULA);
/// the Pi2-Foundation anchoring is algebraic-trivial composition. All three
/// claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F88PopcountCoherenceRegistration.RegisterF88PopcountCoherence"/> +
/// <see cref="F88StaticDyadicAnchorRegistration.RegisterF88StaticDyadicAnchor"/> +
/// <see cref="Pi2OperatorSpaceMirrorRegistration.RegisterPi2OperatorSpaceMirror"/>
/// in the builder pipeline.</para></summary>
public static class F49bCenteredDissipatorPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF49bCenteredDissipatorPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F49bCenteredDissipatorPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var mirror = b.Get<Pi2OperatorSpaceMirrorClaim>();
            return new F49bCenteredDissipatorPi2Inheritance(ladder, mirror);
        });
}
