using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="ThreeLadderHingeClaim"/> (2026-06-15): the three-ladder hinge
/// Q. The disagreement rung k = popcount(i⊕j), the F87 girth ℓ, and the F120 moment j are the two factors
/// of one F87-hardness coefficient on M = A + γQ, hinged by Q (its spectrum is the rung k = N−2k, its
/// action Σ Z_l⊗Z_l projects A = −i[H,·]'s closed walks onto the girth moments t_j = Tr(Z_l H^j)).
/// P_{m,1} = m·Tr(Q·A^{m−1}) = the girth moments at every rung; the rung is essential.
///
/// <para>Tier1Derived (exact, gate-first; self-check battery at N = 2 and N = 3 in the ctor). Typed
/// parents: <see cref="AbsorptionTheoremClaim"/> (the rung k = Q's spectrum) and
/// <see cref="MomentTowerPumpChannelClaim"/> (the moments t_j Q projects onto, the girth ℓ their onset),
/// both registered earlier in the chain. The F87 girth/hardness primitive lives in
/// RCPsiSquared.Diagnostics and is carried in prose. Anchor:
/// <c>docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md</c> §1, §4 +
/// <c>simulations/_three_ladders_bridge.py</c> +
/// <c>compute/RCPsiSquared.Diagnostics/Foundation/LadderHingeWitness.cs</c> (inspect --root ladders).</para></summary>
public static class ThreeLadderHingeClaimRegistration
{
    public static ClaimRegistryBuilder RegisterThreeLadderHingeClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<ThreeLadderHingeClaim>(b =>
            new ThreeLadderHingeClaim(
                b.Get<AbsorptionTheoremClaim>(),
                b.Get<MomentTowerPumpChannelClaim>()));
}
