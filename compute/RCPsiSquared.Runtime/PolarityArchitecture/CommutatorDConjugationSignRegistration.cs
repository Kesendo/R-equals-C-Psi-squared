using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="CommutatorDConjugationSign"/> — F114, the
/// closed-form sign functional ε(σ) = (−1)^{n_Y(σ) + 1} for the action of D-conjugation
/// on the H-commutator superoperator L_σ = −i[σ, ·] in the 4^N Pauli basis.
///
/// <para>Tier1Derived (closed form + bit-exact verification N = 1..4 across 84 single
/// Pauli strings + 18 multi-term cases). Ctor parent
/// <see cref="Pi2KleinV4DephaseSwapGroup"/> (Welle 12) — F114 uses the diagonal D
/// from the Klein-V₄ Claim and reuses its universal-N reduction. Together they
/// characterize D's action on the two main dephase-letter-sensitive structures
/// (Π and L_H).</para>
///
/// <para>F114 published in <c>docs/ANALYTICAL_FORMULAS.md</c>; verifier
/// <c>simulations/_m_level_sign_functional_explore.py</c> (2026-05-27). Surfaced
/// during Welle 15 Task A polish (commit a98fc02) where substantive M_anti
/// equivariance tests at XZ+ZX and YZ+ZY bonds revealed bond-specific sign
/// behavior; this Claim is the closed-form characterization.</para></summary>
public static class CommutatorDConjugationSignRegistration
{
    public static ClaimRegistryBuilder RegisterCommutatorDConjugationSign(
        this ClaimRegistryBuilder builder) =>
        builder.Register<CommutatorDConjugationSign>(b =>
            new CommutatorDConjugationSign(b.Get<Pi2KleinV4DephaseSwapGroup>()));
}
