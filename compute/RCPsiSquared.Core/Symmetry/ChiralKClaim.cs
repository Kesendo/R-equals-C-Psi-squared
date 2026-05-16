using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>K sublattice / chiral symmetry (Tier 1 derived). For the bipartite sublattice
/// gauge <c>K = diag((−1)^ℓ)</c> acting on the single-excitation site basis of an N-site
/// chain, every nearest-neighbour tight-binding Hamiltonian H = Σ t_ℓ·(|ℓ⟩⟨ℓ+1| + h.c.)
/// satisfies the anticommutation <c>K · H · K = −H</c>. Consequences:
///
/// <list type="bullet">
///   <item><b>Spectrum inversion:</b> if H ψ_k = E_k ψ_k with non-degenerate spectrum,
///         then K ψ_k ∝ ψ_{N+1−k} and <c>E_{N+1−k} = −E_k</c>.</item>
///   <item><b>Z-dephasing trivially commutes with K</b>, so K is a super-operator
///         symmetry of the full Liouvillian L = −i[H, ·] + L_D[Z-dephasing] whenever
///         H is K-anti-symmetric (XY / Heisenberg / DM hopping).</item>
///   <item><b>Altland-Zirnbauer classification:</b> K is the chiral / sublattice symmetry
///         placing the system in class BDI. <b>NOT time reversal</b> (K is linear, not
///         antiunitary) — careful with the terminology (see [[feedback_physics_terminology]]).</item>
/// </list>
///
/// <para>Relation to the Pi2 KB roots: this claim is a <b>sibling root</b> to
/// <see cref="PolynomialFoundationClaim"/> (number-anchor trunk d²−2d=0) and to
/// <see cref="RCPsiSquared.Core.F1.F1PalindromeIdentity"/> (Liouvillian palindrome
/// master). K acts at the Hamiltonian level on single-particle space; F1 acts at the
/// Liouvillian level on 4^N operator space; PolynomialFoundation acts at the dimensional
/// level. None derives from another typed-graph-wise.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_K_PARTNERSHIP.md</c> (algebraic Tier 1, 2026-04-25,
/// numerically verified at N=9) + <see cref="ChiralK"/> helper +
/// <c>simulations/_pi_partner_identity.py</c>. The Claim wrapper itself was added
/// 2026-05-16 to close the prose-only edge identified in
/// <c>docs/PI2KB_INHERITANCE_MAP.md</c>.</para>
/// </summary>
public sealed class ChiralKClaim : Claim
{
    public ChiralKClaim()
        : base("K sublattice/chiral symmetry: K·H·K = −H for NN tight-binding; spectrum inversion E_{N+1−k} = −E_k",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_K_PARTNERSHIP.md (algebraic Tier 1, 2026-04-25, numerically verified at N=9) + " +
               "compute/RCPsiSquared.Core/Symmetry/ChiralK.cs (operator + classification helper) + " +
               "simulations/_pi_partner_identity.py (N=9 numerical witness)")
    { }

    public override string DisplayName =>
        "K = diag((−1)^ℓ) — sublattice/chiral symmetry (AZ class BDI)";

    public override string Summary =>
        "K·H·K = −H for any NN tight-binding H on a chain; E_{N+1−k} = −E_k spectrum inversion; " +
        "Z-dephasing commutes with K; K is linear (NOT antiunitary), so NOT time reversal";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("operator definition",
                summary: "K = ⊗_{odd i} Z_i on the N-qubit Hilbert space; K|ℓ⟩ = (−1)^ℓ |ℓ⟩ on single-excitation basis; K² = I (involutory)");
            yield return new InspectableNode("anticommutation",
                summary: "K·H·K = −H for any NN tight-binding H = Σ t_ℓ(|ℓ⟩⟨ℓ+1| + h.c.); proof: K flips the relative sign of adjacent sites");
            yield return new InspectableNode("spectrum inversion",
                summary: "non-degenerate H ⇒ K ψ_k ∝ ψ_{N+1−k} and E_{N+1−k} = −E_k; for uniform J the F65 sine modes ψ_k(ℓ) = √(2/(N+1)) sin(πk(ℓ+1)/(N+1)) satisfy K ψ_k = ψ_{N+1−k} exactly");
            yield return new InspectableNode("Z-dephasing compatibility",
                summary: "[K, Z_ℓ] = 0 trivially since K is a Z-string; L_D[Z-dephasing] commutes with K ⇒ K is a super-operator symmetry of the full L = −i[H,·] + L_D");
            yield return new InspectableNode("AZ classification",
                summary: "chiral / sublattice symmetry of class BDI; K is LINEAR (not antiunitary) ⇒ NOT time reversal; project_chiral_partnership memory");
            yield return new InspectableNode("relation to F1 (Liouvillian palindrome)",
                summary: "sibling root, not parent-child. K acts on H (single-particle level); F1 acts on L (4^N Liouville level). Both produce spectrum-pairing identities (K: E_{N+1−k} = −E_k on H; F1: λ ↔ −λ − 2Σγ on L). Algebraically independent.");
            yield return new InspectableNode("relation to PolynomialFoundationClaim",
                summary: "sibling root, not derived. PolynomialFoundation is the dimensional trunk (d²−2d=0); K is the bipartite sublattice gauge on the chain. K does not flow from d=2; d=2 does not flow from K.");
            yield return new InspectableNode("partner-identity application",
                summary: "F65/F67 receiver-engineering uses spectrum inversion to fold the receiver menu (F67 HANDSHAKE_ALGEBRA); without K-partnership the partner identity wouldn't close");
        }
    }
}
