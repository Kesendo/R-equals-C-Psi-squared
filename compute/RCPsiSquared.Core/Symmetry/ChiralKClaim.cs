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
///         H is K-anti-symmetric (XY / Heisenberg / DM <b>hopping</b>: the hopping part is what
///         anticommutes, and the Heisenberg single-excitation matrix WITH its ZZ diagonal does
///         not, see the F143 paragraph below).</item>
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
/// <para><b>The involution F143 reads its kernel in.</b> K acts on the modes as the reflection
/// k ↦ M − k with M = N + 1 (the bullet above; the map is the permutation of modes, the sign
/// depending on the sine convention), which is exactly the R of
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F143. That entry's seed-rung Gram matrix
/// G = (1/M)·(𝟏𝟏ᵀ + (I + R)/2) has its kernel exactly in the R-odd sector, so "the frozen kernel
/// is the chiral-odd sector" is a statement about THIS K. The reflection needs only K h K = −h and
/// a non-degenerate spectrum, so it survives bond disorder; what needs UNIFORM hopping is F143's
/// closed form, the DST-I modes being where that form is written.</para>
///
/// <para>What is carried here is the involution, and this paragraph is a breadcrumb rather than a
/// claim of ownership: F143's mode-side reading, spec(G) as such, has no typed carrier. Its two
/// numbers do, one block over on the site side, where at M = N + 1 the same matrix appears as B Bᵀ:
/// <see cref="RCPsiSquared.Core.BlockSpectrum.FrozenDivisorClaim"/> types the ⌊N/2⌋ count as
/// <c>FrozenMultiplicity</c>, and <c>MirrorWorld.Divisor.ClockModulus</c> types M for both chains.
/// Neither is read off a matrix; both are adopted as numbers, which is why they carry the count
/// without carrying this reading of it.</para>
///
/// <para>The MODE reading is XY-specific, and the site-index one is not this operator. F143 records
/// that on the Heisenberg chain the modes run k = 0..N−1, where k ↦ N − k is not an involution, so
/// there is no mode-chiral-odd sector to be odd in; what carries is the SITE-index chirality, and
/// there the involution is the site reversal l ↦ N + 1 − l, NOT this K, and NOT M − l either: the
/// site reversal is the same map on both chains, while M is N on the Heisenberg one. On the Heisenberg
/// single-excitation matrix K is not an antisymmetry at all: that spectrum,
/// λ_k = 4cos(kπ/N) + N − 5, is unpaired.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_K_PARTNERSHIP.md</c> (algebraic Tier 1, 2026-04-25,
/// numerically verified at N=9) + <see cref="ChiralK"/> helper +
/// <c>simulations/pi_partner_identity.py</c>. The Claim wrapper itself was added
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
               "simulations/pi_partner_identity.py (N=9 numerical witness)")
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
            yield return new InspectableNode("F143 uses this involution",
                summary: "K acts on the modes as k ↦ M − k, M = N + 1, which is F143's R (the reflection needs only KhK = −h and non-degeneracy, so bond disorder does not break it; uniform hopping is what F143's CLOSED FORM needs). F143's seed-rung Gram matrix G = (1/M)(𝟏𝟏ᵀ + (I + R)/2) has kernel exactly the R-odd sector, dimension ⌊N/2⌋, with gap 1/M above it. Carried here: the involution. Not carried: spec(G) as a mode-side object; its two numbers are typed one block over on the site side, where at M = N + 1 the same matrix appears as B Bᵀ, ⌊N/2⌋ as FrozenDivisorClaim.FrozenMultiplicity and M as MirrorWorld's Divisor.ClockModulus (both adopted as numbers, not read off a matrix). The MODE reading is XY-only, and the site-index chirality that carries to the Heisenberg chain is the site reversal l ↦ N + 1 − l (the same map on both chains; M is N on the Heisenberg one, so M − l is NOT it), not this K.");
            yield return new InspectableNode("partner-identity application",
                summary: "F65/F67 receiver-engineering uses spectrum inversion to fold the receiver menu (F67 HANDSHAKE_ALGEBRA); without K-partnership the partner identity wouldn't close");
        }
    }
}
