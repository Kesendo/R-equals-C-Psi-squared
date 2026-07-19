using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F1;

/// <summary>F1 (Tier 1, proven): the Liouvillian palindrome identity.
///
/// <code>
///     Π · L · Π⁻¹  =  −L − 2Σγ · I
/// </code>
///
/// Equivalently, the per-Pauli-string residual <c>M = Π·L·Π⁻¹ + L + 2Σγ·I</c> vanishes
/// identically for every Heisenberg/XY/Ising/XXZ Hamiltonian under Z-dephasing, uniform
/// or site-dependent.
/// Every Liouvillian eigenvalue λ pairs with −λ − 2Σγ; every decay rate d pairs with
/// 2Σγ − d.
///
/// <para>Verified analytically (<c>docs/proofs/MIRROR_SYMMETRY_PROOF.md</c>) and replaces
/// the brute-force palindrome verification at N=8 (all 65,536 eigenvalues paired, zero
/// exceptions on every topology tested at that size: chain, star, ring, and K₄ plus a
/// disjoint 4-chain (disconnected, 9 bonds, run at γ=0.5); the binary tree is covered at
/// N=4 and N=5).</para>
///
/// <para>Validity: Heisenberg, XY, Ising, XXZ Hamiltonians on any graph, bond terms, PLUS any
/// TRANSVERSE on-site field (any direction in the XY plane against Z-dephasing), so the identity
/// covers strictly more than the bond-only statement. The field's DIRECTION must be the same
/// on every site (magnitudes are free); let it vary between sites and the palindrome itself
/// goes, 4/64 at N=3. Under this fixed Π an X-field reads
/// residual 0 and a Y-field reads 2·max|h|, but that is a gauge artifact: R_z(π/2) per site maps
/// one to the other and is a symmetry of the dephasing, the two Liouvillians are unitarily
/// equivalent with identical spectra, and the correspondingly rotated Π satisfies the identity
/// exactly. The residual under a fixed Π grades the ANGLE to that mirror's preferred axis. A
/// LONGITUDINAL h·ΣZ_i field is the genuine exception and breaks both the identity and the
/// pairing. Z-dephasing,
/// uniform or site-dependent; any N. Mixed dephasing axes across sites also stay palindromic
/// as long as at most two distinct axes appear in any one connected component of the graph.
/// Two Π families exist (P1, P4). The DM bond is palindromic too, but
/// under a site-alternating Π, not this uniform one, so it is NOT in this claim's validity
/// set (see <c>experiments/NON_HEISENBERG_PALINDROME.md</c>). Breaks for depolarizing
/// noise: with the standard convention (total rate γ per site = γ/3 per axis) the residual
/// is (2/3)Σγ, linear in γ and N, measured as the max entry of M IN THE PAULI-STRING BASIS
/// (a max entry is not a norm; the same operator reads (2/9)Σγ in the computational basis,
/// so the basis is part of the constant). Reading depolarizing as rate γ on each of the
/// three axes gives 6Σγ instead. As a function of the per-axis scale s the max residual is
/// piecewise, |−8Nγs + 2Nγ| for s ≥ 1/3 and |−4Nγs + 2Nγ| below; the registered constant sits
/// exactly on the crossing at s = 1/3, so do not extrapolate the first branch downward. Either way the palindrome fails outright, it does not
/// degrade gracefully: the eigenvalue pairing is 0/64 at N=3 in both conventions.</para>
/// </summary>
public sealed class F1PalindromeIdentity : Claim
{
    public F1PalindromeIdentity()
        : base("F1 Liouvillian palindrome: Π·L·Π⁻¹ = −L − 2Σγ·I",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F1 + docs/proofs/MIRROR_SYMMETRY_PROOF.md")
    { }

    public override string DisplayName => "F1: Π·L·Π⁻¹ = −L − 2Σγ·I (Liouvillian palindrome)";

    public override string Summary =>
        "every Liouvillian eigenvalue λ pairs with −λ − 2Σγ; every decay rate d pairs with 2Σγ − d";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("residual",
                summary: "M = Π·L·Π⁻¹ + L + 2Σγ·I vanishes identically for valid (H, D) classes");
            yield return new InspectableNode("Π operator",
                summary: "RCPsiSquared.Core.Symmetry.PiOperator: unitary order-4, signed permutation in Pauli basis");
            yield return new InspectableNode("validity",
                summary: "Heisenberg/XY/Ising/XXZ on any graph, bond terms plus a uniformly directed transverse on-site field; uniform or site-dependent Z-dephasing");
            yield return new InspectableNode("breaks for",
                summary: "depolarizing noise (γ/3 per axis) → max residual (2/3)Σγ in the Pauli basis, linear in γ and N; an on-site field not orthogonal to every dephasing axis in play; T1 alongside CO-AXIAL dephasing (transverse dephasing composes with T1 exactly) or T1 with any on-site field");
            yield return new InspectableNode("verification",
                summary: "replaces brute-force palindrome scan at N=8 (all 65,536 eigenvalues paired, zero exceptions; chain, star, ring, K₄ plus disjoint 4-chain)");
        }
    }
}
