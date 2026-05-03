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
/// identically for every Heisenberg/XY/Ising/XXZ/DM Hamiltonian under uniform Z-dephasing.
/// Every Liouvillian eigenvalue λ pairs with −λ − 2Σγ; every decay rate d pairs with
/// 2Σγ − d.
///
/// <para>Verified analytically (<c>docs/proofs/MIRROR_SYMMETRY_PROOF.md</c>) and replaces
/// the brute-force palindrome verification at N=8 (54,118 eigenvalues).</para>
///
/// <para>Validity: Heisenberg, XY, Ising, XXZ, DM Hamiltonians on any graph; uniform or
/// site-dependent Z-dephasing; any N. Two Π families exist (P1, P4). Breaks for
/// depolarizing noise (residual error = (2/3)Σγ, linear in γ and N).</para>
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
                summary: "Heisenberg/XY/Ising/XXZ/DM on any graph; uniform or site-dependent Z-dephasing");
            yield return new InspectableNode("breaks for",
                summary: "depolarizing noise → residual error (2/3)Σγ, linear in γ and N");
            yield return new InspectableNode("verification",
                summary: "replaces brute-force palindrome scan at N=8 (54,118 eigenvalues, zero exceptions)");
        }
    }
}
