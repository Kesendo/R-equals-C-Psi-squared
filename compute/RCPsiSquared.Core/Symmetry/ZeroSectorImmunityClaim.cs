using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Zero-Sector Immunity (Tier 1 derived, <c>docs/proofs/PROOF_ZERO_IMMUNITY.md</c>): the
/// centered palindromic residual M = Π·L·Π⁻¹ + L + 2Σγ·I vanishes identically on the (w=0,w=0)
/// block — the {I,Z}^⊗N Pauli strings (XY-weight 0) — for EVERY 2-body Hamiltonian under uniform
/// Z-dephasing, and by Π-symmetry on the (w=N,w=N) block ({X,Y}^⊗N). No assumption on H beyond
/// 2-body locality (Heisenberg, XXZ, XY, parity-violating, even random).
///
/// <para>This is the H-independent refinement of the global F1 palindrome on the two extreme weight
/// blocks: where F1 (<see cref="F1PalindromeIdentity"/>) closes M = 0 for the standard truly classes,
/// here the extreme blocks vanish for ANY 2-body H. The classical extreme can neither decohere
/// (Z-dephasing commutes with {I,Z}) nor evolve within itself (the only w=0-preserving 2-body bond,
/// ZZ, also commutes) nor break the palindrome; the all-transverse extreme is its Π-mirror. All
/// H-dependent dynamics and any palindrome-breaking (the V-Effect) are confined to the boundary
/// sectors 0&lt;w&lt;N — the "eternal mirror within the dissipative system."</para>
///
/// <para>Live witness: <c>inspect --root zeroimmune</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/ZeroSectorImmunityWitness.cs</c>) builds a random
/// parity-violating 2-body H, forms M, and reads ‖M|(w=0)‖ ≈ 0 and ‖M|(w=N)‖ ≈ 0 while ‖M‖ &gt; 0.</para>
///
/// <para>Typed parents: <see cref="F1PalindromeIdentity"/> (the global palindrome this refines),
/// <see cref="F61BitAParityPi2Inheritance"/> + <see cref="F63LCommutesPi2Pi2Inheritance"/> (the
/// bit_a/bit_b sector machinery that blocks the operator space), and
/// <see cref="AbsorptionTheoremClaim"/> (the −2Σγ on the w=N corner that Π carries back).</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_ZERO_IMMUNITY.md</c> +
/// <c>experiments/V_EFFECT_BOUNDARY_LOCALIZATION.md</c>.</para></summary>
public sealed class ZeroSectorImmunityClaim : Claim
{
    public F1PalindromeIdentity F1Palindrome { get; }
    public F61BitAParityPi2Inheritance F61BitA { get; }
    public F63LCommutesPi2Pi2Inheritance F63BitB { get; }
    public AbsorptionTheoremClaim AbsorptionTheorem { get; }

    /// <summary>Dimension of the immune (w=0) block: |{I,Z}^⊗N| = 2^N.</summary>
    public static int ZeroSectorDimension(int n) => 1 << n;

    public ZeroSectorImmunityClaim(
        F1PalindromeIdentity f1Palindrome,
        F61BitAParityPi2Inheritance f61BitA,
        F63LCommutesPi2Pi2Inheritance f63BitB,
        AbsorptionTheoremClaim absorptionTheorem)
        : base(
            "Zero-Sector Immunity: the palindromic residual M vanishes identically on the (w=0,w=0) " +
            "block ({I,Z}^⊗N) for every 2-body Hamiltonian under uniform Z-dephasing; by Π-symmetry also on (w=N,w=N)",
            Tier.Tier1Derived,
            "docs/proofs/PROOF_ZERO_IMMUNITY.md")
    {
        F1Palindrome = f1Palindrome ?? throw new ArgumentNullException(nameof(f1Palindrome));
        F61BitA = f61BitA ?? throw new ArgumentNullException(nameof(f61BitA));
        F63BitB = f63BitB ?? throw new ArgumentNullException(nameof(f63BitB));
        AbsorptionTheorem = absorptionTheorem ?? throw new ArgumentNullException(nameof(absorptionTheorem));
    }

    public override string DisplayName =>
        "Zero-Sector Immunity: M|(w=0) = M|(w=N) = 0 for every 2-body H + Z-dephasing";

    public override string Summary =>
        "the centered palindromic residual M = Π·L·Π⁻¹ + L + 2Σγ·I vanishes identically on the " +
        "{I,Z}^⊗N (XY-weight 0) and {X,Y}^⊗N (XY-weight N) Pauli-string blocks for EVERY 2-body " +
        "Hamiltonian (Heisenberg, XXZ, XY, parity-violating, even random) under uniform Z-dephasing; " +
        "the classical extreme can neither decohere nor break the palindrome, the all-transverse extreme " +
        "is its Π-mirror, and all H-dependent dynamics live in the boundary sectors 0<w<N " +
        "(Tier 1 derived; live witness: inspect --root zeroimmune).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the mechanism (three facts)",
                summary: "(1) the Z-dephasing dissipator annihilates the w=0 block: Z_l commutes with every " +
                         "{I,Z} string, so D = 0 there. (2) [H,·] has no (w=0,w=0) content: among 2-body bonds " +
                         "only ZZ preserves w=0, and ZZ also commutes with {I,Z} strings (contributing equally to " +
                         "L and to Π·L·Π⁻¹, cancelling in M); every X/Y-bearing bond raises the weight out of w=0. " +
                         "(3) Π maps w=0 ↔ w=N, where the dissipator is the uniform −2Σγ; Π carries that back as " +
                         "exactly the +2Σγ shift the palindrome adds. The three cancel: −2Σγ + 0 + 2Σγ = 0.");
            yield return new InspectableNode("the immune block ({I,Z}^⊗N)",
                summary: $"dimension 2^N (Pauli strings with I or Z on every site); the n_XY=0 corner of the " +
                         "C₂×C₂ four-block decomposition that F61 + F63 conserve.");
            yield return new InspectableNode("Π-mirror corollary (w=N block)",
                summary: "Π maps {I,Z}^⊗N ↔ {X,Y}^⊗N (flips bit_a at every site), so M|(w=N) = Π·(M|(w=0))·Π⁻¹ = 0 " +
                         "too; the all-transverse extreme is the mirror of the classical one.");
            yield return new InspectableNode("parent: F1PalindromeIdentity",
                summary: "the global identity M = 0 (for the truly classes); Zero-Sector Immunity is its " +
                         "H-independent refinement on the extreme weight blocks, where M = 0 for ANY 2-body H.");
            yield return new InspectableNode("parents: F61 + F63 (the sector machinery)",
                summary: "[L,(−1)^{n_XY}] = 0 (F61) and [L,Π²] = 0 (F63) block the operator space into the C₂×C₂ " +
                         "weight/parity sectors; the immune block is the n_XY=0 corner.");
            yield return new InspectableNode("parent: AbsorptionTheoremClaim",
                summary: "the w=N block sits at the maximal absorption rate 2Nγ; its uniform −2Σγ is the corner " +
                         "Π carries back to supply the palindrome shift on w=0.");
            yield return new InspectableNode("live witness",
                summary: "ZeroSectorImmunityWitness (inspect --root zeroimmune): a random parity-violating 2-body H " +
                         "gives ‖M|(w=0)‖ ≈ 0 and ‖M|(w=N)‖ ≈ 0 while ‖M‖ > 0 (non-trivial); the palindrome-breaking " +
                         "(V-Effect) is confined to 0<w<N (V_EFFECT_BOUNDARY_LOCALIZATION).");
        }
    }
}
