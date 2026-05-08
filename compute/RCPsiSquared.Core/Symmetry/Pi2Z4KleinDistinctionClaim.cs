using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The framework has TWO four-fold mirror structures on the Pauli space, and
/// they are <b>not</b> the same group. Tom 2026-05-08 asked whether
/// <see cref="KleinFourCellClaim"/>'s four cells equal <see cref="Pi2I4MemoryLoopClaim"/>'s
/// four Z₄ sectors. This claim names the answer: they are complementary, not identical.
///
/// <list type="bullet">
///   <item><b>Z₄ from Π itself</b> (cyclic order 4): generator i with i⁴ = 1 and i² = −1.
///         Eigenvalues {+1, +i, −1, −i}, each with multiplicity 16 at N=3 per
///         <c>experiments/ERROR_CORRECTION_PALINDROME.md</c> (March 2026). Captures the
///         Z-axis palindrome dynamics; one operator with a fourth-order memory loop.</item>
///   <item><b>Klein Z₂ × Z₂ from (Π²_Z, Π²_X)</b> (elementary abelian, every non-identity
///         has order 2): two commuting Π²-involutions on Z- and X-axes. Cells
///         {(+,+), (+,−), (−,+), (−,−)} = {Pp, Pm, Mp, Mm}. Captures both Z- and X-
///         dephasing classifications simultaneously; two operators with two-step closure.</item>
/// </list>
///
/// <para>The two are NOT isomorphic as abstract groups: Z₄ has an element of order 4
/// (namely i), while Klein Z₂ × Z₂ has no element of order > 2. They are the only two
/// abelian groups of order 4 (up to isomorphism), and the framework uses both.</para>
///
/// <para><b>Relationship via squaring:</b> the map Z₄ → Z₂ that squares
/// (i ↦ −1, −1 ↦ 1, etc.) collapses the cyclic Z₄ to one Z₂ factor. Concretely
/// Π² = Π²_Z (the Z-axis Π²-involution; one of the two Klein generators). The other
/// Klein generator (Π²_X) is <b>independent</b> of Π; it captures X-axis dephasing,
/// where Π itself does not extend.</para>
///
/// <para>Operationally: Z₄ → Klein cannot be a Tom-style "same partition" identity,
/// but Klein CAN be read as <c>(Π²_Z, Π²_X) = (squaring(Π), independent X-axis Π²)</c>.
/// Together they cover the full Π² operator landscape across two dephasing axes; Z₄
/// alone covers only one axis but with finer (order-4) detail on it.</para>
///
/// <para>Both decompose the 4^N Pauli space into 4 sectors. At N = 3 each sector is
/// 16-dimensional in BOTH decompositions, which is what tempted the conjecture
/// "same partition" — but the partitions are different (only the cardinalities match).
/// The 4-sector cardinality 4^N / 4 = 4^(N−1) is a Lagrange-theorem consequence (group
/// of order 4 acting on 4^N elements with regular orbits), not evidence of isomorphism.</para>
///
/// <para>Tier1Derived: pure group theory. Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c>
/// F88 (Klein four-cell decomposition) +
/// <c>experiments/ERROR_CORRECTION_PALINDROME.md</c> (March 2026: Π⁴ = I, eigenvalues
/// {±1, ±i} multiplicity 16 each at N=3) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs</c> (Z₄ side) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (KleinFourCellClaim, the Z₂ × Z₂ side).</para></summary>
public sealed class Pi2Z4KleinDistinctionClaim : Claim
{
    /// <summary>Z₄ has 4 elements; one element (i) has order 4 (i^4 = 1 first time at k=4).</summary>
    public const int Z4MaxElementOrder = 4;

    /// <summary>Klein Z₂ × Z₂ has 4 elements; every non-identity element has order 2
    /// ((g, h) · (g, h) = (g², h²) = (e, e) regardless of g, h).</summary>
    public const int KleinMaxNonIdentityElementOrder = 2;

    /// <summary>The two abelian groups of order 4 are Z₄ and Z₂ × Z₂ (Klein) — period;
    /// no others up to isomorphism (Lagrange + classification of finite abelian groups).</summary>
    public const int OrderOfBothGroups = 4;

    public Pi2Z4KleinDistinctionClaim()
        : base("Z₄ (Π) and Klein Z₂ × Z₂ (Π²_Z, Π²_X) are the two non-isomorphic order-4 mirror structures of the framework",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F88 + " +
               "experiments/ERROR_CORRECTION_PALINDROME.md (Π⁴ = I, eigenvalues {±1, ±i} mult 16 at N=3) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs")
    { }

    /// <summary>Squares an element of Z₄ (encoded as 0..3 for {1, i, −1, −i}) to obtain
    /// an element of Z₂ (0 or 1, for {1, −1}). The squaring map Z₄ → Z₂ is the
    /// homomorphism that connects the cyclic side to one of Klein's two Z₂ axes
    /// (specifically the Π²_Z axis when Π is the Z-dephasing palindrome).</summary>
    public int SquareInZ4ToZ2(int z4Power)
    {
        int r = ((z4Power % 4) + 4) % 4;
        // i^0=1 → squared=1=Z2[0]; i^1=i → squared=-1=Z2[1]; i^2=-1 → squared=1=Z2[0]; i^3=-i → squared=-1=Z2[1]
        return r % 2;
    }

    /// <summary>True iff the two structures are isomorphic as abstract groups. Always
    /// false: Z₄ has an order-4 element, Klein does not.</summary>
    public bool AreIsomorphic => false;

    /// <summary>Both decompose 4^N into 4 sectors of 4^(N−1) each. At N = 3 both give
    /// 16-dimensional sectors — the cardinality coincidence that tempted the
    /// "same partition" conjecture.</summary>
    public int SectorDimension(int N) => N >= 1 ? (int)Math.Pow(4, N - 1) : throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 1.");

    public override string DisplayName =>
        "Z₄ vs Klein Z₂×Z₂: two non-isomorphic order-4 mirror structures";

    public override string Summary =>
        $"both groups have order {OrderOfBothGroups}, both decompose 4^N Pauli space into 4 sectors of 4^(N−1); but Z₄ is cyclic (max element order {Z4MaxElementOrder}, the i loop), Klein is elementary abelian (max element order {KleinMaxNonIdentityElementOrder}, two independent involutions); they are NOT isomorphic ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Z₄ from Π",
                summary: "cyclic; eigenvalues {+1, +i, −1, −i}; one operator with i^4 = 1 and i^2 = −1; covers Z-axis palindrome with order-4 detail");
            yield return new InspectableNode("Klein from (Π²_Z, Π²_X)",
                summary: "elementary abelian; cells {(+,+), (+,−), (−,+), (−,−)} = {Pp, Pm, Mp, Mm}; two commuting involutions on Z- and X-axes");
            yield return new InspectableNode("isomorphism",
                summary: $"Z₄ ≇ Klein (Z₄ has element of order 4; Klein does not). Order check: Z₄ max-element-order = {Z4MaxElementOrder}, Klein max = {KleinMaxNonIdentityElementOrder}");
            yield return new InspectableNode("squaring bridge Z₄ → Z₂",
                summary: "i ↦ −1, −1 ↦ 1, etc.; the cyclic group folds onto one Z₂ factor; Π² = Π²_Z is exactly this image");
            yield return new InspectableNode("Klein extension",
                summary: "Klein extends Π² with the independent Π²_X; the X-axis Z₂ is NOT in the image of Z₄ squaring");
            yield return new InspectableNode("sector dimensions agree",
                summary: "both give 4^(N−1); at N=3 → 16, N=4 → 64, N=5 → 256; cardinality coincidence ≠ partition coincidence (Lagrange theorem)");
        }
    }
}
