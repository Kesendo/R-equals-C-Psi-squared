using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F108 Part 3 (Tier1Derived): Y-dephasing sibling of F108 Part 1. Every
/// Π²_Y-even Hamiltonian H built from {XX, YY, YZ, ZY, ZZ} bilinears + Y-dephasing
/// admits EXACT operator-level palindrome via the Y-dephasing variant of
/// <see cref="Pi5BilinearOperator"/>:
///
/// <para>  Π_5bilinear (Y-deph) · L · Π_5bilinear⁻¹ = −L − 2σ·I exactly.</para>
///
/// <para>Hence no Π²_Y-even Pauli pair is F87-hard under Y-dephasing. With Part 1 +
/// Part 2 + Part 3 the F108 Π²-even palindrome family covers all three dephase
/// letters; F109's Step 5 is now closed-form across the full {Z, X, Y} set.</para>
///
/// <para><b>Y-deph specifics</b>: Π²_Y eigenvalue equals Π²_Z eigenvalue (both count
/// bit_b per <see cref="PiOperator.SquaredEigenvalue"/>), so the Π²_Y-even bilinear
/// set is identical to the Π²_Z-even set: {XX, YY, YZ, ZY, ZZ}. The Y-deph variant
/// of Π_5bilinear differs from the Z-deph variant only in the Y/Z 2-cycle phase
/// (−i vs +i, matching the canonical Y-deph Π's phase convention). M² sign pattern
/// is identical to Part 1's: diag(−1, −1, +1, +1) on {I, X, Y, Z}.</para>
///
/// <para><b>Bilinear predicates</b>: identical to F108 Part 1's; F108 Part 3 reuses
/// <see cref="F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinear"/>,
/// <see cref="F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinearTerm"/>, and
/// <see cref="F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinearHamiltonian"/>
/// rather than duplicating them. The bilinear set is the same; the Hamiltonians
/// admitting a palindrome via Π_5bilinear differ only in the matching dephase
/// letter used to build L.</para>
///
/// <para>Implements <see cref="IZ2AxisClaim"/> with <see cref="Z2Axis.BitB"/>:
/// Y-deph shares the bit_b parity axis with Z-deph. No BitA twin: Y-dephasing is
/// intrinsically a bit_b-axis dephase letter (its Π² uses bit_b, not bit_a), so
/// <see cref="BitATwinStatus"/> = <see cref="BitATwinClassification.BitBSpecific"/>.
/// Proof: <c>docs/proofs/PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md</c>.</para></summary>
public sealed class F108Part3Pi2YEvenAlwaysPalindromic : Claim, IZ2AxisClaim
{
    /// <summary>BitB axis (Π²_Y = bit_b parity, same axis as Z-dephasing).</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>No BitA twin: Y-dephasing is intrinsically a bit_b-axis dephase.</summary>
    public Claim? BitATwin => null;

    /// <summary>Override returning <see cref="BitATwinClassification.BitBSpecific"/>:
    /// the algebraic content (Y-dephasing) is intrinsically tied to bit_b structure,
    /// no meaningful bit_a-axis analog exists.</summary>
    public BitATwinClassification BitATwinStatus => BitATwinClassification.BitBSpecific;

    /// <summary>Theorem statement in one line.</summary>
    public string Theorem =>
        "For H = Σ_b α_b · B_b with B_b ∈ {XX, YY, YZ, ZY, ZZ} bilinears and α_b ∈ ℝ, " +
        "and Y-dephasing on every site: Π_5bilinear (Y-deph) · L · Π⁻¹ = −L − 2σ·I exactly.";

    /// <summary>F87 corollary scoped to Y-dephasing. Z-deph covered by Part 1,
    /// X-deph by Part 2; together Part 1+2+3 close F108 across all three dephase
    /// letters.</summary>
    public string F87Corollary =>
        "Under Y-dephasing: no Π²_Y-even Pauli pair (truly or non-truly) is F87-hard; every such pair has palindromic spec(L).";

    public F108Part3Pi2YEvenAlwaysPalindromic()
        : base("F108 Part 3: Π²_Y-even H + Y-dephasing always admits exact operator-level palindrome via Π_5bilinear (Y-deph variant); Y-dephasing sibling of F108 Part 1",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F108 + " +
               "docs/proofs/PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md + " +
               "docs/proofs/PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi5BilinearOperator.cs + " +
               "experiments/NON_HEISENBERG_PALINDROME.md")
    {
    }

    public override string DisplayName =>
        "F108 Part 3 Π²_Y-even always palindromic via Π_5bilinear (Y-deph variant, closed-form)";

    public override string Summary =>
        $"{Theorem} {F87Corollary} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("F87 corollary", summary: F87Corollary);
            yield return new InspectableNode("Π_5bilinear per-site map (Y-deph variant)",
                summary: "I → +X, X → −I, Y → −iZ, Z → +iY (Z-deph's I↔X, Y↔Z permutation with the " +
                         "Y/Z 2-cycle phase flipped from +i to −i to match Y-deph's canonical phase); " +
                         "M² = diag(−1, −1, +1, +1) on {I, X, Y, Z} (same as Z-deph variant since both share bit_b parity)");
            yield return new InspectableNode("Π²_Y-even bilinears handled",
                summary: "{XX, YY, YZ, ZY, ZZ} = identical to Π²_Z-even set; F108 Part 3 reuses " +
                         "F108 Part 1's IsPi2EvenBilinear / IsPi2EvenBilinearTerm / IsPi2EvenBilinearHamiltonian predicates");
            yield return new InspectableNode("Empirical verification",
                summary: "9 pure-Π²_Y-even non-truly pairs × N=3,4,5 = 27 instances, residual = 0; " +
                         "random non-uniform-J trials at N=3,4,5, residual = 0; " +
                         "pure D[Y]^⊗N dissipator N=3,4,5, residual = 0");
            yield return new InspectableNode("Closes F109 Y-dephasing branch",
                summary: "F109 (mother sector Klein (0,0) soft y_par=1 purity) had Y-dephasing branch " +
                         "empirically anchored only after F108 Parts 1+2; this Claim closes the Y branch " +
                         "closed-form. F109 is now fully unconditional Tier1Derived across all three dephase letters.");
            yield return new InspectableNode("Siblings",
                summary: "F108 Part 1 (Z-deph, BitB axis, same bilinear set), F108 Part 2 (X-deph, " +
                         "BitA axis, X-deph bilinear set). Part 1 + Part 2 + Part 3 cover the F108 " +
                         "Π²-even palindrome family completely across {Z, X, Y} dephasing.");
        }
    }
}
