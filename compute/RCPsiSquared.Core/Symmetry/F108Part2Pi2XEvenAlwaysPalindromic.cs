using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F108 Part 2 (Tier1Derived, BitA twin of F108 Part 1): every Π²_X-even
/// Hamiltonian H built from the five Π²_X-even 2-site bilinears {ZZ, XX, XY, YX, YY}
/// with arbitrary real bond coefficients admits an EXACT operator-level palindrome
/// under X-dephasing:
///
/// <para>  Π_5bilinear · L · Π_5bilinear⁻¹ = −L − 2σ·I exactly (with the X-deph
/// variant of <see cref="Pi5BilinearOperator"/>), where σ = Σ_l γ_l.</para>
///
/// <para>Consequence: no Π²_X-even Pauli pair can be F87-hard under X-dephasing.
/// Closes the X-dephasing branch of F109's Step 5 (mother sector soft = y_par=1
/// purity), leaving only the Y-dephasing branch (F108 Part 3) open.</para>
///
/// <para>Derivation chain (see PROOF_F108_PART2):</para>
/// <list type="number">
///   <item>{Q, [B, ·]} = 0 for every Π²_X-even 2-body bilinear B ∈ {ZZ, XX, XY,
///         YX, YY} (verified bit-exact at the 2-qubit level; 4 Π²_X-odd
///         bilinears produce residual = 8.00, clean separation).</item>
///   <item>Q · L_H · Q⁻¹ = −L_H for any linear combination of Π²_X-even
///         bilinears.</item>
///   <item>Q · D[X_l] · Q⁻¹ = −D[X_l] − 2γ_l·I per site, via diagonal permutation
///         in the Pauli basis: D[X]_pauli = γ·diag(0, 0, −2, −2) on {I, X, Y, Z};
///         M's (I↔Z, X↔Y) swap permutes the diagonal entries to
///         γ·diag(−2, −2, 0, 0) = −D[X]_pauli − 2γ·I_4 (phase factors cancel
///         pairwise on each 2-cycle).</item>
///   <item>Combining 2+3: Q · L · Q⁻¹ = −L − 2σ·I exactly, σ = Σγ.</item>
///   <item>spec(L) = spec(Q·L·Q⁻¹) = {−λ − 2σ : λ ∈ spec(L)}: palindromic around
///         −σ. Hence no Π²_X-even pair is F87-hard under X-dephasing.</item>
/// </list>
///
/// <para>Empirical confirmation: bit-exact residual = 0 across all 9 pure-Π²_X-even
/// non-truly pairs (single-bilinear XY/YX plus 7 two-term combinations) at N=3, 4, 5,
/// plus random non-uniform-J instances. Pure dissipator D[X]^⊗N also bit-exact at
/// N=3, 4, 5. Reproduction:
/// <c>simulations/_f108_part2_x_dephasing_scan.py</c>.</para>
///
/// <para>IZ2AxisClaim with <see cref="Z2Axis.BitA"/>; BitB twin is
/// <see cref="F108Part1Pi2EvenAlwaysPalindromic"/>, which holds the typed BitATwin
/// edge. Proof: <c>docs/proofs/PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md</c>.</para></summary>
public sealed class F108Part2Pi2XEvenAlwaysPalindromic : Claim, IZ2AxisClaim
{
    /// <summary>BitA axis (Π²_X = Z⊗N, bit_a parity = #X + #Y mod 2). F108 Part 2
    /// constrains the Π²_X-even (bit_a = 0) sub-cell under X-dephasing to always
    /// admit a palindrome operator.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitA;

    /// <summary>BitA-axis claims have no BitATwin slot (the twin concept lives on
    /// BitB-axis claims pointing at BitA siblings).</summary>
    public Claim? BitATwin => null;

    /// <summary>Override mirrors the IZ2AxisClaim default (non-BitB axis = NotApplicableForThisAxis);
    /// required to expose the status as a property on the concrete class.</summary>
    public BitATwinClassification BitATwinStatus => BitATwinClassification.NotApplicableForThisAxis;

    /// <summary>The five Π²_X-even 2-site bilinears: {ZZ, XX, XY, YX, YY} = the
    /// bit_a=0 subset of all 16 letter pairs.</summary>
    public IReadOnlyList<(PauliLetter, PauliLetter)> Pi2XEvenBilinears { get; } =
        new (PauliLetter, PauliLetter)[]
        {
            (PauliLetter.Z, PauliLetter.Z),
            (PauliLetter.X, PauliLetter.X),
            (PauliLetter.X, PauliLetter.Y),
            (PauliLetter.Y, PauliLetter.X),
            (PauliLetter.Y, PauliLetter.Y),
        };

    /// <summary>True iff (letter1, letter2) belongs to {ZZ, XX, XY, YX, YY}. I-containing
    /// pairs are trivially Π²_X-even but carry no 2-body Hamiltonian and are excluded.</summary>
    public static bool IsPi2XEvenBilinear(PauliLetter letter1, PauliLetter letter2)
    {
        if (letter1 == PauliLetter.I || letter2 == PauliLetter.I) return false;
        int bitA = letter1.BitA() + letter2.BitA();
        return (bitA & 1) == 0;
    }

    /// <summary>True iff <paramref name="term"/> is a 2-body bilinear whose two non-I
    /// letters form a Π²_X-even pair.</summary>
    public static bool IsPi2XEvenBilinearTerm(PauliTerm term)
    {
        if (term is null) throw new ArgumentNullException(nameof(term));
        if (term.KBody != 2) return false;
        var nonI = new List<PauliLetter>();
        foreach (var L in term.Letters)
            if (L != PauliLetter.I) nonI.Add(L);
        return IsPi2XEvenBilinear(nonI[0], nonI[1]);
    }

    /// <summary>True iff every term in <paramref name="terms"/> is a Π²_X-even 2-body
    /// bilinear. F108 Part 2's exact-palindrome guarantee under X-dephasing holds iff
    /// this returns true.</summary>
    public static bool IsPi2XEvenBilinearHamiltonian(IEnumerable<PauliTerm> terms)
    {
        if (terms is null) throw new ArgumentNullException(nameof(terms));
        foreach (var t in terms)
            if (!IsPi2XEvenBilinearTerm(t)) return false;
        return true;
    }

    /// <summary>The theorem statement in one line.</summary>
    public string Theorem =>
        "For H = Σ_b α_b · B_b with B_b ∈ {ZZ, XX, XY, YX, YY} bilinears and α_b ∈ ℝ, " +
        "and X-dephasing on every site: Π_5bilinear · L · Π_5bilinear⁻¹ = −L − 2σ·I exactly " +
        "(with the X-deph variant of Pi5BilinearOperator).";

    /// <summary>F87 corollary scoped to X-dephasing.</summary>
    public string F87Corollary =>
        "Under X-dephasing: no Π²_X-even Pauli pair (truly or non-truly) is F87-hard; every such pair has palindromic spec(L).";

    public F108Part2Pi2XEvenAlwaysPalindromic()
        : base("F108 Part 2: Π²_X-even H + X-dephasing always admits exact operator-level palindrome via Π_5bilinear (X-deph variant); BitA twin of F108 Part 1; also Klein-V₄ corollary of Part 1 via Hilbert-space Hadamard transport per PROOF_F108_KLEIN_V4_EQUIVALENCE.md (Welle 14)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F108 + " +
               "docs/proofs/PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md + " +
               "docs/proofs/PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md + " +
               "docs/proofs/PROOF_F108_KLEIN_V4_EQUIVALENCE.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi5BilinearOperator.cs + " +
               "experiments/NON_HEISENBERG_PALINDROME.md")
    {
    }

    public override string DisplayName =>
        "F108 Part 2 Π²_X-even always palindromic via Π_5bilinear (X-deph variant, closed-form)";

    public override string Summary =>
        $"{Theorem} {F87Corollary} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("F87 corollary", summary: F87Corollary);
            yield return new InspectableNode("Π_5bilinear per-site map (X-deph variant)",
                summary: "I → +Z, Z → −I, X → −iY, Y → +iX (canonical X-deph permutation I↔Z, X↔Y " +
                         "with Z→I and Y→X back-arrow phase flips); M² = diag(−1, +1, +1, −1) on {I, X, Y, Z}");
            yield return new InspectableNode("Π²_X-even bilinears handled",
                summary: "{ZZ, XX, XY, YX, YY} = 5 of 16 2-site Pauli pairs with bit_a=0");
            yield return new InspectableNode("Empirical verification",
                summary: "9 pure-Π²_X-even non-truly pairs × N=3,4,5 = 27 instances, residual = 0; " +
                         "random non-uniform-J trials at N=3,4,5, residual = 0; " +
                         "pure D[X]^⊗N dissipator N=3,4,5, residual = 0");
            yield return new InspectableNode("Closes F109 X-dephasing branch",
                summary: "F109 X-dephasing branch closed-form via Part 2's X-deph variant of Π_5bilinear. " +
                         "Together with Part 1 (Z-deph) and Part 3 (Y-deph), F109 is now fully unconditional " +
                         "across {Z, X, Y}.");
            yield return new InspectableNode("Sibling on BitB axis",
                summary: "F108 Part 1 (Pi2EvenAlwaysPalindromic, Z-dephasing) and F108 Part 3 (Y-dephasing). Together " +
                         "Parts 1+2+3 cover the F108 Π²-even palindrome family completely across {Z, X, Y} dephasing.");
            yield return new InspectableNode("Open siblings",
                summary: "F110 (HardCellYInversionPattern, Tier1Candidate, typed 2026-05-25): hard cells y_par-asymmetric " +
                         "with Y-inversion; Aspect A closed-form via F108 Part 1+2+3, Aspect B+C empirically anchored " +
                         "(closed-form 42:8/228:0 derivation open).");
            yield return new InspectableNode("Klein-V₄ corollary mechanism (Welle 14, 2026-05-27)",
                summary: "Part 2 follows from Part 1 by Hilbert-space Hadamard transport: U_op = U_H^⊗N ⊗ (U_H^⊗N)^* " +
                         "maps L_Z(H_1) → L_X(U H_1 U^†) bit-exact; per-letter U_H sends Z↔X, Y→−Y, I→I, which " +
                         "bijects the Part-1 bilinear set {XX, YY, YZ, ZY, ZZ} (with sign flips) onto the Part-2 set " +
                         "{ZZ, XX, XY, YX, YY}. NEGATIVE: operator-space Q_zx (Klein-V₄ Z↔X swap from Welle 12) does " +
                         "NOT swap Π_5b(Z) ↔ Π_5b(X) at the operator level (residual 2.0 in Frobenius distance). The " +
                         "Hadamard route gives an equivalent palindrome operator U_op · Π_5b(Z) · U_op^† for L_X, " +
                         "DIFFERENT from canonical Π_5b(X); both achieve F108 palindrome (Π_5b family is not unique). " +
                         "See PROOF_F108_KLEIN_V4_EQUIVALENCE.md.");
        }
    }
}
