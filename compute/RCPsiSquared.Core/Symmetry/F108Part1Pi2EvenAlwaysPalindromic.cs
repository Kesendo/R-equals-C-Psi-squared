using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F108 Part 1 (Tier1Derived): every Π²_Z-even Hamiltonian H built from the
/// five Π²-even 2-site bilinears {XX, YY, YZ, ZY, ZZ} with arbitrary real bond
/// coefficients admits an EXACT operator-level palindrome under Z-dephasing:
///
/// <para>  Π_5bilinear · L · Π_5bilinear⁻¹ = −L − 2σ·I exactly, where σ = Σ_l γ_l.</para>
///
/// <para>Consequence: no Π²-even Pauli pair can be F87-hard. The 5346+ Π²-even pairs
/// observed soft across F103/F105/F106 are now closed-form, not just empirical.</para>
///
/// <para>Derivation chain (see PROOF_F108_PART1):</para>
/// <list type="number">
///   <item>{Q, [B, ·]} = 0 for every Π²-even 2-body bilinear B ∈ {XX, YY, YZ, ZY, ZZ}
///         (verified bit-exact at the 2-qubit level; 4 Π²-odd bilinears produce
///         residual = 8.00, clean separation).</item>
///   <item>Q · L_H · Q⁻¹ = −L_H for any linear combination of Π²-even bilinears.</item>
///   <item>Q · D[Z_l] · Q⁻¹ = −D[Z_l] − 2γ_l·I per site, via diagonal permutation in
///         the Pauli basis: D[Z]_pauli = γ·diag(0, −2, −2, 0) on {I, X, Y, Z}; M's
///         (I↔X, Y↔Z) swap permutes the diagonal to γ·diag(−2, 0, 0, −2) =
///         −D[Z]_pauli − 2γ·I_4 (phase factors cancel pairwise on each 2-cycle).</item>
///   <item>Combining 2+3: Q · L · Q⁻¹ = −L − 2σ·I exactly, σ = Σγ.</item>
///   <item>spec(L) = spec(Q·L·Q⁻¹) = {−λ − 2σ : λ ∈ spec(L)}: palindromic around −σ.
///         Hence no Π²-even pair can be F87-hard.</item>
/// </list>
///
/// <para>Empirical confirmation: bit-exact residual = 0 across all 9 pure-Π²-even
/// non-truly pairs (single-bilinear YZ/ZY plus 7 two-term combinations) at N=3, 4, 5,
/// plus 15 random non-uniform-J instances + 9 asymmetric J_YZ ≠ J_ZY instances. Pure
/// dissipator D[Z]^⊗N also bit-exact at N=3, 4, 5. Reproduction:
/// <c>simulations/f108_part1_pi_family_scan.py</c>.</para>
///
/// <para>IZ2AxisClaim with <see cref="Z2Axis.BitB"/>; BitA twin
/// <see cref="F108Part2Pi2XEvenAlwaysPalindromic"/> wired as ctor parent (Filled).
/// Y-dephasing analog covered by F108 Part 3 (closed 2026-05-25). Proof:
/// <c>docs/proofs/PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md</c>.</para></summary>
public sealed class F108Part1Pi2EvenAlwaysPalindromic : Claim, IZ2AxisClaim
{
    /// <summary>BitB axis (Π²_Z = X⊗N, bit_b parity = #Y + #Z mod 2). F108 Part 1
    /// constrains the Π²-even (bit_b = 0) sub-cell to always admit a palindrome
    /// operator.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>BitA twin: F108 Part 2 (X-dephasing variant).</summary>
    public Claim? BitATwin => Part2;

    /// <summary>Override mirrors the IZ2AxisClaim default (BitB + non-null twin =
    /// Filled); required to expose the status as a property on the concrete class.</summary>
    public BitATwinClassification BitATwinStatus => BitATwinClassification.Filled;

    /// <summary>Typed BitA-twin parent (F108 Part 2). Records the parent edge in
    /// the inheritance graph and surfaces in the inspector.</summary>
    public F108Part2Pi2XEvenAlwaysPalindromic Part2 { get; }

    /// <summary>The five Π²_Z-even 2-site bilinears that Π_5bilinear handles
    /// with exact palindrome residual. {XX, YY, YZ, ZY, ZZ} = the bit_b=0 subset
    /// of all 16 letter pairs; equivalently the (BitA, BitB) Klein cells (0,0)
    /// {XX, YY, ZZ truly} and (1,0) {YZ, ZY non-truly}.</summary>
    public IReadOnlyList<(PauliLetter, PauliLetter)> Pi2EvenBilinears { get; } =
        new (PauliLetter, PauliLetter)[]
        {
            (PauliLetter.X, PauliLetter.X),
            (PauliLetter.Y, PauliLetter.Y),
            (PauliLetter.Y, PauliLetter.Z),
            (PauliLetter.Z, PauliLetter.Y),
            (PauliLetter.Z, PauliLetter.Z),
        };

    /// <summary>True iff (letter1, letter2) belongs to {XX, YY, YZ, ZY, ZZ}. I-containing
    /// pairs are trivially Π²-even but carry no 2-body Hamiltonian and are excluded.</summary>
    public static bool IsPi2EvenBilinear(PauliLetter letter1, PauliLetter letter2)
    {
        if (letter1 == PauliLetter.I || letter2 == PauliLetter.I) return false;
        // bit_b for each letter: I=0, X=0, Y=1, Z=1. bit_b of pair = sum mod 2.
        int bitB = letter1.BitB() + letter2.BitB();
        return (bitB & 1) == 0;
    }

    /// <summary>True iff <paramref name="term"/> is a 2-body bilinear whose two non-I
    /// letters form a Π²-even pair.</summary>
    public static bool IsPi2EvenBilinearTerm(PauliTerm term)
    {
        if (term is null) throw new ArgumentNullException(nameof(term));
        if (term.KBody != 2) return false;
        var nonI = new List<PauliLetter>();
        foreach (var L in term.Letters)
            if (L != PauliLetter.I) nonI.Add(L);
        if (nonI.Count != 2) return false;
        return IsPi2EvenBilinear(nonI[0], nonI[1]);
    }

    /// <summary>True iff every term in <paramref name="terms"/> is a Π²-even 2-body
    /// bilinear. F108 Part 1's exact-palindrome guarantee holds iff this returns true.</summary>
    public static bool IsPi2EvenBilinearHamiltonian(IEnumerable<PauliTerm> terms)
    {
        if (terms is null) throw new ArgumentNullException(nameof(terms));
        foreach (var t in terms)
            if (!IsPi2EvenBilinearTerm(t)) return false;
        return true;
    }

    /// <summary>The theorem statement in one line.</summary>
    public string Theorem =>
        "For H = Σ_b α_b · B_b with B_b ∈ {XX, YY, YZ, ZY, ZZ} bilinears and α_b ∈ ℝ, " +
        "and Z-dephasing on every site: Π_5bilinear · L · Π_5bilinear⁻¹ = −L − 2σ·I exactly.";

    /// <summary>F87 corollary scoped to Z-dephasing; X-deph covered by F108 Part 2,
    /// Y-deph covered by F108 Part 3.</summary>
    public string F87Corollary =>
        "Under Z-dephasing: no Π²_Z-even Pauli pair (truly or non-truly) is F87-hard; every such pair has palindromic spec(L).";

    public F108Part1Pi2EvenAlwaysPalindromic(F108Part2Pi2XEvenAlwaysPalindromic part2)
        : base("F108 Part 1: Π²-even H + Z-dephasing always admits exact operator-level palindrome via Π_5bilinear (base claim of the Klein-V₄-equivalent F108 family; Parts 2, 3 are Klein-V₄ corollaries per PROOF_F108_KLEIN_V4_EQUIVALENCE.md)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F108 + " +
               "docs/proofs/PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md + " +
               "docs/proofs/PROOF_F108_KLEIN_V4_EQUIVALENCE.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi5BilinearOperator.cs + " +
               "experiments/NON_HEISENBERG_PALINDROME.md")
    {
        Part2 = part2 ?? throw new ArgumentNullException(nameof(part2));
    }

    public override string DisplayName =>
        "F108 Part 1 Π²-even always palindromic via Π_5bilinear (closed-form)";

    public override string Summary =>
        $"{Theorem} {F87Corollary} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("F87 corollary", summary: F87Corollary);
            yield return new InspectableNode("Π_5bilinear per-site map",
                summary: "I → +X, X → −I, Y → +iZ, Z → −iY (P1 permutation with X→I and Z→Y phase flips); M² = diag(−1, −1, +1, +1) on {I, X, Y, Z}");
            yield return new InspectableNode("Π²-even bilinears handled",
                summary: "{XX, YY, YZ, ZY, ZZ} = 5 of 16 2-site Pauli pairs; covers Klein (0,0) truly + (1,0) non-truly");
            yield return new InspectableNode("Empirical verification",
                summary: "9 pure-Π²-even non-truly pairs × N=3,4,5 = 27 instances, residual = 0; " +
                         "15 random non-uniform-J + 9 asymmetric J_YZ≠J_ZY instances, residual = 0; " +
                         "pure D[Z]^⊗N dissipator N=3,4,5, residual = 0");
            yield return new InspectableNode("Closes F109 dependency (Z-dephasing branch)",
                summary: "F109 (mother sector Klein (0,0) soft y_par=1 purity) had Y-dephasing branch empirically " +
                         "anchored only; together with F108 Part 1+2+3 closure on 2026-05-25, F109 is now fully " +
                         "unconditional Tier1Derived across {Z, X, Y}.");
            yield return new InspectableNode("BitA twin (Filled)",
                summary: "F108 Part 2 (Π²_X-even under X-dephasing) typed and wired as ctor parent. " +
                         "BitATwinStatus defaults to Filled per IZ2AxisClaim.");
            yield return new InspectableNode("Open siblings",
                summary: "F110 (HardCellYInversionPattern, Tier1Derived since 2026-06-10, typed 2026-05-25): hard cells y_par-asymmetric " +
                         "with Y-inversion; Aspect A closed-form via F108 Part 1+2+3, Aspect B+C empirically anchored at F103/F105/F106 " +
                         "(closed-form 42:8/228:0 derivation Tier1Derived via F103 §6 counting rule + §7 bipartite mechanism, 2026-06-10).");
            yield return new InspectableNode("Sibling on shared bit_b axis (F112)",
                summary: "F112 (LindbladBitBPiBalance, Tier1Derived for Hermitian H, typed 2026-05-26): orthogonal " +
                         "derived theorem on the same F38/F63 bit_b foundation. F108 Part 1 closes spec(L) palindromy " +
                         "for Π²-even bilinear H + Z-deph; F112 closes Π +i/−i Frobenius balance for arbitrary " +
                         "Hermitian H + bit_b-homogeneous dissipator c. Both Tier1Derived projections of the same " +
                         "bit_b Z₂-grading on the Pauli group; F112 takes F108 Part 1 as typed ctor parent.");
            yield return new InspectableNode("F108 family Klein-V₄ equivalence (Welle 14, 2026-05-27)",
                summary: "F108 Parts 2 and 3 are Klein-V₄ corollaries of Part 1 via complementary mechanisms: " +
                         "Part 3 follows by operator-space D-conjugation (D · Π_5b(Z) · D = Π_5b(Y) bit-exact, " +
                         "bilinear set D-invariant); Part 2 follows by Hilbert-space Hadamard transport " +
                         "(U_op = U_H^⊗N ⊗ (U_H^⊗N)^* maps L_Z to L_X bit-exact, bilinear-set bijection). " +
                         "Operator-level Q_zx and H do NOT swap Π_5b(Z) ↔ Π_5b(X) (negative result). " +
                         "Full proof: PROOF_F108_KLEIN_V4_EQUIVALENCE.md. The three typed Claims remain " +
                         "separate to preserve their independent integration edges, but cross-reference this proof.");
        }
    }
}
