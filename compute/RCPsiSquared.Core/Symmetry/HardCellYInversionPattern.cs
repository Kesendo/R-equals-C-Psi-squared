using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F110 (Tier1Derived, promoted 2026-06-10): F87-hard Klein-homogeneous pairs only appear in the diagonal Klein
/// cell whose Klein index matches the dephase letter's own Klein index, and within
/// that cell the dominant y_par equals y_par(dephase letter). Seventh YParity-axis
/// Claim (after F102/F103/F105/F106/F107/F109).
///
/// <para>Three structural aspects:</para>
/// <list type="bullet">
///   <item><b>Aspect A (closed-form):</b> F87-hard Klein-homogeneous pairs only in the
///         diagonal Klein cell (Z → (0, 1), X → (1, 0), Y → (1, 1)), at every body count,
///         coupling and per-site rate (a pair across the two non-Mother unmatched cells can be
///         hard, TWO_TERM_PALINDROME_KLEIN_ROUTING). Each of the other three cells is coloured: a string of Klein
///         letter K ∉ {I, D} commutes with K^⊗N, a Mother-sector string with every letter
///         string, and a letter string anticommutes with every D jump exactly when its
///         letter is neither I nor D. Right multiplication by such a lit string gives
///         R·L·R⁻¹ = −L† − 2σ, and F158's sufficiency step pairs the spectrum. The F87
///         dissipator-resonance law is the census of this at N = 4, k = 3.</item>
///   <item><b>Aspect B (Y-inversion, derived):</b> Dominant y_par equals
///         y_par(dephase letter): Z/X-deph dominantly y_par=0; Y-deph dominantly
///         y_par=1. Structural reading: the dephase letter's own Y-content
///         determines which y_par value the F87-hard split favors.</item>
///   <item><b>Aspect C (purity at full support, closed-form counting rule and F111):</b> At k=3 N=4: 42:8
///         biased (per F103). At k=3 N=5: identical 42:8, and N-stable from N=4 up
///         (per F105; at N=3 the same cells read 34:0, since rule (b) of the counting
///         rule needs the term placed at two windows). At
///         k=4 N=4: 228:0 fully pure (per F106). These ratios are derived by the
///         diagonal-cell counting rule (F103 §6) and the bipartite-chirality
///         mechanism (§7); the windowed k&lt;N hard-direction converse, once the one
///         remaining edge (F103 §7.3), closed 2026-06-10 (WindowedConverseAllGammaClaim,
///         Pascal-Gram positivity, no residual).</item>
/// </list>
///
/// <para>Tier1Derived (promoted 2026-06-10): Aspect A is closed-form, Aspect B+C are derived
/// by the F103 §6 counting rule + §7 mechanism, and the windowed (k&lt;N) hard-direction
/// converse, the one edge that kept this Tier1Candidate, is now the closed all-γ theorem
/// (WindowedConverseAllGammaClaim, no residual). Proof:
/// <c>docs/proofs/PROOF_F110_HARD_CELL_Y_INVERSION.md</c>.</para></summary>
public sealed class HardCellYInversionPattern : Claim, IZ2AxisClaim
{
    public Z2Axis Z2Axis => Z2Axis.YParity;
    public Claim? BitATwin => null;

    // ============================================================
    // Aspect A: closed-form diagonal Klein cell lookup
    // ============================================================

    /// <summary>F87 dissipator-resonance law: the diagonal Klein cell where F87-hard
    /// Klein-homogeneous pairs can appear for the given dephase letter. By construction, the dephase
    /// letter's own Klein index IS the diagonal cell (per <see cref="PauliLetter"/>
    /// bit_a/bit_b convention: Z = (0, 1), X = (1, 0), Y = (1, 1)). No Klein-homogeneous
    /// pair is hard in any other cell: each of them is coloured by a lit letter string
    /// (PROOF_F110 §2, F103 §8).</summary>
    public static (int BitA, int BitB) DiagonalKleinCellForDephase(PauliLetter dephase)
    {
        if (dephase == PauliLetter.I)
            throw new ArgumentException(
                $"dephase must be X, Y, or Z; got {dephase}", nameof(dephase));
        return (dephase.BitA(), dephase.BitB());
    }

    /// <summary>True iff the (Klein, dephase) pair is on the F87 dissipator-resonance
    /// diagonal, i.e., hard CAN appear in this cell under this dephase letter.
    /// False for all other cells, where hard CANNOT appear (each is coloured,
    /// PROOF_F110 §2).</summary>
    public static bool IsDiagonalCell((int BitA, int BitB) klein, PauliLetter dephase) =>
        klein == DiagonalKleinCellForDephase(dephase);

    // ============================================================
    // Aspect B: Y-inversion structural reading (derived via F103 §6/§7; corollary of F111 at k=N=4)
    // ============================================================

    /// <summary>The dominant y_par value in the F87-hard diagonal cell for the given
    /// dephase letter. Equals y_par(dephase letter) = (#Y in dephase letter) mod 2.
    /// Only Y has #Y = 1, so this is the AND of the letter's bit_a and bit_b (Y is
    /// the only Pauli letter with both bits set). At k = 3 the dominance is biased
    /// (42:8 split per F103/F105); at k = 4 it is fully pure (228:0 per F106).
    /// Y-deph inverts the otherwise-y_par=0-preferred pattern.</summary>
    public static int DominantYParityForDephase(PauliLetter dephase)
    {
        if (dephase == PauliLetter.I)
            throw new ArgumentException(
                $"dephase must be X, Y, or Z; got {dephase}", nameof(dephase));
        return dephase.BitA() & dephase.BitB();
    }

    /// <summary>The F110 theorem statement in one line, covering Aspect A + B + C.</summary>
    public string Theorem =>
        "Aspect A (closed-form): F87-hard Klein-homogeneous pairs only in the diagonal Klein cell matching the dephase letter (Z → (0, 1), X → (1, 0), Y → (1, 1)). " +
        "Aspect B (Y-inversion, derived F103 §6): dominant y_par in hard cell equals y_par(dephase letter); Y-deph inverts to y_par=1. " +
        "Aspect C (purity at full support, derived F103 §6/§7 and F111): pure at N = k for every k (34:0 at k=3, 228:0 at k=4); below it biased at k = 3, 4, 5 (42:8 at k=3, 356:320 at k=4) and pure at k <= 2; Y-inversion preserved.";

    /// <summary>F87 corollary: the dominant y_par value in any F87-hard cell is
    /// determined by y_par(dephase letter), not by the cell's Klein index.</summary>
    public string F87Corollary =>
        "In every F87-hard diagonal Klein cell, the dominant y_par = y_par(dephase letter). At full support (N = k; k = 4 is F106's 228:0, k = 5 reads 2056:0) the dominance is full (purity); below it the dominance is biased (k = 3: 42:8; k = 4: 356:320 at N = 6, 7, 8; k = 5: 6536:4480 at N = 8; Z dephasing; k ≥ 6 below full support not read), read exactly by Theorem 2 of the complement connection (DiagonalCellComplementConnectionTests).";

    /// <summary>Typed Cubic3 parent: <see cref="KleinEightCellClaim"/>. F110
    /// pins the Y-inversion split in F87-hard diagonal cells across the Z₂³
    /// 8-cell decomposition; KleinEightCellClaim is the structural anchor.
    /// Wired 2026-05-26.</summary>
    public KleinEightCellClaim KleinEightParent { get; }

    public HardCellYInversionPattern(KleinEightCellClaim klein8)
        : base("F110 F87-hard Klein-homogeneous pairs only in diagonal Klein cells with Y-inversion (Tier1Derived: Aspect A by the colouring of the three non-diagonal cells, F158; Aspect B+C derived via F103 §6/§7 + the windowed converse closed 2026-06-10); typed Cubic3 parent = KleinEightCellClaim",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F110 + " +
               "docs/proofs/PROOF_F110_HARD_CELL_Y_INVERSION.md + " +
               "docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md + " +
               "docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md + " +
               "experiments/THE_PALINDROME_AS_A_COLOURING.md + " +
               "docs/proofs/PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md + " +
               "docs/proofs/PROOF_F109_MOTHER_SOFT_Y_PARITY_ONE_PURITY.md + " +
               "simulations/f107_f110_route_gate.py (the exact gate of the route) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (KleinEightCellClaim, typed Cubic3 parent)")
    {
        KleinEightParent = klein8 ?? throw new ArgumentNullException(nameof(klein8));
    }

    public override string DisplayName =>
        "F110 hard cells y_par-asymmetric with Y-inversion (Tier1Derived)";

    public override string Summary =>
        $"Aspect A: F87-hard Klein-homogeneous pairs only in the diagonal Klein cell matching the dephase letter (closed-form). " +
        $"Aspect B (derived F103 §6): dominant y_par in hard cell equals y_par(dephase letter); Z/X-deph y_par=0, Y-deph y_par=1. " +
        $"Aspect C (derived F103 §6/§7 and F111): pure at full support (k=4 N=4: 228:0), biased below it (k=3 N=4: 42:8) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("F87 corollary", summary: F87Corollary);
            yield return new InspectableNode("Aspect A: diagonal Klein cell mapping (closed-form)",
                summary: "Z → (0, 1), X → (1, 0), Y → (1, 1). Hard Klein-homogeneous pairs appear only in this cell: the other three are coloured by a lit letter string (F103 §8, F158), at every body count.");
            yield return new InspectableNode("Aspect B: Y-inversion structural reading",
                summary: "Dominant y_par = y_par(dephase letter): Z/X → 0, Y → 1. The Y-letter's y_par=1 inverts the otherwise-y_par=0-preferred pattern. At k = N = 4 closed-form via sibling Claim F111 (HardCellPureDTemplate, 2026-05-25, Tier1Derived since 2026-06-10): hard pairs in diagonal cell contain at least one pure-D template, and pure-D templates have y_par = y_par(D) by construction. At k = 3 the 42:8 dominance follows from the F103 §6 counting rule (see Aspect C).");
            yield return new InspectableNode("Aspect C: purity at full support (§6 closed-form counting rule and F111)",
                summary: "k=3 N=4 (F103): 42:8 biased per diagonal cell. k=3 N=5 (F105): identical 42:8, N-stable from N=4 up (at N=3 the cells read 34:0 and 21:21, the split without rule (b), which needs the term at two windows). k=4 N=4 (F106): 228:0 fully pure with Y-inversion preserved, and biased below full support (356:320 at k=4, N = 6, 7, 8; DiagonalCellComplementConnectionTests).");
            yield return new InspectableNode("Sibling YParity-axis claims",
                summary: "F102 (YParityIndependenceAtK3, Tier1Derived), F103 (F87Z2CubedRefinementN4K3, Tier1Derived), F105 (F87Z2CubedRefinementN5K3, Tier1Derived), F106 (F87Z2CubedRefinementN4K4, Tier1Derived), F107 (TrulyYParityZeroPurity, Tier1Derived), F109 (MotherSoftYParityOnePurity, Tier1Derived), F110 (HardCellYInversionPattern, THIS Claim, Tier1Derived since 2026-06-10), F111 (HardCellPureDTemplate at k=N=4, Tier1Derived since 2026-06-10; sharpens Aspect B). Together the 8 YParity-axis Claims pin the y_par signature of all three F87 trichotomy classes.");
            yield return new InspectableNode("Cross-axis neighbours (BitB and BitA): F108 Parts",
                summary: "F108 Part 1+3 (BitB-axis) and Part 2 (BitA-axis, BitA twin of Part 1) palindromize the Π²-D-even bilinears via Π_5bilinear; their mirrors reach the Π²-D-even strings of even weight. F108 Parts are NOT YParity-axis siblings (per their Z2Axis declarations). Aspect A rests on the colouring, which reaches every body count.");
            yield return new InspectableNode("Promotion record (2026-06-10) + open work",
                summary: "The exact 42:8 (k=3) ratio is derived by the F103 §6 counting rule and the §7 bipartite-chirality mechanism; F111 closes the k=4 228:0 case via the Pure-D Template Rule. The promotion gate, WindowedConverseAllGammaClaim (the all-γ closure of the windowed k<N converse non-bipartite ⟹ hard), closed 2026-06-10 with no residual (girth dichotomy retired R-deg, Pascal-Gram positivity resolved R-sign; PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md), so F110 is Tier1Derived. k = 5 is read exactly by Theorem 2 of the complement connection (DiagonalCellComplementConnectionTests): hard 2056:0 by y_par at N = 5, the dephase letter's y_par dominant at N = 6, 7, 8.");
            yield return new InspectableNode("Cubic3 anchor parent",
                summary: $"KleinEightCellClaim ({KleinEightParent.Tier.Label()}): the Z₂³ 8-cell decomposition anchor for the y_par axis F110 lives on.");
        }
    }
}
