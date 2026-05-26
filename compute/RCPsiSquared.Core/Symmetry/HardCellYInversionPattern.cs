using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F110 (Tier1Candidate): F87-hard pairs only appear in the diagonal Klein
/// cell whose Klein index matches the dephase letter's own Klein index, and within
/// that cell the dominant y_par equals y_par(dephase letter). Seventh YParity-axis
/// Claim (after F102/F103/F105/F106/F107/F109).
///
/// <para>Three structural aspects:</para>
/// <list type="bullet">
///   <item><b>Aspect A (closed-form):</b> F87-hard pairs only in the diagonal Klein
///         cell (Z → (0, 1), X → (1, 0), Y → (1, 1)). Proof chain: F108 Part 1+2+3
///         close Π²-D-even cells (never hard); F107 + F109 close Mother sector Klein
///         (0, 0); F87 dissipator-resonance law selects the single diagonal cell from
///         the two remaining Π²-D-odd non-mother cells (per dephase D).</item>
///   <item><b>Aspect B (Y-inversion, empirical):</b> Dominant y_par equals
///         y_par(dephase letter): Z/X-deph dominantly y_par=0; Y-deph dominantly
///         y_par=1. Structural reading: the dephase letter's own Y-content
///         determines which y_par value the F87-hard split favors.</item>
///   <item><b>Aspect C (k-purity sharpening, empirical):</b> At k=3 N=4: 42:8
///         biased (per F103). At k=3 N=5: identical 42:8 (N-stable per F105). At
///         k=4 N=4: 228:0 fully pure (per F106). Closed-form derivation of these
///         ratios is open (F103 Section 5).</item>
/// </list>
///
/// <para>Tier1Candidate because Aspect A is closed-form but Aspect B+C are
/// empirically anchored only. Proof:
/// <c>docs/proofs/PROOF_F110_HARD_CELL_Y_INVERSION.md</c>.</para></summary>
public sealed class HardCellYInversionPattern : Claim, IZ2AxisClaim
{
    public Z2Axis Z2Axis => Z2Axis.YParity;
    public Claim? BitATwin => null;

    // ============================================================
    // Aspect A: closed-form diagonal Klein cell lookup
    // ============================================================

    /// <summary>F87 dissipator-resonance law: the diagonal Klein cell where F87-hard
    /// pairs can appear for the given dephase letter. By construction, the dephase
    /// letter's own Klein index IS the diagonal cell (per <see cref="PauliLetter"/>
    /// bit_a/bit_b convention: Z = (0, 1), X = (1, 0), Y = (1, 1)). Hard never
    /// appears in any other cell, per F108 Part 1+2+3 (Π²-D-even cells closed) +
    /// F107 + F109 (Mother sector Klein (0, 0) closed).</summary>
    public static (int BitA, int BitB) DiagonalKleinCellForDephase(PauliLetter dephase)
    {
        if (dephase == PauliLetter.I)
            throw new ArgumentException(
                $"dephase must be X, Y, or Z; got {dephase}", nameof(dephase));
        return (dephase.BitA(), dephase.BitB());
    }

    /// <summary>True iff the (Klein, dephase) pair is on the F87 dissipator-resonance
    /// diagonal, i.e., hard CAN appear in this cell under this dephase letter.
    /// False for all other cells, where hard CANNOT appear (closed-form via F108
    /// Part 1+2+3 + F107 + F109).</summary>
    public static bool IsDiagonalCell((int BitA, int BitB) klein, PauliLetter dephase) =>
        klein == DiagonalKleinCellForDephase(dephase);

    // ============================================================
    // Aspect B: Y-inversion structural reading (empirical)
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
        "Aspect A (closed-form): F87-hard pairs only in the diagonal Klein cell matching the dephase letter (Z → (0, 1), X → (1, 0), Y → (1, 1)). " +
        "Aspect B (empirical Y-inversion): dominant y_par in hard cell equals y_par(dephase letter); Y-deph inverts to y_par=1. " +
        "Aspect C (empirical k-sharpening): k=3 (42:8) sharpens to k=4 (228:0) with Y-inversion preserved.";

    /// <summary>F87 corollary: the dominant y_par value in any F87-hard cell is
    /// determined by y_par(dephase letter), not by the cell's Klein index.</summary>
    public string F87Corollary =>
        "In every F87-hard diagonal Klein cell, the dominant y_par = y_par(dephase letter). At k = 4 (F106 N=4) the dominance is full (purity); at k = 3 (F103/F105) it is biased (42:8 split). k ≥ 5 empirical confirmation open.";

    /// <summary>Typed Cubic3 parent: <see cref="KleinEightCellClaim"/>. F110
    /// pins the Y-inversion split in F87-hard diagonal cells across the Z₂³
    /// 8-cell decomposition; KleinEightCellClaim is the structural anchor.
    /// Wired 2026-05-26.</summary>
    public KleinEightCellClaim KleinEightParent { get; }

    public HardCellYInversionPattern(KleinEightCellClaim klein8)
        : base("F110 F87-hard pairs only in diagonal Klein cells with Y-inversion (Tier1Candidate: Aspect A closed-form via F108 Part 1+2+3 + F87 dissipator-resonance; Aspect B+C empirically anchored at F103/F105/F106); typed Cubic3 parent = KleinEightCellClaim",
               Tier.Tier1Candidate,
               "docs/ANALYTICAL_FORMULAS.md F110 + " +
               "docs/proofs/PROOF_F110_HARD_CELL_Y_INVERSION.md + " +
               "docs/proofs/PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md + " +
               "docs/proofs/PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md + " +
               "docs/proofs/PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md + " +
               "docs/proofs/PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md + " +
               "docs/proofs/PROOF_F109_MOTHER_SOFT_Y_PARITY_ONE_PURITY.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (KleinEightCellClaim, typed Cubic3 parent)")
    {
        KleinEightParent = klein8 ?? throw new ArgumentNullException(nameof(klein8));
    }

    public override string DisplayName =>
        "F110 hard cells y_par-asymmetric with Y-inversion (Tier1Candidate)";

    public override string Summary =>
        $"Aspect A: F87-hard pairs only in the diagonal Klein cell matching the dephase letter (closed-form). " +
        $"Aspect B (empirical): dominant y_par in hard cell equals y_par(dephase letter); Z/X-deph y_par=0, Y-deph y_par=1. " +
        $"Aspect C (empirical): k=3 N=4 split 42:8 sharpens to k=4 N=4 split 228:0 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("F87 corollary", summary: F87Corollary);
            yield return new InspectableNode("Aspect A: diagonal Klein cell mapping (closed-form)",
                summary: "Z → (0, 1), X → (1, 0), Y → (1, 1). Hard appears only in this cell per F108 Part 1+2+3 + F107 + F109 exclusion chain.");
            yield return new InspectableNode("Aspect B: Y-inversion structural reading",
                summary: "Dominant y_par = y_par(dephase letter): Z/X → 0, Y → 1. The Y-letter's y_par=1 inverts the otherwise-y_par=0-preferred pattern. At k = N = 4 closed-form (Tier1Candidate) via sibling Claim F111 (HardCellPureDTemplate, 2026-05-25): hard pairs in diagonal cell contain at least one pure-D template, and pure-D templates have y_par = y_par(D) by construction. At k = 3 the 42:8 dominance is empirical (see Aspect C).");
            yield return new InspectableNode("Aspect C: empirical k-purity sharpening",
                summary: "k=3 N=4 (F103): 42:8 biased per diagonal cell. k=3 N=5 (F105): identical 42:8 (N-stable). k=4 N=4 (F106): 228:0 fully pure with Y-inversion preserved.");
            yield return new InspectableNode("Sibling YParity-axis claims",
                summary: "F102 (YParityIndependenceAtK3, Tier1Derived), F103 (F87Z2CubedRefinementN4K3, Tier2Empirical), F105 (F87Z2CubedRefinementN5K3, Tier2Empirical), F106 (F87Z2CubedRefinementN4K4, Tier2Empirical), F107 (TrulyYParityZeroPurity, Tier1Derived), F109 (MotherSoftYParityOnePurity, Tier1Derived), F110 (HardCellYInversionPattern, THIS Claim, Tier1Candidate), F111 (HardCellPureDTemplate at k=N=4, Tier1Candidate; sharpens Aspect B). Together the 8 YParity-axis Claims pin the y_par signature of all three F87 trichotomy classes.");
            yield return new InspectableNode("Cross-axis dependencies (BitB and BitA): F108 Parts",
                summary: "F108 Part 1+3 (BitB-axis): close F107/F109/F110 derivation via Π_5bilinear under Z and Y dephasing. F108 Part 2 (BitA-axis, BitA twin of Part 1): closes the X-deph branch via the Z↔X Π² mirror. F108 Parts are NOT YParity-axis siblings (per their Z2Axis declarations); they are the cross-axis closure mechanism that lets Aspect A be closed-form.");
            yield return new InspectableNode("Open work",
                summary: "Closed-form derivation of the exact 42:8 (k=3) ratio open (F103 Section 5; F111 closes the k=4 228:0 case via Pure-D Template Rule). F111 subclaim (d) Mixed+Mixed=soft closed-form blocks promotion to Tier1Derived. k ≥ 5 empirical confirmation unverified.");
            yield return new InspectableNode("Cubic3 anchor parent",
                summary: $"KleinEightCellClaim ({KleinEightParent.Tier.Label()}): the Z₂³ 8-cell decomposition anchor for the y_par axis F110 lives on.");
        }
    }
}
