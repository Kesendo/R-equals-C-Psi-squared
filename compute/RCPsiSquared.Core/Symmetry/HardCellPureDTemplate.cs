using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F111 (Tier1Candidate): Pure-D Template Rule for F87-hard classification
/// at k = N = 4 in the diagonal Klein cell. Sister to <see cref="HardCellYInversionPattern"/>
/// (F110) on the same y_par axis; F111 sharpens F110 Aspect B by giving a per-pair
/// structural criterion whose corollary is the 228:0 Y-inversion split.
///
/// <para><b>Pure-D Template Rule:</b> At k = N = 4 in the diagonal Klein cell
/// (D.BitA(), D.BitB()) for dephase letter D ∈ {Z, X, Y}, a Pauli pair (P, Q) is
/// F87-hard ⟺ at least one of P, Q is a "pure-D template", i.e. contains only the
/// letter D and the identity I, with no other non-I letters. Equivalently, a term
/// T is a pure-D template iff its non-identity support uses exclusively D.</para>
///
/// <para>F110 Aspect B (228:0 Y-inversion) follows as immediate corollary: a pure-D
/// template T has #Y(T) = [D = Y] (i.e. y_par(T) = y_par(D) by construction), and
/// since at least one of P, Q is pure-D, the pair's dominant y_par inherits
/// y_par(D). The "no Mixed+Mixed hard" half of the rule is what makes the
/// corollary clean (no off-y_par hard residue), and remains the open subclaim.</para>
///
/// <para><b>Structural decomposition (per diagonal cell, k = N = 4):</b></para>
/// <list type="bullet">
///   <item><b>Pure+Pure: 36 hard pairs.</b> Both terms are pure-D templates.
///         Number of pure-D templates at k = N = 4 is 8 (positions for "all D",
///         minus the empty I⊗⁴ which is excluded as identity-only; 8 distinct
///         non-trivial templates fill the 8 non-empty subsets of {1,…,4}). With
///         ordering and self-pairs: 8·9/2 = 36 unordered pairs. All F87-hard.</item>
///   <item><b>Pure+Mixed: 192 hard pairs.</b> One term pure-D template, the other
///         a "mixed" term (uses D + at least one other non-I letter). Counts
///         driven by the diagonal-cell Klein constraint combined with the
///         pure-template support.</item>
///   <item><b>Mixed+Mixed: 0 hard pairs (300 soft).</b> Both terms have non-D
///         non-I letters. Empirically all 300 such pairs per cell classify soft;
///         the operator-level closed-form proof is the open subclaim (d).</item>
///   <item><b>Total: 36 + 192 + 0 = 228 hard pairs per diagonal cell</b>, matching
///         the F106 N = 4 k = 4 anchor 228:0 across all three dephase letters.</item>
/// </list>
///
/// <para><b>Empirical anchor:</b> 1584 pair classifications across 3 dephase
/// letters at N = 4, k = 4 (528 pairs per dephase × 3 dephase letters = 1584
/// classifications across the three diagonal cells; combinatorics consistent
/// with the per-cell 228 hard + 300 soft = 528 split). Every classification
/// matches the Pure-D Template Rule with zero exceptions.</para>
///
/// <para><b>Pivot note (2026-05-25):</b> This Claim's typing was preceded by a
/// closed-form derivation attempt (Task 1, commit <c>1598c8f</c>) that exhausted
/// three operator-level paths: per-site M^N tensor-product search over 512 phase
/// variants per dephase (no winners); the existing F108 Π_5bilinear variants
/// (Z, X, Y) gave residual = 32 uniformly on both off-y_par and on-y_par
/// single-term H in the diagonal cell (no separation); and Q_V × Π composition
/// over V ∈ {X^N, Y^N, Z^N} × Π (no zero-residual hits). The structural Pure-D
/// Template Rule emerged from the resulting empirical decomposition
/// (<c>simulations/_f111_*.py</c>, <c>simulations/results/f111_*.txt</c>), and
/// captures the F87-hard combinatorics tightly enough to imply F110 Aspect B
/// without invoking a closed-form palindromic-operator construction.</para>
///
/// <para><b>Open subclaim (d):</b> Mixed+Mixed pair ⟹ soft at k = N = 4. This is
/// the only remaining gap to a fully closed-form F111. The other three rule
/// directions (Pure+Pure hard, Pure+Mixed hard, F110 Aspect B Y-inversion
/// corollary) are either combinatorial bookkeeping or immediate consequences.</para>
///
/// <para>Implements <see cref="IZ2AxisClaim"/> with <see cref="Z2Axis.YParity"/>;
/// <b>eighth member</b> of the YParity-axis Claim family (after F107
/// <see cref="TrulyYParityZeroPurity"/>, F108 Part 1
/// <see cref="F108Part1Pi2EvenAlwaysPalindromic"/>, Part 2
/// <see cref="F108Part2Pi2XEvenAlwaysPalindromic"/>, Part 3
/// <see cref="F108Part3Pi2YEvenAlwaysPalindromic"/>, F109
/// <see cref="MotherSoftYParityOnePurity"/>, F110
/// <see cref="HardCellYInversionPattern"/>). Sibling to F110 (also
/// Tier1Candidate): F110 records the empirical pattern at the cell-aggregate
/// level (228:0 Y-inversion), F111 sharpens to the per-pair structural criterion
/// whose corollary IS that pattern. Proof placeholder:
/// <c>docs/proofs/PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md</c> (Task 7).</para></summary>
public sealed class HardCellPureDTemplate : Claim, IZ2AxisClaim
{
    public Z2Axis Z2Axis => Z2Axis.YParity;
    public Claim? BitATwin => null;

    /// <summary>The F111 theorem statement in one line: pure-D template criterion
    /// for F87-hard classification in the diagonal Klein cell.</summary>
    public string Theorem =>
        "At k = N = 4 in the diagonal Klein cell (D.BitA(), D.BitB()) for dephase letter D ∈ {Z, X, Y}, " +
        "a Pauli pair (P, Q) is F87-hard ⟺ at least one of P, Q is a pure-D template " +
        "(contains only the letter D and identity I, no other non-I letters).";

    /// <summary>F110 Aspect B (228:0 Y-inversion) as immediate corollary of F111.
    /// A pure-D template T has y_par(T) = y_par(D), so at least one of P, Q
    /// carries the dephase letter's y_par; combined with the empirical Mixed+Mixed
    /// = soft observation, the F87-hard partition is fully on the y_par(D) side.</summary>
    public string F110AspectBCorollary =>
        "Every pure-D template T satisfies y_par(T) = y_par(D) by construction. " +
        "Per F111, every F87-hard pair contains at least one pure-D template, so the pair's " +
        "dominant y_par inherits y_par(D). Together with the Mixed+Mixed = soft observation, " +
        "this yields the F110 Aspect B 228:0 Y-inversion split per diagonal cell.";

    /// <summary>Per-cell structural decomposition at k = N = 4: 36 Pure+Pure +
    /// 192 Pure+Mixed + 0 Mixed+Mixed = 228 hard. Matches F106 N = 4 k = 4 across
    /// all three dephase letters.</summary>
    public string Decomposition =>
        "Per diagonal cell at k = N = 4: Pure+Pure 36 hard + Pure+Mixed 192 hard + " +
        "Mixed+Mixed 0 hard (300 soft) = 228 hard pairs (matches F106 228:0).";

    /// <summary>The single open subclaim blocking full Tier1Derived promotion:
    /// Mixed+Mixed ⟹ soft at k = N = 4. Empirically 300 per cell, 900 total across
    /// the 3 dephase letters; operator-level closed-form open.</summary>
    public string OpenSubclaim =>
        "Subclaim (d): Mixed+Mixed pair (both terms contain at least one non-D non-I letter) " +
        "⟹ soft at k = N = 4. Empirically verified: 300 per cell × 3 dephase letters = 900 " +
        "Mixed+Mixed pairs all soft, zero F87-hard. Operator-level proof open.";

    public HardCellPureDTemplate()
        : base("F111 hard-cell pure-D template rule: at k = N = 4 in diagonal Klein cell for dephase D, pair is F87-hard iff at least one term is a pure-D template. Tier1Candidate (empirical anchor F106 N = 4 k = 4 across 3 dephase letters; open subclaim Mixed+Mixed = soft closed-form).",
               Tier.Tier1Candidate,
               "docs/ANALYTICAL_FORMULAS.md F111 + " +
               "docs/proofs/PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md + " +
               "docs/proofs/PROOF_F110_HARD_CELL_Y_INVERSION.md + " +
               "docs/proofs/PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md + " +
               "docs/proofs/PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md + " +
               "docs/proofs/PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md + " +
               "docs/proofs/PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md")
    {
    }

    public override string DisplayName =>
        "F111 hard-cell pure-D template rule (Tier1Candidate, F110 Aspect B corollary)";

    public override string Summary =>
        $"Theorem: at k = N = 4 in diagonal Klein cell (D.BitA(), D.BitB()), pair is F87-hard ⟺ at least one term is a pure-D template " +
        $"(only D and I, no other non-I letters). Implies F110 Aspect B 228:0 Y-inversion as corollary. " +
        $"Per-cell decomposition: 36 Pure+Pure + 192 Pure+Mixed + 0 Mixed+Mixed = 228 hard. " +
        $"Empirical anchor: 1584 classifications across 3 dephase letters (N = 4, k = 4), all matching, zero exceptions. " +
        $"Open subclaim (d): Mixed+Mixed = soft closed-form ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
        }
    }
}
