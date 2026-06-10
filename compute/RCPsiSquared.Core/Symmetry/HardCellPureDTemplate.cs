using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F111 (Tier1Derived, promoted 2026-06-10): Pure-D Template Rule for F87-hard classification
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
/// corollary clean (no off-y_par hard residue); subclaim (d) closed modulo M via
/// PROOF_F103 §7.4 (see the status paragraph below).</para>
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
///         subclaim (d) is closed modulo M via the chiral-K route (PROOF_F103 §7.4).</item>
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
/// <para><b>Subclaim (d) is CLOSED modulo M</b> (PROOF_F103 §7.4, 2026-05-30):
/// at full support k = N a Mixed+Mixed pair has at most two flip generators,
/// which always admit a linear φ and hence the chiral K, so the hopping graph
/// is bipartite and the pair soft; the earlier operator-level search (3
/// derivation paths, PROOF Section 6) was dissolved by the chiral route. The
/// Tier1Derived promotion gate, the hard-direction converse behind subclaims
/// (a)/(c), CLOSED 2026-06-10: WindowedConverseAllGammaClaim is the all-γ
/// theorem with no residual (girth dichotomy + Pascal-Gram positivity), so the
/// non-bipartite hopping graphs behind (a)/(c) are hard at every γ by theorem
/// rather than by the earlier heuristic dissipator-commute reading; the
/// Pure-Pure pair extension and the F110 Aspect B Y-inversion corollary follow
/// as immediate consequences. The PROOF document is authoritative for the
/// status breakdown.</para>
///
/// <para>Implements <see cref="IZ2AxisClaim"/> with <see cref="Z2Axis.YParity"/>;
/// <b>eighth member</b> of the YParity-axis Claim family (after F107
/// <see cref="TrulyYParityZeroPurity"/>, F108 Part 1
/// <see cref="F108Part1Pi2EvenAlwaysPalindromic"/>, Part 2
/// <see cref="F108Part2Pi2XEvenAlwaysPalindromic"/>, Part 3
/// <see cref="F108Part3Pi2YEvenAlwaysPalindromic"/>, F109
/// <see cref="MotherSoftYParityOnePurity"/>, F110
/// <see cref="HardCellYInversionPattern"/>). Sibling to F110 (both
/// Tier1Derived since 2026-06-10): F110 records the empirical pattern at the cell-aggregate
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
    public string DecompositionPerCell =>
        "Per diagonal cell at k=N=4 (count of pairs by template-membership): " +
        "Pure-Pure pairs: 8*9/2 = 36 (all HARD); " +
        "Pure-Mixed pairs: 8*24 = 192 (all HARD); " +
        "Mixed-Mixed pairs: 24*25/2 = 300 (all SOFT). " +
        "Total hard: 36 + 192 + 0 = 228 (matches F106 anchor exactly, across all 3 dephase letters Z, X, Y).";

    /// <summary>Subclaim (a), originally empirically verified with a heuristic
    /// mechanism (the hard-direction converse closed 2026-06-10,
    /// WindowedConverseAllGammaClaim): pure-D single-term Hamiltonians at k=N=4
    /// in the diagonal cell are F87-HARD.
    /// Mechanism: D[D] commutes with pure-D H (the dissipator letter D commutes
    /// with itself), so L = -i[H, .] + L_D has additive independent spectra;
    /// the combined spectrum is non-palindromic around -σ.</summary>
    public string SubclaimA_PureDSingleTermHard =>
        "Pure-D single-term H at k=N=4 in diagonal cell is F87-HARD. Mechanism: D[D] commutes with pure-D H, so L = L_H + L_D has additive independent spectra; combined spectrum non-palindromic.";

    /// <summary>Subclaim (b), empirically verified; mechanism closed via the
    /// chiral route (PROOF_F103 §7.4): mixed single-term Hamiltonians at k=N=4
    /// in the diagonal cell (contain non-D non-I letters) are F87-SOFT
    /// (palindromic spec, M ≠ 0).</summary>
    public string SubclaimB_MixedSingleTermSoft =>
        "Mixed single-term H at k=N=4 in diagonal cell is F87-SOFT (palindromic spec, M ≠ 0). Mechanism closed via the chiral route: a single Mixed term at full support always admits the chiral K (bipartite hopping graph), PROOF_F103 §7.4.";

    /// <summary>Subclaim (c), originally empirically verified with no closed-form
    /// mechanism (the hard-direction converse closed 2026-06-10,
    /// WindowedConverseAllGammaClaim): pair (Pure-D, Mixed) Hamiltonians at
    /// k=N=4 are F87-HARD.</summary>
    public string SubclaimC_PureMixedPairHard =>
        "Pair (Pure-D, Mixed) H at k=N=4 in diagonal cell is F87-HARD. Closed-form mechanism was open until the hard-direction converse closed 2026-06-10 (WindowedConverseAllGammaClaim).";

    /// <summary>Subclaim (d), CLOSED modulo M (PROOF_F103 §7.4, 2026-05-30):
    /// pair (Mixed, Mixed) Hamiltonians at k=N=4 are F87-SOFT (palindromic spec).
    /// At full support a Mixed+Mixed pair has at most two flip generators, which
    /// always admit the chiral K, hence bipartite hence soft. The F111 promotion
    /// gate is the hard-direction converse behind subclaims (a)/(c), not (d).</summary>
    public string SubclaimD_MixedMixedPairSoft =>
        "Pair (Mixed, Mixed) H at k=N=4 in diagonal cell is F87-SOFT. CLOSED modulo M via PROOF_F103 §7.4 (2026-05-30): at most two flip generators at full support, always admitting the chiral K (bipartite ⟹ soft); the earlier operator-level search (Task 1 Paths 1-3) is dissolved. The F111 promotion gate was the hard-direction converse behind subclaims (a)/(c); it closed 2026-06-10 (WindowedConverseAllGammaClaim, no residual).";

    /// <summary>F87 Y-inversion corollary (parallel to F110's F87Corollary): every
    /// pure-D template T has y_par(T) = y_par(D) by construction (templates contain
    /// only D and I; #Y(template) = #Y(D) since I has #Y=0; only Y has #Y=1 mod 2
    /// of itself). Combined with the F111 rule, every F87-hard pair at k=N=4 in the
    /// diagonal cell satisfies y_par(pair) = y_par(D); the F106 N=4 k=4 228:0 split
    /// across {Z, X, Y} is a Pure-D Template Rule corollary.</summary>
    public string YInversionCorollary =>
        "Pure-D templates have y_par = y_par(D) by construction (templates contain only D and I; #Y(template) = #Y(D) since I has #Y=0; only Y has #Y=1 mod 2 of itself). Therefore at k=N=4 in diagonal cell, every F87-hard pair has y_par(pair) = y_par(D). Equivalently: F106 N=4 k=4 228:0 split across {Z, X, Y} is a Pure-D Template Rule corollary.";

    // ============================================================
    // Pure-D Template Rule helpers (static)
    // ============================================================

    /// <summary>Returns true iff <paramref name="term"/> is a pure-<paramref name="dephase"/>
    /// template at the non-trivial level: contains at least one <paramref name="dephase"/>
    /// letter and no other non-I letter. The 8 non-trivial pure-D templates at length 4
    /// correspond to the 8 non-empty subsets of positions {1,2,3,4} carrying D (everything
    /// else I). The all-I term is excluded as identity-only (mirrors the F111 decomposition
    /// docstring: "8 non-trivial templates").</summary>
    public static bool IsPureDTemplate(PauliTerm term, PauliLetter dephase)
    {
        if (term is null) throw new ArgumentNullException(nameof(term));
        if (dephase == PauliLetter.I)
            throw new ArgumentException(
                $"dephase must be X, Y, or Z; got {dephase}", nameof(dephase));
        bool sawD = false;
        foreach (var letter in term.Letters)
        {
            if (letter == PauliLetter.I) continue;
            if (letter != dephase) return false;
            sawD = true;
        }
        return sawD;
    }

    /// <summary>Returns true iff the pair (<paramref name="p"/>, <paramref name="q"/>)
    /// lies in the F111 scope: k = N = 4 strings (both terms have length 4) in the diagonal
    /// Klein cell (D.BitA(), D.BitB()) for the given <paramref name="dephase"/> letter.
    ///
    /// <para>Note: "K4N4" in this method name means string length N = 4 only. The k_body
    /// (popcount of non-I letters) is not constrained by this scope check; the F111 rule's
    /// scope is k_body ≤ 4 within the diagonal cell at N = 4 (the diagonal-cell Klein
    /// constraint already excludes most low-k_body terms; e.g., all-I has Klein (0, 0) and
    /// is never in the diagonal cell for any D ∈ {X, Y, Z}).</para></summary>
    public static bool IsInDiagonalCellAtK4N4(PauliTerm p, PauliTerm q, PauliLetter dephase)
    {
        if (p is null) throw new ArgumentNullException(nameof(p));
        if (q is null) throw new ArgumentNullException(nameof(q));
        if (dephase == PauliLetter.I)
            throw new ArgumentException(
                $"dephase must be X, Y, or Z; got {dephase}", nameof(dephase));
        if (p.N != 4 || q.N != 4) return false;
        var expectedKlein = (dephase.BitA(), dephase.BitB());
        return p.KleinIndex == expectedKlein && q.KleinIndex == expectedKlein;
    }

    /// <summary>Pure-D Template Rule prediction: returns the F111-predicted F87-hardness
    /// for the pair (<paramref name="p"/>, <paramref name="q"/>) under <paramref name="dephase"/>.
    /// Returns false if the pair is outside F111 scope (not in the diagonal cell at k = N = 4).
    /// Returns true iff at least one of P, Q is a pure-D template.</summary>
    public static bool IsPredictedHardAtK4N4(PauliTerm p, PauliTerm q, PauliLetter dephase)
    {
        if (p is null) throw new ArgumentNullException(nameof(p));
        if (q is null) throw new ArgumentNullException(nameof(q));
        if (dephase == PauliLetter.I)
            throw new ArgumentException(
                $"dephase must be X, Y, or Z; got {dephase}", nameof(dephase));
        if (!IsInDiagonalCellAtK4N4(p, q, dephase)) return false;
        return IsPureDTemplate(p, dephase) || IsPureDTemplate(q, dephase);
    }

    /// <summary>F110 Aspect B corollary check: if the Pure-D Template Rule predicts
    /// the pair (<paramref name="p"/>, <paramref name="q"/>) to be F87-hard under
    /// <paramref name="dephase"/>, then both terms must carry y_par = y_par(D).
    ///
    /// <para>Semantics:</para>
    /// <list type="bullet">
    ///   <item>Throws <see cref="ArgumentException"/> when the pair is y_par-inhomogeneous
    ///         (<c>p.YParity != q.YParity</c>). The F111 scope is the Z₂³-homogeneous
    ///         partition (Klein + y_par both shared); y_par-inhomogeneous pairs are out of
    ///         scope and silently returning true would mask real corollary failures.
    ///         <c>Z2HomogeneousKBodyEnumeration</c> (RCPsiSquared.Diagnostics.Tests.F87)
    ///         filters non-homogeneous pairs upstream, so this throw never triggers on
    ///         enumeration-driven call sites.</item>
    ///   <item>Returns true (vacuous-pass) when the rule does NOT predict hard for this pair
    ///         (the corollary speaks only about predicted-hard pairs).</item>
    ///   <item>Returns true (pass) when the rule predicts hard and y_par(p) = y_par(q) =
    ///         y_par(D).</item>
    ///   <item>Returns false (fail) only when the rule predicts hard, the pair is
    ///         y_par-homogeneous, and the shared y_par differs from y_par(D).</item>
    /// </list></summary>
    public static bool VerifyYInversionCorollaryAtK4N4(PauliTerm p, PauliTerm q, PauliLetter dephase)
    {
        if (p is null) throw new ArgumentNullException(nameof(p));
        if (q is null) throw new ArgumentNullException(nameof(q));
        if (dephase == PauliLetter.I)
            throw new ArgumentException(
                $"dephase must be X, Y, or Z; got {dephase}", nameof(dephase));
        if (p.YParity != q.YParity)
            throw new ArgumentException(
                $"Pair must be y_par-homogeneous (got p.YParity={p.YParity}, q.YParity={q.YParity}); " +
                "F111 scope is the Z₂³-homogeneous partition. Filter via Z2HomogeneousKBodyEnumeration upstream.",
                nameof(p));
        if (!IsPredictedHardAtK4N4(p, q, dephase)) return true;
        var expectedYpar = dephase == PauliLetter.Y ? 1 : 0;
        return p.YParity == expectedYpar;
    }

    /// <summary>Typed Cubic3 parent: <see cref="KleinEightCellClaim"/>. F111
    /// sharpens F110 Aspect B via the per-pair Pure-D Template Rule across the
    /// Z₂³ 8-cell decomposition; KleinEightCellClaim is the structural anchor.
    /// Wired 2026-05-26.</summary>
    public KleinEightCellClaim KleinEightParent { get; }

    public HardCellPureDTemplate(KleinEightCellClaim klein8)
        : base("F111 hard-cell pure-D template rule: at k = N = 4 in diagonal Klein cell for dephase D, pair is F87-hard iff at least one term is a pure-D template. Tier1Derived (promoted 2026-06-10: the hard-direction converse gate closed via WindowedConverseAllGammaClaim, Pascal-Gram positivity; Mixed+Mixed = soft closed modulo M via PROOF_F103 §7.4); typed Cubic3 parent = KleinEightCellClaim.",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F111 + " +
               "docs/proofs/PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md + " +
               "docs/proofs/PROOF_F110_HARD_CELL_Y_INVERSION.md + " +
               "docs/proofs/PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md + " +
               "docs/proofs/PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md + " +
               "docs/proofs/PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md + " +
               "docs/proofs/PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (KleinEightCellClaim, typed Cubic3 parent)")
    {
        KleinEightParent = klein8 ?? throw new ArgumentNullException(nameof(klein8));
    }

    public override string DisplayName =>
        "F111 hard-cell pure-D template rule (Tier1Derived, F110 Aspect B corollary)";

    public override string Summary =>
        $"Theorem: at k = N = 4 in diagonal Klein cell (D.BitA(), D.BitB()), pair is F87-hard ⟺ at least one term is a pure-D template " +
        $"(only D and I, no other non-I letters). Implies F110 Aspect B 228:0 Y-inversion as corollary. " +
        $"Per-cell decomposition: 36 Pure+Pure + 192 Pure+Mixed + 0 Mixed+Mixed = 228 hard. " +
        $"Empirical anchor: 1584 classifications across 3 dephase letters (N = 4, k = 4), all matching, zero exceptions. " +
        $"Promoted 2026-06-10: the hard-direction converse gate closed (WindowedConverseAllGammaClaim, no residual) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("Y-inversion corollary", summary: YInversionCorollary);
            yield return new InspectableNode("F110 Aspect B corollary", summary: F110AspectBCorollary);
            yield return new InspectableNode("Per-cell decomposition", summary: DecompositionPerCell);
            yield return new InspectableNode("Subclaim (a): pure-D single-term H is HARD",
                summary: SubclaimA_PureDSingleTermHard);
            yield return new InspectableNode("Subclaim (b): mixed single-term H is SOFT",
                summary: SubclaimB_MixedSingleTermSoft);
            yield return new InspectableNode("Subclaim (c): (Pure-D, Mixed) pair is HARD",
                summary: SubclaimC_PureMixedPairHard);
            yield return new InspectableNode("Subclaim (d) closed modulo M: (Mixed, Mixed) pair is SOFT",
                summary: SubclaimD_MixedMixedPairSoft);
            yield return new InspectableNode("Sister claims on YParity axis",
                summary: "F102 (YParityIndependenceAtK3, Tier1Derived). F103 (F87Z2CubedRefinementN4K3, Tier1Derived). " +
                         "F105 (F87Z2CubedRefinementN5K3, Tier1Derived). F106 (F87Z2CubedRefinementN4K4, Tier1Derived). " +
                         "F107 (TrulyYParityZeroPurity, Tier1Derived). F109 (MotherSoftYParityOnePurity, Tier1Derived unconditional). " +
                         "F110 (HardCellYInversionPattern, Tier1Derived since 2026-06-10): Y-inversion across k=3 and k=4. " +
                         "F111 (THIS, Tier1Derived since 2026-06-10): sharper per-pair structural rule than F110 Aspect B; F111 implies F110 Aspect B at k=N=4 as immediate corollary.");
            yield return new InspectableNode("Cross-axis dependencies (BitB and BitA): F108 Parts",
                summary: "F108 Part 1+3 (BitB-axis): close F107/F109/F110 derivation via Π_5bilinear under Z and Y dephasing. F108 Part 2 (BitA-axis, BitA twin of Part 1): closes the X-deph branch via the Z↔X Π² mirror. F108 Parts are NOT YParity-axis sisters (per their Z2Axis declarations); they are the cross-axis closure mechanism that grounds F107/F109/F110/F111's diagonal-cell scope.");
            yield return new InspectableNode("Promotion record (2026-06-10) + open siblings",
                summary: "Hard-direction converse behind subclaims (a)/(c), the F111 promotion gate: it reduced to the windowed converse typed as WindowedConverseAllGammaClaim (RCPsiSquared.Diagnostics.F87), CLOSED 2026-06-10 with no residual (girth dichotomy retired R-deg, Pascal-Gram positivity resolved R-sign) ⟹ F111 Tier1Derived. Subclaim (d) Mixed+Mixed = soft is CLOSED modulo M via PROOF_F103 §7.4. " +
                         "F110 Aspect C (k=3 ratio 42:8): derived by the F103 §6 counting rule + §7 bipartite mechanism. " +
                         "Pure-D Template Rule at k > 4 or N > 4: empirically unverified. " +
                         "Hardware QPU confirmation at k ≥ 3: open (no F87 QPU confirmations beyond Marrakesh k = 2).");
            yield return new InspectableNode("Cubic3 anchor parent",
                summary: $"KleinEightCellClaim ({KleinEightParent.Tier.Label()}): the Z₂³ 8-cell decomposition anchor for the y_par axis F111 lives on.");
        }
    }
}
