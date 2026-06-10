using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 surface for F115, the windowed F87 hardness combinatorial theory, as a single
/// registry-queryable Claim. Wraps the static <see cref="WindowedObstructionScan"/> helper: a k-body
/// diagonal-cell Mixed term's X/Y positions read as a GF(2)[x] polynomial p(x) (bit j ↦ x^j), even
/// popcount so (1 + x) | p. The whole §7 diagonal-cell hard/soft line collapses to one integer test
/// (PROOF_F103 §7.7): a pair is HARD iff its two X/Y window-masks have different (1 + x)-adic
/// valuations, SOFT iff equal. The obstruction (minimal odd 𝔽₂-relation) has a derived size law
/// max = min(2W − 1, 2k − 3), W = N − k + 1 windows.
///
/// <para>F115 is the mask-combinatorial reading of the same diagonal-cell hard/soft split that
/// <see cref="F87DiagonalCellBipartiteWitnessSet"/> reads through the chiral-K bipartite criterion,
/// so it is wired as that set's child (single typed parent). No <c>ChainSystem</c> is involved: the
/// criterion is pure GF(2) bit arithmetic on the term masks (no Hamiltonian, no Liouvillian), which
/// is why the helper scales far past the eigendecomposition range.</para>
///
/// <para>Tier: Tier1Derived. The hard/soft VERDICT rests on the §7.5/§7.6 converse
/// (non-bipartite ⟹ hard), derived modulo standard perturbation theory (the first-order-block
/// premise closed via degenerate-PT + analyticity, §7.6); the obstruction size-law and counts are
/// pure derived combinatorics, and per §7.10 the obstruction size is spectrally inert (only the
/// (1 + x)-valuation crosses into the spectrum). Promoted from Tier1Candidate in the formal
/// promotion pass (2026-06-08), together with the parent F87DiagonalCellBipartiteWitnessSet.</para>
///
/// <para>Scope of the Tier1Derived label: the criterion itself (valuations-differ ⟺ non-bipartite is
/// exact GF(2)) plus the genericity result (non-bipartite ⟹ hard for all but finitely many γ). The
/// all-γ closure including the physical operating point is isolated as
/// <see cref="WindowedConverseAllGammaClaim"/> (Tier1Derived since 2026-06-10, no residual: Pascal-Gram
/// positivity); it is a downstream strengthening, not a parent of this claim, so the tiers stay consistent.</para>
///
/// <para>The per-size shape of the obstruction (the MacWilliams kernel) splits into a closed floor
/// and a number-theoretic middle: the size-3 floor closes: c(D,3)=3D−1 per reduced max-degree,
/// d=0 base T0(k)=(k−1)²(k−2)/2 (<see cref="WindowedObstructionScan.TriangleHardCountBaseD0"/>),
/// d-layered by the same 2^(d-1) reduction; the monomial column (one generator a monomial x^j) is
/// polynomial of degree β−1 in the other's weight β; the ceiling (max size 2D+1) is the repunit
/// pair (count → 2). The middle (size ≥ 5) is genuinely non-polynomial, it is the distribution of
/// weighted coprime polynomial pairs in GF(2)[x] (the closed-form frontier is sharp: a size cell is
/// polynomial iff some popcount split carries a monomial; first irregular cell (3,2) at size 5).
/// See experiments/F115_OBSTRUCTION_DISTRIBUTION.md Finding 6.</para>
///
/// <para>Anchor: <c>docs/ANALYTICAL_FORMULAS.md</c> F115 +
/// <c>docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md</c> §7.7-7.9 +
/// <see cref="WindowedObstructionScan"/>.</para></summary>
public sealed class WindowedHardnessClaim : Claim
{
    /// <summary>One self-check in the witness battery: a named assertion about
    /// <see cref="WindowedObstructionScan"/> with its expected vs actual verdict.</summary>
    public readonly record struct BatteryCase(
        string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    public IReadOnlyList<BatteryCase> Cases { get; }

    public WindowedHardnessClaim()
        : base("F115 windowed hardness: one-number (1+x)-valuation criterion + obstruction size law",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F115 + " +
               "docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md §7.7-7.9 + " +
               "compute/RCPsiSquared.Diagnostics/F87/WindowedObstructionScan.cs")
    {
        Cases = BuildBattery();
    }

    public int PassCount => Cases.Count(c => c.Passes);

    public override string DisplayName =>
        $"F115 windowed hardness (GF(2)[x] (1+x)-valuation criterion, {Cases.Count} witnesses)";

    public override string Summary =>
        "hard iff the two X/Y masks have different (1+x)-adic valuations; obstruction size " +
        $"min(2W-1, 2k-3); {PassCount}/{Cases.Count} PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("criterion",
                summary: "IsHardPair(p1, p2) = ValuationAtOnePlusX(p1) != ValuationAtOnePlusX(p2)");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " +
                             (c.Passes ? "PASS" : "FAIL"));
        }
    }

    /// <summary>The witness battery, each case self-checked against <see cref="WindowedObstructionScan"/>:
    /// a hard pair with DIFFERING (1+x)-valuations, a soft pair with EQUAL (1+x)-valuations, and a
    /// size-law case (the extremal body-bound family saturating 2k-3, plus the windowed
    /// min(2W-1, 2k-3) maximum at k=3, N=4). Mask values are the ones the
    /// <see cref="WindowedObstructionScanTests"/> family pins.</summary>
    private static IReadOnlyList<BatteryCase> BuildBattery()
    {
        var cases = new List<BatteryCase>();

        // Hard: differing (1+x)-valuations (extremal body-bound family at k=4):
        // p1 = x + x^3, p2 = 1 + x^3. The terms I·X·I·Y and X·I·I·Y.
        const int hardK = 4;
        ulong hard1 = (1UL << 1) | (1UL << (hardK - 1)); // x + x^{k-1}
        ulong hard2 = 1UL | (1UL << (hardK - 1));        // 1 + x^{k-1}
        bool hardVerdict = WindowedObstructionScan.IsHardPair(hard1, hard2);
        int hv1 = WindowedObstructionScan.ValuationAtOnePlusX(hard1);
        int hv2 = WindowedObstructionScan.ValuationAtOnePlusX(hard2);
        cases.Add(new BatteryCase(
            Name: "hard: differing (1+x)-valuations",
            Detail: $"k={hardK}, p1={Bin(hard1)} (v={hv1}), p2={Bin(hard2)} (v={hv2})",
            Expected: "hard",
            Actual: (hardVerdict && hv1 != hv2) ? "hard" : "soft"));

        // Soft: equal (1+x)-valuations, two DISTINCT masks (1+x = 0b11 and x+x^2 = 0b110, both v=1).
        ulong soft1 = 0b11;  // 1 + x      (v = 1)
        ulong soft2 = 0b110; // x + x^2    (v = 1)
        bool softVerdict = WindowedObstructionScan.IsHardPair(soft1, soft2);
        int sv1 = WindowedObstructionScan.ValuationAtOnePlusX(soft1);
        int sv2 = WindowedObstructionScan.ValuationAtOnePlusX(soft2);
        cases.Add(new BatteryCase(
            Name: "soft: equal (1+x)-valuations",
            Detail: $"p1={Bin(soft1)} (v={sv1}), p2={Bin(soft2)} (v={sv2})",
            Expected: "soft",
            Actual: (!softVerdict && sv1 == sv2) ? "soft" : "hard"));

        // Size law (body bound): the extremal family saturates 2k-3 via GcdFormulaSize.
        int bodyBound = 2 * hardK - 3;
        int gcdSize = WindowedObstructionScan.GcdFormulaSize(hard1, hard2);
        cases.Add(new BatteryCase(
            Name: "size law: extremal family saturates 2k-3",
            Detail: $"k={hardK}, GcdFormulaSize(p1, p2)={gcdSize}, 2k-3={bodyBound}",
            Expected: bodyBound.ToString(CultureInfo.InvariantCulture),
            Actual: gcdSize.ToString(CultureInfo.InvariantCulture)));

        // Size law (windowed maximum): max minimal-odd-cycle == min(2W-1, 2k-3) at k=3, N=4.
        const int sk = 3, sn = 4;
        int sw = sn - sk + 1;
        int expectedMax = Math.Min(2 * sw - 1, 2 * sk - 3);
        var scan = WindowedObstructionScan.Scan(sk, sn);
        int actualMax = scan.MinOddCycleSizes.Count == 0 ? -1 : scan.MinOddCycleSizes.Keys.Max();
        cases.Add(new BatteryCase(
            Name: "size law: windowed max = min(2W-1, 2k-3)",
            Detail: $"k={sk}, N={sn}, W={sw}, max-obstruction={actualMax}, min(2W-1,2k-3)={expectedMax}",
            Expected: expectedMax.ToString(CultureInfo.InvariantCulture),
            Actual: actualMax.ToString(CultureInfo.InvariantCulture)));

        // Size-3 floor (MacWilliams-kernel floor): the d-layered triangle counts sum to the all-d
        // TriangleHardMaskCount. T0(k)=(k-1)^2(k-2)/2 is the d=0 base; c(D,3)=3D-1 the per-max-degree law.
        const int fk = 5;
        long triSum = 0;
        for (int d = 0; d <= fk - 3; d++) triSum += WindowedObstructionScan.TriangleHardCountByGRestDegree(fk, d);
        cases.Add(new BatteryCase(
            Name: "size-3 floor: d-layered triangle counts sum to TriangleHardMaskCount",
            Detail: $"k={fk}, Σ_d TriangleHardCountByGRestDegree = {triSum}, " +
                    $"TriangleHardMaskCount = {WindowedObstructionScan.TriangleHardMaskCount(fk)}, " +
                    $"d=0 base T0 = {WindowedObstructionScan.TriangleHardCountBaseD0(fk)}",
            Expected: WindowedObstructionScan.TriangleHardMaskCount(fk).ToString(CultureInfo.InvariantCulture),
            Actual: triSum.ToString(CultureInfo.InvariantCulture)));

        return cases;
    }

    private static string Bin(ulong mask) =>
        "0b" + Convert.ToString((long)mask, 2);
}
