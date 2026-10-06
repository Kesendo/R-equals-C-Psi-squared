using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F106 (Tier1Derived): Z₂³ refinement of the F87 trichotomy anchored at
/// N=4, k=4 (4248 Klein-homogeneous + Y-par-homogeneous Pauli pairs, new enumeration
/// vs F103's 294 at k=3). Third concrete derived class of
/// <see cref="F87Z2CubedRefinementBase"/>; siblings: <see cref="F87Z2CubedRefinementN4K3"/>
/// (F103, N=4 k=3) and <see cref="F87Z2CubedRefinementN5K3"/> (F105, N=5 k=3).
///
/// <para>F105 confirmed F85's N-stability prediction bit-exactly at the y_par axis
/// for k=3 (N=4 to N=5 counts identical). F106 tests k-stability: do F103's
/// structural patterns survive at larger k? F85 does NOT predict k-stability
/// of the y_par sub-refinement, and the Klein (0,0) enum balance shifts from
/// 45/21 at k=3 to 780/300 at k=4. The outcome (PROOF_F106): the two purity
/// statements and the off-diagonal B/C partition hold, the diagonal splits reshape.</para>
///
/// <para>Regenerate via:
/// <c>dotnet test "compute\RCPsiSquared.Diagnostics.Tests" --filter "Category=SLOW_F106_BATCH"</c>
/// (~2-3min PLINQ on 24 cores, ~40min sequential; writes
/// <c>simulations/results/f87_z2cubed_split_n4_k4_counts.json</c>).
/// Proof: <c>docs/proofs/PROOF_F106_F87_Z2_CUBED_REFINEMENT_N4K4.md</c>.</para></summary>
public sealed class F87Z2CubedRefinementN4K4 : F87Z2CubedRefinementBase
{
    public override int N => 4;
    public override int K => 4;
    public override int TotalPairs => 4248;

    public override TrulyCounts TrulyPurity { get; }
    public override HardDiagonalSplit HardDiagonal { get; }
    public override DiagonalSoftSplit DiagonalSoft { get; }
    public override MotherSoftCounts MotherSoft { get; }
    public override OffDiagonalSoftPatterns OffDiagonalSoft { get; }

    public F87Z2CubedRefinementN4K4(KleinEightCellClaim klein8)
        : base("F106 Z₂³ refinement of the F87 trichotomy at N=4 k=4 (4248 pairs, 5 sub-statements); typed Cubic3 parent = KleinEightCellClaim",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F106 + " +
               "docs/proofs/PROOF_F106_F87_Z2_CUBED_REFINEMENT_N4K4.md + " +
               "simulations/results/f87_z2cubed_split_n4_k4_counts.json + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (KleinEightCellClaim, typed Cubic3 parent)",
               klein8)
    {
        // Frozen counts derived from the SLOW_F106_BATCH PLINQ run
        // (simulations/results/f87_z2cubed_split_n4_k4_counts.json, 12744 classifications).
        TrulyPurity = new TrulyCounts(
            TotalTrulyClassifications: 3924,
            YParityOneCount: 0);

        HardDiagonal = new HardDiagonalSplit(
            ZDephKlein01: (228, 0),
            XDephKlein10: (228, 0),
            YDephKlein11: (0, 228));

        DiagonalSoft = new DiagonalSoftSplit(
            ZDephKlein01: (300, 528),
            XDephKlein10: (300, 528),
            YDephKlein11: (528, 300));

        MotherSoft = new MotherSoftCounts(
            ZDephCounts: (0, 300),
            XDephCounts: (0, 300),
            YDephCounts: (0, 300));

        OffDiagonalSoft = new OffDiagonalSoftPatterns(
            new Dictionary<(int KleinA, int KleinB, char Dephase), (int YPar0, int YPar1)>
            {
                // 6 off-diagonal cells at k=4: 3 Pattern C (y_par=1-pure) and 3 Pattern B,
                // the whole cell soft; at k=4 every non-mother cell holds 32 strings of each
                // y_par, so Pattern B reads (528, 528) where k=3 had (55, 21)/(21, 55)
                // (PROOF_F103 section 8).
                { (0, 1, 'X'), (0, 528) },
                { (0, 1, 'Y'), (528, 528) },
                { (1, 0, 'Y'), (0, 528) },
                { (1, 0, 'Z'), (0, 528) },
                { (1, 1, 'X'), (528, 528) },
                { (1, 1, 'Z'), (528, 528) },
            });
    }

    public override string DisplayName =>
        "F106 F87 Z₂³ refinement (N=4 k=4, 4248 pairs, 5 sub-statements)";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            foreach (var child in base.ExtraChildren)
                yield return child;
            yield return new InspectableNode("Scope (out-of-scope items)",
                summary: "The diagonal cell beyond N = 4 is read exactly by Theorem 2 of the complement connection " +
                         "(DiagonalCellComplementConnectionTests): hard (292, 128) by y_par at k = 4, N = 5, " +
                         "(356, 320) at N = 6, 7, 8. The full batches of every cell at N = 5, k = 4 and N = 6, " +
                         "k = 3 are not run. The off-diagonal B/C " +
                         "partition is derived at every k (PROOF_F103 section 8).");
        }
    }
}
