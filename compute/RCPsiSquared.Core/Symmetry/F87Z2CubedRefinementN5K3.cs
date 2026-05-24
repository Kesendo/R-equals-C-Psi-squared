using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F105 (Tier1Derived): Z₂³ refinement of the F87 trichotomy anchored at
/// N=5, k=3 (same 294 Pauli-pair enumeration as F103, larger chain). Sibling derived
/// class of <see cref="F87Z2CubedRefinementBase"/>; parallel to
/// <see cref="F87Z2CubedRefinementN4K3"/> (F103) at the next chain length.
///
/// <para>F85 (k-body generalization) predicts the Π²-class trichotomy is N-stable
/// for any k. If F85's N-stability lifts to the y_par sub-refinement, F103's frozen
/// counts (300 truly, 42:8 hard with Y-inversion, 13:13 diagonal soft, 0:21 mother
/// soft, Pattern B+C off-diagonal) survive bit-exactly at N=5. The actual counts are
/// captured by re-running F104's C# classifier on all 294 pairs × 3 dephase letters
/// at N=5 (~3h batch); see <see cref="F87Z2CubedEnumerationN5K3Tool"/>.</para>
///
/// <para>Regenerate via:
/// <c>dotnet test "compute\RCPsiSquared.Diagnostics.Tests" --filter "Category=SLOW_F105_BATCH"</c>
/// (~3h, writes <c>simulations/results/f87_z2cubed_split_n5_k3_counts.json</c>).
/// Proof: <c>docs/proofs/PROOF_F105_F87_Z2_CUBED_REFINEMENT_N5K3.md</c>.</para></summary>
public sealed class F87Z2CubedRefinementN5K3 : F87Z2CubedRefinementBase
{
    public override int N => 5;
    public override int K => 3;
    public override int TotalPairs => 294;

    public override TrulyCounts TrulyPurity { get; }
    public override HardDiagonalSplit HardDiagonal { get; }
    public override DiagonalSoftSplit DiagonalSoft { get; }
    public override MotherSoftCounts MotherSoft { get; }
    public override OffDiagonalSoftPatterns OffDiagonalSoft { get; }

    public F87Z2CubedRefinementN5K3()
        : base("F105 Z₂³ refinement of the F87 trichotomy at N=5 k=3 (294 pairs, 5 sub-statements)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F105 + " +
               "docs/proofs/PROOF_F105_F87_Z2_CUBED_REFINEMENT_N5K3.md + " +
               "simulations/results/f87_z2cubed_split_n5_k3_counts.json")
    {
        TrulyPurity = new TrulyCounts(
            TotalTrulyClassifications: 300,
            YParityOneCount: 0);

        HardDiagonal = new HardDiagonalSplit(
            ZDephKlein01: (42, 8),
            XDephKlein10: (42, 8),
            YDephKlein11: (8, 42));

        DiagonalSoft = new DiagonalSoftSplit(
            ZDephKlein01: (13, 13),
            XDephKlein10: (13, 13),
            YDephKlein11: (13, 13));

        MotherSoft = new MotherSoftCounts(
            ZDephCounts: (0, 21),
            XDephCounts: (0, 21),
            YDephCounts: (0, 21));

        OffDiagonalSoft = new OffDiagonalSoftPatterns(
            new Dictionary<(int KleinA, int KleinB, char Dephase), (int YPar0, int YPar1)>
            {
                { (0, 1, 'Y'), (55, 21) },
                { (1, 1, 'Z'), (21, 55) },
                { (1, 1, 'X'), (21, 55) },
                { (0, 1, 'X'), (0, 21) },
                { (1, 0, 'Z'), (0, 21) },
                { (1, 0, 'Y'), (0, 21) },
            });
    }

    public override string DisplayName =>
        "F105 F87 Z₂³ refinement (N=5 k=3, 294 pairs, 5 sub-statements)";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            foreach (var child in base.ExtraChildren)
                yield return child;
            yield return new InspectableNode("Scope (out-of-scope items)",
                summary: "F106 covers N=4 k=4 (separate spec). F107+ open: N=6 k=3 (requires " +
                         "block-spectrum Classify; N=5 k=4 batch ~42h). Closed-form derivation " +
                         "of 42:8 (if N-stable) or N=5 split (if not) remains open.");
        }
    }
}
