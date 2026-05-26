using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F103 (Tier1Derived): Z₂³ refinement of the F87 trichotomy anchored at
/// N=4, k=3 (294 Pauli pairs). First concrete derived class of
/// <see cref="F87Z2CubedRefinementBase"/>; sibling: <c>F87Z2CubedRefinementN5K3</c>
/// (F105, planned).
///
/// <para>Empirical anchor: 294 Z₂³-homogeneous + Y-parity-homogeneous k_body=3
/// Pauli pairs at N=4, classified under each of Z, X, Y single-letter dephasing
/// via the Python framework's <c>classify_pauli_pair</c> (and re-verified in C#
/// by F104), bucketed by (Klein × dephase letter × y_par × trichotomy class).
/// Five structural patterns:</para>
///
/// <list type="number">
///   <item>truly is y_par=0-pure across all 300 truly classifications.</item>
///   <item>hard in diagonal Klein cells (Klein matches dephase letter) splits
///         42:8 with Y-inversion: Z-deph (0,1) = (42, 8); X-deph (1,0) = (42, 8);
///         Y-deph (1,1) = (8, 42) inverted because Y carries y_par=1.</item>
///   <item>Same diagonal Klein cells contain a soft 13:13 y_par-symmetric split.</item>
///   <item>Mother sector (0,0) soft is y_par=1-pure: all 3 letter cells = (0, 21).</item>
///   <item>Off-diagonal soft cells (6 cells) split into Pattern B (proportional,
///         3 cells) and Pattern C (y_par=1-pure, 3 cells).</item>
/// </list>
///
/// <para>Regenerate via: <c>simulations/f87_z2cubed_split_n4_k3.py</c> (~60s),
/// or C# re-verification via <c>F104KBodyTrichotomyVerificationTests</c>.
/// Proof: <c>docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md</c>.</para></summary>
public sealed class F87Z2CubedRefinementN4K3 : F87Z2CubedRefinementBase
{
    public override int N => 4;
    public override int K => 3;
    public override int TotalPairs => 294;

    public override TrulyCounts TrulyPurity { get; }
    public override HardDiagonalSplit HardDiagonal { get; }
    public override DiagonalSoftSplit DiagonalSoft { get; }
    public override MotherSoftCounts MotherSoft { get; }
    public override OffDiagonalSoftPatterns OffDiagonalSoft { get; }

    public F87Z2CubedRefinementN4K3(KleinEightCellClaim klein8)
        : base("F103 Z₂³ refinement of the F87 trichotomy at N=4 k=3 (294 pairs, 5 sub-statements); typed Cubic3 parent = KleinEightCellClaim",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F103 + " +
               "docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md + " +
               "simulations/f87_z2cubed_split_n4_k3.py + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (KleinEightCellClaim, typed Cubic3 parent)",
               klein8)
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
        "F103 F87 Z₂³ refinement (N=4 k=3, 294 pairs, 5 sub-statements)";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            foreach (var child in base.ExtraChildren)
                yield return child;
            yield return new InspectableNode("Scope (closed and out-of-scope items)",
                summary: "Closed by F104 (C# k_body=3 classifier lift verified bit-exact). " +
                         "Out of scope: N>4 or k>3 generalization (F105+ series), closed-form " +
                         "derivation of 42:8, hardware confirmation of k>=3 F87.");
        }
    }
}
