using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F103 (Tier1Derived): Z₂³ refinement of the F87 trichotomy at k_body=3
/// (N=4 empirical anchor). After F102 established that Y-parity is independent
/// of the Klein (bit_a, bit_b) signature at k_body≥3, this Claim captures the
/// observed structure of how Y-parity refines the existing F87 trichotomy
/// (truly / soft / hard) cells.
///
/// <para>Empirical anchor: 294 Z₂³-homogeneous + Y-parity-homogeneous k_body=3
/// Pauli pairs at N=4, classified under each of Z, X, Y single-letter dephasing
/// via the Python framework's <c>classify_pauli_pair</c>, bucketed by
/// (Klein × dephase letter × y_par × trichotomy class). Five structural patterns
/// crystallized:</para>
///
/// <list type="number">
///   <item>truly is y_par=0-pure across all 300 truly classifications.</item>
///   <item>hard in diagonal Klein cells (Klein matches dephase letter) splits
///         42:8 with Y-inversion: Z-deph (0,1) = (42, 8); X-deph (1,0) = (42, 8);
///         Y-deph (1,1) = (8, 42) inverted because Y carries y_par=1.</item>
///   <item>Same diagonal Klein cells contain a soft 13:13 y_par-symmetric split.</item>
///   <item>Mother sector (0,0) soft is y_par=1-pure: all 3 letter cells = (0, 21).</item>
///   <item>Off-diagonal soft cells (6 cells, Klein non-mother and Klein ≠ dephase
///         Klein) split into Pattern B (proportional to enumeration breakdown,
///         3 cells) and Pattern C (y_par=1-pure, 3 cells).</item>
/// </list>
///
/// <para>Implements <see cref="IZ2AxisClaim"/> with <see cref="Z2Axis.YParity"/>;
/// <see cref="BitATwin"/> is null (Y-axis Claims do not carry a BitA-twin slot).
/// The five frozen-count records are exposed as properties and as inspectable
/// children for <c>rcpsi inspect --claim F87Z2CubedRefinement</c>.</para>
///
/// <para>Regenerate via: <c>simulations/f87_z2cubed_split_n4_k3.py</c> (~60s).
/// Proof / empirical anchor: <c>docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md</c>.</para></summary>
public sealed class F87Z2CubedRefinement : Claim, IZ2AxisClaim
{
    public Z2Axis Z2Axis => Z2Axis.YParity;
    public Claim? BitATwin => null;

    public TrulyYParityZeroPurity TrulyPurity { get; }
    public HardDiagonalSplit42To8 HardDiagonal { get; }
    public DiagonalSoftSplit13To13 DiagonalSoft { get; }
    public MotherSoftYParityOnePurity MotherSoft { get; }
    public OffDiagonalSoftPatterns OffDiagonalSoft { get; }

    public F87Z2CubedRefinement()
        : base("F103 Z₂³ refinement of the F87 trichotomy at k_body=3 (N=4 anchor, 294 pairs, 5 sub-statements)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F103 + " +
               "docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md + " +
               "simulations/f87_z2cubed_split_n4_k3.py")
    {
        TrulyPurity = new TrulyYParityZeroPurity(
            TotalTrulyClassifications: 300,
            YParityOneCount: 0);

        HardDiagonal = new HardDiagonalSplit42To8(
            ZDephKlein01: (42, 8),
            XDephKlein10: (42, 8),
            YDephKlein11: (8, 42));

        DiagonalSoft = new DiagonalSoftSplit13To13(
            ZDephKlein01: (13, 13),
            XDephKlein10: (13, 13),
            YDephKlein11: (13, 13));

        MotherSoft = new MotherSoftYParityOnePurity(
            ZDephCounts: (0, 21),
            XDephCounts: (0, 21),
            YDephCounts: (0, 21));

        OffDiagonalSoft = new OffDiagonalSoftPatterns(
            new Dictionary<(int KleinA, int KleinB, char Dephase), (int YPar0, int YPar1)>
            {
                // Pattern B (proportional to enum breakdown)
                { (0, 1, 'Y'), (55, 21) },
                { (1, 1, 'Z'), (21, 55) },
                { (1, 1, 'X'), (21, 55) },
                // Pattern C (y_par=1-pure)
                { (0, 1, 'X'), (0, 21) },
                { (1, 0, 'Z'), (0, 21) },
                { (1, 0, 'Y'), (0, 21) },
            });
    }

    public override string DisplayName =>
        "F103 F87 Z₂³ refinement (N=4 k=3, 294 pairs, 5 sub-statements)";

    public override string Summary =>
        $"Y-parity refinement of F87 trichotomy: truly y_par=0-pure (300 classifications), " +
        $"hard diagonal 42:8 with Y-inversion, diagonal soft 13:13 symmetric, mother soft y_par=1-pure, " +
        $"off-diagonal soft Pattern B (proportional, 3) + Pattern C (y_par=1-pure, 3) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Truly y_par=0 purity",
                summary: $"Across all 12 (Klein × dephase) cells, y_par=1 truly count = {TrulyPurity.YParityOneCount}; " +
                         $"total truly classifications = {TrulyPurity.TotalTrulyClassifications}");
            yield return new InspectableNode("Hard diagonal split 42:8 with Y-inversion",
                summary: $"Z-deph Klein (0,1) hard = {HardDiagonal.ZDephKlein01}; " +
                         $"X-deph Klein (1,0) hard = {HardDiagonal.XDephKlein10}; " +
                         $"Y-deph Klein (1,1) hard = {HardDiagonal.YDephKlein11} (inverted)");
            yield return new InspectableNode("Diagonal soft 13:13 y_par-symmetric",
                summary: $"Z-deph (0,1) soft = {DiagonalSoft.ZDephKlein01}; " +
                         $"X-deph (1,0) soft = {DiagonalSoft.XDephKlein10}; " +
                         $"Y-deph (1,1) soft = {DiagonalSoft.YDephKlein11}");
            yield return new InspectableNode("Mother (0,0) soft y_par=1-pure",
                summary: $"Z-deph = {MotherSoft.ZDephCounts}; X-deph = {MotherSoft.XDephCounts}; Y-deph = {MotherSoft.YDephCounts}");
            yield return new InspectableNode("Off-diagonal soft patterns (6 cells)",
                summary: $"3 Pattern B (proportional to enum) + 3 Pattern C (y_par=1-pure); " +
                         $"keys: (KleinA, KleinB, Dephase) -> (y0, y1); see Cells property");
            yield return new InspectableNode("Scope (out of scope)",
                summary: "C# k_body≥3 classifier lift, N>4 or k>3 generalization, closed-form derivation " +
                         "of 42:8, hardware confirmation of k≥3 F87. See spec out-of-scope section.");
        }
    }
}

public sealed record TrulyYParityZeroPurity(int TotalTrulyClassifications, int YParityOneCount);

public sealed record HardDiagonalSplit42To8(
    (int YPar0, int YPar1) ZDephKlein01,
    (int YPar0, int YPar1) XDephKlein10,
    (int YPar0, int YPar1) YDephKlein11);

public sealed record DiagonalSoftSplit13To13(
    (int YPar0, int YPar1) ZDephKlein01,
    (int YPar0, int YPar1) XDephKlein10,
    (int YPar0, int YPar1) YDephKlein11);

public sealed record MotherSoftYParityOnePurity(
    (int YPar0, int YPar1) ZDephCounts,
    (int YPar0, int YPar1) XDephCounts,
    (int YPar0, int YPar1) YDephCounts);

public sealed record OffDiagonalSoftPatterns(
    IReadOnlyDictionary<(int KleinA, int KleinB, char Dephase), (int YPar0, int YPar1)> Cells);
