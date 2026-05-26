using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Abstract base for the F87 Z₂³ refinement family of typed Claims.
/// Each concrete subclass anchors a specific (N, k) regime with frozen empirical
/// counts derived from re-classifying Pauli pairs under Z, X, Y dephasing.
///
/// <para>Family members (planned series):</para>
/// <list type="bullet">
///   <item>F87Z2CubedRefinementN4K3 (F103): N=4 k=3, 294 pairs</item>
///   <item>F87Z2CubedRefinementN5K3 (F105): N=5 k=3, 294 pairs (same enumeration, larger chain)</item>
///   <item>F87Z2CubedRefinementN4K4 (F106): N=4 k=4, 4248 pairs (new enumeration)</item>
/// </list>
///
/// <para>All members share the same five-sub-statement structure (Truly purity,
/// Hard diagonal split, Diagonal soft split, Mother soft purity, Off-diagonal
/// soft patterns). The 5 record TYPES carry historical F103-numerics-bearing
/// names; the values per (N, k) are set by each subclass.</para>
///
/// <para>Implements <see cref="IZ2AxisClaim"/> with <see cref="Z2Axis.YParity"/>;
/// <see cref="BitATwin"/> is null. Each concrete subclass registers separately
/// into PolarityCubeMap, so YParityClaims.Count grows by 1 per regime.</para></summary>
public abstract class F87Z2CubedRefinementBase : Claim, IZ2AxisClaim
{
    public Z2Axis Z2Axis => Z2Axis.YParity;
    public Claim? BitATwin => null;

    /// <summary>Chain length N (number of qubits) for this anchor.</summary>
    public abstract int N { get; }

    /// <summary>k-body order (number of letters per Pauli term) for this anchor.</summary>
    public abstract int K { get; }

    /// <summary>Total Klein-homogeneous + Y-par-homogeneous pair count at this (N, k).
    /// Independent of N for k=3 (294); k=4 yields 4248; etc.</summary>
    public abstract int TotalPairs { get; }

    public abstract TrulyCounts TrulyPurity { get; }
    public abstract HardDiagonalSplit HardDiagonal { get; }
    public abstract DiagonalSoftSplit DiagonalSoft { get; }
    public abstract MotherSoftCounts MotherSoft { get; }
    public abstract OffDiagonalSoftPatterns OffDiagonalSoft { get; }

    /// <summary>Typed Cubic3 parent: <see cref="KleinEightCellClaim"/>. The F87
    /// Z₂³ refinement family enumerates pairs across the 8-cell decomposition
    /// (Klein × dephase letter × y_par × trichotomy class); KleinEightCellClaim
    /// is the structural anchor. Wired 2026-05-26.</summary>
    public KleinEightCellClaim KleinEightParent { get; }

    protected F87Z2CubedRefinementBase(
        string displayName,
        Tier tier,
        string anchor,
        KleinEightCellClaim klein8)
        : base(displayName, tier, anchor)
    {
        KleinEightParent = klein8 ?? throw new ArgumentNullException(nameof(klein8));
    }

    public override string Summary =>
        $"Y-parity refinement of F87 trichotomy at N={N} k={K} ({TotalPairs} pairs): " +
        $"truly y_par=0-pure ({TrulyPurity.TotalTrulyClassifications} classifications, " +
        $"{TrulyPurity.YParityOneCount} y_par=1), " +
        $"hard diagonal {HardDiagonal.ZDephKlein01} with Y-inversion, " +
        $"diagonal soft {DiagonalSoft.ZDephKlein01} " +
        $"({(DiagonalSoft.ZDephKlein01.Item1 == DiagonalSoft.ZDephKlein01.Item2 ? "symmetric" : "asymmetric")}), " +
        $"mother soft {MotherSoft.ZDephCounts}, " +
        $"off-diagonal soft {OffDiagonalSoft.Cells.Count} cells ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("N", N);
            yield return InspectableNode.RealScalar("k", K);
            yield return InspectableNode.RealScalar("Total pairs", TotalPairs);
            yield return new InspectableNode("Truly y_par=0 purity",
                summary: $"Across all 12 (Klein × dephase) cells, y_par=1 truly count = {TrulyPurity.YParityOneCount}; " +
                         $"total truly classifications = {TrulyPurity.TotalTrulyClassifications}");
            yield return new InspectableNode("Hard diagonal split (matches dephase letter)",
                summary: $"Z-deph Klein (0,1) hard = {HardDiagonal.ZDephKlein01}; " +
                         $"X-deph Klein (1,0) hard = {HardDiagonal.XDephKlein10}; " +
                         $"Y-deph Klein (1,1) hard = {HardDiagonal.YDephKlein11}");
            yield return new InspectableNode("Diagonal soft split (same 3 cells as hard)",
                summary: $"Z-deph (0,1) soft = {DiagonalSoft.ZDephKlein01}; " +
                         $"X-deph (1,0) soft = {DiagonalSoft.XDephKlein10}; " +
                         $"Y-deph (1,1) soft = {DiagonalSoft.YDephKlein11}");
            yield return new InspectableNode("Mother (0,0) soft",
                summary: $"Z-deph = {MotherSoft.ZDephCounts}; X-deph = {MotherSoft.XDephCounts}; Y-deph = {MotherSoft.YDephCounts}");
            yield return new InspectableNode("Off-diagonal soft patterns",
                summary: $"{OffDiagonalSoft.Cells.Count} cells; keys: (KleinA, KleinB, Dephase) -> (y0, y1); see Cells property");
            yield return new InspectableNode("Cubic3 anchor parent",
                summary: $"KleinEightCellClaim ({KleinEightParent.Tier.Label()}): the (bit_a, bit_b, y_par) 8-cell decomposition that this F87 Z₂³ refinement enumerates pairs across.");
        }
    }
}

public sealed record TrulyCounts(int TotalTrulyClassifications, int YParityOneCount);

public sealed record HardDiagonalSplit(
    (int YPar0, int YPar1) ZDephKlein01,
    (int YPar0, int YPar1) XDephKlein10,
    (int YPar0, int YPar1) YDephKlein11);

public sealed record DiagonalSoftSplit(
    (int YPar0, int YPar1) ZDephKlein01,
    (int YPar0, int YPar1) XDephKlein10,
    (int YPar0, int YPar1) YDephKlein11);

public sealed record MotherSoftCounts(
    (int YPar0, int YPar1) ZDephCounts,
    (int YPar0, int YPar1) XDephCounts,
    (int YPar0, int YPar1) YDephCounts);

public sealed record OffDiagonalSoftPatterns(
    IReadOnlyDictionary<(int KleinA, int KleinB, char Dephase), (int YPar0, int YPar1)> Cells);
