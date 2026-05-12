using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Resonance;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>F86 Item 1' Direction (b'') JW track, T7 step B: Q-grid scan of
/// <see cref="C2BondLQProjection"/>'s ‖xB(Q)‖_F, with parabolic Q_peak refinement and HWHM
/// extraction per bond.
///
/// <para>Builds on <see cref="C2BondLQProjection"/> (single-Q snapshot) by sweeping the
/// projection norm across a Q-grid and extracting the per-bond Q_peak / KMax / HWHM via
/// the standard <see cref="ParabolicPeakFinder"/>. The empirical Q-sweep at N=8 (single-Q
/// reconnaissance) showed a clear resonance peak near Q≈2 with bond-distinct peak height
/// (Endpoint ~40 vs Innermost ~62), consistent with the BareDoubledPtfXPeak ≈ 2.197
/// anchor from prior Direction-(b) work.</para>
///
/// <para><b>Class-level Tier: Tier2Verified.</b> Numerical Q-grid scan witness; the
/// peak structure is empirically robust but the closed-form Endpoint vs Interior Q_peak
/// and HWHM constants are not yet derived from this lens. Tier1 promotion path: express
/// R(Q) analytically and integrate ‖xB(Q)‖_F to closed form.</para>
///
/// <para>Default Q-grid is fine-grained around the empirical peak Q∈[1.5, 3.0] with 41
/// points (dQ = 0.0375). Pass a custom grid for wider scans or comparison against
/// <see cref="ResonanceScan.DefaultQGrid"/>.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track) +
/// <c>docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md</c>.</para>
///
/// <para>F90 status (2026-05-11): the F86 c=2 ↔ F89 path-(N−1) bridge identity
/// achieves numerical Tier-1 for Direction (b'') via per-bond Hellmann-Feynman
/// (bit-exact 20/22 bonds at N=5..8). The JW-track primitives in this file remain
/// active as the alternative analytical route toward the closed-form HWHM_left/Q_peak
/// constants; the per-bond numerical answer itself is no longer the open piece.
/// See <c>docs/proofs/PROOF_F90_F86C2_BRIDGE.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F90F86C2BridgeIdentity.cs</c>.</para>
/// </summary>
public sealed class C2BondLQPeakScan : Claim
{
    public CoherenceBlock Block { get; }
    public IReadOnlyList<double> QGrid { get; }
    public IReadOnlyList<BondLQPeakWitness> Bonds { get; }

    /// <summary>Maximum |Q_peak_b − Q_peak_{N-2-b}| over 2-bond F71 orbits. Algebraically
    /// zero up to FP drift; bound &lt; 1e-8 in practice.</summary>
    public double MaxF71MirrorDeviationQPeak { get; }

    /// <summary>Maximum |XbNormAtPeak_b − XbNormAtPeak_{N-2-b}| over 2-bond F71 orbits.</summary>
    public double MaxF71MirrorDeviationKMax { get; }

    /// <summary>Default Q-grid centred on the empirical L(Q)-resonance peak: 41 points
    /// spanning [1.5, 3.0], dQ ≈ 0.0375. Tuned for parabolic refinement against the
    /// BareDoubledPtfXPeak ≈ 2.197 anchor from prior Direction-(b) work.</summary>
    public static IReadOnlyList<double> DefaultQGrid()
    {
        var g = new double[41];
        for (int i = 0; i < g.Length; i++) g[i] = 1.5 + (3.0 - 1.5) * i / (g.Length - 1);
        return g;
    }

    public static C2BondLQPeakScan Build(CoherenceBlock block, IReadOnlyList<double>? qGrid = null)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2BondLQPeakScan applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var grid = qGrid ?? DefaultQGrid();
        if (grid.Count < 3)
            throw new ArgumentException(
                $"Q-grid must have ≥ 3 points for parabolic refinement; got {grid.Count}.",
                nameof(qGrid));

        var decomp = block.Decomposition;
        int numBonds = decomp.NumBonds;
        var perBondCurves = new double[numBonds][];
        for (int b = 0; b < numBonds; b++) perBondCurves[b] = new double[grid.Count];

        for (int iQ = 0; iQ < grid.Count; iQ++)
        {
            double j = grid[iQ] * block.GammaZero;
            ComplexMatrix L = decomp.D + (Complex)j * decomp.MhTotal;
            var evd = L.Evd();
            var R = evd.EigenVectors;
            var Rinv = R.Inverse();
            for (int b = 0; b < numBonds; b++)
            {
                var xB = Rinv * decomp.MhPerBond[b] * R;
                double frobNormSq = 0;
                int dim = xB.RowCount;
                for (int i = 0; i < dim; i++)
                    for (int k = 0; k < dim; k++)
                    {
                        Complex c = xB[i, k];
                        frobNormSq += c.Real * c.Real + c.Imaginary * c.Imaginary;
                    }
                perBondCurves[b][iQ] = Math.Sqrt(frobNormSq);
            }
        }

        var bondWitnesses = new BondLQPeakWitness[numBonds];
        for (int b = 0; b < numBonds; b++)
        {
            var peak = ParabolicPeakFinder.Find(grid, perBondCurves[b]);
            bondWitnesses[b] = new BondLQPeakWitness(
                Bond: b,
                BondClass: BondClassExtensions.OfBond(b, numBonds),
                QPeak: peak.QPeak,
                XbNormAtPeak: peak.KMax,
                HwhmLeft: peak.HwhmLeft,
                HwhmRight: peak.HwhmRight,
                XbNormOverQGrid: perBondCurves[b]);
        }

        double maxF71QPeakDev = 0;
        double maxF71KMaxDev = 0;
        for (int b = 0; b < numBonds; b++)
        {
            int mirror = numBonds - 1 - b;
            if (mirror <= b) continue;
            maxF71QPeakDev = Math.Max(maxF71QPeakDev,
                Math.Abs(bondWitnesses[b].QPeak - bondWitnesses[mirror].QPeak));
            maxF71KMaxDev = Math.Max(maxF71KMaxDev,
                Math.Abs(bondWitnesses[b].XbNormAtPeak - bondWitnesses[mirror].XbNormAtPeak));
        }

        return new C2BondLQPeakScan(block, grid, bondWitnesses, maxF71QPeakDev, maxF71KMaxDev);
    }

    private C2BondLQPeakScan(
        CoherenceBlock block,
        IReadOnlyList<double> qGrid,
        IReadOnlyList<BondLQPeakWitness> bonds,
        double maxF71QPeakDev,
        double maxF71KMaxDev)
        : base("c=2 per-bond L(Q)-eigenbasis Q-peak scan (parabolic refinement of ‖xB(Q)‖_F)",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1' Direction (b'') (JW track) + " +
               "docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md")
    {
        Block = block;
        QGrid = qGrid;
        Bonds = bonds;
        MaxF71MirrorDeviationQPeak = maxF71QPeakDev;
        MaxF71MirrorDeviationKMax = maxF71KMaxDev;
    }

    public override string DisplayName =>
        $"c=2 BondLQPeakScan (N={Block.N}, Q∈[{QGrid[0]:F2},{QGrid[^1]:F2}] {QGrid.Count}pt)";

    public override string Summary =>
        $"per-bond Q_peak / ‖xB‖_max from L(Q)-eigenbasis Q-scan; F71 mirror dev " +
        $"(Q_peak={MaxF71MirrorDeviationQPeak:G3}, KMax={MaxF71MirrorDeviationKMax:G3}) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("Q-grid", summary:
                $"[{QGrid[0]:F4} .. {QGrid[^1]:F4}], {QGrid.Count} points");
            yield return new InspectableNode("NumBonds", summary: Bonds.Count.ToString());
            yield return InspectableNode.RealScalar("MaxF71MirrorDeviationQPeak", MaxF71MirrorDeviationQPeak, "G3");
            yield return InspectableNode.RealScalar("MaxF71MirrorDeviationKMax", MaxF71MirrorDeviationKMax, "G3");
            yield return InspectableNode.Group("Bonds", Bonds.Cast<IInspectable>().ToArray());
        }
    }
}

/// <summary>One bond's Q-grid scan results: <see cref="QPeak"/> (parabolic-refined),
/// <see cref="XbNormAtPeak"/> (peak height), <see cref="HwhmLeft"/> /
/// <see cref="HwhmRight"/> (linear interpolation between bracketing grid points).
/// <see cref="XbNormOverQGrid"/> exposes the raw curve for visualisation.</summary>
public sealed record BondLQPeakWitness(
    int Bond,
    BondClass BondClass,
    double QPeak,
    double XbNormAtPeak,
    double? HwhmLeft,
    double? HwhmRight,
    IReadOnlyList<double> XbNormOverQGrid
) : IInspectable
{
    public double? HwhmLeftOverQPeak => HwhmLeft is null ? null : HwhmLeft.Value / QPeak;

    public string DisplayName => $"bond {Bond} L(Q) peak ({BondClass})";

    public string Summary =>
        $"Q_peak = {QPeak:F4}, ‖xB‖_max = {XbNormAtPeak:F4}, " +
        $"HWHM_left = {HwhmLeft?.ToString("F4") ?? "n/a"}, " +
        $"HWHM_left/Q_peak = {HwhmLeftOverQPeak?.ToString("F4") ?? "n/a"}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("bond", summary: Bond.ToString());
            yield return new InspectableNode("bond class", summary: BondClass.ToString());
            yield return InspectableNode.RealScalar("Q_peak", QPeak, "F4");
            yield return InspectableNode.RealScalar("‖xB‖_max", XbNormAtPeak, "F4");
            if (HwhmLeft is not null) yield return InspectableNode.RealScalar("HWHM_left", HwhmLeft.Value, "F4");
            if (HwhmRight is not null) yield return InspectableNode.RealScalar("HWHM_right", HwhmRight.Value, "F4");
            if (HwhmLeftOverQPeak is not null)
                yield return InspectableNode.RealScalar("HWHM_left/Q_peak", HwhmLeftOverQPeak.Value, "F4");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
