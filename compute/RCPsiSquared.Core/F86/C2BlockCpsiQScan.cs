using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 c=2 block CΨ Q-sweep with one bond perturbed: tests whether the 1/4
/// Mandelbrot boundary inherited at the c=2 level (per <see cref="C2BlockCpsiTrajectory"/>)
/// shows a Q_EP signature in the non-uniform-J regime.
///
/// <para><b>Setup:</b> uniform background J = Q·γ₀ across all bonds, plus an additional
/// perturbation Δ on a single bond at <see cref="PerturbedBondIndex"/>. Sweep Q across
/// the grid; at each Q snapshot CΨ_block at <see cref="SnapshotTime"/>. The expected
/// signature: at Q = Q_EP for the orbit containing the perturbed bond, the 2-level
/// effective EP coalesces; if F86 Q_EP ↔ R=CΨ² 1/4 inheritance is real (open Question A),
/// CΨ_block(Q) should show a feature (peak/dip/inflection) at Q = Q_EP.</para>
///
/// <para><b>Tier:</b> Tier2Verified — empirical witness sweep, no analytic prediction.
/// Q_EP per orbit values (N=5): Endpoint Q_EP_E ≈ 1.138, mid-chain Q_EP_I ≈ 0.674,
/// computed as Q_peak / BareDoubledPtfXPeak with x_peak = 2.196910.</para>
/// </summary>
public sealed class C2BlockCpsiQScan : Claim
{
    public CoherenceBlock Block { get; }
    public int PerturbedBondIndex { get; }
    public double PerturbationDelta { get; }
    public double SnapshotTime { get; }
    public IReadOnlyList<double> QGrid { get; }
    public IReadOnlyList<double> CBlockAtSnapshot { get; }
    public IReadOnlyList<double> PsiBlockAtSnapshot { get; }
    public IReadOnlyList<double> CPsiBlockAtSnapshot { get; }

    public C2BlockCpsiQScan(
        CoherenceBlock block, int perturbedBondIndex, double perturbationDelta,
        double snapshotTime,
        IReadOnlyList<double> qGrid,
        IReadOnlyList<double> cBlockAtSnapshot,
        IReadOnlyList<double> psiBlockAtSnapshot,
        IReadOnlyList<double> cPsiBlockAtSnapshot)
        : base("c=2 block CΨ Q-sweep with single-bond perturbation",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md Layer 2")
    {
        Block = block;
        PerturbedBondIndex = perturbedBondIndex;
        PerturbationDelta = perturbationDelta;
        SnapshotTime = snapshotTime;
        QGrid = qGrid;
        CBlockAtSnapshot = cBlockAtSnapshot;
        PsiBlockAtSnapshot = psiBlockAtSnapshot;
        CPsiBlockAtSnapshot = cPsiBlockAtSnapshot;
    }

    public static C2BlockCpsiQScan Build(
        CoherenceBlock block, int perturbedBondIndex, double perturbationDelta,
        double snapshotTime, IReadOnlyList<double> qGrid)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2BlockCpsiQScan applies only to c=2; got c={block.C}.", nameof(block));
        if (perturbedBondIndex < 0 || perturbedBondIndex >= block.NumBonds)
            throw new ArgumentOutOfRangeException(nameof(perturbedBondIndex),
                $"bond index must be in [0, {block.NumBonds - 1}]; got {perturbedBondIndex}.");

        var cTraj = new List<double>(qGrid.Count);
        var psiTraj = new List<double>(qGrid.Count);
        var cPsiTraj = new List<double>(qGrid.Count);

        var couplings = new double[block.NumBonds];
        var snapshotTimes = new[] { snapshotTime };

        foreach (double q in qGrid)
        {
            double background = q * block.GammaZero;
            for (int i = 0; i < block.NumBonds; i++) couplings[i] = background;
            couplings[perturbedBondIndex] = background + perturbationDelta;

            var trajectory = C2BlockCpsiTrajectory.BuildPerBond(
                block, q, couplings, snapshotTimes);

            cTraj.Add(trajectory.CBlockTrajectory[0]);
            psiTraj.Add(trajectory.PsiBlockTrajectory[0]);
            cPsiTraj.Add(trajectory.CPsiBlockTrajectory[0]);
        }

        return new C2BlockCpsiQScan(
            block, perturbedBondIndex, perturbationDelta, snapshotTime,
            qGrid.ToArray(), cTraj, psiTraj, cPsiTraj);
    }

    /// <summary>Index of the Q-grid entry minimizing CΨ_block — the "deepest crossing"
    /// in the sweep. If Q_EP-inheritance is real, this index should align with one of
    /// the orbit Q_EP values (Endpoint or mid-chain).</summary>
    public int ArgMinIndex
    {
        get
        {
            int idx = 0;
            for (int i = 1; i < CPsiBlockAtSnapshot.Count; i++)
                if (CPsiBlockAtSnapshot[i] < CPsiBlockAtSnapshot[idx]) idx = i;
            return idx;
        }
    }

    public double ArgMinQ => QGrid[ArgMinIndex];
    public double MinCPsi => CPsiBlockAtSnapshot[ArgMinIndex];

    public override string DisplayName =>
        $"C2BlockCpsiQScan (N={Block.N}, bond={PerturbedBondIndex}, Δ={PerturbationDelta:F3}, t_snap={SnapshotTime:F2})";

    public override string Summary =>
        $"{QGrid.Count} Q points; min CΨ={MinCPsi:F4} at Q={ArgMinQ:F3} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("N", Block.N);
            yield return InspectableNode.RealScalar("perturbed bond", PerturbedBondIndex);
            yield return InspectableNode.RealScalar("Δ (perturbation)", PerturbationDelta, "F4");
            yield return InspectableNode.RealScalar("snapshot t", SnapshotTime, "F4");
            yield return InspectableNode.RealScalar("argmin Q", ArgMinQ, "F4");
            yield return InspectableNode.RealScalar("min CΨ_block", MinCPsi, "F6");

            var rows = new List<IInspectable>();
            for (int i = 0; i < QGrid.Count; i++)
            {
                rows.Add(new InspectableNode(
                    $"Q={QGrid[i]:F3}",
                    summary: $"C={CBlockAtSnapshot[i]:F4}, Ψ={PsiBlockAtSnapshot[i]:F4}, CΨ={CPsiBlockAtSnapshot[i]:F4}"));
            }
            yield return InspectableNode.Group("sweep (Q, C, Ψ, CΨ)", rows.ToArray());
        }
    }
}
