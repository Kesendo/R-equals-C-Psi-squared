using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>Dissipator D's Frobenius² decomposition in the JW dispersion-cluster basis.
/// Composes <see cref="JwBlockBasis"/> + <see cref="JwDispersionStructure"/> and splits
/// <c>D' = U⁻¹·D·U</c> into diagonal, within-cluster off-diagonal, and between-cluster
/// (cluster-mixing) Frobenius² pieces.
///
/// <para>The between-cluster piece dominates the within-cluster (~10-13% vs ~3-7% across
/// N=4..6) — D mixes JW dispersion-clusters, which is the EP-creating mechanism at finite Q.
/// Cluster-cardinality-symmetry: same-size clusters have bit-identical intra-cluster
/// D-Frobenius² (verified via <see cref="MaxIntraClusterFrobeniusDeviationForSameSize"/>),
/// rolling F71-mirror plus δ ↔ −δ cosine-identity symmetry into one empirical fact.</para>
///
/// <para>Tier2Verified runtime decomposition; the underlying JW-mixing + cluster-cardinality-
/// symmetry is Tier1 algebraic from sine-mode parity. Anchors:
/// <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track) +
/// <c>docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md</c>.</para>
/// </summary>
public sealed class JwDispersionDProjection : Claim
{
    public CoherenceBlock Block { get; }
    public JwBlockBasis JwBasis { get; }
    public JwDispersionStructure DispersionStructure { get; }

    /// <summary>‖diag(D')‖_F² = Σ_α |D'[α, α]|².</summary>
    public double DiagonalFrobeniusSquared { get; }

    /// <summary>Σ_α≠β, same-cluster |D'[α, β]|² (intra-cluster off-diagonal).</summary>
    public double WithinClusterFrobeniusSquared { get; }

    /// <summary>Σ_α≠β, different-cluster |D'[α, β]|² (cluster-mixing — EP-relevant).</summary>
    public double BetweenClusterFrobeniusSquared { get; }

    public double TotalFrobeniusSquared => DiagonalFrobeniusSquared + WithinClusterFrobeniusSquared + BetweenClusterFrobeniusSquared;

    /// <summary>Per-cluster intra-cluster diagonal/off-diagonal Frobenius² breakdown.</summary>
    public IReadOnlyList<JwClusterDProjection> ClusterProjections { get; }

    /// <summary>Maximum |Δ‖D'_cluster‖²| between any two clusters of the same size. Algebraically
    /// zero (F71-mirror + δ ↔ −δ symmetry); FP residual bounds drift below 1e-12.</summary>
    public double MaxIntraClusterFrobeniusDeviationForSameSize { get; }

    public static JwDispersionDProjection Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"JwDispersionDProjection applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var jw = JwBlockBasis.Build(block);
        var disp = JwDispersionStructure.Build(block.N);

        var DJw = jw.Uinv * block.Decomposition.D * jw.U;

        var tripleToCluster = new Dictionary<(int, int, int), int>();
        for (int c = 0; c < disp.Clusters.Count; c++)
            foreach (var t in disp.Clusters[c].Triples)
                tripleToCluster[(t.K, t.K1, t.K2)] = c;

        int Mtot = jw.Triples.Count;
        var alphaToCluster = new int[Mtot];
        for (int alpha = 0; alpha < Mtot; alpha++)
        {
            var t = jw.Triples[alpha];
            alphaToCluster[alpha] = tripleToCluster[(t.K, t.K1, t.K2)];
        }

        double diag = 0, within = 0, between = 0;
        var perClusterDiag = new double[disp.Clusters.Count];
        var perClusterOff = new double[disp.Clusters.Count];
        for (int i = 0; i < Mtot; i++)
            for (int j = 0; j < Mtot; j++)
            {
                Complex z = DJw[i, j];
                double mag2 = z.Real * z.Real + z.Imaginary * z.Imaginary;
                if (i == j)
                {
                    diag += mag2;
                    perClusterDiag[alphaToCluster[i]] += mag2;
                }
                else if (alphaToCluster[i] == alphaToCluster[j])
                {
                    within += mag2;
                    perClusterOff[alphaToCluster[i]] += mag2;
                }
                else
                {
                    between += mag2;
                }
            }

        var clusterProjections = new JwClusterDProjection[disp.Clusters.Count];
        for (int c = 0; c < disp.Clusters.Count; c++)
        {
            clusterProjections[c] = new JwClusterDProjection(
                Cluster: disp.Clusters[c],
                DiagonalFrobeniusSquared: perClusterDiag[c],
                OffDiagonalFrobeniusSquared: perClusterOff[c]);
        }

        double maxDev = 0;
        foreach (var group in clusterProjections.GroupBy(c => c.Cluster.Triples.Count).Where(g => g.Count() > 1))
        {
            var arr = group.ToArray();
            for (int i = 0; i < arr.Length; i++)
                for (int k = i + 1; k < arr.Length; k++)
                {
                    maxDev = Math.Max(maxDev, Math.Abs(arr[i].DiagonalFrobeniusSquared - arr[k].DiagonalFrobeniusSquared));
                    maxDev = Math.Max(maxDev, Math.Abs(arr[i].OffDiagonalFrobeniusSquared - arr[k].OffDiagonalFrobeniusSquared));
                }
        }

        return new JwDispersionDProjection(block, jw, disp, diag, within, between, clusterProjections, maxDev);
    }

    private JwDispersionDProjection(
        CoherenceBlock block,
        JwBlockBasis jwBasis,
        JwDispersionStructure dispersionStructure,
        double diag,
        double within,
        double between,
        IReadOnlyList<JwClusterDProjection> clusterProjections,
        double maxDevSameSize)
        : base("c=2 dissipator D' = U⁻¹·D·U Frobenius² decomposition by JW dispersion clusters",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1' Direction (b'') (JW track) + " +
               "docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md")
    {
        Block = block;
        JwBasis = jwBasis;
        DispersionStructure = dispersionStructure;
        DiagonalFrobeniusSquared = diag;
        WithinClusterFrobeniusSquared = within;
        BetweenClusterFrobeniusSquared = between;
        ClusterProjections = clusterProjections;
        MaxIntraClusterFrobeniusDeviationForSameSize = maxDevSameSize;
    }

    public override string DisplayName =>
        $"c=2 D' Frobenius² decomposition (N={Block.N}, {DispersionStructure.Clusters.Count} clusters)";

    public override string Summary =>
        $"diag = {DiagonalFrobeniusSquared:F4} ({100.0 * DiagonalFrobeniusSquared / TotalFrobeniusSquared:F1}%), " +
        $"within = {WithinClusterFrobeniusSquared:F4} ({100.0 * WithinClusterFrobeniusSquared / TotalFrobeniusSquared:F1}%), " +
        $"between = {BetweenClusterFrobeniusSquared:F4} ({100.0 * BetweenClusterFrobeniusSquared / TotalFrobeniusSquared:F1}%) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return InspectableNode.RealScalar("Total ‖D'‖_F²", TotalFrobeniusSquared, "F4");
            yield return InspectableNode.RealScalar("Diagonal ‖D'‖_F²", DiagonalFrobeniusSquared, "F4");
            yield return InspectableNode.RealScalar("Within-cluster ‖D'‖_F²", WithinClusterFrobeniusSquared, "F4");
            yield return InspectableNode.RealScalar("Between-cluster ‖D'‖_F²", BetweenClusterFrobeniusSquared, "F4");
            yield return InspectableNode.RealScalar("MaxIntraClusterFrobeniusDeviationForSameSize",
                MaxIntraClusterFrobeniusDeviationForSameSize, "G3");
            yield return InspectableNode.Group("ClusterProjections",
                ClusterProjections.Cast<IInspectable>().ToArray());
        }
    }
}

/// <summary>One JW-dispersion-cluster's D-Frobenius² breakdown.</summary>
public sealed record JwClusterDProjection(
    JwDispersionCluster Cluster,
    double DiagonalFrobeniusSquared,
    double OffDiagonalFrobeniusSquared
) : IInspectable
{
    public double TotalFrobeniusSquared => DiagonalFrobeniusSquared + OffDiagonalFrobeniusSquared;
    public double OffOverDiagRatio =>
        DiagonalFrobeniusSquared > 0 ? OffDiagonalFrobeniusSquared / DiagonalFrobeniusSquared : 0.0;

    public string DisplayName => $"cluster δ={Cluster.Delta:F4} (size {Cluster.Triples.Count})";

    public string Summary =>
        $"diag={DiagonalFrobeniusSquared:F4}, off={OffDiagonalFrobeniusSquared:F4}, ratio={OffOverDiagRatio:F4}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return Cluster;
            yield return InspectableNode.RealScalar("DiagonalFrobeniusSquared", DiagonalFrobeniusSquared, "F4");
            yield return InspectableNode.RealScalar("OffDiagonalFrobeniusSquared", OffDiagonalFrobeniusSquared, "F4");
            yield return InspectableNode.RealScalar("Off/Diag ratio", OffOverDiagRatio, "F4");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
