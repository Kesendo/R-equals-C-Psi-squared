using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>F86 Item 1' Direction (b'') JW track, T11 (Pfad A step 3): the dissipator D's
/// Frobenius² decomposition in the JW dispersion-cluster basis.
///
/// <para>Composes <see cref="JwBlockBasis"/> (T9) and <see cref="JwDispersionStructure"/> (T10).
/// Computes <c>D' = U⁻¹·D·U</c> and decomposes its Frobenius² into three orthogonal pieces:
/// diagonal (within each JW basis vector), within-cluster off-diagonal (intra-dispersion-
/// degenerate-cluster), and between-cluster off-diagonal (cluster-mixing).</para>
///
/// <para>Empirical structure (verified bit-exact N=4..6):</para>
/// <list type="bullet">
///   <item>~80-86% diagonal</item>
///   <item>~3-7% within-cluster (intra-cluster off-diagonal)</item>
///   <item>~10-13% between-cluster (cluster-mixing — EP-relevant)</item>
/// </list>
///
/// <para><b>EP-relevance:</b> the between-cluster Frobenius² dominates the within-cluster
/// — D mixes JW dispersion-clusters, NOT just intra-cluster. The EP-bildung at finite Q
/// comes from this cluster-mixing: two eigenvalues from different δ-clusters (different
/// MhTotal-eigenvalues) collide as Q grows because D couples them. The dominant
/// between-cluster D-coupling pair determines the leading EP per bond — the foundation for
/// the closed-form 2-block effective Hamiltonian (Pfad A step 4).</para>
///
/// <para><b>Cluster-cardinality-symmetry (verified bit-exact N=4..6):</b> clusters with the
/// same size have <i>bit-identical</i> intra-cluster D-Frobenius² (diagonal and off-diagonal).
/// E.g. at N=6, all 6 size-7 clusters share diag=0.3729 and off=0.0269. This is F71-mirror
/// invariance plus the δ ↔ −δ cosine-identity symmetry rolled into a single empirical fact.
/// Verified per <see cref="MaxIntraClusterFrobeniusDeviationForSameSize"/>.</para>
///
/// <para><b>Class-level Tier: Tier2Verified.</b> Numerical decomposition; the underlying
/// algebraic structure (JW-mixing of D + cluster-cardinality-symmetry) is Tier1 algebraic
/// content from sine-mode parity, but this primitive ships the runtime Frobenius witnesses
/// at machine precision.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track) +
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

        // Cluster-cardinality-symmetry witness: clusters with the same size should have
        // bit-identical intra-cluster diag and off Frobenius².
        double maxDev = 0;
        var bySize = clusterProjections.GroupBy(c => c.Cluster.Triples.Count).Where(g => g.Count() > 1);
        foreach (var group in bySize)
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
