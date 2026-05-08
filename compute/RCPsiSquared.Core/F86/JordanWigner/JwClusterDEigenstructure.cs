using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>Per-JW-cluster intra-cluster D-eigenvalue spectrum: for each cluster c of size k,
/// diagonalises the k×k Hermitian sub-matrix <c>W_c = D'_{cluster c}</c>. The k eigenvalues
/// are the L-eigenvalues for this cluster at Q=0.
///
/// <para>Cluster-cardinality-symmetry at the eigenvalue level: same-size clusters have
/// identical W_c spectra. Stronger than T11's Frobenius²-level symmetry — same-size W_c
/// matrices are unitarily equivalent. At N=4 γ=0.05 all 4 size-4 clusters share
/// {−6γ, −4.4γ, −3.6γ, −2γ}; verified via
/// <see cref="MaxIntraClusterEigenvalueDeviationForSameSize"/>
/// below <see cref="SymmetryTolerance"/> bit-exact.</para>
///
/// <para>Tier1Derived (F71-mirror invariance + cosine-identity δ ↔ −δ symmetry + W_c
/// hermiticity). Anchors: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'')
/// (JW track) + <c>docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md</c>.</para>
/// </summary>
public sealed class JwClusterDEigenstructure : Claim
{
    public const double SymmetryTolerance = 1e-10;

    public CoherenceBlock Block { get; }
    public JwBlockBasis JwBasis { get; }
    public JwDispersionStructure DispersionStructure { get; }

    /// <summary>Per-cluster eigenvalue spectrum of W_c = D'_{cluster c}. Sorted ascending
    /// per cluster.</summary>
    public IReadOnlyList<JwClusterEigenstructureWitness> Clusters { get; }

    /// <summary>Maximum |Δλ| between any two clusters of the same size (cluster-cardinality-
    /// symmetry at eigenvalue level). Algebraically zero (W_c matrices are unitarily
    /// equivalent across same-size clusters); the witness bounds FP drift.</summary>
    public double MaxIntraClusterEigenvalueDeviationForSameSize { get; }

    public static JwClusterDEigenstructure Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"JwClusterDEigenstructure applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var jw = JwBlockBasis.Build(block);
        var disp = JwDispersionStructure.Build(block.N);

        var DJw = jw.Uinv * block.Decomposition.D * jw.U;

        var tripleToAlpha = new Dictionary<(int, int, int), int>();
        for (int alpha = 0; alpha < jw.Triples.Count; alpha++)
        {
            var t = jw.Triples[alpha];
            tripleToAlpha[(t.K, t.K1, t.K2)] = alpha;
        }

        var witnesses = new JwClusterEigenstructureWitness[disp.Clusters.Count];
        for (int c = 0; c < disp.Clusters.Count; c++)
        {
            var cluster = disp.Clusters[c];
            int k = cluster.Triples.Count;
            int[] indices = cluster.Triples.Select(t => tripleToAlpha[(t.K, t.K1, t.K2)]).ToArray();

            var W = ComplexMatrix.Build.Dense(k, k);
            for (int i = 0; i < k; i++)
                for (int j = 0; j < k; j++)
                    W[i, j] = DJw[indices[i], indices[j]];

            var evd = W.Evd();
            var eigvals = evd.EigenValues.Select(z => z.Real).OrderBy(x => x).ToArray();
            witnesses[c] = new JwClusterEigenstructureWitness(cluster, eigvals);
        }

        double maxDev = 0;
        foreach (var group in witnesses.GroupBy(w => w.Cluster.Triples.Count).Where(g => g.Count() > 1))
        {
            var arr = group.ToArray();
            int k = arr[0].Eigenvalues.Count;
            for (int i = 0; i < arr.Length; i++)
                for (int j = i + 1; j < arr.Length; j++)
                    for (int e = 0; e < k; e++)
                        maxDev = Math.Max(maxDev, Math.Abs(arr[i].Eigenvalues[e] - arr[j].Eigenvalues[e]));
        }

        if (maxDev > SymmetryTolerance)
            throw new InvalidOperationException(
                $"Cluster-cardinality-symmetry violated: max same-size eigenvalue deviation = {maxDev:E2} " +
                $"exceeds tolerance {SymmetryTolerance:G3}");

        return new JwClusterDEigenstructure(block, jw, disp, witnesses, maxDev);
    }

    private JwClusterDEigenstructure(
        CoherenceBlock block,
        JwBlockBasis jwBasis,
        JwDispersionStructure dispersionStructure,
        IReadOnlyList<JwClusterEigenstructureWitness> witnesses,
        double maxDev)
        : base("c=2 per-JW-cluster D-eigenvalue spectrum (cluster-cardinality-symmetry at eigenvalue level)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1' Direction (b'') (JW track) + " +
               "docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md")
    {
        Block = block;
        JwBasis = jwBasis;
        DispersionStructure = dispersionStructure;
        Clusters = witnesses;
        MaxIntraClusterEigenvalueDeviationForSameSize = maxDev;
    }

    public override string DisplayName =>
        $"c=2 JW cluster D-eigenstructure (N={Block.N}, {Clusters.Count} clusters)";

    public override string Summary =>
        $"per-cluster W_c eigenvalues; cluster-cardinality-symmetry at eigenvalue level " +
        $"(max-Δλ same-size = {MaxIntraClusterEigenvalueDeviationForSameSize:G3}) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("ClusterCount", summary: Clusters.Count.ToString());
            yield return InspectableNode.RealScalar("MaxIntraClusterEigenvalueDeviationForSameSize",
                MaxIntraClusterEigenvalueDeviationForSameSize, "G3");
            yield return InspectableNode.Group("Clusters", Clusters.Cast<IInspectable>().ToArray());
        }
    }
}

/// <summary>One cluster's intra-cluster D-eigenvalue spectrum.</summary>
public sealed record JwClusterEigenstructureWitness(
    JwDispersionCluster Cluster,
    IReadOnlyList<double> Eigenvalues
) : IInspectable
{
    public string DisplayName => $"cluster δ={Cluster.Delta:F4} (size {Cluster.Triples.Count})";

    public string Summary =>
        $"eigenvalues = [{string.Join(", ", Eigenvalues.Select(e => e.ToString("F4")))}]";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return Cluster;
            yield return new InspectableNode("eigenvalues",
                summary: $"[{string.Join(", ", Eigenvalues.Select(e => e.ToString("F6")))}]");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
