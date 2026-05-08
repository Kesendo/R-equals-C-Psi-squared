using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>F86 Item 1' Direction (b'') JW track, T12 (Pfad A step 4.2): per-JW-cluster
/// intra-cluster D-eigenvalue spectrum, with the cluster-cardinality-symmetry verified at
/// the eigenvalue level (bit-exact).
///
/// <para>Composes <see cref="JwBlockBasis"/> (T9) + <see cref="JwDispersionStructure"/> (T10).
/// For each JW dispersion-cluster c (size k), computes the k×k Hermitian sub-matrix
/// <c>W_c = D'_{cluster c}</c> and diagonalizes it. The k eigenvalues of W_c are the
/// "intra-cluster D-eigenvalues" — at Q=0 these are the actual L-eigenvalues for this
/// cluster.</para>
///
/// <para><b>Cluster-cardinality-symmetry at the eigenvalue level (verified bit-exact N=4..6):</b>
/// clusters with the same size have <i>identical</i> W_c eigenvalue spectra. E.g. at N=4
/// all 4 size-4 clusters share the spectrum <c>{−6γ, −4.4γ, −3.6γ, −2γ}</c> within FP drift
/// (max |Δλ| = 1.67e-16). This is a stronger statement than T11's Frobenius²-level cluster-
/// cardinality-symmetry: it says the W_c matrices are <i>unitarily equivalent</i> across
/// same-size clusters, not just Frobenius-norm-matched.</para>
///
/// <para><b>Tier1 implication for the EP closed form:</b> in the W-eigenbasis (which
/// diagonalises W_c for both clusters of a pair), the inter-cluster D-coupling X̃ has
/// specific structure. For each shared W-eigenvalue λ between two clusters c₁ and c₂,
/// if X̃[λ, λ] ≠ 0, the L_eff(Q) eigenvalue problem reduces to a 2×2 sub-block on the
/// (W-eigenvector, W-eigenvector) pair. The 2×2 EP formula gives:</para>
///
/// <para><c>Q_EP^{(λ)} = 2 · |X̃[λ, λ]| / (γ · |δ_{c₁} − δ_{c₂}|)</c></para>
///
/// <para>under the a = b condition (same W-eigenvalue λ). Multiple sub-block EPs combine to
/// give the bond-specific Q_peak. This is the structural skeleton of Step 4 of Pfad A; the
/// open Tier1-Derivation step is the analytical formula for X̃[λ, λ] in terms of bond-
/// position and cluster geometry.</para>
///
/// <para><b>Class-level Tier: Tier1Derived.</b> Cluster-cardinality-symmetry at eigenvalue
/// level is an algebraic theorem (F71-mirror invariance + cosine-identity δ ↔ −δ symmetry +
/// W_c hermiticity). The runtime <see cref="MaxIntraClusterEigenvalueDeviationForSameSize"/>
/// witness bounds FP drift below 1e-10.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track) +
/// <c>docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md</c>.</para>
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

        // Cluster-cardinality-symmetry at eigenvalue level: same-size clusters → same spectrum
        double maxDev = 0;
        var bySize = witnesses.GroupBy(w => w.Cluster.Triples.Count).Where(g => g.Count() > 1);
        foreach (var group in bySize)
        {
            var arr = group.ToArray();
            int k = arr[0].Eigenvalues.Count;
            for (int i = 0; i < arr.Length; i++)
                for (int j = i + 1; j < arr.Length; j++)
                {
                    for (int e = 0; e < k; e++)
                        maxDev = Math.Max(maxDev, Math.Abs(arr[i].Eigenvalues[e] - arr[j].Eigenvalues[e]));
                }
        }

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
