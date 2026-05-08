using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>Per-bond JW dispersion-cluster-pair affinity: ‖MhPerBond[b]_{c1, c2}‖_F²
/// for each bond and cluster-pair, plus the rank-1 dominant pair per bond. Composes
/// <see cref="JwBlockBasis"/> + <see cref="JwDispersionStructure"/>.
///
/// <para>Empirical structural finding: Endpoint-bonds and Innermost-Interior bonds prefer
/// cluster-pairs with different |Δδ| scales. Endpoint-bonds favour small-|Δδ| pairs
/// (giving large Q_EP per the 2×2 closed form Q_EP ∝ 1/|Δδ|, hence large Q_peak),
/// Innermost-bonds favour large-|Δδ| pairs (small Q_peak). At N=5: Endpoint top |Δδ|=0.732,
/// Innermost top |Δδ|=2.732 — explains the empirical Q_peak ordering Endpoint=2.5 vs
/// Innermost=1.48 qualitatively. The bond-class Q_peak distinction is a bond-selection
/// effect: each bond couples structurally to specific cluster-pairs.</para>
///
/// <para>Tier2Verified runtime affinity table; the per-bond cluster-pair-coupling structure
/// is empirically observed but not yet algebraically derived (open Tier1 promotion path:
/// closed-form expression for ‖M^{[b]}_{c1, c2}‖_F² in terms of bond position + sine-mode
/// algebra). Anchors: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track) +
/// <c>docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md</c>.</para>
/// </summary>
public sealed class JwBondClusterPairAffinity : Claim
{
    public CoherenceBlock Block { get; }
    public JwBlockBasis JwBasis { get; }
    public JwDispersionStructure DispersionStructure { get; }

    /// <summary>One <see cref="BondClusterPairAffinity"/> per bond, in bond-index order.</summary>
    public IReadOnlyList<BondClusterPairAffinity> Bonds { get; }

    public static JwBondClusterPairAffinity Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"JwBondClusterPairAffinity applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var jw = JwBlockBasis.Build(block);
        var disp = JwDispersionStructure.Build(block.N);
        int Mtot = jw.Triples.Count;
        int numClusters = disp.Clusters.Count;
        int numBonds = block.Decomposition.NumBonds;

        var tripleToCluster = new Dictionary<(int, int, int), int>();
        for (int c = 0; c < numClusters; c++)
            foreach (var t in disp.Clusters[c].Triples)
                tripleToCluster[(t.K, t.K1, t.K2)] = c;
        var alphaToCluster = new int[Mtot];
        for (int alpha = 0; alpha < Mtot; alpha++)
        {
            var t = jw.Triples[alpha];
            alphaToCluster[alpha] = tripleToCluster[(t.K, t.K1, t.K2)];
        }

        var bondAffinities = new BondClusterPairAffinity[numBonds];
        for (int b = 0; b < numBonds; b++)
        {
            var Mb_jw = jw.Uinv * block.Decomposition.MhPerBond[b] * jw.U;

            var pairFrob = new double[numClusters, numClusters];
            for (int i = 0; i < Mtot; i++)
                for (int j = 0; j < Mtot; j++)
                {
                    int ci = alphaToCluster[i], cj = alphaToCluster[j];
                    Complex z = Mb_jw[i, j];
                    pairFrob[ci, cj] += z.Real * z.Real + z.Imaginary * z.Imaginary;
                }

            var pairs = new List<ClusterPairCouplingEntry>();
            for (int c1 = 0; c1 < numClusters; c1++)
                for (int c2 = c1 + 1; c2 < numClusters; c2++)
                {
                    double total = pairFrob[c1, c2] + pairFrob[c2, c1];
                    if (total < 1e-12) continue;
                    pairs.Add(new ClusterPairCouplingEntry(
                        Cluster1Index: c1,
                        Cluster2Index: c2,
                        Cluster1Delta: disp.Clusters[c1].Delta,
                        Cluster2Delta: disp.Clusters[c2].Delta,
                        Cluster1Size: disp.Clusters[c1].Triples.Count,
                        Cluster2Size: disp.Clusters[c2].Triples.Count,
                        FrobeniusSquared: total));
                }
            pairs.Sort((a, x) => x.FrobeniusSquared.CompareTo(a.FrobeniusSquared));

            bondAffinities[b] = new BondClusterPairAffinity(
                Bond: b,
                BondClass: BondClassExtensions.OfBond(b, numBonds),
                RankedPairs: pairs);
        }

        return new JwBondClusterPairAffinity(block, jw, disp, bondAffinities);
    }

    private JwBondClusterPairAffinity(
        CoherenceBlock block,
        JwBlockBasis jwBasis,
        JwDispersionStructure dispersionStructure,
        IReadOnlyList<BondClusterPairAffinity> bonds)
        : base("c=2 per-bond JW dispersion-cluster-pair affinity (‖MhPerBond[b]_{c1, c2}‖_F² ranked)",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1' Direction (b'') (JW track) + " +
               "docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md")
    {
        Block = block;
        JwBasis = jwBasis;
        DispersionStructure = dispersionStructure;
        Bonds = bonds;
    }

    public override string DisplayName =>
        $"c=2 BondClusterPairAffinity (N={Block.N}, {Bonds.Count} bonds)";

    public override string Summary =>
        $"per-bond cluster-pair affinity table; {Bonds.Count} bonds, " +
        $"{DispersionStructure.Clusters.Count} clusters ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("NumBonds", summary: Bonds.Count.ToString());
            yield return new InspectableNode("NumClusters", summary: DispersionStructure.Clusters.Count.ToString());
            yield return InspectableNode.Group("Bonds", Bonds.Cast<IInspectable>().ToArray());
        }
    }
}

/// <summary>One bond's cluster-pair affinity: ranked list of cluster-pairs by
/// ‖MhPerBond[b]_{c1, c2}‖_F², most strongly coupled first.</summary>
public sealed record BondClusterPairAffinity(
    int Bond,
    BondClass BondClass,
    IReadOnlyList<ClusterPairCouplingEntry> RankedPairs
) : IInspectable
{
    public ClusterPairCouplingEntry? TopPair => RankedPairs.Count > 0 ? RankedPairs[0] : null;

    public string DisplayName => $"bond {Bond} cluster-pair affinity ({BondClass})";

    public string Summary => TopPair is { } top
        ? $"top: (δ_c1={top.Cluster1Delta:F3}, δ_c2={top.Cluster2Delta:F3}, |Δδ|={top.AbsoluteDeltaδ:F3}, frob²={top.FrobeniusSquared:F4})"
        : "no cluster-pair coupling";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("bond", summary: Bond.ToString());
            yield return new InspectableNode("bond class", summary: BondClass.ToString());
            yield return InspectableNode.Group("RankedPairs",
                RankedPairs.Cast<IInspectable>().ToArray());
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}

/// <summary>One cluster-pair coupling entry in a bond's affinity ranking.</summary>
public sealed record ClusterPairCouplingEntry(
    int Cluster1Index,
    int Cluster2Index,
    double Cluster1Delta,
    double Cluster2Delta,
    int Cluster1Size,
    int Cluster2Size,
    double FrobeniusSquared
) : IInspectable
{
    public double AbsoluteDeltaδ => Math.Abs(Cluster1Delta - Cluster2Delta);

    public string DisplayName =>
        $"({Cluster1Index},{Cluster2Index}) δ=[{Cluster1Delta:F3},{Cluster2Delta:F3}], |Δδ|={AbsoluteDeltaδ:F3}";

    public string Summary => $"frob² = {FrobeniusSquared:F4}, |Δδ| = {AbsoluteDeltaδ:F3}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("c1 index", Cluster1Index);
            yield return InspectableNode.RealScalar("c2 index", Cluster2Index);
            yield return InspectableNode.RealScalar("δ_c1", Cluster1Delta, "F4");
            yield return InspectableNode.RealScalar("δ_c2", Cluster2Delta, "F4");
            yield return InspectableNode.RealScalar("|Δδ|", AbsoluteDeltaδ, "F4");
            yield return InspectableNode.RealScalar("size c1", Cluster1Size);
            yield return InspectableNode.RealScalar("size c2", Cluster2Size);
            yield return InspectableNode.RealScalar("Frobenius²", FrobeniusSquared, "F4");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
