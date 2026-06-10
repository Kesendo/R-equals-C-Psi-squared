using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>Dispersion-degenerate clusters of L_H on the c=2 (n=1, n+1=2) coherence block:
/// JW triples (k, k₁, k₂) grouped by δ = ε_k − ε_{k₁} − ε_{k₂} where ε_k is the OBC
/// sine-mode dispersion (T1 <see cref="XyJordanWignerModes"/>).
///
/// <para>Tier1Derived (combinatorial + δ-grouping is exact algebra; <see cref="Tolerance"/>
/// bounds FP drift). Total triple count = N · C(N, 2) verified at construction.
/// Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track).</para>
///
/// <para>F90 status (2026-05-11): the F86 c=2 ↔ F89 path-(N−1) bridge identity
/// achieves numerical Tier-1 for Direction (b'') via per-bond Hellmann-Feynman
/// (bit-exact 20/22 bonds at N=5..8). The JW-track primitives in this file remain
/// active as the alternative analytical route toward the closed-form HWHM_left/Q_peak
/// constants; the per-bond numerical answer itself is no longer the open piece.
/// See <c>docs/proofs/PROOF_F90_F86C2_BRIDGE.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F90F86C2BridgeIdentity.cs</c>.</para>
///
/// <para><b>Chiral ± pairing named 2026-06-10:</b> the global k-reflection k → N+1−k maps
/// the triple (k, k₁, k₂) to (N+1−k, N+1−k₂, N+1−k₁) and flips δ → −δ because
/// ε_{N+1−k} = −ε_k (the <see cref="RCPsiSquared.Core.Symmetry.ChiralKClaim"/> BDI spectrum
/// inversion carried by <see cref="XyJordanWignerModes"/>, this claim's typed ancestor), so
/// every cluster has a mirror cluster of equal size at −δ; verified at construction via
/// <see cref="ChiralPartnerIndex"/> + <see cref="ChiralClusterPairingResidual"/>.</para>
/// </summary>
public sealed class JwDispersionStructure : Claim
{
    /// <summary>Tolerance for grouping δ values into clusters. Pure FP drift; clusters are
    /// algebraically exact.</summary>
    public const double Tolerance = 1e-10;

    public int N { get; }
    public XyJordanWignerModes Modes { get; }

    /// <summary>One <see cref="JwDispersionCluster"/> per distinct δ value, sorted ascending.</summary>
    public IReadOnlyList<JwDispersionCluster> Clusters { get; }

    /// <summary>Index of each cluster's chiral mirror cluster under the global k-reflection
    /// k → N+1−k (self-index for the δ ≈ 0 self-paired clusters). Equal sizes and
    /// δ-negation are verified at construction; the pairing is an involution.</summary>
    public IReadOnlyList<int> ChiralPartnerIndex { get; }

    /// <summary>max over clusters of |δ + δ_mirror|: pure FP drift (≲ 1.4e-15 observed at
    /// N ≤ 10); the ± pairing is algebraically exact.</summary>
    public double ChiralClusterPairingResidual { get; }

    public int SelfPairedClusterCount => ChiralPartnerIndex.Where((p, i) => p == i).Count();

    public int TotalTriples => Clusters.Sum(c => c.Triples.Count);
    public int LargestClusterSize => Clusters.Max(c => c.Triples.Count);
    public int TriplesInDegenerateClusters => Clusters.Where(c => c.Triples.Count > 1).Sum(c => c.Triples.Count);
    public int DegenerateClusterCount => Clusters.Count(c => c.Triples.Count > 1);
    public int SingletonCount => Clusters.Count(c => c.Triples.Count == 1);

    public static JwDispersionStructure Build(int N)
    {
        if (N < 2)
            throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 2 for c=2 (n=1, n+1=2) block; got {N}.");

        var modes = XyJordanWignerModes.Build(N, J: 1.0);

        var entries = new List<(JwTriple triple, double δ)>(N * N * (N - 1) / 2);
        for (int k = 1; k <= N; k++)
            for (int k1 = 1; k1 <= N; k1++)
                for (int k2 = k1 + 1; k2 <= N; k2++)
                {
                    double δ = modes.Dispersion[k - 1] - modes.Dispersion[k1 - 1] - modes.Dispersion[k2 - 1];
                    entries.Add((new JwTriple(k, k1, k2), δ));
                }

        int expectedCount = N * N * (N - 1) / 2;
        if (entries.Count != expectedCount)
            throw new InvalidOperationException(
                $"Triple-count mismatch: enumerated {entries.Count}, expected {expectedCount} for N={N}.");

        entries.Sort((a, b) => a.δ.CompareTo(b.δ));

        var clusters = new List<JwDispersionCluster>();
        var current = new List<JwTriple>();
        double currentδ = double.NaN;
        foreach (var (triple, δ) in entries)
        {
            if (current.Count == 0 || Math.Abs(δ - currentδ) > Tolerance)
            {
                if (current.Count > 0) clusters.Add(new JwDispersionCluster(currentδ, current));
                current = new List<JwTriple> { triple };
                currentδ = δ;
            }
            else
            {
                current.Add(triple);
            }
        }
        if (current.Count > 0) clusters.Add(new JwDispersionCluster(currentδ, current));

        // Chiral ± pairing witness (named 2026-06-10): k → N+1−k maps (k, k₁, k₂) to
        // (N+1−k, N+1−k₂, N+1−k₁) and flips δ → −δ (ε_{N+1−k} = −ε_k, ChiralKClaim's BDI
        // spectrum inversion), so every cluster mirrors onto exactly one cluster of equal
        // size at −δ. Verified structurally here; FP drift only (clusters are exact algebra).
        var tripleToCluster = new Dictionary<(int, int, int), int>();
        for (int c = 0; c < clusters.Count; c++)
            foreach (var t in clusters[c].Triples)
                tripleToCluster[(t.K, t.K1, t.K2)] = c;

        var chiralPartner = new int[clusters.Count];
        double chiralResidual = 0.0;
        for (int c = 0; c < clusters.Count; c++)
        {
            int partner = -1;
            foreach (var t in clusters[c].Triples)
            {
                int pc = tripleToCluster[(N + 1 - t.K, N + 1 - t.K2, N + 1 - t.K1)];
                if (partner < 0) partner = pc;
                else if (partner != pc)
                    throw new InvalidOperationException(
                        $"Chiral pairing witness failed at N={N}: cluster {c} (δ={clusters[c].Delta:E6}) " +
                        $"reflects into two clusters ({partner} and {pc}).");
            }
            chiralPartner[c] = partner;
            if (clusters[partner].Triples.Count != clusters[c].Triples.Count)
                throw new InvalidOperationException(
                    $"Chiral pairing witness failed at N={N}: cluster {c} (size {clusters[c].Triples.Count}) " +
                    $"mirrors cluster {partner} (size {clusters[partner].Triples.Count}); sizes must match.");
            double residual = Math.Abs(clusters[partner].Delta + clusters[c].Delta);
            if (residual > 2.0 * Tolerance)
                throw new InvalidOperationException(
                    $"Chiral pairing witness failed at N={N}: |δ_{c} + δ_{partner}| = {residual:E3} " +
                    $"exceeds 2·Tolerance; the cluster table must be antisymmetric under k-reflection.");
            chiralResidual = Math.Max(chiralResidual, residual);
        }

        return new JwDispersionStructure(N, modes, clusters, chiralPartner, chiralResidual);
    }

    private JwDispersionStructure(int n, XyJordanWignerModes modes, IReadOnlyList<JwDispersionCluster> clusters,
        IReadOnlyList<int> chiralPartnerIndex, double chiralClusterPairingResidual)
        : base("c=2 JW-mode dispersion-degenerate clusters of L_H on (n=1, n+1=2)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1' Direction (b'') (JW track) + " +
               "textbook XY OBC sine-mode dispersion")
    {
        N = n;
        Modes = modes;
        Clusters = clusters;
        ChiralPartnerIndex = chiralPartnerIndex;
        ChiralClusterPairingResidual = chiralClusterPairingResidual;
    }

    public override string DisplayName =>
        $"c=2 JW dispersion structure (N={N}, {Clusters.Count} δ-clusters, max size {LargestClusterSize})";

    public override string Summary =>
        $"{TotalTriples} triples → {Clusters.Count} clusters, " +
        $"{TriplesInDegenerateClusters}/{TotalTriples} ({100.0 * TriplesInDegenerateClusters / TotalTriples:F1}%) " +
        $"in size-≥2 clusters, max size {LargestClusterSize} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("N", N);
            yield return InspectableNode.RealScalar("TotalTriples", TotalTriples);
            yield return InspectableNode.RealScalar("Cluster count", Clusters.Count);
            yield return InspectableNode.RealScalar("LargestClusterSize", LargestClusterSize);
            yield return InspectableNode.RealScalar("TriplesInDegenerateClusters", TriplesInDegenerateClusters);
            yield return InspectableNode.RealScalar("SingletonCount", SingletonCount);
            yield return new InspectableNode("chiral ± pairing (ChiralKClaim)",
                summary: "global k-reflection k → N+1−k maps (k, k₁, k₂) → (N+1−k, N+1−k₂, N+1−k₁) " +
                         "and flips δ → −δ (ε_{N+1−k} = −ε_k, ChiralKClaim's BDI spectrum inversion): " +
                         $"all {Clusters.Count} clusters pair off ± with equal sizes " +
                         $"({SelfPairedClusterCount} self-paired at δ ≈ 0), " +
                         $"residual {ChiralClusterPairingResidual:E2}");
            yield return InspectableNode.Group("Clusters", Clusters.Cast<IInspectable>().ToArray());
        }
    }
}

/// <summary>One dispersion-degenerate cluster: all JW triples (k, k₁, k₂) sharing the same
/// L_H eigenvalue δ = ε_k − ε_{k₁} − ε_{k₂} (within <see cref="JwDispersionStructure.Tolerance"/>).
/// </summary>
public sealed record JwDispersionCluster(
    double Delta,
    IReadOnlyList<JwTriple> Triples
) : IInspectable
{
    public string DisplayName => $"δ = {Delta:F4} (size {Triples.Count})";

    public string Summary =>
        $"δ = {Delta:F4}, {Triples.Count} triple(s): " +
        string.Join(", ", Triples.Select(t => $"({t.K},{t.K1},{t.K2})"));

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("δ", Delta, "F4");
            yield return InspectableNode.RealScalar("size", Triples.Count);
            foreach (var t in Triples) yield return t;
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
