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

        return new JwDispersionStructure(N, modes, clusters);
    }

    private JwDispersionStructure(int n, XyJordanWignerModes modes, IReadOnlyList<JwDispersionCluster> clusters)
        : base("c=2 JW-mode dispersion-degenerate clusters of L_H on (n=1, n+1=2)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1' Direction (b'') (JW track) + " +
               "textbook XY OBC sine-mode dispersion")
    {
        N = n;
        Modes = modes;
        Clusters = clusters;
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
