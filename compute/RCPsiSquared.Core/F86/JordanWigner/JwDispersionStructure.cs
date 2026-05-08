using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>F86 Item 1' Direction (b'') JW track, T10 (Pfad A step 2): enumerate the
/// dispersion-degenerate sub-spaces of L_H on the c=2 (n=1, n+1=2) coherence block.
///
/// <para>For each JW triple α = (k, k₁, k₂) the unitary part of L(Q) has eigenvalue
/// δ_α = ε_k − ε_{k₁} − ε_{k₂} (at unit J), where ε_k = 2·cos(π·k/(N+1)) is the OBC
/// sine-mode dispersion (T1 <see cref="XyJordanWignerModes"/>). Triples with the same δ
/// form dispersion-degenerate clusters that, together with the dissipator's off-diagonal
/// JW-mixing, host the EP-creating level repulsion at finite Q.</para>
///
/// <para>Empirical structure (verified bit-exact N=4..8):</para>
/// <list type="bullet">
///   <item>67%–84% of triples sit in degenerate clusters (size ≥ 2)</item>
///   <item>Max cluster size grows linearly with N: 4 (N=4) → 10 (N=8)</item>
///   <item>Cluster structure is symmetric under δ ↔ −δ via the cosine identity
///   ε_k = −ε_{N+1−k}</item>
/// </list>
///
/// <para>Three algebraic patterns generate degenerate triples at δ = −ε_m:</para>
/// <list type="bullet">
///   <item><b>Pattern A (trivial inversion):</b> (k, m, k) for any k ≠ m gives δ = −ε_m
///   independent of k</item>
///   <item><b>Pattern B (complement pair):</b> (N+1−m, k₁, k₂) with k₁+k₂ = N+1 gives
///   δ = ε_{N+1−m} − 0 = −ε_m via ε_{k₁} + ε_{k₂} = 0</item>
///   <item><b>Pattern C (residual):</b> rare specific solutions ε_{k₁} + ε_{k₂} = 2·ε_m</item>
/// </list>
///
/// <para>The primitive exposes the empirical clustering as data; algebraic-pattern
/// interpretation is an open Tier1-Derivation step (next: bond-specific D-coupling
/// matrix in cluster basis, then closed-form K-resonance). This primitive is the
/// foundation for that work.</para>
///
/// <para><b>Class-level Tier: Tier1Derived.</b> Combinatorial enumeration + δ-value
/// grouping is exact algebra; the runtime <see cref="Tolerance"/> bounds FP drift.
/// Total triple count = N · C(N, 2) (verified at construction).</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track) +
/// textbook XY OBC sine-mode dispersion.</para>
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
