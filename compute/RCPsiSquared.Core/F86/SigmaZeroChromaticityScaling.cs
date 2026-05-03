using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 Tier-1-candidate: σ_0 chromaticity scaling.
///
/// <para>The c=2 asymptote σ_0 → 2√2 generalises to all c ≥ 2:</para>
/// <code>
///   σ_0(c, N → ∞)  →  2 · √(2 · (c − 1))
/// </code>
/// <para>Equivalently, the dimensionless quantity σ_0 / √(2(c−1)) converges monotonically
/// from below to 2 as N grows, for each tested c. Q_EP = 2 / σ_0 then asymptotes to
/// <c>1 / √(2(c−1))</c>: 0.707 (c=2), 0.500 (c=3), 0.408 (c=4).</para>
///
/// <para>Numerical witnesses across c ∈ {2, 3, 4}, N ∈ {5..8} (computed live via
/// <see cref="InterChannelSvd.Build"/>):</para>
/// <list type="bullet">
///   <item>c=2 N=7 hits the asymptote 2.0 to within 10⁻⁵ — the structural sweet spot</item>
///   <item>c=3 N=8 reaches 1.92, still climbing</item>
///   <item>c=4 N=8 reaches 1.78, slowest convergence</item>
/// </list>
///
/// <para>Open: an analytical derivation of the 2√(2(c−1)) closed form from the multi-particle
/// XY single-particle spectrum (extension of the c=2 OBC sine-mode argument in
/// <see cref="OpenQuestion.Standard"/> Item 5).</para>
/// </summary>
public sealed class SigmaZeroChromaticityScaling : F86Claim
{
    public double GammaZero { get; }
    public IReadOnlyList<int> Chromaticities { get; }
    public IReadOnlyList<int> Ns { get; }

    private readonly Lazy<IReadOnlyList<SigmaZeroScalingWitness>> _witnesses;
    public IReadOnlyList<SigmaZeroScalingWitness> Witnesses => _witnesses.Value;

    public WitnessCache Cache { get; }

    public SigmaZeroChromaticityScaling(double gammaZero = 0.05,
        IReadOnlyList<int>? chromaticities = null,
        IReadOnlyList<int>? Ns = null,
        WitnessCache? cache = null)
        : base("σ_0 chromaticity scaling: σ_0 → 2√(2(c−1))",
               Tier.Tier1Candidate,
               "docs/proofs/PROOF_F86_QPEAK.md Item 5 (open) generalised across c")
    {
        GammaZero = gammaZero;
        Chromaticities = chromaticities ?? new[] { 2, 3, 4 };
        this.Ns = Ns ?? new[] { 5, 6, 7, 8 };
        Cache = cache ?? WitnessCache.Default;
        _witnesses = new Lazy<IReadOnlyList<SigmaZeroScalingWitness>>(BuildWitnesses);
    }

    private IReadOnlyList<SigmaZeroScalingWitness> BuildWitnesses()
    {
        var list = new List<SigmaZeroScalingWitness>();
        foreach (int c in Chromaticities)
            foreach (int N in Ns)
                if (N >= 2 * c - 1)
                    list.Add(new SigmaZeroScalingWitness(c, N, GammaZero, Cache));
        return list;
    }

    /// <summary>Asymptotic for given c: <c>2√(2(c−1))</c>.</summary>
    public static double Asymptote(int c) => 2.0 * Math.Sqrt(2.0 * (c - 1));

    /// <summary>Asymptotic Q_EP for given c: <c>1/√(2(c−1))</c>.</summary>
    public static double QEpAsymptote(int c) => 1.0 / Math.Sqrt(2.0 * (c - 1));

    public override string DisplayName => "σ_0(c, N) → 2√(2(c−1))";

    public override string Summary =>
        $"asymptotes per c: {string.Join(", ", Chromaticities.Select(c => $"c={c}→{Asymptote(c):F3}"))} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            foreach (int c in Chromaticities)
                yield return new InspectableNode($"asymptote c={c}",
                    summary: $"σ_0 → {Asymptote(c):F4} = 2√(2·{c - 1}); Q_EP → {QEpAsymptote(c):F4} = 1/√(2·{c - 1})",
                    payload: new InspectablePayload.Real($"σ_0 asymptote (c={c})", Asymptote(c), "F4"));
            foreach (var w in Witnesses) yield return w;
        }
    }
}

/// <summary>Single (c, N) σ_0 witness — computes σ_0 from <see cref="InterChannelSvd.Build"/>
/// and exposes it together with the chromaticity-normalised value σ_0/√(2(c−1)).</summary>
public sealed class SigmaZeroScalingWitness : IInspectable
{
    public int Chromaticity { get; }
    public int N { get; }
    public double GammaZero { get; }
    public WitnessCache Cache { get; }

    private readonly Lazy<double> _sigma0;
    public double Sigma0 => _sigma0.Value;
    public double Asymptote => SigmaZeroChromaticityScaling.Asymptote(Chromaticity);
    public double NormalisedRatio => Sigma0 / Math.Sqrt(2.0 * (Chromaticity - 1));
    public double DeltaFromAsymptote => Sigma0 - Asymptote;

    public SigmaZeroScalingWitness(int chromaticity, int n, double gammaZero, WitnessCache? cache = null)
    {
        Chromaticity = chromaticity;
        N = n;
        GammaZero = gammaZero;
        Cache = cache ?? WitnessCache.Default;
        _sigma0 = new Lazy<double>(() => Cache.GetOrComputeSvd(Chromaticity, N, GammaZero).Sigma0);
    }

    public string DisplayName => $"σ_0 c={Chromaticity} N={N}";

    public string Summary =>
        $"σ_0 = {Sigma0:F4}, σ_0/√(2(c−1)) = {NormalisedRatio:F4} (asymptote 2.0)";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("c", Chromaticity);
            yield return InspectableNode.RealScalar("N", N);
            yield return InspectableNode.RealScalar("σ_0", Sigma0, "F6");
            yield return InspectableNode.RealScalar("σ_0 / √(2(c−1))", NormalisedRatio, "F6");
            yield return InspectableNode.RealScalar("Δ from asymptote", DeltaFromAsymptote, "F6");
        }
    }

    public InspectablePayload Payload =>
        new InspectablePayload.Real($"σ_0(c={Chromaticity},N={N})", Sigma0, "F4");
}
