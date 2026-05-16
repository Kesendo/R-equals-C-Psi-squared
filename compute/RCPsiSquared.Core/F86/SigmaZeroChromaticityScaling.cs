using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 Tier-2-empirical: σ_0 chromaticity scaling. The 2√(2(c−1)) value is a
/// trajectory crossing point, not an asymptote — σ_0(c=2, N) keeps growing monotonically
/// past 2√2 (N=7: 2.8284 → N=11: 2.8561); Aitken extrapolation suggests true N→∞ limit
/// ~2.85–2.89. Same retraction lesson as Q_peak Endpoint = csc(π/(N+1)) (see
/// <see cref="RetractedClaim"/>): N=7 coincidence, not asymptote.
///
/// <para>Numerical witnesses across c ∈ {2, 3, 4}, N ∈ {5..11} (computed live via
/// <see cref="InterChannelSvd.Build"/>):</para>
/// <list type="bullet">
///   <item>c=2: N=5 → 1.955, N=6 → 1.981, <b>N=7 → 2.000</b> (crossing), N=8 → 2.008,
///         N=9 → 2.014, N=10 → 2.017, N=11 → 2.020</item>
///   <item>c=3 N=8 reaches 1.92, still climbing; sweet spot N_c* not yet reached at N=8</item>
///   <item>c=4 N=9 reaches 1.84, slowest convergence; sweet spot N_c* far above N=9</item>
/// </list>
///
/// <para>Bridge to F86 c=2 g_eff_Endpoint: σ_0(c=2, N) · √(3/8) approximates the empirical
/// g_eff_E (= 4.39382/Q_peak_Endpoint via the C2HwhmRatio composition) with Δ ≤ 0.01 for
/// N ≥ 6 and Δ = 0.005 at N=7 (sweet-spot crossing again); Δ = 0.064 at N=5. Higher-c
/// bridge √(3/8) untestable without an empirical g_eff table at c=3, c=4. Open: the true
/// σ_0(c, N → ∞) closed form (anchored in <see cref="F86OpenQuestions.Standard"/> Item 5).</para>
/// </summary>
public sealed class SigmaZeroChromaticityScaling : Claim
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
        : base("σ_0 chromaticity scaling: 2√(2(c−1)) trajectory crossing (asymptote retracted)",
               Tier.Tier2Empirical,
               "docs/proofs/PROOF_F86_QPEAK.md Item 3 (open: closed-form asymptote) + " +
               "docs/superpowers/syntheses/2026-05-07-sigma0-bridge-sweep.md (2√(2(c−1)) refuted as asymptote, preserved as N=7 crossing)")
    {
        GammaZero = gammaZero;
        Chromaticities = chromaticities ?? new[] { 2, 3, 4 };
        // Default Ns extended 2026-05-08 to include N=9: this is the witness that observes
        // c=2 ratio 2.014 > 2.0, falsifying the "monotone from below" reading of the asymptote.
        this.Ns = Ns ?? new[] { 5, 6, 7, 8, 9 };
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

    /// <summary>Trajectory-crossing value for given c: <c>2√(2(c−1))</c>. NOT an asymptote
    /// (refuted 2026-05-08; see class summary). At c=2 this equals 2√2 = 2.8284 and is
    /// crossed bit-exactly at N=7. Higher-c crossing N is unknown from current data.</summary>
    public static double Asymptote(int c) => 2.0 * Math.Sqrt(2.0 * (c - 1));

    /// <summary>Q_EP value at the trajectory-crossing N: <c>1/√(2(c−1))</c>. NOT an
    /// asymptotic Q_EP (refuted 2026-05-08). Preserved for symmetry with
    /// <see cref="Asymptote"/>; both are sweet-spot crossings, not limits.</summary>
    public static double QEpAsymptote(int c) => 1.0 / Math.Sqrt(2.0 * (c - 1));

    /// <summary><b>Sweet-spot identity:</b> σ_0(c=2, N=7) = 2√2 bit-exactly (verified to 9·10⁻¹⁶
    /// independently, reproducible across full block-L SVD and the framework primitive
    /// computation). This is a Tier-1-candidate <b>finite-N structural identity</b> at a
    /// specific (c=2, N=7) point — distinct from the refuted asymptotic claim
    /// σ_0 → 2√(2(c−1)). The c=3 trajectory does NOT cross its analogous c-asymptote 4.0
    /// within tested N≤10, suggesting the N=7 c=2 identity is c=2-specific, not a generic
    /// pattern. Tier-1 candidate awaiting analytical proof via OBC sine-mode + JW reduction
    /// (Item 5 / "Item 3" in <see cref="F86.F86OpenQuestions.Standard"/> /
    /// <c>docs/proofs/PROOF_F86_QPEAK.md</c>). See
    /// <see cref="VerifySweetSpotIdentity_C2N7"/> for the empirical witness check.</summary>
    public static readonly double SigmaZeroSweetSpotIdentity_C2N7 = 2.0 * Math.Sqrt(2.0);

    /// <summary>Empirical witness: σ_0(c=2, N=7) computed via <see cref="InterChannelSvd"/>
    /// matches <see cref="SigmaZeroSweetSpotIdentity_C2N7"/> bit-exactly (to within
    /// machine precision). Returns <c>true</c> iff |σ_0 − 2√2| &lt; <paramref name="tolerance"/>.</summary>
    public static bool VerifySweetSpotIdentity_C2N7(double tolerance = 1e-14)
    {
        var block = new CoherenceBlock(N: 7, n: 1, gammaZero: 0.05);
        var svd = InterChannelSvd.Build(block, hd1: 1, hd2: 3);
        return Math.Abs(svd.Sigma0 - SigmaZeroSweetSpotIdentity_C2N7) < tolerance;
    }

    public override string DisplayName => "σ_0(c, N) crosses 2√(2(c−1)) at sweet-spot N";

    public override string Summary =>
        $"crossings per c: {string.Join(", ", Chromaticities.Select(c => $"c={c}→{Asymptote(c):F3}"))} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            foreach (int c in Chromaticities)
                yield return new InspectableNode($"trajectory-crossing value c={c}",
                    summary: $"σ_0 = {Asymptote(c):F4} = 2√(2·{c - 1}) crossed at sweet-spot N (=7 at c=2); " +
                             $"Q_EP = {QEpAsymptote(c):F4} = 1/√(2·{c - 1}) at that N. NOT an asymptote.",
                    payload: new InspectablePayload.Real($"σ_0 crossing value (c={c})", Asymptote(c), "F4"));
            foreach (var w in Witnesses) yield return w;
        }
    }
}

/// <summary>Single (c, N) σ_0 witness; computes σ_0 from <see cref="InterChannelSvd.Build"/>
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
        $"σ_0 = {Sigma0:F4}, σ_0/√(2(c−1)) = {NormalisedRatio:F4} (crosses 2.0 at sweet-spot N; not the asymptote)";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("c", Chromaticity);
            yield return InspectableNode.RealScalar("N", N);
            yield return InspectableNode.RealScalar("σ_0", Sigma0, "F6");
            yield return InspectableNode.RealScalar("σ_0 / √(2(c−1))", NormalisedRatio, "F6");
            yield return InspectableNode.RealScalar("Δ from 2√(2(c−1)) crossing value", DeltaFromAsymptote, "F6");
        }
    }

    public InspectablePayload Payload =>
        new InspectablePayload.Real($"σ_0(c={Chromaticity},N={N})", Sigma0, "F4");
}
