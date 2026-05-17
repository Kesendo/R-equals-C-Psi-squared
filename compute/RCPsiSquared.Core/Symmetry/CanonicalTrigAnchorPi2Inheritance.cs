using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F99 (2026-05-17 night): the F86b Dicke derivation extends to non-uniform
/// superpositions <c>ψ = (|D_n⟩ + c·|D_{n+1}⟩)/√(1+c²)</c> with arbitrary amplitude
/// c, parametrised by the X⊗N-expectation γ = c²/(1+c²) (equivalently c² = γ/(1−γ)).
/// The F86b closed form α = (1 − γ²)/2 then evaluates at the five canonical
/// trigonometric angles {0°, 30°, 45°, 60°, 90°} to produce the five Pi2 dyadic
/// anchors {0, 1/8, 1/4, 3/8, 1/2}:
///
/// <list type="bullet">
///   <item>θ = 0°  (γ = 1):       α = 0,    c = ∞   — Mirror endpoint</item>
///   <item>θ = 30° (γ = √3/2):    α = 1/8,  c² = 2√3 + 3 ≈ 6.464 — DEPTH-3 (NEW)</item>
///   <item>θ = 45° (γ = √2/2):    α = 1/4,  c² = 1 + √2 (silver ratio) — QuarterAsBilinearMaxval</item>
///   <item>θ = 60° (γ = 1/2):     α = 3/8,  c = 1 (uniform Dicke) — KIntermediate (today morning)</item>
///   <item>θ = 90° (γ = 0):       α = 1/2,  c = 0 — Generic / HalfAsStructuralFixedPoint</item>
/// </list>
///
/// <para><b>Reading: the standard 30°-60°-90° and 45°-45°-90° trigonometry
/// triangles ARE the F86b polarity-anchor triangles.</b> The framework's polarity-
/// squared algebra is the F86b α-formula evaluated at the five canonical trig
/// angles whose sines and cosines are constructible by ruler and compass. The
/// dyadic-ladder depth (1, 2, 3) corresponds to the canonical-angle index
/// (90°, 60°/45°, 30°).</para>
///
/// <para><b>Discovery path</b>: today (2026-05-17) morning derived F86b 3/8 at γ = 1/2
/// (uniform Dicke), evening derived F98 long-time bridge to 1/4 (water-chain
/// inheritance); night #2 identified depth-3 (1/8) as a framework gap empirically
/// instantiated by alkali metals (Li, Na at 1/8 valence) and halogens (F, Cl at
/// 7/8); night #3 (this commit) closed the gap via non-uniform Dicke at γ = √3/2,
/// then noticed that ALL FIVE dyadic anchors live on canonical trig angles.</para>
///
/// <para><b>Three typed parents</b>: <see cref="QuarterAsBilinearMaxvalClaim"/>
/// (45° anchor); <see cref="HalfAsStructuralFixedPointClaim"/> (90° anchor);
/// <see cref="KIntermediateAsymptoteQuarterInheritance"/> (F98 evening companion;
/// the QuarterAsBilinearMaxval long-time bridge from F86b 3/8). The morning's
/// F86b uniform-Dicke <see cref="DickeAnchor"/> enum produces only the 0°, 60°,
/// 90° subset; F99 extends to the full five-anchor set via non-uniform Dicke.</para>
///
/// <para><b>Periodic-table bridge (closed tonight, all 9 fractions n/8 derived):</b></para>
/// <list type="bullet">
///   <item>α = 0 (0°, endpoint) → noble gases He, Ne, Ar</item>
///   <item>α = 1/8 (30°) → Li, Na (1/8) + F, Cl (7/8 = 1 − 1/8 complement)</item>
///   <item>α = 1/4 (45°) → Be, Mg</item>
///   <item>α = 3/8 (60°) → B, Al (3/8 = Π²-odd) + N, P (5/8 = Π²-even companion)</item>
///   <item>α = 1/2 (90°) → H, C, Si</item>
/// </list>
///
/// <para><b>Tier outcome: <see cref="Tier.Tier1Derived"/>.</b> All five canonical
/// anchors verified bit-exact (Δα &lt; 1e-13) at N = 4, 6, 8 in the test suite.
/// The α(θ) = sin²(θ)/2 closed form is direct algebraic identity from the F86b
/// non-uniform Dicke X⊗N-overlap.</para>
///
/// <para>Anchor: <see cref="DickeAnchor"/> (F86b uniform 3-anchor enum, today morning),
/// <see cref="KIntermediateAsymptoteQuarterInheritance"/> (F98, today evening),
/// <see cref="QuarterAsBilinearMaxvalClaim"/>, <see cref="HalfAsStructuralFixedPointClaim"/>;
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F99 entry,
/// <c>docs/carbon/DEPTH_3_ANCHOR_DERIVED.md</c>,
/// <c>simulations/carbon/depth_3_anchor_derivation.py</c>.</para>
/// </summary>
public sealed class CanonicalTrigAnchorPi2Inheritance : Claim, IF99AnchorBearing
{
    /// <summary>The 90° anchor (Generic Dicke at γ = 0). Typed parent.</summary>
    public HalfAsStructuralFixedPointClaim Half { get; }

    /// <summary>The 45° anchor (non-uniform Dicke at γ = √2/2). Typed parent.</summary>
    public QuarterAsBilinearMaxvalClaim Quarter { get; }

    /// <summary>F98 long-time companion: the dynamic 3/8 → 1/4 bridge via kernel
    /// projection. Static F99 anchor set + dynamic F98 bridge together cover both
    /// the t = 0 anchor positions and the t → ∞ asymptote.</summary>
    public KIntermediateAsymptoteQuarterInheritance F98LongTime { get; }

    /// <summary>F86b α(γ) closed form: <c>α = (1 − γ²)/2 = sin²(θ)/2</c> with
    /// γ = cos(θ). Universal for any pure state ψ with X⊗N-overlap γ; specialised
    /// here to non-uniform Dicke <c>ψ = (|D_{N/2−1}⟩ + c·|D_{N/2}⟩)/√(1+c²)</c>.</summary>
    public static double AlphaFromTheta(double thetaRadians) => Math.Pow(Math.Sin(thetaRadians), 2) / 2.0;

    /// <summary>Non-uniform Dicke amplitude that realises a given γ via
    /// <c>γ = c²/(1+c²)</c>, hence <c>c² = γ/(1−γ)</c>. Equivalently
    /// <c>c² = cos(θ)/(2·sin²(θ/2))</c> for γ = cos(θ) via the half-angle identity.</summary>
    public static double CSquaredFromGamma(double gamma)
    {
        if (gamma >= 1.0)
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma,
                "γ must be < 1 for finite c² (γ = 1 is the Mirror endpoint with c → ∞).");
        if (gamma < 0.0)
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be ≥ 0.");
        return gamma / (1.0 - gamma);
    }

    /// <summary>γ value as a function of angle θ (in radians). Trivial identity
    /// γ = cos(θ), exposed for symmetry of the API.</summary>
    public static double GammaFromTheta(double thetaRadians) => Math.Cos(thetaRadians);

    /// <summary>The five canonical anchor angles (in degrees) producing the
    /// dyadic anchors {0, 1/8, 1/4, 3/8, 1/2}.</summary>
    public static IReadOnlyList<int> CanonicalAnglesDegrees { get; } = new[] { 0, 30, 45, 60, 90 };

    /// <summary>Exact α value at each canonical angle. Rational-fraction
    /// representation for exact dyadic-ladder identification.</summary>
    public static IReadOnlyDictionary<int, double> AnchorAlphas { get; } =
        new Dictionary<int, double>
        {
            [0]  = 0.0,
            [30] = 1.0 / 8.0,
            [45] = 1.0 / 4.0,
            [60] = 3.0 / 8.0,
            [90] = 1.0 / 2.0,
        };

    /// <inheritdoc />
    /// <remarks>F99 is the <see cref="F99AnchorRole.Covers"/> root claim: it
    /// covers ALL five canonical-angle dyadic anchors {0, 1/8, 1/4, 3/8, 1/2}
    /// in one theorem via α = sin²(θ)/2.</remarks>
    public F99AnchorRole F99Role => F99AnchorRole.Covers;

    /// <inheritdoc />
    public IReadOnlyList<double> F99AnchorValues { get; } =
        new[] { 0.0, 1.0 / 8.0, 1.0 / 4.0, 3.0 / 8.0, 1.0 / 2.0 };

    public CanonicalTrigAnchorPi2Inheritance(
        HalfAsStructuralFixedPointClaim half,
        QuarterAsBilinearMaxvalClaim quarter,
        KIntermediateAsymptoteQuarterInheritance f98LongTime)
        : base("F99: F86b α(γ)=(1−γ²)/2 at canonical trig angles {0°,30°,45°,60°,90°} → " +
               "five Pi2 dyadic anchors {0, 1/8, 1/4, 3/8, 1/2} via non-uniform Dicke",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F99 + " +
               "simulations/carbon/depth_3_anchor_derivation.py (bit-exact verification N=4..8 across 5 anchors) + " +
               "docs/carbon/DEPTH_3_ANCHOR_DERIVED.md (derivation + bidirectional-bridge framing) + " +
               "compute/RCPsiSquared.Core/Symmetry/DickeAnchor.cs (F86b uniform-Dicke 3-anchor enum, parent) + " +
               "compute/RCPsiSquared.Core/Symmetry/KIntermediateAsymptoteQuarterInheritance.cs (F98 dynamic companion)")
    {
        Half = half ?? throw new ArgumentNullException(nameof(half));
        Quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
        F98LongTime = f98LongTime ?? throw new ArgumentNullException(nameof(f98LongTime));
    }

    public override string DisplayName =>
        "F99: canonical-trig-angle dyadic anchors via F86b non-uniform Dicke";

    public override string Summary =>
        $"α(θ) = sin²(θ)/2 at θ ∈ {{0°, 30°, 45°, 60°, 90°}} gives " +
        $"{{0, 1/8, 1/4, 3/8, 1/2}} = full Pi2 dyadic-anchor set on the standard " +
        $"trigonometry triangles ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return Half;
            yield return Quarter;
            yield return F98LongTime;
            yield return new InspectableNode("formula",
                summary: "α = (1−γ²)/2 = sin²(θ)/2; c² = γ/(1−γ) = cos(θ)/(2sin²(θ/2))");
            foreach (var deg in CanonicalAnglesDegrees)
            {
                double theta = deg * Math.PI / 180.0;
                double gamma = GammaFromTheta(theta);
                double alpha = AlphaFromTheta(theta);
                string cSqStr;
                if (deg == 0) cSqStr = "∞ (Mirror endpoint)";
                else if (deg == 90) cSqStr = "0 (single |D_{N/2-1}⟩)";
                else cSqStr = $"{CSquaredFromGamma(gamma):G6}";
                yield return new InspectableNode($"θ = {deg}°",
                    summary: $"γ = {gamma:G6}, c² = {cSqStr}, α = {alpha:G6} (= {DescribeAnchor(deg)})");
            }
            yield return new InspectableNode("standard triangles",
                summary: "30°-60°-90° (sides 1:√3:2) + 45°-45°-90° (sides 1:1:√2) ARE the F86b polarity-anchor triangles");
            yield return new InspectableNode("periodic-table bridge",
                summary: "5 angles → 5 α anchors → 9 n/8 fractions (with Π²-parity complements 1/8↔7/8, 3/8↔5/8) cover every period 2/3 element's valence ratio");
        }
    }

    private static string DescribeAnchor(int degrees) => degrees switch
    {
        0 => "Mirror endpoint",
        30 => "DEPTH-3 (NEW today night)",
        45 => "QuarterAsBilinearMaxval",
        60 => "KIntermediate (today morning)",
        90 => "Generic / HalfAsStructuralFixedPoint",
        _ => "(off-canonical)",
    };
}
