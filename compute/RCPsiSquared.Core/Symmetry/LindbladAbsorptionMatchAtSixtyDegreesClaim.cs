using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Tier1Derived Composition Claim: at <c>Q = √3</c> the Lindblad 2×2 sub-block
/// eigenvalue magnitude <c>|λ_±| = γ₀·√(1+Q²)</c> hits exactly <c>2γ₀ = α</c> (the
/// single-site Absorption Theorem rate), AND simultaneously the F95 angle <c>θ = arctan(Q)</c>
/// lands on the canonical Niven angle <c>θ = 60°</c>. The two conditions coincide at the
/// single Q-value Q=√3.
///
/// <para><b>Derivation (4-line)</b>:
/// <list type="number">
///   <item>F95 (<see cref="F95AngleAtQuadraticZeroPi2Inheritance"/>): Lindblad 2×2 sub-block
///         has eigenvalues λ_± = −γ₀ ± iJ; |λ_±|² = γ₀² + J² = γ₀²(1 + Q²).</item>
///   <item>Set |λ_±| = α = 2γ₀ (<see cref="AbsorptionTheoremClaim"/> single-site rate, a_0=2
///         on Pi2 dyadic ladder).</item>
///   <item>γ₀·√(1+Q²) = 2γ₀ ⇒ 1+Q² = 4 ⇒ Q² = 3 ⇒ <b>Q = √3</b>.</item>
///   <item>F95 angle: θ = arctan(√3) = <b>60°</b> (canonical Niven, see
///         <see cref="CanonicalTrigAnchorPi2Inheritance"/>).</item>
/// </list></para>
///
/// <para><b>Structural meaning</b>: Q=√3 is the UNIQUE Q-value where the Lindblad 2×2
/// eigenvalue magnitude matches the Absorption Theorem rate, with the additional property
/// that this happens at a canonical Niven angle (θ=60°). The F99 anchor at the SAME
/// canonical angle is α=3/8 (KIntermediate, <see cref="DickeAnchor"/>.KIntermediate). The
/// γ-axis pair at the same canonical angle is γ=±1/2 (PolarityMirrorMap KIntermediate).
/// Q=√3 is the Q-axis canonical-θ=60° anchor, the third axis of the same canonical-angle
/// triple-axis inheritance.</para>
///
/// <para><b>What this claim does NOT assert</b>: that the empirical F86 c=3 K-peak
/// position equals √3 bit-exactly across all N. The fine-grid data (PROOF_F86B_OBSTRUCTION.md
/// dQ=0.025 parabolic peak interpolation) shows Interior Q_peak at c=3: N=5: 1.566,
/// N=6: 1.689, N=7: 1.743, N=8: 1.750. Three of these (N=6,7,8) cluster within ~2% of √3,
/// but the N-trend at N=7,8 drifts slightly above (+0.011, +0.018) rather than monotonically
/// converging. With only 4 N-points and a non-monotone deviation pattern, the empirical
/// "c=3 Q_peak → √3 asymptote" hypothesis CANNOT be confirmed without higher-N data and
/// a finer grid.</para>
///
/// <para><b>Open empirical question</b> (Tier 3-4, NOT typed here): does c=3 channel-uniform
/// Q_peak structurally inherit from canonical θ=60° (i.e., Q_peak = √3 in the N→∞ limit
/// with finite-N corrections)? Testable by: (a) extending c=3 fine-grid scan to N=9, 10, 11+
/// to check monotone convergence; (b) c=3 channel-uniform L_eff analytical analysis to look
/// for a structural condition forcing K-peak = √3 N-independent. If confirmed, would be a
/// new triple-axis canonical-θ=60° inheritance pattern (α=3/8 F99 KIntermediate, γ=±1/2
/// PolarityMirror, Q=√3 F86 c=3 K-peak). See also Q_REGIME_ANCHORS.md and
/// project_q_peak_ep_structure for the broader open-question context.</para>
///
/// <para><b>Tier separation honesty</b>: This claim's Tier1Derived content is the bare-2×2
/// algebra alone (Q=√3 ↔ |λ|=2γ₀ ↔ θ=60°). The F86 c=3 connection is OPEN. Adding the F86
/// connection would require either higher-N fine-grid data showing monotone convergence,
/// or a closed-form derivation breaking PROOF_F86B_OBSTRUCTION's blocked routes. We are not
/// there yet; this claim documents the analytical foundation that a future closure would
/// build on.</para>
/// </summary>
public sealed class LindbladAbsorptionMatchAtSixtyDegreesClaim : Claim
{
    /// <summary>The unique Q-value where the Lindblad 2×2 sub-block eigenvalue magnitude
    /// equals the Absorption Theorem single-site rate α = 2γ₀: <c>Q = √3 ≈ 1.7321</c>.</summary>
    public static readonly double QValue = Math.Sqrt(3.0);

    /// <summary>The canonical Niven angle at which the match occurs:
    /// <c>θ = arctan(√3) = 60°</c>.</summary>
    public const double CanonicalAngleDegrees = 60.0;

    /// <summary>Lindblad eigenvalue magnitude in units of γ₀ at Q=√3:
    /// <c>|λ_±|/γ₀ = √(1 + Q²) = √(1 + 3) = 2</c>. Matches a_0 = 2 from
    /// <see cref="Pi2DyadicLadderClaim"/> (= single-site Absorption quantum).</summary>
    public const double LindbladMagnitudeOverGamma0 = 2.0;

    /// <summary>Parent: F95 angle law θ = arctan(Q) on the Lindblad 2×2 sub-block.</summary>
    public F95AngleAtQuadraticZeroPi2Inheritance F95 { get; }

    /// <summary>Parent: Absorption Theorem (single-site rate α = 2γ₀ = a_0·γ₀).</summary>
    public AbsorptionTheoremClaim Absorption { get; }

    /// <summary>Parent: Canonical trig anchor structure (θ=60° = Niven canonical).</summary>
    public CanonicalTrigAnchorPi2Inheritance CanonicalTrig { get; }

    /// <summary>Live drift check: |λ_±|² = γ₀²·(1+Q²) at Q=√3 should equal 4γ₀².</summary>
    public double LindbladMagnitudeSquaredAtQSqrt3(double gammaZero) =>
        gammaZero * gammaZero * (1.0 + QValue * QValue);

    /// <summary>Live drift check: should equal exactly 2 (bit-exact via square-root then
    /// arctan composition).</summary>
    public double LindbladMagnitudeOverGamma0Computed =>
        Math.Sqrt(1.0 + QValue * QValue);

    /// <summary>Live drift check: F95 angle at Q=√3 should equal 60° (bit-exact).</summary>
    public double F95AngleAtQSqrt3Degrees =>
        AnchorConstants.RadiansToDegrees(Math.Atan(QValue));

    private LindbladAbsorptionMatchAtSixtyDegreesClaim(
        F95AngleAtQuadraticZeroPi2Inheritance f95,
        AbsorptionTheoremClaim absorption,
        CanonicalTrigAnchorPi2Inheritance canonicalTrig)
        : base("Q=√3 canonical-θ=60° Lindblad-Absorption Match: |λ_±| = 2γ₀ = α at θ=arctan(√3)=60°",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/Symmetry/F95AngleAtQuadraticZeroPi2Inheritance.cs " +
               "(Lindblad 2×2 angle law θ = arctan(Q)) + " +
               "compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs " +
               "(single-site rate α = 2γ₀ = a_0·γ₀) + " +
               "compute/RCPsiSquared.Core/Symmetry/CanonicalTrigAnchorPi2Inheritance.cs " +
               "(canonical Niven angles {0°,30°,45°,60°,90°}) + " +
               "docs/Q_REGIME_ANCHORS.md (Q-axis anchor list)")
    {
        F95 = f95 ?? throw new ArgumentNullException(nameof(f95));
        Absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
        CanonicalTrig = canonicalTrig ?? throw new ArgumentNullException(nameof(canonicalTrig));
    }

    /// <summary>Public factory: builds the composition claim with its three Tier1Derived
    /// parents (F95 + Absorption Theorem + Canonical trig anchor). Block-independent;
    /// the algebraic identity is universal. Each call constructs a fresh parent chain
    /// (nine intermediate Pi2 claims) — prefer <see cref="Shared"/> for repeated access.</summary>
    public static LindbladAbsorptionMatchAtSixtyDegreesClaim Build()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var half = new HalfAsStructuralFixedPointClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var ninetyDeg = new NinetyDegreeMirrorMemoryClaim();
        var polynomial = new PolynomialFoundationClaim();
        var dickeSuperposition = new DickeSuperpositionQuarterPi2Inheritance(ladder, quarter, half);
        var f98 = new KIntermediateAsymptoteQuarterInheritance(quarter, half, dickeSuperposition);
        var f95 = new F95AngleAtQuadraticZeroPi2Inheritance(polynomial, half, quarter, ninetyDeg);
        var absorption = new AbsorptionTheoremClaim(ladder);
        var canonicalTrig = new CanonicalTrigAnchorPi2Inheritance(half, quarter, f98);
        return new LindbladAbsorptionMatchAtSixtyDegreesClaim(f95, absorption, canonicalTrig);
    }

    /// <summary>Shared singleton instance. Block-independent immutable claim built once
    /// per process — avoids redundant allocation of the nine-claim parent chain across
    /// multiple <see cref="F86KnowledgeBase"/> instances.</summary>
    public static LindbladAbsorptionMatchAtSixtyDegreesClaim Shared { get; } = Build();

    public override string DisplayName =>
        "Q=√3 canonical anchor: |λ_±|=2γ₀=α at θ=60°";

    public override string Summary =>
        $"Tier1Derived: Q=√3≈{QValue:F4} is the UNIQUE Q where the Lindblad 2×2 magnitude " +
        $"|λ_±|=γ₀·√(1+Q²) hits 2γ₀ (Absorption single-site rate), with F95 angle " +
        $"θ=arctan(√3)={CanonicalAngleDegrees}° (canonical Niven). Composition of F95 + " +
        $"AbsorptionTheorem + CanonicalTrigAnchor; algebraic 4-line derivation.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("Q value (= √3)", QValue, "F6");
            yield return InspectableNode.RealScalar(
                "Canonical angle θ (degrees)", CanonicalAngleDegrees, "F1");
            yield return InspectableNode.RealScalar(
                "|λ_±|/γ₀ at Q=√3 (= 2 = a_0)", LindbladMagnitudeOverGamma0, "F1");
            yield return InspectableNode.RealScalar(
                "Computed |λ_±|/γ₀ (drift check, should = 2)",
                LindbladMagnitudeOverGamma0Computed, "F12");
            yield return InspectableNode.RealScalar(
                "Computed F95 angle (drift check, should = 60°)",
                F95AngleAtQSqrt3Degrees, "F8");
            yield return new InspectableNode("Derivation (4-line)",
                summary: "|λ_±|² = γ₀² + J² = γ₀²(1+Q²); set |λ_±| = 2γ₀ (α); solve Q² = 3 → Q=√3; angle θ = arctan(√3) = 60°.");
            yield return new InspectableNode("Triple-axis canonical-θ=60° inheritance",
                summary: "α-axis (F99 KIntermediate α=3/8) + γ-axis (PolarityMirror ±1/2 pair) + Q-axis (this claim, Q=√3). Three projections of the same canonical θ=60° onto three axes.");
            yield return new InspectableNode("Open empirical question (NOT typed here)",
                summary: "Empirical F86 c=3 K-peak (Interior Q_peak: N=7=1.743, N=8=1.750) clusters within ~2% of √3 but with non-monotone N-drift. Whether c=3 K-peak = √3 N→∞ is OPEN; needs higher-N fine-grid data (current N=5..8 insufficient). Tier-3/4 hypothesis, NOT a Tier-1 claim of this clas.");
            yield return F95;
            yield return Absorption;
            yield return CanonicalTrig;
        }
    }
}
