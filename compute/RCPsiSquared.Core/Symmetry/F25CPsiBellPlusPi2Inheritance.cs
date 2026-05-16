using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F25 closed form (Tier 1, proven monotonicity in PROOF_MONOTONICITY_CPSI):
///
/// <code>
///   CΨ(t) = f · (1 + f²) / 6,        f = e^{−4γt}
///
///   dCΨ/dt = −2γ · f · (1 + 3f²) / 3
///
///   Crossing at f* = 0.8612  (from f*(1 + f*²) = 3/2)
///   K = γ · t_cross = 0.0374
/// </code>
///
/// <para>F25 is the mother claim of <see cref="F57DwellTimeQuarterPi2Inheritance"/>:
/// F57's Bell+ K_dwell prefactor 1.080088 = 2/1.851701 derives directly from
/// F25's <c>|dCΨ/dt|</c> at the crossing. F57 → F25 is the typed
/// mother-claim inheritance edge (parallel pattern to F77 → F75).</para>
///
/// <para>Three Pi2-Foundation anchors:</para>
///
/// <list type="bullet">
///   <item><b>DecayRateCoefficient = 4</b>: in <c>f = e^{−4γt}</c>. Live from
///         <see cref="Pi2DyadicLadderClaim.Term"/>(−1) = <c>a_{−1}</c>.
///         Same anchor as F76 t-decay (e^{−4γ₀t}), F61/F63 4-block per parity,
///         F66 (multiplicity N+1 with a_{−1}), F77 correction denominator.</item>
///   <item><b>Coefficient2 = 2 = a_0</b>: in <c>dCΨ/dt = −2γ f(1+3f²)/3</c>.
///         <see cref="Pi2DyadicLadderClaim.Term"/>(0) = polynomial root d.
///         Same anchor as F1, F66, F86 Q_EP, F70 pair-bound, F81/F82.</item>
///   <item><b>CrossingThreshold = 1/4 = a_3</b>: F25 monotone strictly
///         decreasing crosses 1/4 once (Bell+ crosses fold). Same anchor as
///         F57, Dicke, F60 fold, F62 fold.</item>
/// </list>
///
/// <para>State-specific (NOT Pi2-anchored):</para>
///
/// <list type="bullet">
///   <item>Bell+ crossing <c>f* = 0.8612</c> from <c>f*(1+f*²) = 3/2</c>
///         (cubic root in f).</item>
///   <item>Bell+ <c>K = γ · t_cross = 0.0374</c>.</item>
///   <item>Bell+ <c>|dCΨ/dt|_{t_cross} = 1.851701</c>; F57's prefactor
///         <c>1.080088 = 2 / 1.851701</c>.</item>
///   <item>Denominator "6" in CΨ(t) and "3" in dCΨ/dt: combinatorial
///         normalization (Bell+ specific).</item>
/// </list>
///
/// <para>F25 + F57 mother-corollary chain:</para>
///
/// <code>
///   F25: CΨ(t) = f(1+f²)/6 closed form
///    ↓ "F57's Bell+ prefactor 1.080088 = 2 / |dCΨ/dt|_{t_cross}"
///   F57: t_dwell = 2δ / |dCΨ/dt|_{t_cross}
/// </code>
///
/// <para>Tier1Derived: F25 is Tier 1 proven (PROOF_MONOTONICITY_CPSI);
/// O(1) evaluation instead of ODE solver. The Pi2-Foundation anchoring is
/// algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F25 (line 436) +
/// <c>docs/proofs/PROOF_MONOTONICITY_CPSI.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F25CPsiBellPlusPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;

    /// <summary>1/4 = (1/2)² bilinear-apex maxval — the typed parent that
    /// grounds F25's <c>CrossingThreshold</c> on the Quarter axis (same anchor
    /// as F57 + Dicke + F60 fold + F62 fold). Added 2026-05-16 as a typed
    /// ctor parent.</summary>
    public QuarterAsBilinearMaxvalClaim Quarter { get; }

    /// <summary>The "4" decay-rate coefficient in <c>f = e^{−4γt}</c>. Live
    /// from <see cref="Pi2DyadicLadderClaim.Term"/>(−1) = <c>a_{−1}</c>.
    /// Same anchor as F76 (e^{−4γ₀t} mirror-pair coherence decay), F61/F63
    /// 4-block per parity, F66 multiplicity, F77 correction.</summary>
    public double DecayRateCoefficient => _ladder.Term(-1);

    /// <summary>The "2" coefficient in <c>dCΨ/dt = −2γ f(1+3f²)/3</c>.
    /// Live from <see cref="Pi2DyadicLadderClaim.Term"/>(0) = <c>a_0</c>
    /// = polynomial root d.</summary>
    public double Coefficient2 => _ladder.Term(0);

    /// <summary>The CΨ crossing threshold: <c>1/4 = a_3</c>. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(3); same anchor as F57's
    /// CrossingThreshold and Dicke + F60 + F62 fold.</summary>
    public double CrossingThreshold => _ladder.Term(3);

    /// <summary>The Bell+ crossing parameter: <c>f* = 0.8612</c> (state-specific,
    /// cubic root of <c>f(1+f²) = 3/2</c>; NOT Pi2-anchored).</summary>
    public double BellPlusFCross => 0.8612;

    /// <summary>The Bell+ K-invariant: <c>K = γ · t_cross = 0.0374</c>
    /// (state-specific, NOT Pi2-anchored).</summary>
    public double BellPlusKInvariant => 0.0374;

    /// <summary>F57's Bell+ prefactor: <c>1.080088 = 2 / 1.851701</c>; the
    /// "2" is <see cref="Coefficient2"/> (Pi2-anchored), the "1.851701"
    /// is <c>|dCΨ/dt|_{t_cross}</c> (state-specific, NOT Pi2).</summary>
    public double BellPlusF57Prefactor => 1.080088;

    /// <summary>Live closed form: <c>CΨ(t) = f(1+f²)/6, f = e^{−4γt}</c>.
    /// Throws for negative γ or t.</summary>
    public double CPsiAtTime(double gamma, double t)
    {
        if (gamma < 0.0) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be ≥ 0.");
        if (t < 0.0) throw new ArgumentOutOfRangeException(nameof(t), t, "t must be ≥ 0.");
        double f = Math.Exp(-DecayRateCoefficient * gamma * t);
        return f * (1.0 + f * f) / 6.0;
    }

    /// <summary>Live closed form: <c>dCΨ/dt = −2γ f(1+3f²)/3</c>.</summary>
    public double DCPsiDtAtTime(double gamma, double t)
    {
        if (gamma < 0.0) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be ≥ 0.");
        if (t < 0.0) throw new ArgumentOutOfRangeException(nameof(t), t, "t must be ≥ 0.");
        double f = Math.Exp(-DecayRateCoefficient * gamma * t);
        return -Coefficient2 * gamma * f * (1.0 + 3.0 * f * f) / 3.0;
    }

    /// <summary>Live drift check: at t = 0, CΨ(0) = 1·(1+1)/6 = 1/3 (Bell+
    /// initial value, above fold).</summary>
    public bool BellPlusInitialIsOneThird() =>
        Math.Abs(CPsiAtTime(0.05, 0.0) - 1.0 / 3.0) < 1e-12;

    /// <summary>Live drift check: at f = 0.8612, CΨ ≈ 0.25 (= a_3).</summary>
    public bool CrossingFConsistency()
    {
        // f * (1 + f²) / 6 should equal 0.25 at f = 0.8612
        double f = BellPlusFCross;
        double cpsi = f * (1.0 + f * f) / 6.0;
        return Math.Abs(cpsi - CrossingThreshold) < 1e-3;   // 4-decimal Bell+ f* tabulated
    }

    public F25CPsiBellPlusPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        QuarterAsBilinearMaxvalClaim quarter)
        : base("F25 CΨ(t) = f(1+f²)/6 inherits from Pi2-Foundation: 4 = a_{-1} decay rate, 2 = a_0 coefficient, 1/4 = a_3 crossing (Quarter)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F25 + " +
               "docs/proofs/PROOF_MONOTONICITY_CPSI.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (QuarterAsBilinearMaxval)")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        Quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
    }

    public override string DisplayName =>
        "F25 CΨ Bell+ Z-dephasing closed form as Pi2-Foundation a_{-1} + a_0 + a_3 inheritance";

    public override string Summary =>
        $"CΨ(t) = f(1+f²)/6, f = e^{{-4γt}}: decay rate 4 = a_{{-1}}; |dCΨ/dt| coefficient 2 = a_0; crossing 1/4 = a_3; " +
        $"Bell+ specific: f* = 0.8612, K = 0.0374, F57 prefactor 1.080088 = 2 / 1.851701 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F25 closed form",
                summary: "CΨ(t) = f(1+f²)/6, f = e^{-4γt}; dCΨ/dt = -2γf(1+3f²)/3; Tier 1 proven monotonic; O(1) evaluation");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "DecayRateCoefficient = a_{-1} = 4 (same as F76 mirror-pair decay); Coefficient2 = a_0 = 2 (polynomial root); CrossingThreshold = a_3 = 1/4 (fold)");
            yield return InspectableNode.RealScalar("DecayRateCoefficient (= a_{-1} = 4)", DecayRateCoefficient);
            yield return InspectableNode.RealScalar("Coefficient2 (= a_0 = 2)", Coefficient2);
            yield return InspectableNode.RealScalar("CrossingThreshold (= a_3 = 1/4)", CrossingThreshold);
            yield return InspectableNode.RealScalar("BellPlusFCross (state-specific)", BellPlusFCross);
            yield return InspectableNode.RealScalar("BellPlusKInvariant (state-specific)", BellPlusKInvariant);
            yield return InspectableNode.RealScalar("BellPlusF57Prefactor (= 2/1.851701)", BellPlusF57Prefactor);
            yield return new InspectableNode("F25 ↔ F57 mother-corollary chain",
                summary: "F25 closed form gives |dCΨ/dt|_{t_cross} = 1.851701 for Bell+; F57's prefactor 1.080088 = 2 / 1.851701 (the 2 IS Coefficient2 = a_0). Pattern parallel to F75 → F77 mother-claim chain.");
            yield return new InspectableNode("Verifications",
                summary: $"CΨ(0) = 1·(1+1)/6 = 1/3 (Bell+ initial, above fold; drift check: {BellPlusInitialIsOneThird()}); CΨ at f* = 0.8612 ≈ 0.25 (drift check: {CrossingFConsistency()})");
            yield return new InspectableNode("State-specific values (NOT Pi2-anchored)",
                summary: "f* = 0.8612 cubic root of f(1+f²)=3/2; K=0.0374 hardware gauge; 1.851701 = |dCΨ/dt|_{t_cross}; 6 in denominator = combinatorial Bell+ normalization");
        }
    }
}
