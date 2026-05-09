using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F55 closed form (Tier 1 proven from D6; ANALYTICAL_FORMULAS line 1064):
///
/// <code>
///   K_death = ln(10) ≈ 2.303      dose K = γ·t for 99% absorption of slowest mortal mode
///   K_death / K_fold ~ 2.3         ratio to CΨ = 1/4 crossing dose
///   Immortal modes = N + 1         zero-absorption-rate modes (pure {I, Z} content)
/// </code>
///
/// <para>F55 is the structural absorption-dose result for a Heisenberg chain
/// under Z-dephasing. Derivation: 99% absorption of the slowest mortal mode
/// means e^{−rate_min·t} = 0.01, i.e. rate_min·t = ln(100) = 2·ln(10).
/// By F50 (D6), rate_min = 2γ, so γ·t = ln(10) = K_death. The "2" in
/// rate_min = 2γ and the "2" in 2·ln(10) cancel exactly, leaving the
/// universal K_death = ln(10) independent of N, γ, topology.</para>
///
/// <para>The N+1 immortal modes have ⟨n_XY⟩ = 0 (purely {I, Z} Pauli content,
/// invisible to the Z-dephasing dissipator). Complete absorption is therefore
/// impossible while the F1 palindrome holds: the cavity always retains light.</para>
///
/// <para>Pi2-Foundation anchor:</para>
/// <list type="bullet">
///   <item><b>RateMinCoefficient = 2 = a_0</b>: in rate_min = 2γ. Live transitively
///         from <see cref="Pi2DyadicLadderClaim.Term"/>(0) via
///         <see cref="F50WeightOneDegeneracyPi2Inheritance.DecayRateFactor"/>.
///         The cancellation 2γ·t = 2·ln(10) → γ·t = ln(10) is exact because
///         rate_min and the "2" in ln(100) = 2·ln(10) share the same a_0.</item>
/// </list>
///
/// <para>The "ln(10)" itself is irrational and has no Pi2 anchor — it comes
/// from the 99% absorption convention. The structural Pi2 content is the
/// rate_min = 2γ inheritance from F50.</para>
///
/// <para>Tier1Derived: F55 follows directly from F50 (rate_min = 2γ) plus the
/// 99%-absorption convention; valid for any Heisenberg chain, Z-dephasing,
/// all N. Pi2-Foundation anchoring is composition through F50 + Pi2DyadicLadder.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F55 (line 1064) +
/// <c>experiments/TRAPPED_LIGHT_LOCALIZATION.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F50WeightOneDegeneracyPi2Inheritance.cs</c>
/// (rate_min = 2γ universal anchor) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F55UniversalAbsorptionDosePi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F50WeightOneDegeneracyPi2Inheritance _f50;

    /// <summary>The "2" coefficient in rate_min = 2γ (transitively from F50).</summary>
    public double RateMinCoefficient => _ladder.Term(0);

    /// <summary>K_death = ln(10) ≈ 2.303 — universal absorption dose for 99%
    /// absorption of the slowest mortal mode. Independent of N, γ, topology.</summary>
    public static readonly double KDeath = Math.Log(10.0);

    /// <summary>Absorption fraction at K_death: 0.99 (99% by construction).</summary>
    public const double AbsorptionFraction = 0.99;

    /// <summary>Approximate ratio K_death / K_fold ≈ 2.3 from
    /// ANALYTICAL_FORMULAS F55 line 1067. K_fold is the dose at the CΨ = 1/4
    /// crossing.</summary>
    public const double KDeathOverKFoldApprox = 2.3;

    /// <summary>Number of immortal (zero-absorption-rate) modes: N + 1. These
    /// are the purely-{I, Z} Pauli operators that commute with Z-dephasing.</summary>
    public int ImmortalModeCount(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F55 requires N ≥ 1.");
        return N + 1;
    }

    /// <summary>Time at which K_death dose is reached: t_death = K_death / γ
    /// = ln(10) / γ. The slowest mortal mode is at 99% absorption at this time.</summary>
    public double AbsorptionTime(double gammaZero)
    {
        if (gammaZero <= 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be > 0.");
        return KDeath / gammaZero;
    }

    /// <summary>Drift check: 2γ · t_death = 2·ln(10) = ln(100).</summary>
    public bool DerivationConsistencyHolds(double gammaZero, double tolerance = 1e-12)
    {
        double t = AbsorptionTime(gammaZero);
        double rateMin = RateMinCoefficient * gammaZero;
        return Math.Abs(rateMin * t - Math.Log(100.0)) < tolerance;
    }

    /// <summary>Drift check: F55's rate_min IS F50's universal weight-1 eigenvalue position
    /// |Re(λ)| = 2γ. Verifies the F50 → F55 typed inheritance is live.</summary>
    public bool RateMinMatchesF50(double gammaZero)
    {
        return Math.Abs(RateMinCoefficient * gammaZero - Math.Abs(_f50.EigenvaluePosition(gammaZero))) < 1e-12;
    }

    public F55UniversalAbsorptionDosePi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F50WeightOneDegeneracyPi2Inheritance f50)
        : base("F55 universal absorption dose K_death = ln(10) ≈ 2.303 from F50 rate_min = 2γ; 99% absorption time t = ln(10)/γ; immortal modes = N+1",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F55 + " +
               "experiments/TRAPPED_LIGHT_LOCALIZATION.md + " +
               "compute/RCPsiSquared.Core/Symmetry/F50WeightOneDegeneracyPi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f50 = f50 ?? throw new ArgumentNullException(nameof(f50));
    }

    public override string DisplayName =>
        "F55 universal absorption dose K_death as Pi2-Foundation a_0 + F50 inheritance";

    public override string Summary =>
        $"K_death = ln(10) ≈ {KDeath:F4} (99% absorption); rate_min = 2γ (= a_0·γ via F50); immortal modes = N+1; cavity always retains light ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F55 closed form",
                summary: $"K_death = ln(10) ≈ {KDeath:F6}; absorption time t = ln(10)/γ; derivation: 99% absorption ⇒ e^{{-rate_min·t}} = 0.01 ⇒ rate_min·t = ln(100); F50: rate_min = 2γ ⇒ γ·t = ln(10)");
            yield return InspectableNode.RealScalar("RateMinCoefficient (= a_0 = 2 via F50)", RateMinCoefficient);
            yield return InspectableNode.RealScalar("KDeath (= ln(10))", KDeath);
            yield return new InspectableNode("F50 inheritance",
                summary: $"F55's rate_min = 2γ IS F50's universal weight-1 eigenvalue position; F50.DecayRateFactor (= {_f50.DecayRateFactor}) is the same '2' as F55's RateMinCoefficient");
            yield return new InspectableNode("immortal modes",
                summary: "N+1 modes have rate = 0 (purely {I, Z} Pauli content, invisible to Z-dephasing); complete absorption impossible while F1 palindrome holds; the cavity always retains light");
            yield return new InspectableNode("verified at γ=0.05",
                summary: $"AbsorptionTime(0.05) = {AbsorptionTime(0.05):G6}; rate_min·t = {RateMinCoefficient * 0.05 * AbsorptionTime(0.05):G6} (= ln(100) = {Math.Log(100):G6})");
            yield return new InspectableNode("topology universality",
                summary: "K_death is N-, γ-, and topology-independent (chain, star, ring, complete, tree all give same dose); inherits F50's topology universality");
        }
    }
}
