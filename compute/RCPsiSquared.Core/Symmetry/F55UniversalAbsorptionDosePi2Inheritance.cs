using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F55 closed form (Tier 1 proven from D6; ANALYTICAL_FORMULAS line 1064):
///
/// <code>
///   K_death = ln(10) ≈ 2.303      dose K = γ·t for 99% absorption of slowest mortal mode
///   K_death / K_fold ~ 62          ratio to CΨ = 1/4 crossing dose
///                                  (K_fold = 0.03735 from F25; 2.3 was K_death
///                                   itself, i.e. K_fold taken as 1)
///   Immortal modes = N + 1         zero-absorption-rate modes (pure {I, Z} content)
/// </code>
///
/// <para>F55 is the structural absorption-dose result for a Heisenberg chain
/// under Z-dephasing. Derivation: 99% absorption of the slowest mortal mode
/// means e^{−rate_min·t} = 0.01, i.e. rate_min·t = ln(100) = 2·ln(10).
/// Where rate_min = 2γ, γ·t = ln(10) = K_death: the "2" in rate_min and the
/// "2" in 2·ln(10) cancel exactly, leaving K_death = ln(10) independent of
/// N, γ and topology.</para>
///
/// <para><b>What supplies rate_min = 2γ, and when.</b> Not F50: F50 counts the
/// modes sitting AT 2γ, which is existence, not minimality. Minimality is D6,
/// and D6 holds only above an N-dependent coupling threshold Q*_gap(N) in Q = J/γ
/// (0.5000, 0.8002, 1.3422, 1.8194 at N = 2, 3, 4, 5 on the Heisenberg chain,
/// in PAULI-J units; the C# engines use spin J/4, so convert before comparing).
/// Distinct from CoherenceHorizonClaim's Q*(N) ladder, which is a different
/// threshold sharing the symbol.
/// Below it the slowest mortal mode is Zeno-suppressed, rate_min is far smaller
/// and the dose is correspondingly larger. K_death = ln(10) is therefore the
/// strong-coupling value, not a value at every Q.</para>
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
/// <para>The "ln(10)" itself is irrational and has no Pi2 anchor; it comes
/// from the 99% absorption convention. The structural Pi2 content is the
/// rate_min = 2γ inheritance from F50.</para>
///
/// <para>Tier1Derived above Q*_gap(N): F55 follows from D6 (rate_min = 2γ) plus the
/// 99%-absorption convention; valid for any Heisenberg chain, Z-dephasing,
/// all N, in the strong-coupling regime. Pi2-Foundation anchoring is composition through F50 + Pi2DyadicLadder.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F55 +
/// <c>experiments/TRAPPED_LIGHT_LOCALIZATION.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F50WeightOneDegeneracyPi2Inheritance.cs</c>
/// (the 2γ level: existence anchor, not minimality) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F55UniversalAbsorptionDosePi2Inheritance : Claim, IZ2AxisClaim
{

    /// <summary>The F1² / Π²_Z axis (bit_b parity, n_Y + n_Z mod 2). The
    /// canonical Pi²-Inheritance axis. The bit_a-twin (Π²_X / F61 axis) is
    /// currently not typed for this Claim.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>The typed bit_a-twin sibling, if one exists. Currently null
    /// (no bit_a twin is typed for this Claim; this is an open slot in the
    /// cubic-architecture coverage).</summary>
    public Claim? BitATwin => null;

    // Absorption-Theorem descendant; bit_a twin holds by the Hadamard X↔Z duality (PROOF_BIT_A_TWIN_VIA_HADAMARD.md).
    public BitATwinClassification BitATwinStatus => BitATwinClassification.CoveredByHadamardDuality;
    public Pi2DyadicLadderClaim Ladder { get; }
    public F50WeightOneDegeneracyPi2Inheritance F50 { get; }
    /// <summary>The "2" coefficient in rate_min = 2γ (transitively from F50, which
    /// anchors the coefficient and the existence of the level, not its minimality).</summary>
    public double RateMinCoefficient => Ladder.Term(0);

    /// <summary>K_death = ln(10) ≈ 2.303; absorption dose for 99% absorption of
    /// the slowest mortal mode. Independent of N, γ and topology, but only where
    /// rate_min = 2γ, i.e. above the coupling threshold Q*_gap(N) of D6.</summary>
    public static readonly double KDeath = Math.Log(10.0);

    /// <summary>Absorption fraction at K_death: 0.99 (99% by construction).</summary>
    public const double AbsorptionFraction = 0.99;

    /// <summary>Ratio K_death / K_fold ≈ 61.7. K_fold is the dose at the
    /// CΨ = 1/4 crossing, given in closed form by F25: f*(1 + f*²) = 3/2 puts
    /// CΨ at exactly 1/4 with K_fold = 0.03735, so 2.302585 / 0.03735 = 61.65.
    /// (This read 2.3 until 2026-07, which is K_death itself, i.e. K_fold
    /// taken as 1.)</summary>
    public const double KDeathOverKFoldApprox = 61.65;

    /// <summary>The fold dose itself: K = gamma*t at the CΨ = 1/4 crossing,
    /// from F25's closed form. Pinned in Symphony's chain-collapse lens too.</summary>
    public const double KFold = 0.03735;

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

    /// <summary>Drift check: F55's rate_min sits at F50's weight-1 eigenvalue position
    /// |Re(λ)| = 2γ. This checks the typed F50 → F55 inheritance is live; it does NOT
    /// establish that 2γ is the MINIMUM, which is D6's and holds only above Q*_gap(N).</summary>
    public bool RateMinMatchesF50(double gammaZero)
    {
        return Math.Abs(RateMinCoefficient * gammaZero - Math.Abs(F50.EigenvaluePosition(gammaZero))) < 1e-12;
    }

    public F55UniversalAbsorptionDosePi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F50WeightOneDegeneracyPi2Inheritance f50)
        : base("F55 absorption dose K_death = ln(10) ≈ 2.303 from D6 rate_min = 2γ (holds above the coupling threshold Q*_gap(N); F50 gives existence at 2γ, not minimality); 99% absorption time t = ln(10)/γ; immortal modes = N+1",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F55 + " +
               "experiments/TRAPPED_LIGHT_LOCALIZATION.md + " +
               "compute/RCPsiSquared.Core/Symmetry/F50WeightOneDegeneracyPi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        Ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        F50 = f50 ?? throw new ArgumentNullException(nameof(f50));
    }

    public override string DisplayName =>
        "F55 absorption dose K_death as Pi2-Foundation a_0 + F50 inheritance (strong-coupling regime)";

    public override string Summary =>
        $"K_death = ln(10) ≈ {KDeath:F4} (99% absorption) where rate_min = 2γ, i.e. above the D6 coupling threshold Q*_gap(N); immortal modes = N+1; cavity always retains light ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F55 closed form",
                summary: $"K_death = ln(10) ≈ {KDeath:F6}; absorption time t = ln(10)/γ; derivation: 99% absorption ⇒ e^{{-rate_min·t}} = 0.01 ⇒ rate_min·t = ln(100); F50: rate_min = 2γ ⇒ γ·t = ln(10)");
            yield return InspectableNode.RealScalar("RateMinCoefficient (= a_0 = 2 via F50)", RateMinCoefficient);
            yield return InspectableNode.RealScalar("KDeath (= ln(10))", KDeath);
            yield return new InspectableNode("F50 inheritance",
                summary: $"F55's rate_min sits at F50's weight-1 eigenvalue position (existence; minimality is D6, above Q*_gap(N)); F50.DecayRateFactor (= {F50.DecayRateFactor}) is the same '2' as F55's RateMinCoefficient");
            yield return new InspectableNode("immortal modes",
                summary: "N+1 modes have rate = 0 (purely {I, Z} Pauli content, invisible to Z-dephasing); complete absorption impossible while F1 palindrome holds; the cavity always retains light");
            yield return new InspectableNode("verified at γ=0.05",
                summary: $"AbsorptionTime(0.05) = {AbsorptionTime(0.05):G6}; rate_min·t = {RateMinCoefficient * 0.05 * AbsorptionTime(0.05):G6} (= ln(100) = {Math.Log(100):G6})");
            yield return new InspectableNode("topology universality (above Q*_gap)",
                summary: "above Q*_gap(N): K_death is N-, γ- and topology-independent (chain, star, ring, complete, tree all give the same dose). Below Q*_gap(N) the slowest mortal mode is Zeno-suppressed and the dose is larger; the threshold, not the dose, is what carries N and the Hamiltonian");
        }
    }
}
