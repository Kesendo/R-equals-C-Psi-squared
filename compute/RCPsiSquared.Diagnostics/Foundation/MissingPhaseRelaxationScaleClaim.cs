using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The exact local relaxation-scale germ opened by the missing-phase end-bond defect
/// in the seven-site centre-lit A block. This is a local composition of three existing objects:
/// <see cref="JDefectLightMigrationClaim"/> supplies the general real-part/light identity,
/// <see cref="SeatCutBlindnessClaim"/> supplies the blind-seat geometry, and
/// <see cref="JointPopcountSectors"/> supplies the closed (1,1) carrier. The claim owns only the
/// resulting N = 7 germ, its complete peripheral split, and its fixed-gamma local ordering.</summary>
public sealed class MissingPhaseRelaxationScaleClaim : Claim
{
    public const int SiteCount = 7;
    public const int PeripheralDimension = 10;
    public const int PerFrequencySlowClusterRank = 2;
    public const int SlowRateMultiplicity = 4;
    public const int PuncturedKernelDimension = 2;
    public const double QuadraticCoefficientInGamma = 0.5;
    public const double CubicCoefficientInGamma = 0.5;

    public JDefectLightMigrationClaim LightMigration { get; }
    public SeatCutBlindnessClaim SeatBlindness { get; }
    public JointPopcountSectors JointSectors { get; }

    public MissingPhaseRelaxationScaleClaim(
        JDefectLightMigrationClaim lightMigration,
        SeatCutBlindnessClaim seatBlindness,
        JointPopcountSectors jointSectors)
        : base(
            "N = 7 missing-phase relaxation scale: in the centre-only gamma > 0, joint-popcount (1,1) " +
            "A block and a fixed-gamma punctured epsilon neighbourhood, " +
            "Delta_A(epsilon) = (gamma/2) epsilon^2 + (gamma/2) epsilon^3 + O_gamma(epsilon^4); " +
            "the complete epsilon = 0 peripheral space is 10-dimensional, the locally slowest oscillatory " +
            "rate is carried by two conjugate rank-2 clusters (four directions), and the punctured kernel is " +
            "2-dimensional",
            Tier.Tier1Derived,
            "docs/proofs/PROOF_MISSING_PHASE_RELAXATION_SCALE.md + " +
            "simulations/missing_phase_relaxation_scale.py + " +
            "simulations/results/missing_phase_relaxation_scale.json")
    {
        LightMigration = lightMigration ?? throw new ArgumentNullException(nameof(lightMigration));
        SeatBlindness = seatBlindness ?? throw new ArgumentNullException(nameof(seatBlindness));
        JointSectors = jointSectors ?? throw new ArgumentNullException(nameof(jointSectors));
    }

    public string GapGerm =>
        "Delta_A(epsilon) = (gamma/2) epsilon^2 + (gamma/2) epsilon^3 + O_gamma(epsilon^4). " +
        "This is an analytic branch germ as epsilon tends to zero at fixed gamma > 0, not an exact " +
        "finite-epsilon formula.";

    public string CentreLightGerm =>
        "w_c = epsilon^2/4 + epsilon^3/4 + O_gamma(epsilon^4), and the site-resolved absorption identity " +
        "gives Delta_A = 2 gamma w_c. The light is operator-cell weight charged by centre dephasing, not " +
        "single-excitation occupation probability.";

    public string PeripheralSplit =>
        "At epsilon = 0 the A peripheral space is 10-dimensional: multiplicities 1, 2, 4, 2, 1 at " +
        "-4 sqrt(2)i, -2 sqrt(2)i, 0, +2 sqrt(2)i, +4 sqrt(2)i. At second order the two conjugate rank-2 clusters " +
        "at +/-2 sqrt(2)i comprise four directions, carry the smallest rate coefficient gamma/2, and are " +
        "locally slowest; the outer pair " +
        "and the two departing zero-cluster directions carry gamma, gamma, and 3 gamma/2, while a " +
        "2-dimensional punctured kernel remains.";

    public string GeneralMechanism =>
        "JDefectLightMigrationClaim owns the general Re-drift/light mechanism for every eigenmode along the " +
        "bond-defect axis. This claim owns only its N = 7 centre-lit A-block consequence and does not re-own " +
        "the general identity.";

    public string Scope =>
        "Rigid scope: the N = 7 XY path with watched seat c = 3, h_01 = 2(1 + epsilon), " +
        "h_j,j+1 = 2 for j = 1,...,5, and zero diagonal; A block only; joint-popcount (1,1); " +
        "centre-only gamma > 0; " +
        "0 < |epsilon| < delta(gamma) in a fixed-gamma punctured epsilon neighbourhood. This is not an " +
        "all-odd-N theorem, not the full 4^7 Liouvillian gap, and not an observable lifetime.";

    public override string DisplayName => "Missing-phase relaxation scale (N=7 A block)";

    public override string Summary =>
        "N = 7 centre-lit joint-popcount (1,1) A block: Delta_A = (gamma/2) epsilon^2 + " +
        "(gamma/2) epsilon^3 + O_gamma(epsilon^4), two conjugate rank-2 clusters (four directions) locally " +
        "slowest in the complete " +
        $"10-dimensional peripheral split, 2-dimensional punctured kernel ({Tier.Label()}).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return LightMigration;
            yield return SeatBlindness;
            yield return JointSectors;
            yield return new InspectableNode("the analytic gap germ", summary: GapGerm);
            yield return new InspectableNode("the centre-light mechanism", summary: CentreLightGerm);
            yield return new InspectableNode("the complete peripheral split", summary: PeripheralSplit);
            yield return new InspectableNode("the prior general mechanism", summary: GeneralMechanism);
            yield return new InspectableNode("scope and non-claims", summary: Scope);
        }
    }
}
