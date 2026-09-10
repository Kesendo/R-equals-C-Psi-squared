using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The certified local remainder at the N=6 Route-B crossing N6-E-A2-T-007.
/// F163 proves the crossing semisimple and states in the same breath that the algebraic
/// corollary supplies no remainder bound and no explicit ε radius. This claim owns exactly
/// that missing piece, for one crossing: a complex-ball certificate that both EP2 branches
/// are holomorphic on an explicit η disk and that the corrected Taylor tail is bounded there.
/// One crossing, conservative enclosures, no radius uniform over the 266-locus family.</summary>
public sealed class RouteBN6RemainderBoundClaim : Claim
{
    /// <summary>N and the block dimension belong to F163, the typed parent, and are not restated
    /// here as constants of this claim.</summary>
    public const string Crossing = "N6-E-A2-T-007";

    /// <summary>Working precision of the ball producer, in bits.</summary>
    public const int PrecisionBits = 512;

    /// <summary>The certified eta radii as negative powers of two, per end profile.
    /// Both branches of a profile are covered by the smaller of the profile's two radii.</summary>
    public const int OneCommonRadiusExponent = 20;
    public const int EvenCommonRadiusExponent = 18;
    public const int OddCommonRadiusExponent = 20;

    /// <summary>eta = epsilon for one/even and eta = epsilon^2 for odd, so the odd profile's
    /// eta radius 2^-20 is an epsilon radius of 2^-10.</summary>
    public const int OddEpsilonRadiusExponent = 10;

    /// <summary>The Newton image is required strictly inside this radius; every CERTIFIED image
    /// radius is smaller, the largest being 0.0038654854588457965, which the proof reports
    /// upper-rounded as 0.003865486. Failed candidate attempts and the displaced-center control
    /// record far larger image bounds; they are what rejection looks like, not counterexamples.</summary>
    public const int ImageRadiusExponent = 8;

    /// <summary>The largest certified Newton image, exactly as the certificate records it.
    /// The proof reports it upper-rounded to 0.003865486; this is the value, not the rounding.</summary>
    public const double LargestCertifiedNewtonImage = 0.0038654854588457965;

    /// <summary>The worst-branch relative remainder at half the profile's COMMON eta radius,
    /// in percent, upper-rounded to five decimals. Which branch is worst is decided by the
    /// product rho/(|z_center|-rho) times s^2/(1-s) and by neither factor alone: on one and even
    /// the wider branch carries a 3.26 times larger prefactor and still loses, because sitting at
    /// s = 1/8 instead of s = 1/2 costs it a factor 28. It is read off the certificate, not argued.</summary>
    public const double OneWorstBranchHalfRadiusPercent = 0.03630;
    public const double EvenWorstBranchHalfRadiusPercent = 0.03630;
    public const double OddWorstBranchHalfRadiusPercent = 0.04582;

    /// <summary>F163 supplies the exact premises this computation starts from: the rank, the
    /// semisimplicity, the characteristic factorization, the unique repeated root, and the count
    /// and order assignment of the local branches. All of it is premise, not result rederived here.
    /// What this computation adds from below is that |z_center| stays far above rho, so z(0) is
    /// nonzero and the response really is exactly of order eta.</summary>
    public RouteBN6A2UnfoldingClaim Unfolding { get; }

    public RouteBN6RemainderBoundClaim(RouteBN6A2UnfoldingClaim unfolding)
        : base("At the N=6 Route-B A2 crossing N6-E-A2-T-007 of the 90-dimensional (1,2) coherence " +
               "block (XY, Delta=0, uniform gamma=1, interior bonds one), each of the three end " +
               "profiles one/even/odd opens two EP2 branches that are holomorphic in eta on an " +
               "explicitly certified disk, with eta = epsilon for one and even and eta = epsilon^2 " +
               "for odd. A 512-bit complex-ball Banach contraction with a fixed preconditioner " +
               "certifies common eta radii 2^-20, 2^-18 and 2^-20, the last being an epsilon radius " +
               "of 2^-10. On each disk q(eta) = q0 + eta z(eta) and lambda(eta) = lambda0 + eta h(eta) " +
               "with z and h holomorphic and enclosed within rho = 2^-8 of their centers, every Newton " +
               "image radius lying strictly inside rho. The corrected Taylor tail, corrected by the exact " +
               "Taylor coefficients c = z(0) and d = z'(0), obeys the relative bound " +
               "rho/(|z_center|-rho) * s^2/(1-s) with s = |eta|/R and R the certified BRANCH-SPECIFIC " +
               "eta radius, measured against the actual displacement q(eta)-q0 rather than the leading " +
               "term. At half the common eta radius the WORST-BRANCH bounds are 0.03630 percent, " +
               "0.03630 percent and 0.04582 percent; rounding c and d adds their rounding errors on top. " +
               "At eta = 0 the crossing stays " +
               "semisimple, which is F163's theorem and a premise here. These are sufficient " +
               "neighborhoods: a failed larger candidate means the enclosure failed and does not " +
               "locate a singularity or a Taylor radius. It supplies no physical-real-q operating point, " +
               "no different-N result and no all-locus radius. One crossing, one N, no uniform family radius.",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_ROUTE_B_N6_REMAINDER_BOUND.md + " +
               "simulations/route_b_n6_remainder_ball.py + simulations/route_b_n6_ball_base.py + " +
               "simulations/results/route_b_n6_remainder_ball.json")
    {
        Unfolding = unfolding ?? throw new ArgumentNullException(nameof(unfolding));
    }

    public override string DisplayName =>
        "Certified local remainder at the N=6 crossing N6-E-A2-T-007 (complex-ball, one crossing)";

    public override string Summary =>
        "Both EP2 branches holomorphic on certified eta disks 2^-20 / 2^-18 / 2^-20 for the " +
        "one/even/odd end profiles (2^-10 in epsilon for odd), corrected Taylor tail bounded " +
        "relative to the actual displacement; no uniform radius over the 266-locus family (Tier 1, derived)";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("what F163 leaves open, and this closes",
                summary: "F163 proves all 266 loci semisimple and states that the algebraic corollary alone " +
                         "supplies no remainder bound and no explicit epsilon radius. This claim supplies both, " +
                         "for the single crossing " + Crossing + ". The family theorem still carries no uniform radius.");
            yield return new InspectableNode("the certificate, not a rounded model",
                summary: PrecisionBits + "-bit complex ball arithmetic in which every enclosure is outward-rounded; " +
                         "floating point only chooses starting points, the pivot minor, the preconditioner and the " +
                         "branch and lambda centers, and a successful ball solve, not the float step, establishes " +
                         "invertibility. The contraction fixes that preconditioner Y, " +
                         "bounds the row sums of I-Y*dF/dx below one and keeps " +
                         "the Newton image strictly inside rho = 2^-" + ImageRadiusExponent + ". The tabulated decimals are " +
                         "reported values; the enclosures, not the rounded entries, carry the proof.");
            yield return new InspectableNode("one word list, three houses",
                summary: "The end profiles are named one/even/odd here, after the proof, one/equal/opposite in " +
                         "F163 and in the arc, and One-sided/Symmetric/Antisymmetric in the local-validity note. " +
                         "All three lists name the same three profiles in the same order: a single detuned end, " +
                         "both ends detuned alike, both ends detuned oppositely.");
            yield return new InspectableNode("the odd profile is analytic without a cut",
                summary: "The reflection similarity T_epsilon = E + epsilon*O turns the opposite-end profile into a " +
                         "problem analytic in eta = epsilon^2, so no square-root branch cut is crossed and the same " +
                         "contraction applies. Its eta radius 2^-" + OddCommonRadiusExponent + " reads as |epsilon| <= 2^-" +
                         OddEpsilonRadiusExponent + ", and half the eta radius means |epsilon| <= 1/(1024*sqrt(2)), not half of 1/1024.");
            yield return new InspectableNode("the fences",
                summary: "One crossing of one N, complex-q analytic continuation, sufficient rather than maximal " +
                         "neighborhoods. Not an all-N statement, not a uniform family radius, not a physical " +
                         "measurement, and not a location of any singularity.");
            yield return Unfolding;
        }
    }
}
