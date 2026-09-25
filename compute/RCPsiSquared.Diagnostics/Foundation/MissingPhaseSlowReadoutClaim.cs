using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The prepared physical readout of the existing N=7 slow spectral cluster.
/// The spectral theorem owns the rate and local minimality; F70 removes the intercopy B
/// contribution to pair marginals; the surviving blind ray's even-site support, the zero mode of the
/// bipartite path living on its majority sublattice, removes the slow A residue from the three odd-site
/// pairs; the chiral pairing E_+ = C conj(E_-) C gives the mixed even/odd antisymmetry and the leakage
/// contraction. Exact evidence is recomputed by MissingPhaseSlowReadoutWitness.</summary>
public sealed class MissingPhaseSlowReadoutClaim : Claim
{
    public MissingPhaseRelaxationScaleClaim RelaxationScale { get; }
    public F70DeltaNSelectionRulePi2Inheritance SelectionRule { get; }
    public ChiralKClaim Chirality { get; }

    public MissingPhaseSlowReadoutClaim(
        MissingPhaseRelaxationScaleClaim relaxationScale,
        F70DeltaNSelectionRulePi2Inheritance selectionRule,
        ChiralKClaim chirality)
        : base(
            "N=7 missing-phase slow readout: the physical end-odd preparation has a complete rank-two " +
            "slow residue with ZZ_06 coefficient tending to -1/2, hence paired oscillation amplitude " +
            "tending to 1; four even pairs couple at epsilon, the leakage at epsilon^2, three odd pairs " +
            "never; every pair loses the uniform sine quadrature, which a Z-tagged three-site current reads",
            Tier.Tier1Derived,
            "docs/proofs/PROOF_MISSING_PHASE_SLOW_READOUT.md + " +
            "compute/RCPsiSquared.Diagnostics/Foundation/MissingPhaseSlowReadoutWitness.cs " +
            "(inspect --root missingphasereadout)")
    {
        RelaxationScale = relaxationScale ?? throw new ArgumentNullException(nameof(relaxationScale));
        SelectionRule = selectionRule ?? throw new ArgumentNullException(nameof(selectionRule));
        Chirality = chirality ?? throw new ArgumentNullException(nameof(chirality));
    }

    public string Preparation =>
        "N=7 XY path, centre Z dephasing, h_01=2(1+epsilon), other hoppings 2. " +
        "d=(|0>-|6>)/sqrt(2), and the physical state is (d+X^7 d)/sqrt(2); A_init=B_init=dd-dagger. " +
        "The intercopy blocks have |DeltaN|=5, so F70 makes their two-site partial traces vanish.";

    public string Projection =>
        "With the surviving blind-ray projector P_v and simple K=-ih-2gamma|3><3| projectors E_-/E_+, " +
        "the complete A Riesz projector is P_-^A(X)=E_- X P_v + P_v X E_+^dagger. " +
        "Both members of the rank-two cluster are required. This is not an orthogonal light projector.";

    public string PreparedSignal =>
        "For M_-=P_-^A(dd-dagger), the slow ZZ_06 contribution is " +
        "2 Re[f(epsilon) exp(lambda_-(epsilon)t)], f=-2(M_00+M_66), f(0)=-1/2. " +
        "Its envelope amplitude is 2|f| exp(-Delta_A t), with 2|f| tending to 1. " +
        "The uniform XX_06 residue is -1/4 and ||M_-(0)||_F^2=1/4. Analytic continuity " +
        "gives nonzero coupling in a sufficiently small neighbourhood at each fixed gamma>0.";

    public string LeakageGerm =>
        "For the original reflection-even projector Q_R, " +
        "Tr(Q_R M_-)=(-1+i sqrt(2)gamma)epsilon^2/16+O_gamma(epsilon^3). " +
        "The exact live witness derives this coefficient from the normalized blind-vector derivative " +
        "and an 8x8 bordered eigenvector-derivative solve; it does not sample finite-epsilon roots.";

    public string ReadoutZeros =>
        "Pairs (1,3), (1,5), (3,5), with zero-based seats, lose the entire slow contribution for this " +
        "preparation at every epsilon of the local branch: the surviving blind ray v_r has no odd-site " +
        "entry, and each entry (a,b) of E_- X P_v + P_v X E_+^dagger is a term carrying (v_r)_b plus a " +
        "term carrying (v_r)_a. At epsilon=0 the " +
        "even pairs (0,2), (0,4), (2,6), (4,6) miss the residue as well, since e+ and e- agree on the even " +
        "sites; they switch on at first order, with pair-map Frobenius norms sqrt(22)/16, sqrt(6)/16 and " +
        "sqrt(12 gamma^2+6)/16 per |epsilon| for (0,2), (4,6) and each of (0,4), (2,6), so for small " +
        "epsilon != 0 the three odd pairs are the only null pairs. The chiral pairing E_+ = C conj(E_-) C " +
        "makes the mixed even/odd cells antisymmetric. At epsilon=0 the sine quadrature vanishes on all 21 " +
        "physical pair maps, while the isolated paired decoder contribution has half trace norm 1/2 at every " +
        "phase. A complex residue requires (M_ab+M_ba)/2, not Re(M_ab), before conjugate reconstruction.";

    public string ThreeSpinReadout =>
        "At epsilon=0, Y0 X1 Z_k with any site k in {2,...,6} reads the sine quadrature every pair " +
        "loses: i sqrt(2)/8 on M_-(0), sqrt(2)/4 on the sine and 0 on the cosine quadrature, because Z_k is " +
        "+1 on one copy of the (0,1) coherence and -1 on the other. Three sites are the minimum before " +
        "decoding, and the flipped copy is what forces it: on one copy the pair Y0 X1 reads the same " +
        "component, while ZZ, XX and the leakage read only the A block and are unchanged on one copy. On " +
        "the uniform trajectory d(theta) = (v + cos theta u - i sin theta b1)/sqrt(2), theta = 2 sqrt(2) t, " +
        "<Y0X1Z3> = sin theta (1 + cos theta)/(2 sqrt(2)) and <Z0Z6> = 1 - (1 + cos theta)^2/2 at every " +
        "gamma >= 0, and d<Z0Z6>/dt = 8 <Y0X1Z3> because i[H, Z0Z6] = 2(Y0X1 - X0Y1)Z6 + 2(X5Y6 - Y5X6)Z0 " +
        "and the trajectory is reflection-symmetric. At t+ = pi/(4 sqrt(2)) and t- = 3 pi/(4 sqrt(2)) all " +
        "21 pair pages agree, while 24 of the 35 three-site pages, those holding site 1 or 5 and an even " +
        "site, tell the two apart.";

    public string Scope =>
        "Fixed N=7, XY and centre-only gamma>0; local in epsilon at each fixed gamma. " +
        "The complete signal also contains a stationary term and faster A contributions. " +
        "The three-site minimum is claimed at epsilon=0 only. No gamma-uniform neighbourhood, full 4^7 Liouvillian " +
        "gap, all-odd-N theorem, hardware qualification, or lifetime of the nonlinear d_out/d_2, which " +
        "subtract the matched gamma=0 run as a moving unitary reference, is claimed.";

    public override string DisplayName => "Missing-phase slow physical readout (N=7)";
    public override string Summary =>
        "The existing end preparation exposes the N=7 slow rate through ZZ_06 with amplitude tending to 1; " +
        "leakage residue starts at epsilon^2; a Z-tagged three-site current reads what every pair loses. " +
        "Exact live root: missingphasereadout. Local fixed-gamma coupling, not a d_out/d_2 lifetime.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return RelaxationScale;
            yield return SelectionRule;
            yield return Chirality;
            yield return new InspectableNode("physical preparation", Preparation);
            yield return new InspectableNode("complete Riesz projection", Projection);
            yield return new InspectableNode("prepared endpoint signal", PreparedSignal);
            yield return new InspectableNode("leakage germ", LeakageGerm);
            yield return new InspectableNode("physical output zeros and quadratures", ReadoutZeros);
            yield return new InspectableNode("three-spin reading (proof section 7.1)", ThreeSpinReadout);
            yield return new InspectableNode("live exact evidence",
                "inspect --root missingphasereadout; MissingPhaseSlowReadoutWitness recomputes in Q(sqrt(2),i) " +
                "the uniform residue and the cluster's completeness, the physical pair maps with their null " +
                "tiers and first-order coefficients, the decoder quadratures, the three-spin reading in the " +
                "full 2^7 space and the leakage derivative.");
            yield return new InspectableNode("scope", Scope);
        }
    }
}
