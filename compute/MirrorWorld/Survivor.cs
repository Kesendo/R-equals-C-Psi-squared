namespace MirrorWorld;

// The survivor: the slowest non-stationary Liouvillian mode, the thing that lasts longest. T1 at the Zeno end
// (strong dephasing, Q = J/gamma -> 0; PROOF_WEIGHT1_DEGENERACY, Theorem E (c) and (d)): in every filling
// p = 1..N-1 the slowest mode is the same density standing wave n(j) ~ cos(pi*(j-1/2)/N) (j = 1..N), dark
// (<n_XY> ~ Q^2/N^2), and at leading order every filling holds it at the same rate. On the world's own XY chain
// (zz = 0) the fillings stay tied at every order, the ladder W = sum_l c_l^dag (.) c_l carrying each filling's
// slowest mode onto the next, so there is no single survivor. A ZZ term (zz != 0; Restless's zz = 1 with hopping
// J = 2 is the isotropic Heisenberg chain)
// picks the filling at order J^4/gamma^3: a hop that changes the ZZ energy is a slower hop, and on the chain a hop
// changes it when the two sites beyond it disagree, which happens most often at half filling. So with ZZ the
// survivor is the HALF-FILLING mode at even N, R-odd and X-odd (the odd partner of its own mirror), and at odd N
// the two central fillings hold it together, X^N exchanging them (SURVIVOR_FLIP_AND_REFLECTION_ODD measures it at
// finite Q, where the filling is measured and not derived, and an uneven profile of rates can reorder the diagonal
// blocks, as GAMMA_AS_BINDING's does while the (0,1) band edge holds the survivor). The full-L survivor of the XY chain hands over to the (0,1) single-excitation band edge at
// Q_h, which is distinct from the single-excitation EP Q* from N=4 (CoherenceHorizonClaim); this object carries
// those two for the XY chain only. Only the slope of Q*(N) is asymptotically 2/pi.
public sealed class Survivor : GameObject
{
    public int N { get; }
    public double Zz { get; }                            // the ZZ term's coefficient, as Restless's zz; 0 = the XY chain
    public Survivor(World world, int n, double zz = 0.0) : base(world) { N = n; Zz = zz; }

    public bool FillingDegenerate => N > 2 && (Zz == 0.0 || N == 3);     // XY: every filling ties, at every order;
                                                                         // at N = 3 the central pair is every filling
    public bool HasHalfFillingSurvivor => N % 2 == 0 && (Zz != 0.0 || N == 2); // with ZZ: p = N/2 an integer only at even N;
                                                                               // at N = 2 the half is the only filling
    public bool HasCentralPairSurvivor => Zz != 0.0 && N % 2 == 1;       // with ZZ at odd N: p = (N-1)/2 and (N+1)/2, tied

    // SE EP: tabulated through N=24, asymptotic estimate after; the XY chain's, so none is carried with ZZ.
    public double? Qstar => Zz == 0.0 ? Formulas.Qstar(N) : null;

    // Full-L floor crossing of the XY chain, adopted from CoherenceHorizonClaim. N=2,3 are exact;
    // N=4,5 are rounded numerical readings. No later Q_h census and no census with ZZ is carried here.
    public double? HandoverQ => Zz != 0.0 ? null : N switch
    {
        2 => 1.0,
        3 => Math.Sqrt(2.0),
        4 => 1.87854,
        5 => 2.37217,
        _ => null,
    };

    public override IReadOnlyList<string> Own => new[] { "regime", "filling", "parity", "Qstar", "Qh" };
}
