namespace MirrorWorld;

// The survivor: the slowest non-stationary Liouvillian mode, the thing that lasts longest. T1:
// in the strong-dephasing (low-Q) regime it is the HALF-FILLING mode (k=N/2, even N only), an
// antisymmetric density standing wave n(j) ~ cos(pi*(j-1/2)/N), R-odd and X-odd (the odd partner of
// its own mirror), dark (<n_XY> ~ Q^2/N^2). The full-L survivor hands over to the (0,1)
// single-excitation band edge at Q_h, which is distinct from the single-excitation EP Q* from N=4.
// Odd N has no (N/2,N/2) sector, so no half-filling survivor (SURVIVOR_FLIP_AND_REFLECTION_ODD,
// CoherenceHorizonClaim). Only the slope of Q*(N) is asymptotically 2/pi.
public sealed class Survivor : GameObject
{
    public int N { get; }
    public Survivor(World world, int n) : base(world) => N = n;

    public bool HasHalfFillingSurvivor => N % 2 == 0;    // k=N/2 integer only at even N
    public double Qstar => Formulas.Qstar(N);            // SE EP: tabulated through N=24, asymptotic estimate after

    // Full-L floor crossing, adopted from CoherenceHorizonClaim. N=2,3 are exact;
    // N=4,5 are rounded numerical readings. No later Q_h census is carried here.
    public double? HandoverQ => N switch
    {
        2 => 1.0,
        3 => Math.Sqrt(2.0),
        4 => 1.87854,
        5 => 2.37217,
        _ => null,
    };

    public override IReadOnlyList<string> Own => new[] { "regime", "parity", "Qstar", "Qh" };
}
