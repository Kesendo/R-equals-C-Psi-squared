namespace MirrorWorld;

// The survivor: the slowest non-stationary Liouvillian mode, the thing that lasts longest. T1:
// in the strong-dephasing (low-Q) regime it is the HALF-FILLING mode (k=N/2, even N only), an
// antisymmetric density standing wave n(j) ~ cos(pi*(j-1/2)/N), R-odd and X-odd (the odd partner of
// its own mirror), dark (<n_XY> ~ Q^2/N^2). Above the coherence horizon Q*(N) ~ 2N/pi the survivor
// hands over to the (0,1) single-excitation band edge (<n_XY>=1, rate -2*gamma, oscillating). Odd N
// has no (N/2,N/2) sector, so no half-filling survivor (SURVIVOR_FLIP_AND_REFLECTION_ODD,
// ClockHandLadder Q*(N)=2N/pi).
public sealed class Survivor : GameObject
{
    public int N { get; }
    public Survivor(World world, int n) : base(world) => N = n;

    public bool HasHalfFillingSurvivor => N % 2 == 0;    // k=N/2 integer only at even N
    public double Qstar => Formulas.Qstar(N);            // exact T1 values (2N/pi only asymptotically)

    public override IReadOnlyList<string> Own => new[] { "regime", "parity", "Qstar" };
}
