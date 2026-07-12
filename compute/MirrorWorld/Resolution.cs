namespace MirrorWorld;

// The individuation reading: the watching does not only price a coherence, it draws the boundary that
// makes it an object. The Pair's rate -2*gamma*k, read the other way, is a spectral LINEWIDTH (the
// width a line has): k=0 (agreement, the diagonal) has width 0, infinitely sharp -- the immortal,
// perfectly individuated object; larger k is broader, shorter-lived, less an object. Pair carries the
// width, Clock carries Q = J/gamma, Survivor carries the horizon Q*(N); but no object yet asks the one
// question this one owns: ONE object or two? Two damped modes coalesce at the EP Q = Q*(N) (Survivor's
// handover): below it, over-damped, one merged object; above it, oscillating, two distinct objects. So
// gamma (via Q) sets what counts as a thing -- the rate is the object-maker, not only the price. Sober
// here (the meaning lives in the docs); the two Own outputs are the width and the one-or-two verdict.
public sealed class Resolution : GameObject
{
    public int N { get; }
    public double J { get; }
    public double Gamma { get; }

    public Resolution(World world, int n, double j, double gamma) : base(world)
    {
        N = n;
        J = j;
        Gamma = gamma;
    }

    // left: the pair's rate read as a boundary. |Re lambda| = 2*gamma*k, the same number Pair.Rate holds.
    public double Linewidth(int k) => 2.0 * Gamma * k;

    // the adopted inputs to the verdict: Q from the Clock, Q*(N) the coherence horizon from Survivor.
    public double Q => J / Gamma;
    public double Qstar => Formulas.Qstar(N);

    // left: the one question no other object asks. Above the horizon the mode pair oscillates (two
    // distinct objects); below, it is over-damped (one merged object). The rate draws the boundary.
    public bool Resolved => Q > Qstar;

    public override IReadOnlyList<string> Own => new[] { "linewidth", "resolved" };
}
