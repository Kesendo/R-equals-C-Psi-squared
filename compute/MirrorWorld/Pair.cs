namespace MirrorWorld;

// A pair of versions |i><j|: the first real object in the world. It inherits the frame (x,y,z) from
// the World (right) and produces its own two outputs (left): the disagreement count k = popcount(i^j),
// the only thing the watching reads, and the rate Re lambda = -2*gamma*k it decays at. Empty world:
// no Hamiltonian, only the watching, so the pair sits still at its rate.
public sealed class Pair : GameObject
{
    public int I { get; }
    public int J { get; }
    public double Gamma { get; }

    public Pair(World world, int i, int j, double gamma) : base(world)
    {
        I = i;
        J = j;
        Gamma = gamma;
    }

    // left: the pair's own.
    public int Disagreement
    {
        get
        {
            int x = I ^ J, c = 0;
            while (x != 0) { c += x & 1; x >>= 1; }
            return c;
        }
    }

    public double Rate => -2.0 * Gamma * Disagreement;   // Re lambda

    public override IReadOnlyList<string> Own => new[] { "disagreement", "rate" };
}
