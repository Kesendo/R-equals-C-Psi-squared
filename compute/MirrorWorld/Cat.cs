namespace MirrorWorld;

// Schrodinger's cat, composed from the atoms -- the k=N sighting of the same law the double slit is the
// k=1 sighting of ("k=1 is the pattern, k=2 the pair, k=N the cat", SCHRODINGERS_CAT_TRANSLATED.md). It IS
// Field at N with just the two definite branches |0..0> (dead) and |1..1> (alive) on the diagonal (k=0,
// immortal poles) and the one coherence between them |0..0><1..1> (disagreement k=N, the maximal, the
// "both at once"), paying the deepest rate -2*gamma*N = -2*Sigma_gamma. Nothing new is computed; it names
// and reads what Field and Pair already do, so the cat is recognizable. The point it makes: the bigger the
// superposition (larger N, k=N), the FASTER the between decoheres -- N times faster than the slit's k=1 --
// which is why we never catch a macroscopic thing "both at once". Meaning in the doc, not here.
public sealed class Cat : GameObject
{
    public int N { get; }
    public double Gamma { get; }
    readonly Field cloud;
    readonly Pair coherence;

    public Cat(World world, int n, double gamma) : base(world)
    {
        N = n;
        Gamma = gamma;
        int alive = (1 << n) - 1;               // |1..1>: every bit set (the other pole)
        cloud = new Field(world, n, gamma);
        cloud[0, 0] = 1.0;                       // the "dead" branch |0..0><0..0| (a definite pole, immortal)
        cloud[alive, alive] = 1.0;              // the "alive" branch |1..1><1..1| (the other pole, immortal)
        cloud[0, alive] = 1.0;                  // the cat coherence |0..0><1..1|: disagreement k=N
        coherence = new Pair(world, 0, alive, gamma);
    }

    public void Watch(double dt) => cloud.Step(dt);    // one tick of the watching
    public double T => cloud.T;

    // left: the phenomenon read off the composed atoms.
    public double Branches => cloud.Structure;         // the two definite poles (dead, alive): immortal
    public double CatCoherence => cloud.Novelty;       // the "both at once", the maximal-disagreement between
    public double CoherenceRate => coherence.Rate;     // -2*gamma*N = -2*Sigma_gamma, the deepest rate

    public override IReadOnlyList<string> Own => new[] { "branches", "catCoherence" };
}
