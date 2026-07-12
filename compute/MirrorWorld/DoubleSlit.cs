namespace MirrorWorld;

// The double slit, composed from the atoms -- the access layer the world was missing: the phenomenon
// assembled under its own name, so it is recognizable where Pair + Field alone were not. It IS Field at
// N=1: two "places" |L>, |R>, seeded with the two humps (the diagonal) and the one between |L><R| (the
// off-diagonal). Nothing new is computed; this only NAMES and reads what Pair and Field already do:
//   humps  = Field.Structure  -- the immortal diagonal (k=0), the two plain humps (the "particle face"),
//   fringe = Field.Novelty    -- the between |L><R|, the k=1 coherence paying -2gamma (the "wave face"),
//   the between's rate -2*gamma*k = -2gamma is the generator of |rho_LR(t)| = |rho_LR(0)| e^{-2gamma t}.
// The MEANING (why the pattern is what not-being-watched looks like) lives in the docs, not here:
// docs/quantum/DOUBLE_SLIT_TRANSLATED.md, which already breadcrumbs back to Pair.cs.
public sealed class DoubleSlit : GameObject
{
    public double Gamma { get; }
    readonly Field cloud;
    readonly Pair between;

    public DoubleSlit(World world, double gamma) : base(world)
    {
        Gamma = gamma;
        cloud = new Field(world, 1, gamma);           // N=1: the two-place route register
        cloud.SeedUniform();                          // the two humps (diagonal) + the one between (off-diagonal)
        between = new Pair(world, 0b0, 0b1, gamma);    // |L><R|: disagreement k=1
    }

    public void Watch(double dt) => cloud.Step(dt);    // one tick of the watching
    public double T => cloud.T;

    // left: the phenomenon read off the composed atoms.
    public double Humps => cloud.Structure;            // the immortal diagonal: the two plain humps
    public double Fringe => cloud.Novelty;             // the between |L><R|: the pattern that pays
    public double Visibility => cloud.Novelty / cloud.Structure;   // V = 2|rho_LR|/(rho_LL+rho_RR): 1 at the balanced seed, in [0,1]
    public double BetweenRate => between.Rate;         // -2*gamma*k at k=1 = -2gamma (the e^{-2gamma t} law)

    public override IReadOnlyList<string> Own => new[] { "humps", "fringe", "visibility" };
}
