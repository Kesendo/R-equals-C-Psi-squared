namespace MirrorWorld;

// The field of possibilities (rule 1 of ClaudeTasks/DIAGONAL_PROTOCOL_GAME.md): a weight on every pair
// |i><j| over N units. It inherits the frame from the World (right) and produces its own split (left):
// structure = the diagonal (k=0) that stays, novelty = the off-diagonal that fades. One Step applies the
// one question (rule 2): each weight keeps (1 + Re lambda * dt) of itself, Re lambda the pair's own rate
// -2*gamma*k; structure costs nothing and stays, novelty fades faster the more it disagrees. No Hamiltonian
// yet, only the watching -- the bedrock split (rules 1-3) running in time, nothing more. (Step 2 adds the
// inner restlessness that is BORN from structure; for now novelty can only fade.)
public sealed class Field : GameObject
{
    public int N { get; }
    public double Gamma { get; }
    public double T { get; private set; }

    readonly int dim;
    readonly double[,] w;              // weight on pair (i,j), i,j in [0, 2^N)
    readonly double[,] rate;           // Re lambda per pair, read once from Pair (time-independent)
    readonly int[,] dis;               // disagreement k per pair, read once from Pair
    readonly (int i, int j)[] alive;   // the cells that can change (k>0); the immortal diagonal is skipped

    public Field(World world, int n, double gamma) : base(world)
    {
        N = n;
        Gamma = gamma;
        dim = 1 << n;
        w = new double[dim, dim];
        rate = new double[dim, dim];
        dis = new int[dim, dim];
        var aliveList = new List<(int, int)>();
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
            {
                var p = new Pair(world, i, j, gamma);   // the atom: its own disagreement and its own rate
                rate[i, j] = p.Rate;
                dis[i, j] = p.Disagreement;
                if (i < j && p.Disagreement != 0) aliveList.Add((i, j));   // step only the upper off-diagonal:
                                                                          // k=0 is immortal (held), the lower is the mirror
            }
        alive = aliveList.ToArray();
    }

    public int Dim => dim;

    // the lower triangle is the mirror of the upper (rho = rho-dagger): only the upper is stored, set or read.
    public double this[int i, int j]
    {
        get => i <= j ? w[i, j] : w[j, i];
        set { if (i <= j) w[i, j] = value; else w[j, i] = value; }
    }

    // the knower's census: cells the runner steps (the upper off-diagonal) vs the immortal diagonal (held).
    public int AliveCount => alive.Length;
    public int ImmortalCount => dim;

    // seed the cloud (rule 1): equal weight on every pair. only the upper triangle is stored; the lower is
    // the mirror of it (rho = rho-dagger). tests overwrite a single pair instead.
    public void SeedUniform()
    {
        for (int i = 0; i < dim; i++)
            for (int j = i; j < dim; j++)
                w[i, j] = 1.0;
    }

    // one tick of the one question (rule 2): w *= 1 + Re lambda * dt. only the alive (k>0) are stepped;
    // the immortal diagonal (k=0, rate 0) is held -- 100% known waste, never recomputed (the knower's cut).
    public void Step(double dt)
    {
        foreach (var (i, j) in alive)
            w[i, j] *= 1.0 + rate[i, j] * dt;
        T += dt;
    }

    // the magnitude carried at each disagreement k = 0..N (k=0 is the structure, k>0 the novelty). Reads
    // the upper triangle + diagonal; each off-diagonal cell stands for its mirror twin, so counts double.
    public double[] WeightByDisagreement()
    {
        var byK = new double[N + 1];
        for (int i = 0; i < dim; i++)
            for (int j = i; j < dim; j++)
                byK[dis[i, j]] += (i == j ? 1.0 : 2.0) * Math.Abs(w[i, j]);
        return byK;
    }

    // left: structure (the diagonal that stays) and novelty (the off-diagonal that fades).
    public double Structure
    {
        get { double s = 0; for (int i = 0; i < dim; i++) s += w[i, i]; return s; }
    }

    public double Novelty => WeightByDisagreement().Skip(1).Sum();

    public override IReadOnlyList<string> Own => new[] { "structure", "novelty" };
}
