using System.Numerics;

namespace MirrorWorld;

// The living world (step 2 of ClaudeTasks/DIAGONAL_PROTOCOL_GAME.md): the empty world's watching plus the
// inner restlessness, a Hamiltonian on the geometry. The full Lindblad loop rho-dot = -i[H,rho] + D[rho]:
// the handshake H (a flip-flop on each bond) takes population (the diagonal, structure) and makes coherence
// (off-diagonal, novelty) -- novelty BORN from structure (rule 4) -- while D (the watching, the Pair rate
// -2*gamma*k) culls it by disagreement. So novelty no longer only fades; it is born and culled, settling
// into a living balance. The empty Field (H off) was step 1. Stepped with a small RK4 (H couples the cells,
// no per-cell closed form). Full rho here; the knower's cuts on the living world (the mirror, the forbidden
// blocks) come next.
public sealed class Restless : GameObject
{
    public int N { get; }
    public double J { get; }
    public double Gamma { get; }
    public double T { get; private set; }

    readonly int dim;
    readonly Complex[,] rho;        // the living weights (rho = rho-dagger)
    readonly double[,] mask;        // the watching: D[rho]_ij = -2*gamma*k * rho_ij, read from Pair
    readonly int[,] dis;            // disagreement k per cell, read from Pair
    readonly double[,] h;           // the handshake H: flip-flop on the bonds (real, symmetric)
    readonly (int a, int b)[] bonds; // the geometry the handshake rides (default: a chain)
    readonly int[] pc;              // popcount per basis index (the excitation number H conserves)
    readonly HashSet<(int p, int q)> occupied = new();   // the joint-popcount blocks the seed lives in
    (int i, int j)[]? alive;        // the cells inside the occupied blocks (F63); the rest is forbidden

    // antiWatching = the rules turned around (2026-07-03, the mirror's rho-level face): the watching
    // reads AGREEMENT instead of disagreement, rate -2*gamma*(N-k). The mirror guarantees this is the
    // SAME world read through the bra complement: anti-rho(t) = rho(t) * X^N exactly (H commutes with
    // X^N, and k(i, complement j) = N - k(i,j)), so nothing was taken -- the conservation law just
    // moves from the trace (the diagonal dies fastest here) to the ANTI-trace, sum of rho(i, ~i).
    public Restless(World world, int n, double j, double gamma, (int a, int b)[]? bonds = null,
        bool antiWatching = false) : base(world)
    {
        N = n;
        J = j;
        Gamma = gamma;
        dim = 1 << n;
        this.bonds = bonds ?? Topology.Chain(n);
        rho = new Complex[dim, dim];
        mask = new double[dim, dim];
        dis = new int[dim, dim];
        for (int i = 0; i < dim; i++)
            for (int jj = 0; jj < dim; jj++)
            {
                var p = new Pair(world, i, jj, gamma);   // the atom: its own disagreement and its own rate
                mask[i, jj] = antiWatching ? -2.0 * gamma * (n - p.Disagreement) : p.Rate;
                dis[i, jj] = p.Disagreement;
            }
        pc = new int[dim];
        for (int i = 0; i < dim; i++) pc[i] = BitOperations.PopCount((uint)i);
        h = BuildHandshake();
    }

    public int Dim => dim;
    public Complex this[int i, int j] => rho[i, j];

    // seed a population |s><s| with weight amp (structure: one possibility held, no novelty).
    public void Seed(int s, double amp = 1.0) { rho[s, s] = amp; occupied.Add((pc[s], pc[s])); alive = null; }

    // seed one coherence |i><j| (+ its Hermitian twin), real amplitude.
    public void SeedCoherence(int i, int j, double amp)
    {
        rho[i, j] = amp; rho[j, i] = amp;
        occupied.Add((pc[i], pc[j])); occupied.Add((pc[j], pc[i])); alive = null;
    }

    // cut (c): the loop stays in the seed's joint-popcount blocks (F63, [L,Pi^2]=0); the rest is forbidden,
    // never run. Alive = the cells whose (popcount row, popcount col) block the seed occupies.
    (int i, int j)[] Alive => alive ??= BuildAlive();
    (int i, int j)[] BuildAlive()
    {
        var list = new List<(int, int)>();
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                if (occupied.Contains((pc[i], pc[j]))) list.Add((i, j));
        return list.ToArray();
    }
    public int AliveCount => Alive.Length;
    public int ForbiddenCount => dim * dim - AliveCount;

    // the handshake: H = J * sum over the bonds (sigma+_a sigma-_b + h.c.), the flip-flop that hops one
    // excitation between joined sites -- the smallest inner motion that turns a population into a coherence.
    double[,] BuildHandshake()
    {
        var H = new double[dim, dim];
        foreach (var (a, b) in bonds)
            for (int s = 0; s < dim; s++)
                if (((s >> a) & 1) == 1 && ((s >> b) & 1) == 0)
                {
                    int s2 = (s & ~(1 << a)) | (1 << b);   // move the excitation a -> b
                    H[s, s2] += J;
                    H[s2, s] += J;
                }
        return H;
    }

    // rho-dot = -i[H,rho] + D[rho]: the restlessness (commutator) plus the watching (element-wise mask).
    Complex[,] Rhs(Complex[,] x)
    {
        var r = new Complex[dim, dim];
        foreach (var (i, j) in Alive)               // cut (c): only the occupied blocks, the rest forbidden
        {
            Complex hx = Complex.Zero, xh = Complex.Zero;
            for (int m = 0; m < dim; m++)
            {
                hx += h[i, m] * x[m, j];
                xh += x[i, m] * h[m, j];
            }
            r[i, j] = -Complex.ImaginaryOne * (hx - xh) + mask[i, j] * x[i, j];
        }
        return r;
    }

    // one RK4 tick of the living loop.
    public void Step(double dt)
    {
        var k1 = Rhs(rho);
        var k2 = Rhs(Axpy(rho, k1, dt / 2));
        var k3 = Rhs(Axpy(rho, k2, dt / 2));
        var k4 = Rhs(Axpy(rho, k3, dt));
        foreach (var (i, j) in Alive)
            rho[i, j] += (dt / 6) * (k1[i, j] + 2 * k2[i, j] + 2 * k3[i, j] + k4[i, j]);
        T += dt;
    }

    Complex[,] Axpy(Complex[,] a, Complex[,] b, double s)
    {
        var r = new Complex[dim, dim];
        foreach (var (i, j) in Alive) r[i, j] = a[i, j] + s * b[i, j];
        return r;
    }

    // the magnitude carried at each disagreement k = 0..N (k=0 the structure/populations, k>0 the novelty).
    public double[] WeightByDisagreement()
    {
        var byK = new double[N + 1];
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                byK[dis[i, j]] += rho[i, j].Magnitude;
        return byK;
    }

    // left: structure (the trace, the populations that persist) and novelty (the off-diagonal coherence).
    public double Structure { get { double s = 0; for (int i = 0; i < dim; i++) s += rho[i, i].Real; return s; } }
    public double Novelty => WeightByDisagreement().Skip(1).Sum();

    // the anti-trace: the anti-diagonal sum rho(i, ~i), k = N everywhere. In the normal world it is the
    // fastest-dying content; in the anti-watched world it is the conserved one (the trace's mirror twin).
    public double AntiStructure
    {
        get
        {
            double s = 0;
            for (int i = 0; i < dim; i++) s += rho[i, dim - 1 - i].Real;
            return s;
        }
    }

    public override IReadOnlyList<string> Own => new[] { "structure", "novelty", "antistructure" };
}
