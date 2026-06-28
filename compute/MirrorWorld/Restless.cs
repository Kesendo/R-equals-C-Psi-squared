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

    public Restless(World world, int n, double j, double gamma) : base(world)
    {
        N = n;
        J = j;
        Gamma = gamma;
        dim = 1 << n;
        rho = new Complex[dim, dim];
        mask = new double[dim, dim];
        dis = new int[dim, dim];
        for (int i = 0; i < dim; i++)
            for (int jj = 0; jj < dim; jj++)
            {
                var p = new Pair(world, i, jj, gamma);   // the atom: its own disagreement and its own rate
                mask[i, jj] = p.Rate;
                dis[i, jj] = p.Disagreement;
            }
        h = BuildHandshake();
    }

    public int Dim => dim;
    public Complex this[int i, int j] => rho[i, j];

    // seed a pure population |s><s| (structure: one possibility held, no novelty).
    public void Seed(int s) => rho[s, s] = Complex.One;

    // seed one coherence |i><j| (+ its Hermitian twin), real amplitude.
    public void SeedCoherence(int i, int j, double amp) { rho[i, j] = amp; rho[j, i] = amp; }

    // the handshake: H = J * sum over chain bonds (sigma+_a sigma-_b + h.c.), the flip-flop that hops one
    // excitation between neighbours -- the smallest inner motion that turns a population into a coherence.
    double[,] BuildHandshake()
    {
        var H = new double[dim, dim];
        for (int a = 0; a + 1 < N; a++)
        {
            int b = a + 1;
            for (int s = 0; s < dim; s++)
                if (((s >> a) & 1) == 1 && ((s >> b) & 1) == 0)
                {
                    int s2 = (s & ~(1 << a)) | (1 << b);   // move the excitation a -> b
                    H[s, s2] += J;
                    H[s2, s] += J;
                }
        }
        return H;
    }

    // rho-dot = -i[H,rho] + D[rho]: the restlessness (commutator) plus the watching (element-wise mask).
    Complex[,] Rhs(Complex[,] x)
    {
        var r = new Complex[dim, dim];
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
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
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                rho[i, j] += (dt / 6) * (k1[i, j] + 2 * k2[i, j] + 2 * k3[i, j] + k4[i, j]);
        T += dt;
    }

    static Complex[,] Axpy(Complex[,] a, Complex[,] b, double s)
    {
        int d = a.GetLength(0);
        var r = new Complex[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                r[i, j] = a[i, j] + s * b[i, j];
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

    public override IReadOnlyList<string> Own => new[] { "structure", "novelty" };
}
