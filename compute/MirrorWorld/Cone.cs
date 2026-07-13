using System.Numerics;

namespace MirrorWorld;

// The memory cut: the living world stored on only its block, not the full 4^N. For a single excitation the
// block is the N singly-occupied sites, so rho is N x N (not 2^N x 2^N) -- the same dynamics Restless runs
// in block (1,1), but in O(N^2) memory instead of O(4^N). This is what unblocks large N: at N=100 a single
// excitation lives in 10^4 cells, where the full Liouvillian (and even Restless's dense rho) is impossible.
// rho-dot = -i[H,rho] + D[rho]: H the tight-binding hop along the bonds; D the uniform -4*gamma on every
// coherence (two distinct single-excitation states always disagree in exactly 2 bits). The light-cone made
// reachable -- a dynamical, spatial picture the static spectrum cannot show, at an N the spectrum cannot hold.
public sealed class Cone : GameObject
{
    public int N { get; }
    public double J { get; }
    public double Gamma { get; }
    public double T { get; private set; }

    readonly Complex[,] rho;        // N x N: rho[a,b] = coherence between excitation-at-a and excitation-at-b
    readonly double[,] hop;         // the tight-binding H on the sites (J on a bond, else 0)
    int[][]? nbr;                   // per site: the bonded sites (rebuilt from hop when a bond changes)
    double[][]? nbrJ;               // per site: the matching couplings

    public Cone(World world, int n, double j, double gamma, (int a, int b)[]? bonds = null) : base(world)
    {
        N = n;
        J = j;
        Gamma = gamma;
        rho = new Complex[n, n];
        hop = new double[n, n];
        foreach (var (a, b) in bonds ?? Topology.Chain(n)) { hop[a, b] = j; hop[b, a] = j; }
    }

    // re-tune one bond's coupling (the defect knob of the walk-time step); symmetric. Configuration: set it
    // before seeding and running, the way the constructor's bonds are.
    public void SetBond(int a, int b, double j) { hop[a, b] = j; hop[b, a] = j; nbr = null; }

    // H is sparse (a chain has two neighbours); walk only the bonds instead of all N columns.
    void RebuildNeighbours()
    {
        nbr = new int[N][]; nbrJ = new double[N][];
        for (int a = 0; a < N; a++)
        {
            var sites = new List<int>();
            for (int m = 0; m < N; m++) if (hop[a, m] != 0.0) sites.Add(m);
            nbr[a] = sites.ToArray();
            nbrJ[a] = sites.Select(m => hop[a, m]).ToArray();
        }
    }

    public int Sites => N;
    public Complex this[int a, int b] => rho[a, b];

    // the excitation starts at one site (a population, no novelty).
    public void Seed(int site) => rho[site, site] = Complex.One;

    // rho-dot = -i[H,rho] + D[rho]; D = -4 gamma on every off-diagonal (disagreement 2), 0 on the diagonal.
    Complex[,] Rhs(Complex[,] x)
    {
        if (nbr is null) RebuildNeighbours();
        var r = new Complex[N, N];
        for (int a = 0; a < N; a++)
            for (int b = 0; b < N; b++)
            {
                Complex hx = Complex.Zero, xh = Complex.Zero;
                var na = nbr![a]; var ja = nbrJ![a];
                for (int i = 0; i < na.Length; i++) hx += ja[i] * x[na[i], b];
                var nb = nbr[b]; var jb = nbrJ[b];
                for (int i = 0; i < nb.Length; i++) xh += x[a, nb[i]] * jb[i];   // hop symmetric: hop[m,b] = hop[b,m]
                double deph = (a == b) ? 0.0 : -4.0 * Gamma;
                r[a, b] = -Complex.ImaginaryOne * (hx - xh) + deph * x[a, b];
            }
        return r;
    }

    // one RK4 tick.
    public void Step(double dt)
    {
        var k1 = Rhs(rho);
        var k2 = Rhs(Axpy(rho, k1, dt / 2));
        var k3 = Rhs(Axpy(rho, k2, dt / 2));
        var k4 = Rhs(Axpy(rho, k3, dt));
        for (int a = 0; a < N; a++)
            for (int b = 0; b < N; b++)
                rho[a, b] += (dt / 6) * (k1[a, b] + 2 * k2[a, b] + 2 * k3[a, b] + k4[a, b]);
        T += dt;
    }

    Complex[,] Axpy(Complex[,] u, Complex[,] v, double s)
    {
        var r = new Complex[N, N];
        for (int a = 0; a < N; a++)
            for (int b = 0; b < N; b++)
                r[a, b] = u[a, b] + s * v[a, b];
        return r;
    }

    public double Population(int site) => rho[site, site].Real;          // the excitation density at a site

    // left: structure (the trace, the one excitation that persists) and novelty (the spread-out coherence).
    public double Structure { get { double s = 0; for (int a = 0; a < N; a++) s += rho[a, a].Real; return s; } }
    public double Novelty
    {
        get
        {
            double s = 0;
            for (int a = 0; a < N; a++)
                for (int b = 0; b < N; b++)
                    if (a != b) s += rho[a, b].Magnitude;
            return s;
        }
    }

    public override IReadOnlyList<string> Own => new[] { "structure", "novelty" };
}
