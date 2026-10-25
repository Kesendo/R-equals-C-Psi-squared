using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The decoder's α-dictionary read in the single-excitation block. The canonical painter protocol
/// (open XY chain, H = (J/2)·Σ(XX + YY), Z dephasing γ(ZρZ − ρ), the bonding mode ψ₁ as carrier, a bond
/// defect δJ·(XX + YY)/2) never leaves the block |i⟩⟨j| of one excitation, and the painters read only the
/// per-site purity P_i = ½(1 + ⟨Z_i⟩²) = ½(1 + (1 − 2p_i)²) from its diagonal. So the N × N block carries
/// the whole reading: hopping J on the bonds, every coherence |i⟩⟨j| (i ≠ j) charged 4γ, the populations
/// p_i = r_ii. On that block this class runs the painters' own steps: the clean, defected and guard
/// trajectories on <see cref="PaintersMovement"/>'s time grid (t ∈ [0, 1/γ], 60 points), the same
/// <see cref="PaintersMovement.FitAlpha"/>, f_i = (α_i − 1)/δJ, and the same reliability rule (|f| ≤ 10 and
/// the guard at δJ/2 agreeing). It reaches N far past the 4^N pipeline's ceiling (N ≤ 6); its Python
/// predecessor, which first read the painter on this block, is simulations/handshake_rk_block.py.
///
/// <para>Readings: the per-bond letters F[b]; the chain mirror's odd weight w(b) of each letter (the share
/// of ‖f^b‖² that is odd under i ↦ N − 1 − i); the worst signed cosine over all sites and over the sites
/// reliable for every bond. Since the chain, the watching and ψ₁ are mirror-symmetric, the letter of bond
/// N − 2 − b is the mirror image of the letter of bond b (exactly for the protocol, to the fit's stopping
/// width in the code), so cos(f^b, f^{N−2−b}) = 1 − 2·w(b) for every chain-mirror pair, and at even N the
/// centre bond's letter is mirror-even, w = 0.</para></summary>
public sealed class SingleExcitationPainterDictionary
{
    public const double DefaultDeltaJ = 0.02;
    public const int DefaultTimePoints = 60;

    public int N { get; }
    public double J { get; }
    public double Gamma { get; }
    public double DeltaJ { get; }
    public IReadOnlyList<double> TimeGrid { get; }

    /// <summary>F[b][i] = (α_i − 1)/δJ for a defect on bond b.</summary>
    public double[][] F { get; }
    /// <summary>The painters' reliability flag per bond and site.</summary>
    public bool[][] Reliable { get; }

    public SingleExcitationPainterDictionary(int n, double j = 1.0, double gamma = 0.05,
                                             double deltaJ = DefaultDeltaJ, int timePoints = DefaultTimePoints)
    {
        // the generator is a dense N² × N² matrix with one eigendecomposition per trajectory; N = 24 is 576²
        if (n < 3 || n > 24) throw new ArgumentOutOfRangeException(nameof(n), n, "N must lie in 3..24");
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be positive");
        if (deltaJ <= 0) throw new ArgumentOutOfRangeException(nameof(deltaJ), deltaJ, "δJ must be positive");
        if (timePoints < 2) throw new ArgumentOutOfRangeException(nameof(timePoints), timePoints, "at least two time points");
        N = n; J = j; Gamma = gamma; DeltaJ = deltaJ;
        var t = new double[timePoints];
        for (int s = 0; s < timePoints; s++) t[s] = (1.0 / gamma) * s / (timePoints - 1);
        TimeGrid = t;

        var psi = new double[n];
        double norm = 0;
        for (int l = 0; l < n; l++) { psi[l] = Math.Sin(Math.PI * (l + 1) / (n + 1)); norm += psi[l] * psi[l]; }
        var r0 = new Complex[n * n];
        for (int a = 0; a < n; a++)
            for (int b = 0; b < n; b++) r0[a * n + b] = psi[a] * psi[b] / norm;

        _r0 = r0;
        _clean = new Lazy<CleanBasis>(BuildClean);
        _cleanOnGrid = new Lazy<(double[][] P, double[][] Rate)>(() => CleanPopulations(TimeGrid));
        _carrierLetters = new Lazy<double[][]>(BuildCarrierLetters);
        _linearF = new Lazy<double[][]>(() =>
            LinearLetters(Enumerable.Range(0, n - 1).Select(b => BondOperator(n, b)).ToArray()));
        var pA = Purities(-1, 0.0, r0, t);
        F = new double[n - 1][];
        Reliable = new bool[n - 1][];
        for (int b = 0; b < n - 1; b++)
        {
            var pB = Purities(b, deltaJ, r0, t);
            var pG = Purities(b, deltaJ / 2, r0, t);
            F[b] = new double[n];
            Reliable[b] = new bool[n];
            for (int i = 0; i < n; i++)
            {
                double f = (PaintersMovement.FitAlpha(t, pA[i], pB[i]) - 1.0) / deltaJ;
                double fg = (PaintersMovement.FitAlpha(t, pA[i], pG[i]) - 1.0) / (deltaJ / 2);
                F[b][i] = f;
                Reliable[b][i] = Math.Abs(f) <= PaintersMovement.FMax
                              && Math.Abs(f - fg) <= PaintersMovement.ConsistencyTol * (Math.Abs(fg) + PaintersMovement.FFloor);
            }
        }
    }

    private readonly Complex[] _r0;
    private readonly Lazy<double[][]> _linearF;
    private readonly Lazy<CleanBasis> _clean;
    private readonly Lazy<(double[][] P, double[][] Rate)> _cleanOnGrid;
    private readonly Lazy<double[][]> _carrierLetters;

    /// <summary>The clean generator's eigenbasis: L = V·diag(λ)·V⁻¹, and c = V⁻¹·vec(r₀).</summary>
    private sealed record CleanBasis(Complex[] Lambda, ComplexMatrix V,
                                     MathNet.Numerics.LinearAlgebra.Factorization.LU<Complex> Lu, Complex[] C);

    /// <summary>−i[v, ·] on row-major vec(r) for a real symmetric site operator v.</summary>
    private static ComplexMatrix Commutator(double[,] v)
    {
        int n = v.GetLength(0), d = n * n;
        var L = ComplexMatrix.Build.Dense(d, d);
        for (int x = 0; x < n; x++)
            for (int y = 0; y < n; y++)
            {
                if (v[x, y] == 0.0) continue;
                for (int k = 0; k < n; k++)
                {
                    L[x * n + k, y * n + k] += -Complex.ImaginaryOne * v[x, y];   // −i·v_xy·r_yk lands on r_xk
                    L[k * n + y, k * n + x] += Complex.ImaginaryOne * v[x, y];    // +i·r_kx·v_xy lands on r_ky
                }
            }
        return L;
    }

    /// <summary>The hopping (XX + YY)/2 of one bond as a site operator: 1 on (b, b + 1) and (b + 1, b).</summary>
    public static double[,] BondOperator(int n, int bond)
    {
        if (bond < 0 || bond > n - 2) throw new ArgumentOutOfRangeException(nameof(bond), bond, "bond must lie in 0..N−2");
        var v = new double[n, n];
        v[bond, bond + 1] = v[bond + 1, bond] = 1.0;
        return v;
    }

    /// <summary>The angle π·m/q with the integer m reduced modulo 2q first, so it lies in [0, 2π) at every N and a
    /// sine or cosine of it carries a few ε.</summary>
    private static double PiFraction(long m, long q) => Math.PI * (((m % (2 * q)) + 2 * q) % (2 * q)) / q;

    /// <summary>The clean chain's mode ψ_k(l) = √(2/(N+1))·sin(πk(l+1)/(N+1)), k = 1..N; ψ₁ is the carrier.</summary>
    public static double[] Mode(int n, int k)
    {
        if (k < 1 || k > n) throw new ArgumentOutOfRangeException(nameof(k), k, "k must lie in 1..N");
        var m = new double[n];
        for (int l = 0; l < n; l++) m[l] = Math.Sqrt(2.0 / (n + 1)) * Math.Sin(PiFraction((long)k * (l + 1), n + 1));
        return m;
    }

    /// <summary>The carrier coupling c_k(b) = ⟨ψ_k|V_b|ψ₁⟩ = ψ_k(b)ψ₁(b+1) + ψ_k(b+1)ψ₁(b) (k = 2..N F124's location
    /// dictionary M[b, k], k = 1 the strength), in closed form with x = π/(N+1):
    /// (2/(N+1))·[cos((k+1)x/2)·cos((k−1)(b+3/2)x) − cos((k−1)x/2)·cos((k+1)(b+3/2)x)]. At k = N both terms
    /// vanish (cos((N+1)x/2) = 0 and cos((N+1)(b+3/2)x) = 0), the K-partner rule.</summary>
    public static double CarrierCoupling(int n, int k, int bond)
    {
        if (k < 1 || k > n) throw new ArgumentOutOfRangeException(nameof(k), k, "k must lie in 1..N");
        if (bond < 0 || bond > n - 2) throw new ArgumentOutOfRangeException(nameof(bond), bond, "bond must lie in 0..N−2");
        long q = 2L * (n + 1);   // x/2 = π/q and (b + 3/2)·x = π·(2b + 3)/q
        return 2.0 / (n + 1) * (Math.Cos(PiFraction(k + 1, q)) * Math.Cos(PiFraction((long)(k - 1) * (2 * bond + 3), q))
                              - Math.Cos(PiFraction(k - 1, q)) * Math.Cos(PiFraction((long)(k + 1) * (2 * bond + 3), q)));
    }

    private ComplexMatrix Generator(int bond, double dj)
    {
        int n = N;
        var h = new double[n, n];
        for (int b = 0; b < n - 1; b++)
            h[b, b + 1] = h[b + 1, b] = J + (b == bond ? dj : 0.0);
        var L = Commutator(h);
        for (int a = 0; a < n; a++)
            for (int c = 0; c < n; c++)
                if (a != c) L[a * n + c, a * n + c] += -4.0 * Gamma;
        return L;
    }

    private CleanBasis BuildClean()
    {
        var evd = Generator(-1, 0.0).Evd();
        var lu = evd.EigenVectors.LU();
        return new CleanBasis(evd.EigenValues.ToArray(), evd.EigenVectors, lu,
                              lu.Solve(ComplexVector.Build.DenseOfArray(_r0)).ToArray());
    }

    /// <summary>The clean populations p_i(t) and their rates ṗ_i(t) at arbitrary times, from the eigenvectors.</summary>
    internal (double[][] P, double[][] Rate) CleanPopulations(IReadOnlyList<double> times)
    {
        var cb = _clean.Value;
        int n = N, d = n * n;
        var p = new double[n][]; var rate = new double[n][];
        for (int i = 0; i < n; i++) { p[i] = new double[times.Count]; rate[i] = new double[times.Count]; }
        for (int s = 0; s < times.Count; s++)
            for (int i = 0; i < n; i++)
            {
                Complex v = 0, dv = 0;
                for (int k = 0; k < d; k++)
                {
                    var e = cb.V[i * n + i, k] * Complex.Exp(cb.Lambda[k] * times[s]) * cb.C[k];
                    v += e; dv += cb.Lambda[k] * e;
                }
                p[i][s] = v.Real; rate[i][s] = dv.Real;
            }
        return (p, rate);
    }

    /// <summary>φ₁(z) = (e^z − 1)/z: its series through z⁴/120 where |z| &lt; 1e-3 (truncation about |z|⁵/720,
    /// 1e-18), the quotient above (cancellation error about ε/|z| ≤ 2e-13).</summary>
    private static Complex Phi1(Complex z) =>
        Complex.Abs(z) < 1e-3 ? 1 + z * (1.0 / 2 + z * (1.0 / 6 + z * (1.0 / 24 + z / 120))) : (Complex.Exp(z) - 1) / z;

    /// <summary>The linear letter of a perturbation v (a real symmetric site operator, h ↦ h + ε·v), with no fit:
    /// f_i = ⟨g_i, ∂P_i⟩/⟨g_i, g_i⟩ over the time grid, where g_i(t) = t·Ṗ_i(t) is the clean purity's response
    /// to a time stretch and ∂P_i = dP_i/dε at ε = 0, exact through the eigenbasis (the integral
    /// ∫₀ᵗ e^{λ_k(t−s)} e^{λ_l s} ds = t·e^{λ_l t}·φ₁((λ_k − λ_l)t)). It is the letter the α fit returns as δJ → 0
    /// when the fit stretches the exact clean purity. The painters' fit is not that fit: it stretches a spline clamped
    /// to its end value beyond the grid, so for α &gt; 1 the last, heaviest point stops moving; its small-δJ letters meet
    /// this one to about 1 %, plus the spline's derivative error.</summary>
    public double[] LinearLetter(double[,] v) => LinearLetters(new[] { v })[0];

    /// <summary>The linear letters of several perturbations, sharing the divided differences of each time point.</summary>
    public double[][] LinearLetters(IReadOnlyList<double[,]> vs)
    {
        foreach (var v in vs)
        {
            if (v.GetLength(0) != N || v.GetLength(1) != N) throw new ArgumentException("v must be N × N", nameof(vs));
            for (int x = 0; x < N; x++)
                for (int w = 0; w < x; w++)
                    if (v[x, w] != v[w, x]) throw new ArgumentException("v must be symmetric", nameof(vs));
        }
        var cb = _clean.Value;
        var lambda = cb.Lambda;
        int n = N, d = n * n, m = vs.Count;
        var bc = new Complex[m][,];                        // (V⁻¹ B V)_kl · c_l per perturbation
        for (int q = 0; q < m; q++)
        {
            var btm = cb.Lu.Solve(Commutator(vs[q]) * cb.V);
            bc[q] = new Complex[d, d];
            for (int k = 0; k < d; k++)
                for (int l = 0; l < d; l++) bc[q][k, l] = btm[k, l] * cb.C[l];
        }
        var diagRows = new Complex[n, d];
        for (int i = 0; i < n; i++)
            for (int k = 0; k < d; k++) diagRows[i, k] = cb.V[i * n + i, k];
        var (p, rate) = _cleanOnGrid.Value;
        var num = new double[m, n]; var den = new double[n];
        var e = new Complex[d];
        var dd = new Complex[d, d];
        var y = new Complex[d];
        var g = new double[n]; var slope = new double[n];
        for (int s = 0; s < TimeGrid.Count; s++)
        {
            double t = TimeGrid[s];
            for (int k = 0; k < d; k++) e[k] = Complex.Exp(lambda[k] * t);
            // the divided differences (e^{λ_k t} − e^{λ_l t})/(λ_k − λ_l), shared by every perturbation
            for (int k = 0; k < d; k++)
                for (int l = 0; l < d; l++)
                {
                    var z = (lambda[k] - lambda[l]) * t;
                    dd[k, l] = Complex.Abs(z) < 1e-3 ? t * e[l] * Phi1(z) : (e[k] - e[l]) / (lambda[k] - lambda[l]);
                }
            for (int i = 0; i < n; i++)
            {
                slope[i] = -2.0 * (1.0 - 2.0 * p[i][s]);   // dP/dp
                g[i] = t * slope[i] * rate[i][s];
                den[i] += g[i] * g[i];
            }
            for (int q = 0; q < m; q++)
            {
                var b = bc[q];
                for (int k = 0; k < d; k++)
                {
                    Complex acc = 0;
                    for (int l = 0; l < d; l++) acc += dd[k, l] * b[k, l];
                    y[k] = acc;
                }
                for (int i = 0; i < n; i++)
                {
                    Complex dp = 0;
                    for (int k = 0; k < d; k++) dp += diagRows[i, k] * y[k];
                    num[q, i] += g[i] * slope[i] * dp.Real;
                }
            }
        }
        var f = new double[m][];
        for (int q = 0; q < m; q++)
        {
            f[q] = new double[n];
            for (int i = 0; i < n; i++) f[q][i] = num[q, i] / den[i];
        }
        return f;
    }

    /// <summary>LinearF[b] = the linear letter of bond b, the δJ → 0 limit of the exactly stretched fit.</summary>
    public double[][] LinearF => _linearF.Value;

    /// <summary>The letter u_k of the carrier channel (1, k): the linear letter of |ψ₁⟩⟨ψ_k| + |ψ_k⟩⟨ψ₁| for
    /// k ≥ 2 and of |ψ₁⟩⟨ψ₁| for k = 1. A bond's letter is Σ_k c_k(b)·u_k plus the channels that do not touch
    /// the carrier.</summary>
    public double[] CarrierChannelLetter(int k)
    {
        if (k < 1 || k > N) throw new ArgumentOutOfRangeException(nameof(k), k, "k must lie in 1..N");
        return _carrierLetters.Value[k - 1];
    }

    private double[][] BuildCarrierLetters()
    {
        var a = Mode(N, 1);
        var vs = new List<double[,]>();
        for (int k = 1; k <= N; k++)
        {
            var b = Mode(N, k);
            var v = new double[N, N];
            for (int x = 0; x < N; x++)
                for (int z = 0; z < N; z++)
                    v[x, z] = k == 1 ? a[x] * a[z] : a[x] * b[z] + b[x] * a[z];
            vs.Add(v);
        }
        return LinearLetters(vs);
    }

    /// <summary>The per-site populations on the time grid with bond's coupling J + δJ (bond −1: the clean chain).</summary>
    internal double[][] DefectPopulations(int bond, double dj)
    {
        if (bond < -1 || bond > N - 2) throw new ArgumentOutOfRangeException(nameof(bond), bond, "bond must lie in −1..N−2");
        return Populations(bond, dj, _r0, TimeGrid.ToArray());
    }

    private double[][] Purities(int bond, double dj, Complex[] r0, double[] t) =>
        Populations(bond, dj, r0, t)
            .Select(row => row.Select(pop => 0.5 * (1.0 + (1.0 - 2.0 * pop) * (1.0 - 2.0 * pop))).ToArray()).ToArray();

    private double[][] Populations(int bond, double dj, Complex[] r0, double[] t)
    {
        int n = N, d = n * n;
        var L = Generator(bond, dj);
        var evd = L.Evd();
        var V = evd.EigenVectors;
        var c0 = V.Solve(ComplexVector.Build.DenseOfArray(r0));
        var lambda = evd.EigenValues;
        var p = new double[n][];
        for (int i = 0; i < n; i++) p[i] = new double[t.Length];
        var e = ComplexVector.Build.Dense(d);
        for (int s = 0; s < t.Length; s++)
        {
            for (int k = 0; k < d; k++) e[k] = Complex.Exp(lambda[k] * t[s]) * c0[k];
            var r = V * e;
            for (int i = 0; i < n; i++) p[i][s] = r[i * n + i].Real;
        }
        return p;
    }

    /// <summary>The mirror-odd weight w(b) = ‖(f − Rf)/2‖² / ‖f‖² of bond b's letter, R: i ↦ N − 1 − i.</summary>
    public double MirrorOddWeight(int bond)
    {
        var f = F[bond];
        double odd = 0, all = 0;
        for (int i = 0; i < N; i++) { double o = (f[i] - f[N - 1 - i]) / 2; odd += o * o; all += f[i] * f[i]; }
        return odd / all;
    }

    public static double Cos(double[] a, double[] b, bool[]? mask = null)
    {
        double dot = 0, na = 0, nb = 0;
        for (int i = 0; i < a.Length; i++)
        {
            if (mask != null && !mask[i]) continue;
            dot += a[i] * b[i]; na += a[i] * a[i]; nb += b[i] * b[i];
        }
        return dot / Math.Sqrt(na * nb);
    }

    /// <summary>The sites reliable for every bond.</summary>
    public bool[] SitesReliableForEveryBond()
    {
        var m = new bool[N];
        for (int i = 0; i < N; i++) m[i] = Enumerable.Range(0, N - 1).All(b => Reliable[b][i]);
        return m;
    }

    /// <summary>The most anti-collinear pair: the minimum signed cos over all sites (reliableOnly false) or
    /// over the sites reliable for every bond.</summary>
    public (double Cos, int A, int B) WorstAntiCollinearPair(bool reliableOnly)
    {
        var mask = reliableOnly ? SitesReliableForEveryBond() : null;
        (double Cos, int A, int B) worst = (double.PositiveInfinity, -1, -1);
        for (int a = 0; a < N - 1; a++)
            for (int b = a + 1; b < N - 1; b++)
            {
                double c = Cos(F[a], F[b], mask);
                if (c < worst.Cos) worst = (c, a, b);
            }
        return worst;
    }
}
