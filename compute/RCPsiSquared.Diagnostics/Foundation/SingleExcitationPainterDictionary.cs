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

    private double[][] Purities(int bond, double dj, Complex[] r0, double[] t)
    {
        int n = N, d = n * n;
        var L = ComplexMatrix.Build.Dense(d, d);
        // −i[h, r] with h the hopping matrix (amplitude J, J + δJ on the defect bond); row-major vec(r)
        for (int b = 0; b < n - 1; b++)
        {
            double hop = J + (b == bond ? dj : 0.0);
            foreach (var (x, y) in new[] { (b, b + 1), (b + 1, b) })
                for (int k = 0; k < n; k++)
                {
                    L[x * n + k, y * n + k] += -Complex.ImaginaryOne * hop;   // −i·h_xy·r_yk lands on r_xk
                    L[k * n + y, k * n + x] += Complex.ImaginaryOne * hop;    // +i·r_kx·h_xy lands on r_ky
                }
        }
        for (int a = 0; a < n; a++)
            for (int c = 0; c < n; c++)
                if (a != c) L[a * n + c, a * n + c] += -4.0 * Gamma;
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
            for (int i = 0; i < n; i++)
            {
                double pop = r[i * n + i].Real;
                p[i][s] = 0.5 * (1.0 + (1.0 - 2.0 * pop) * (1.0 - 2.0 * pop));
            }
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
