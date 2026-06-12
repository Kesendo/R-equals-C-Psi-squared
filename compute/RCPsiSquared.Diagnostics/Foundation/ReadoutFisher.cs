using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Lindblad;
using CVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

public enum ReadoutBasis { Z, X, Y }

/// <summary>Classical Fisher information of a basis readout along a fixed-K trajectory,
/// by forward difference in δJ. Pure and self-contained: builds its own trajectories from
/// the same Core pieces the Symphony's painters use; never reaches into Symphony internals.</summary>
public static class ReadoutFisher
{
    static readonly Matrix<Complex> PX = Matrix<Complex>.Build.DenseOfArray(
        new Complex[,] { { 0, 1 }, { 1, 0 } });
    static readonly Matrix<Complex> PY = Matrix<Complex>.Build.DenseOfArray(
        new Complex[,] { { 0, new Complex(0, -1) }, { new Complex(0, 1), 0 } });
    static readonly Matrix<Complex> I2 = Matrix<Complex>.Build.DenseIdentity(2);
    static readonly Matrix<Complex> Had = Matrix<Complex>.Build.DenseOfArray(
        new Complex[,] { { 1, 1 }, { 1, -1 } }) / Math.Sqrt(2.0);
    static readonly Matrix<Complex> Ymeas = Had * Matrix<Complex>.Build.DenseOfDiagonalArray(
        new Complex[] { 1, new Complex(0, -1) });

    static Matrix<Complex> Site(Matrix<Complex> p, int site, int n)
    {
        var m = Matrix<Complex>.Build.DenseIdentity(1);
        for (int s = 0; s < n; s++) m = m.KroneckerProduct(s == site ? p : I2);
        return m;
    }

    /// <summary>The k=1 sine mode Σ_l sin(π(l+1)/(N+1))|1_l⟩, site 0 = MSB (the F67
    /// bonding receiver, the canonical decoder carrier).</summary>
    public static CVector BondingState(int n)
    {
        int d = 1 << n;
        var psi = CVector.Build.Dense(d);
        double norm = 0;
        for (int l = 0; l < n; l++)
        { double a = Math.Sin(Math.PI * (l + 1) / (n + 1)); psi[1 << (n - 1 - l)] = a; norm += a * a; }
        return psi / Math.Sqrt(norm);
    }

    /// <summary>XY chain H = Σ_b J_b/2 (X_bX_{b+1}+Y_bY_{b+1}), one optional defect bond at J+δJ.</summary>
    public static Matrix<Complex> XyChain(int n, double j, int? defectBond, double deltaJ)
    {
        int d = 1 << n;
        var h = Matrix<Complex>.Build.Dense(d, d);
        for (int b = 0; b < n - 1; b++)
        {
            double jb = j + (defectBond == b ? deltaJ : 0.0);
            h += (jb / 2.0) * (Site(PX, b, n) * Site(PX, b + 1, n)
                             + Site(PY, b, n) * Site(PY, b + 1, n));
        }
        return h;
    }

    public static double[] KGrid(double gamma, double kMax, int points)
    {
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma));
        if (kMax <= 0) throw new ArgumentOutOfRangeException(nameof(kMax));
        if (points < 1) throw new ArgumentOutOfRangeException(nameof(points));
        var ts = new double[points];
        for (int i = 0; i < points; i++)
            ts[i] = (kMax / gamma) * (i + 1) / points;   // (0, kMax/γ], no t=0 point
        return ts;
    }

    /// <summary>ρ(t) on the grid: ONE eigendecomposition of the Z-dephasing Liouvillian,
    /// unrolled. Uniform γ per site.</summary>
    public static IReadOnlyList<Matrix<Complex>> Trajectory(
        int n, double j, double gamma, int? defectBond, double deltaJ, double[] times)
    {
        int d = 1 << n;
        var h = XyChain(n, j, defectBond, deltaJ);
        var gammas = Enumerable.Repeat(gamma, n).ToArray();
        var L = PauliDephasingDissipator.BuildZ(h, gammas);   // d²×d², vec[a·d+b] = ρ[a,b]
        var evd = L.Evd();
        var R = evd.EigenVectors; var Rinv = R.Inverse();
        var psi = BondingState(n);
        var rho0 = CVector.Build.Dense(d * d);
        for (int a = 0; a < d; a++) for (int b = 0; b < d; b++)
            rho0[a * d + b] = psi[a] * Complex.Conjugate(psi[b]);
        var c0 = Rinv * rho0;
        var outRhos = new List<Matrix<Complex>>(times.Length);
        foreach (var t in times)
        {
            var ct = CVector.Build.Dense(d * d,
                i => c0[i] * Complex.Exp(evd.EigenValues[i] * t));
            var vec = R * ct;
            var rho = Matrix<Complex>.Build.Dense(d, d, (a, b) => vec[a * d + b]);
            outRhos.Add(rho);
        }
        return outRhos;
    }

    /// <summary>The n-fold per-site rotation mapping the given Pauli basis to Z, or
    /// null for Z (identity, hoisted out of the time loop by the FI core).</summary>
    static Matrix<Complex>? BasisRotation(int n, ReadoutBasis basis)
    {
        if (basis == ReadoutBasis.Z) return null;
        Matrix<Complex> u1 = basis == ReadoutBasis.X ? Had : Ymeas;
        var u = Matrix<Complex>.Build.DenseIdentity(1);
        for (int s = 0; s < n; s++) u = u.KroneckerProduct(u1);
        return u;
    }

    /// <summary>Per-outcome probabilities given a precomputed basis rotation (null = Z):
    /// p = diag(U ρ U†), normalized.</summary>
    static double[] ProbsRotated(Matrix<Complex> rho, Matrix<Complex>? u)
    {
        int d = rho.RowCount;
        var rr = u is null ? rho : u * rho * u.ConjugateTranspose();
        var p = new double[d]; double sum = 0;
        // 1e-12 floor: structurally-zero outcomes cancel in the forward difference; FI is
        // invariant to this value across 1e-10..1e-14 (verified 2026-06-12) - do not tune.
        for (int i = 0; i < d; i++) { p[i] = Math.Max(rr[i, i].Real, 1e-12); sum += p[i]; }
        for (int i = 0; i < d; i++) p[i] /= sum;
        return p;
    }

    /// <summary>Per-outcome probabilities of measuring every site in the given Pauli
    /// basis: p = diag(U ρ U†), U the per-site rotation mapping the basis to Z.</summary>
    public static double[] Probs(Matrix<Complex> rho, ReadoutBasis basis)
    {
        int d = rho.RowCount; int n = (int)Math.Round(Math.Log2(d));
        return ProbsRotated(rho, BasisRotation(n, basis));
    }

    /// <summary>FI from precomputed trajectories (share the clean run across bases/bonds).</summary>
    public static double FiMax(IReadOnlyList<Matrix<Complex>> clean,
        IReadOnlyList<Matrix<Complex>> perturbed, double deltaJ, ReadoutBasis basis)
    {
        int n = (int)Math.Round(Math.Log2(clean[0].RowCount));
        var u = BasisRotation(n, basis);
        double best = 0;
        for (int k = 0; k < clean.Count; k++)
        {
            var p0 = ProbsRotated(clean[k], u); var p1 = ProbsRotated(perturbed[k], u);
            double fi = 0;
            for (int i = 0; i < p0.Length; i++)
            { double dp = (p1[i] - p0[i]) / deltaJ; fi += dp * dp / p0[i]; }
            best = Math.Max(best, fi);
        }
        return best;
    }

    /// <summary>Location discrimination from precomputed trajectories.</summary>
    public static double DiscriminationMax(IReadOnlyList<Matrix<Complex>> clean,
        IReadOnlyList<Matrix<Complex>> a, IReadOnlyList<Matrix<Complex>> b, ReadoutBasis basis)
    {
        int n = (int)Math.Round(Math.Log2(clean[0].RowCount));
        var u = BasisRotation(n, basis);
        double best = 0;
        for (int k = 0; k < clean.Count; k++)
        {
            var p0 = ProbsRotated(clean[k], u);
            var pa = ProbsRotated(a[k], u); var pb = ProbsRotated(b[k], u);
            double dsc = 0;
            for (int i = 0; i < p0.Length; i++)
            { double dd = pa[i] - pb[i]; dsc += dd * dd / p0[i]; }
            best = Math.Max(best, dsc);
        }
        return best;
    }

    /// <summary>max over the K-grid of the classical FI of the readout w.r.t. δJ at the
    /// given bond (forward difference: clean vs defected trajectory). Forward-difference
    /// surrogate at the given δJ (O(δJ) bias, ~1% at δJ=0.02); reference values in the
    /// tests carry that convention.</summary>
    public static double StrengthFiMax(int n, double j, double gamma, int defectBond,
        double deltaJ, ReadoutBasis basis, double kMax, int points)
    {
        var ts = KGrid(gamma, kMax, points);
        var clean = Trajectory(n, j, gamma, null, 0.0, ts);
        var pert = Trajectory(n, j, gamma, defectBond, deltaJ, ts);
        return FiMax(clean, pert, deltaJ, basis);
    }

    /// <summary>max over the K-grid of the location-discrimination surrogate between a
    /// defect at bondA vs bondB (same δJ): Σ_i (p_A−p_B)²/p_clean.</summary>
    public static double LocationDMax(int n, double j, double gamma, int bondA, int bondB,
        double deltaJ, ReadoutBasis basis, double kMax, int points)
    {
        var ts = KGrid(gamma, kMax, points);
        var clean = Trajectory(n, j, gamma, null, 0.0, ts);
        var a = Trajectory(n, j, gamma, bondA, deltaJ, ts);
        var b = Trajectory(n, j, gamma, bondB, deltaJ, ts);
        return DiscriminationMax(clean, a, b, basis);
    }
}
