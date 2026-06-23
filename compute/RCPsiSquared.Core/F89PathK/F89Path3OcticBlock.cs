using System.Collections.Generic;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;  // disambiguate from System.Numerics.Vector<T> (SIMD): bare Vector<Complex> is CS0104 here (both namespaces imported) — repo convention, see InspectablePayload.cs:3

namespace RCPsiSquared.Core.F89PathK;

/// <summary>The path-3 (SE,DE) coherence sub-block of the N=4 Liouvillian:
/// SE = single-excitation ket (4 sites), DE = double-excitation bra (6 pairs).
/// L = −i·M_SE(ket index) + i·M_DE(bra index) + diag(−2γ overlap / −6γ no-overlap),
/// with M_SE/M_DE the 2J nearest-neighbour hops. The R=+1 (site reflection s→3−s)
/// symmetric sector (12-dim) carries the octic whose (3q⁴+q²−1)² discriminant
/// factor locates the DIABOLIC degeneracy at q_EP. Convention matches
/// <see cref="F89BlockLiouvillian"/> (H = J·(XX+YY), SE hop 2J).</summary>
public static class F89Path3OcticBlock
{
    private static readonly (int A, int B)[] DePairs =
        { (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3) };

    private static int Idx(int seSite, int dePair) => seSite * 6 + dePair;

    private static int PairIndex(int a, int b)
    {
        int lo = a < b ? a : b, hi = a < b ? b : a;
        for (int p = 0; p < 6; p++)
            if (DePairs[p].A == lo && DePairs[p].B == hi) return p;
        return -1;
    }

    /// <summary>The full 24×24 (SE,DE) coherence block at (J, γ).</summary>
    public static ComplexMatrix BuildSeDeBlock(double j, double gamma)
    {
        if (gamma < 0)
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "gamma must be ≥ 0.");

        // M_SE: nearest-neighbour hop on the 4-chain, amplitude 2J
        var mse = new double[4, 4];
        for (int a = 0; a < 4; a++)
            for (int b = 0; b < 4; b++)
                if (a - b == 1 || b - a == 1) mse[a, b] = 2 * j;

        // M_DE: move one excitation of the pair to an adjacent in-range empty site
        var mde = new double[6, 6];
        for (int col = 0; col < 6; col++)
        {
            var (pj, pk) = DePairs[col];
            foreach (int nj in new[] { pj - 1, pj + 1 })
                if (nj >= 0 && nj <= 3 && nj != pk)
                {
                    int row = PairIndex(nj, pk);
                    if (row >= 0) mde[row, col] += 2 * j;
                }
            foreach (int nk in new[] { pk - 1, pk + 1 })
                if (nk >= 0 && nk <= 3 && nk != pj)
                {
                    int row = PairIndex(pj, nk);
                    if (row >= 0) mde[row, col] += 2 * j;
                }
        }

        var l = ComplexMatrix.Build.Dense(24, 24);
        for (int i = 0; i < 4; i++)
            for (int p = 0; p < 6; p++)
            {
                int col = Idx(i, p);
                for (int i2 = 0; i2 < 4; i2++)                       // SE side (ket): −i·M_SE (opposite signs = ket vs bra side of the commutator)
                    if (mse[i2, i] != 0)
                        l[Idx(i2, p), col] += new Complex(0, -mse[i2, i]);
                for (int p2 = 0; p2 < 6; p2++)                       // DE side (bra): +i·M_DE
                    if (mde[p, p2] != 0)
                        l[Idx(i, p2), col] += new Complex(0, mde[p, p2]);
                var (dj, dk) = DePairs[p];                           // dissipative diagonal
                double diag = (i == dj || i == dk) ? -2 * gamma : -6 * gamma;
                l[col, col] += new Complex(diag, 0);
            }
        return l;
    }

    /// <summary>The 12×12 R=+1 (site reflection s→3−s) symmetric sector.</summary>
    public static ComplexMatrix BuildSeDeSymBlock(double j, double gamma)
    {
        var full = BuildSeDeBlock(j, gamma);
        int[] perm = { 3, 2, 1, 0 };
        int ReflPair(int p)
        {
            var (a, b) = DePairs[p];
            return PairIndex(perm[a], perm[b]);
        }

        var cols = new List<ComplexVector>();
        var handled = new HashSet<int>();
        for (int i = 0; i < 4; i++)
            for (int p = 0; p < 6; p++)
            {
                int n = Idx(i, p);
                if (handled.Contains(n)) continue;
                int m = Idx(perm[i], ReflPair(p));
                var v = ComplexVector.Build.Dense(24);
                // perm {3,2,1,0} is fixed-point-free, so n==m never fires here; kept for parity with the general Python reflection_projector.
                if (n == m) { v[n] = Complex.One; }
                else
                {
                    double s = 1.0 / System.Math.Sqrt(2);
                    v[n] = s; v[m] = s; handled.Add(m);
                }
                handled.Add(n);
                cols.Add(v);
            }
        var pPlus = ComplexMatrix.Build.DenseOfColumnVectors(cols);  // 24×12
        return pPlus.ConjugateTranspose() * full * pPlus;            // 12×12
    }

    /// <summary>The Gaussian-integer coefficients of F_8(λ, q=2) — the certifying
    /// specialization for the Galois group. Index = λ-power (0..8); coefficient k is
    /// re[k] + i·im[k]. Monic over Z[i]. This is the q=2 evaluation of the Gate-0-verified
    /// octic literal (build_octic in simulations/f89_path3_octic_galois.py, reproduced
    /// bit-exact from the 12×12 symbolic charpoly). Modulo the split prime 𝔭|5 (F_5, i↦2)
    /// it factors to cycle type (5,2,1), which certifies Gal(F_8) = S_8
    /// (see <see cref="RCPsiSquared.Core.Numerics.OcticGaloisCertificate"/>).</summary>
    public static (long[] Re, long[] Im) OcticCoefficientsAtQ2() => (
        Re: new long[] { 8407296, 7833600, 2976000, 671232, 99808, 10112, 720, 32, 1 },
        Im: new long[] { -5275648, -2760704, -671744, -110592, -10240, -512, 0, 0, 0 });
}
