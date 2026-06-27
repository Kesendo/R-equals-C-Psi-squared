using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace RCPsiSquared.Core.Spectrum;

/// <summary>The free-fermion-additivity probe for the F89 Door-C CSR sequel
/// (docs/superpowers/plans/2026-06-27-f89-door-c-csr-integrability-sweep.md, review round 1):
/// the sample-size-INDEPENDENT mechanism tracker for "does breaking the (SE,DE) Liouvillian block's
/// free-fermion additivity drive the fixed-q CSR toward Ginibre?".
///
/// <para>The (SE,DE) coherence bra is a two-excitation (two-magnon) state of the open XXZ chain
/// H(Δ) = J·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}) + J·Δ·Σ_b Z_b Z_{b+1}, reference = all-down |↓…↓⟩,
/// a magnon = one up-spin. At Δ=0 the chain is free fermions (Jordan-Wigner is boundary-term-free on an
/// OPEN chain), so the two-magnon excitation energies are EXACT sums of single-magnon excitation
/// energies — the additivity residual is bit-exact zero. The ZZ term is a contact interaction between
/// adjacent magnons; as Δ turns on it breaks the additivity EVEN THOUGH the XXZ Hamiltonian stays
/// Bethe-integrable. That distinction (free-fermion additivity ≠ Bethe integrability) is exactly what
/// the Door-C sequel must separate, and this residual is its honest, sample-size-independent guard.</para>
///
/// <para>Energies are measured relative to the vacuum (all-down) energy E_vac = J·Δ·(N−1), so the
/// residual tracks the genuine interaction, not the trivial doubly-counted ZZ reference constant.</para></summary>
public static class TwoMagnonAdditivity
{
    /// <summary>All w-subsets of {0,…,n−1}, each sorted ascending, in lexicographic order.</summary>
    private static List<int[]> Configs(int n, int w)
    {
        var res = new List<int[]>();
        void Rec(int start, List<int> cur)
        {
            if (cur.Count == w) { res.Add(cur.ToArray()); return; }
            for (int s = start; s < n; s++) { cur.Add(s); Rec(s + 1, cur); cur.RemoveAt(cur.Count - 1); }
        }
        Rec(0, new List<int>());
        return res;
    }

    /// <summary>ZZ diagonal energy J·Δ·Σ_b z_b z_{b+1} of a configuration (z = −1 on up-sites, +1 on
    /// down-sites), open chain bonds b = 0..n−2.</summary>
    private static double ZzEnergy(int n, double j, double delta, HashSet<int> up)
    {
        int s = 0;
        for (int b = 0; b < n - 1; b++)
        {
            int zb = up.Contains(b) ? -1 : 1;
            int zb1 = up.Contains(b + 1) ? -1 : 1;
            s += zb * zb1;
        }
        return j * delta * s;
    }

    /// <summary>The single-magnon (weight-1) Hamiltonian, N×N over the basis {up-site i : i=0..N−1}.
    /// Off-diagonal hopping 2J between adjacent sites (open chain); diagonal = the ZZ energy.</summary>
    public static double[,] SingleMagnonHamiltonian(int n, double j, double delta)
    {
        var h = new double[n, n];
        for (int i = 0; i < n; i++)
        {
            h[i, i] = ZzEnergy(n, j, delta, new HashSet<int> { i });
            foreach (int s2 in new[] { i - 1, i + 1 })
                if (s2 >= 0 && s2 < n) h[s2, i] += 2.0 * j;
        }
        return h;
    }

    /// <summary>The two-magnon (weight-2) Hamiltonian, C(N,2)×C(N,2) over the 2-subset basis. Each
    /// magnon hops 2J to an adjacent EMPTY site (hard core); diagonal = the ZZ energy. On an open chain
    /// this is exactly the free-fermion two-particle problem at Δ=0.</summary>
    public static double[,] TwoMagnonHamiltonian(int n, double j, double delta)
    {
        var configs = Configs(n, 2);
        int d = configs.Count;
        var index = new Dictionary<(int, int), int>();
        for (int t = 0; t < d; t++) index[(configs[t][0], configs[t][1])] = t;

        var h = new double[d, d];
        for (int col = 0; col < d; col++)
        {
            var c = configs[col];
            var up = new HashSet<int> { c[0], c[1] };
            h[col, col] = ZzEnergy(n, j, delta, up);
            foreach (int site in c)
                foreach (int s2 in new[] { site - 1, site + 1 })
                {
                    if (s2 < 0 || s2 >= n || up.Contains(s2)) continue;      // off-chain or hard-core blocked
                    int other = site == c[0] ? c[1] : c[0];
                    int a = Math.Min(s2, other), b = Math.Max(s2, other);
                    h[index[(a, b)], col] += 2.0 * j;
                }
        }
        return h;
    }

    /// <summary>The additivity residual ‖sort(E_DE) − sort({ε_a+ε_b : a&lt;b})‖₂, both spectra measured
    /// as excitation energies above the vacuum E_vac = J·Δ·(N−1). Bit-exact zero at Δ=0; grows as the
    /// ZZ contact interaction breaks free-fermion additivity.</summary>
    public static double Residual(int n, double j, double delta)
    {
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(n), n, "N must be ≥ 2.");
        double eVac = j * delta * (n - 1);

        var single = SymEigs(SingleMagnonHamiltonian(n, j, delta)).Select(x => x - eVac).ToArray();
        var dbl = SymEigs(TwoMagnonHamiltonian(n, j, delta)).Select(x => x - eVac).ToArray();

        var sums = new List<double>();
        for (int a = 0; a < single.Length; a++)
            for (int b = a + 1; b < single.Length; b++)
                sums.Add(single[a] + single[b]);

        var sa = sums.OrderBy(x => x).ToArray();
        var da = dbl.OrderBy(x => x).ToArray();
        double r = 0;
        for (int t = 0; t < sa.Length; t++) { double diff = da[t] - sa[t]; r += diff * diff; }
        return Math.Sqrt(r);
    }

    private static double[] SymEigs(double[,] m)
        => Matrix<double>.Build.DenseOfArray(m).Evd().EigenValues.Select(c => c.Real).ToArray();
}
