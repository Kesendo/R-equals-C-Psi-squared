using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The (SE,DE) = (weight-1 ket, weight-2 bra) coherence block of the XXZ-chain Liouvillian,
/// parameterized by q = J/γ AND the ZZ-anisotropy Δ, for the F89 path-4 (N=5) integrability test
/// (docs/superpowers/plans/2026-06-27-f89-path4-delta-test.md). Convention (bare Pauli, matching
/// F89Path3OcticBlock / f89_zz_break_gate.py): H(Δ) = J·Σ_b(X_bX_{b+1}+Y_bY_{b+1}) + J·Δ·Σ_b Z_bZ_{b+1},
/// γ=1, q=J. L = −i[H,·] + Z-dephasing:
/// <list type="bullet">
/// <item>dephasing diagonal −2·n_diff (Absorption Theorem 2γ·HammingDistance(ket,bra), ∈{−2,−6} here);</item>
/// <item>ket excitation hops −2qi, bra excitation hops +2qi (the −i(H⊗I − I⊗Hᵀ) split), NN, exclusion;</item>
/// <item>the Δ·ZZ term is DIAGONAL in the computational basis and Hermitian, so the AT rate is UNTOUCHED
/// (Re λ = −2γ⟨n_XY⟩ preserved); it adds only the frequency −i·qΔ·(zz(ket) − zz(bra)), with
/// zz(c) = Σ_bonds(+1 if the two sites are equal, −1 if they differ) = Σ_b ⟨c|Z_bZ_{b+1}|c⟩.</item>
/// </list>
/// The S₂ site-reflection (s → N−1−s) R=+1 symmetric sector carries the diabolics (same sector as
/// F89Path3OcticBlock.BuildSeDeSymBlock at N=4 / the path-k residual at N≥5). At Δ=0 this reproduces the
/// XY (SE,DE) block exactly (the trusted anchor). The Δ-test object: does a path-4 diabolic flip defective
/// or lift as Δ turns on, generalizing DIABOLIC_BY_INTEGRABILITY's N=4 gate off N=4.</summary>
public static class XxzCoherenceBlock
{
    // all bitmasks on n sites with exactly w set bits (excitations), ascending.
    private static List<int> Weight(int n, int w)
    {
        var res = new List<int>();
        for (int m = 0; m < (1 << n); m++)
            if (System.Numerics.BitOperations.PopCount((uint)m) == w) res.Add(m);
        return res;
    }

    // zz(c) = Σ_{bond (b,b+1)} ⟨c|Z_bZ_{b+1}|c⟩ = Σ_b (+1 if bits b,b+1 equal, −1 if differ).
    private static int Zz(int n, int c)
    {
        int s = 0;
        for (int b = 0; b < n - 1; b++)
            s += (((c >> b) & 1) == ((c >> (b + 1)) & 1)) ? 1 : -1;
        return s;
    }

    // reverse the n-bit string (the site reflection s → n−1−s).
    private static int Reflect(int n, int m)
    {
        int r = 0;
        for (int s = 0; s < n; s++)
            if ((m & (1 << s)) != 0) r |= 1 << (n - 1 - s);
        return r;
    }

    private static (List<(int ket, int bra)> basis, Dictionary<(int, int), int> index) Basis(int n)
    {
        var kets = Weight(n, 1);
        var bras = Weight(n, 2);
        var basis = new List<(int, int)>();
        var index = new Dictionary<(int, int), int>();
        foreach (var k in kets)
            foreach (var b in bras) { index[(k, b)] = basis.Count; basis.Add((k, b)); }
        return (basis, index);
    }

    /// <summary>The full (SE,DE) coherence block at complex q and real Δ.</summary>
    public static Complex[,] BuildFull(int n, Complex q, double delta)
    {
        var (basis, index) = Basis(n);
        int d = basis.Count;
        var l = new Complex[d, d];
        for (int col = 0; col < d; col++)
        {
            var (kc, bc) = basis[col];
            int nDiff = System.Numerics.BitOperations.PopCount((uint)(kc ^ bc));
            // AT dephasing rate (−2·n_diff) + the Δ·ZZ frequency (−i·qΔ·(zz_ket − zz_bra)), q complex.
            l[col, col] += new Complex(-2.0 * nDiff, 0)
                         + (-Complex.ImaginaryOne) * q * (delta * (Zz(n, kc) - Zz(n, bc)));
            for (int s = 0; s < n; s++)                                  // ket excitation hops −2qi
                if ((kc & (1 << s)) != 0)
                    foreach (int s2 in new[] { s - 1, s + 1 })
                        if (s2 >= 0 && s2 < n && (kc & (1 << s2)) == 0)
                            l[index[((kc & ~(1 << s)) | (1 << s2), bc)], col] += new Complex(0, -2) * q;
            for (int s = 0; s < n; s++)                                  // bra excitation hops +2qi
                if ((bc & (1 << s)) != 0)
                    foreach (int s2 in new[] { s - 1, s + 1 })
                        if (s2 >= 0 && s2 < n && (bc & (1 << s2)) == 0)
                            l[index[(kc, (bc & ~(1 << s)) | (1 << s2))], col] += new Complex(0, 2) * q;
        }
        return l;
    }

    /// <summary>The R=+1 (site reflection s→N−1−s) symmetric sector of the full block. Same construction as
    /// F89Path3OcticBlock.BuildSeDeSymBlock, generalized to any N (the reflection commutes with XX+YY, ZZ and
    /// uniform dephasing, so the sector is invariant including the Δ·ZZ term).</summary>
    public static Matrix<Complex> BuildSym(int n, Complex q, double delta)
    {
        var (basis, index) = Basis(n);
        int d = basis.Count;
        var full = Matrix<Complex>.Build.DenseOfArray(BuildFull(n, q, delta));

        var cols = new List<Complex[]>();
        var handled = new HashSet<int>();
        for (int col = 0; col < d; col++)
        {
            if (handled.Contains(col)) continue;
            var (kc, bc) = basis[col];
            int mcol = index[(Reflect(n, kc), Reflect(n, bc))];
            var v = new Complex[d];
            if (mcol == col) v[col] = Complex.One;                      // reflection-fixed coherence
            else { double s = 1.0 / Math.Sqrt(2); v[col] = s; v[mcol] = s; handled.Add(mcol); }
            handled.Add(col);
            cols.Add(v);
        }
        var p = Matrix<Complex>.Build.Dense(d, cols.Count, (r, c) => cols[c][r]);
        return p.ConjugateTranspose() * full * p;
    }

    /// <summary>The spectrum of the R=+1 symmetric (SE,DE) sector at (q, Δ). At Δ=0, N=4 this reproduces
    /// F89Path3OcticBlock.BuildSeDeSymBlock(q, 1)'s spectrum (the trusted XY anchor).</summary>
    public static Complex[] SeDeSymSpectrum(int n, Complex q, double delta)
        => BuildSym(n, q, delta).Evd().EigenValues.ToArray();
}
