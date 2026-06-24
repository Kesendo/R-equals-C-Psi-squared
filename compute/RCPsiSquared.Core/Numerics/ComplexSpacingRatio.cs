using System;
using System.Collections.Generic;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace RCPsiSquared.Core.Numerics;

/// <summary>The complex spacing ratio (Sá, Ribeiro, Prosen, PRX 2020), the non-Hermitian RMT
/// diagnostic for a spectrum in the complex plane: for each eigenvalue λ_k, z_k = (NN−λ_k)/(NNN−λ_k)
/// where NN/NNN are its nearest and next-nearest neighbours. ⟨|z|⟩ and ⟨cos arg z⟩ classify the
/// spectrum: a 2D-Poisson cloud (integrable / symmetry-fragmented) reads ⟨|z|⟩≈0.658, ⟨cos⟩≈0; a
/// GinUE spectrum (dissipative quantum chaos) reads ⟨|z|⟩≈0.738, ⟨cos⟩≈−0.241 (level repulsion).
/// The reference classes are produced by the same diagnostic (no hardcoded comparison values).</summary>
public static class ComplexSpacingRatio
{
    /// <summary>⟨|z|⟩, ⟨cos arg z⟩, and the distinct-point count over a complex spectrum. Exact
    /// degeneracies are collapsed first (round to 1e-9): the dephased Liouvillian is massively
    /// degenerate, and a coincident NN would give a spurious z=0, so the CSR is meaningful only on
    /// the distinct spectrum. Returns NaN if fewer than 10 distinct points remain.</summary>
    public static (double meanAbs, double meanCos, int count) Of(IReadOnlyList<Complex> points)
    {
        var seen = new HashSet<(long, long)>();
        var pts = new List<Complex>(points.Count);
        foreach (var p in points)
            if (seen.Add(((long)Math.Round(p.Real * 1e9), (long)Math.Round(p.Imaginary * 1e9))))
                pts.Add(p);

        int n = pts.Count;
        if (n < 10) return (double.NaN, double.NaN, n);

        double sumAbs = 0, sumCos = 0;
        int used = 0;
        for (int k = 0; k < n; k++)
        {
            double d1 = double.PositiveInfinity, d2 = double.PositiveInfinity;
            int i1 = -1, i2 = -1;
            for (int j = 0; j < n; j++)
            {
                if (j == k) continue;
                double d = (pts[j] - pts[k]).Magnitude;
                if (d < d1) { d2 = d1; i2 = i1; d1 = d; i1 = j; }
                else if (d < d2) { d2 = d; i2 = j; }
            }
            if (i2 < 0) continue;
            var denom = pts[i2] - pts[k];
            if (denom.Magnitude > 1e-12)
            {
                var z = (pts[i1] - pts[k]) / denom;
                sumAbs += z.Magnitude;
                sumCos += Math.Cos(z.Phase);
                used++;
            }
        }
        return used == 0 ? (double.NaN, double.NaN, n) : (sumAbs / used, sumCos / used, n);
    }

    /// <summary>The integrable reference: ⟨|z|⟩, ⟨cos⟩ of a 2D-Poisson cloud (uniform points in the
    /// unit disk), computed live with a seeded RNG so it is reproducible and never hardcoded.</summary>
    public static (double meanAbs, double meanCos) PoissonDiskReference(int count, int seed)
    {
        var r = new Random(seed);
        var pts = new List<Complex>(count);
        while (pts.Count < count)
        {
            double x = 2 * r.NextDouble() - 1, y = 2 * r.NextDouble() - 1;
            if (x * x + y * y <= 1.0) pts.Add(new Complex(x, y));
        }
        var (a, c, _) = Of(pts);
        return (a, c);
    }

    /// <summary>The dissipative-chaos reference: ⟨|z|⟩, ⟨cos⟩ of a GinUE spectrum (eigenvalues of an
    /// n×n complex Gaussian matrix), via the managed MathNet EVD (no native dependency, as the
    /// live witnesses use). CSR is scale-invariant, so the entries need no 1/√n normalisation.</summary>
    public static (double meanAbs, double meanCos) GinueReference(int n, int seed)
    {
        var r = new Random(seed);
        var m = Matrix<Complex>.Build.Dense(n, n, (_, _) => new Complex(Gauss(r), Gauss(r)));
        var vals = m.Evd().EigenValues;
        var pts = new List<Complex>(n);
        for (int i = 0; i < n; i++) pts.Add(vals[i]);
        var (a, c, _) = Of(pts);
        return (a, c);
    }

    private static double Gauss(Random r)
    {
        double u1 = 1.0 - r.NextDouble(), u2 = 1.0 - r.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }
}
