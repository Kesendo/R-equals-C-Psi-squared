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
/// The reference classes are produced by the same diagnostic (no hardcoded comparison values).
///
/// <para><see cref="ZValues"/> exposes the per-point z's so they can be POOLED across many spectra
/// (the methodologically correct way to build a CSR statistic over a q-sweep or a disorder ensemble:
/// compute z within each fixed-parameter spectrum, then pool the dimensionless z's — never concatenate
/// raw eigenvalues across spectra, which superimposes independent point processes and fakes Poisson).
/// <see cref="PoissonDiskZValues"/> / <see cref="GinueZValues"/> give the per-point z's of one
/// finite-size reference draw, so a pool of small draws is a finite-size-matched reference rather than
/// the asymptotic 0.658 / 0.738.</para></summary>
public static class ComplexSpacingRatio
{
    /// <summary>Dedup (round to 1e-9) then the per-point z's: for each distinct point, z = (NN−λ)/(NNN−λ).
    /// Returns the z list (one per point with a non-degenerate NNN denominator) and the distinct-point
    /// count. The dephased Liouvillian is massively degenerate, and a coincident NN would give a spurious
    /// z=0, so the CSR is meaningful only on the distinct spectrum.</summary>
    private static (List<Complex> zs, int distinct) Compute(IReadOnlyList<Complex> points)
    {
        var seen = new HashSet<(long, long)>();
        var pts = new List<Complex>(points.Count);
        foreach (var p in points)
            if (seen.Add(((long)Math.Round(p.Real * 1e9), (long)Math.Round(p.Imaginary * 1e9))))
                pts.Add(p);

        int n = pts.Count;
        var zs = new List<Complex>(n);
        if (n < 10) return (zs, n);

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
                zs.Add((pts[i1] - pts[k]) / denom);
        }
        return (zs, n);
    }

    /// <summary>⟨|z|⟩, ⟨cos arg z⟩, and the distinct-point count over a complex spectrum. Returns NaN if
    /// fewer than 10 distinct points remain (or none has a valid neighbour ratio).</summary>
    public static (double meanAbs, double meanCos, int count) Of(IReadOnlyList<Complex> points)
    {
        var (zs, distinct) = Compute(points);
        if (distinct < 10 || zs.Count == 0) return (double.NaN, double.NaN, distinct);

        double sumAbs = 0, sumCos = 0;
        foreach (var z in zs) { sumAbs += z.Magnitude; sumCos += Math.Cos(z.Phase); }
        return (sumAbs / zs.Count, sumCos / zs.Count, distinct);
    }

    /// <summary>The per-point complex spacing ratios z_k over one spectrum (dedup at 1e-9). Empty if
    /// fewer than 10 distinct points. Pool these across spectra to build a CSR statistic over a q-sweep
    /// or a disorder ensemble, then bootstrap the pooled list.</summary>
    public static IReadOnlyList<Complex> ZValues(IReadOnlyList<Complex> points) => Compute(points).zs;

    /// <summary>A 2D-Poisson cloud: uniform points in the unit disk, seeded RNG (reproducible).</summary>
    private static List<Complex> PoissonDiskCloud(int count, int seed)
    {
        var r = new Random(seed);
        var pts = new List<Complex>(count);
        while (pts.Count < count)
        {
            double x = 2 * r.NextDouble() - 1, y = 2 * r.NextDouble() - 1;
            if (x * x + y * y <= 1.0) pts.Add(new Complex(x, y));
        }
        return pts;
    }

    /// <summary>A GinUE spectrum: eigenvalues of an n×n complex Gaussian matrix (managed MathNet EVD,
    /// no native dependency). CSR is scale-invariant, so the entries need no 1/√n normalisation.</summary>
    private static List<Complex> GinueCloud(int n, int seed)
    {
        var r = new Random(seed);
        var m = Matrix<Complex>.Build.Dense(n, n, (_, _) => new Complex(Gauss(r), Gauss(r)));
        var vals = m.Evd().EigenValues;
        var pts = new List<Complex>(n);
        for (int i = 0; i < n; i++) pts.Add(vals[i]);
        return pts;
    }

    /// <summary>The integrable reference: ⟨|z|⟩, ⟨cos⟩ of a 2D-Poisson cloud, computed live (never
    /// hardcoded).</summary>
    public static (double meanAbs, double meanCos) PoissonDiskReference(int count, int seed)
    {
        var (a, c, _) = Of(PoissonDiskCloud(count, seed));
        return (a, c);
    }

    /// <summary>The dissipative-chaos reference: ⟨|z|⟩, ⟨cos⟩ of a GinUE spectrum, computed live.</summary>
    public static (double meanAbs, double meanCos) GinueReference(int n, int seed)
    {
        var (a, c, _) = Of(GinueCloud(n, seed));
        return (a, c);
    }

    /// <summary>The per-point z's of one finite-size 2D-Poisson draw, for pooling a finite-size-matched
    /// reference at the measurement's per-spectrum size.</summary>
    public static IReadOnlyList<Complex> PoissonDiskZValues(int count, int seed) => ZValues(PoissonDiskCloud(count, seed));

    /// <summary>The per-point z's of one finite-size GinUE draw, for pooling a finite-size-matched
    /// reference at the measurement's per-spectrum size.</summary>
    public static IReadOnlyList<Complex> GinueZValues(int n, int seed) => ZValues(GinueCloud(n, seed));

    private static double Gauss(Random r)
    {
        double u1 = 1.0 - r.NextDouble(), u2 = 1.0 - r.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }
}
