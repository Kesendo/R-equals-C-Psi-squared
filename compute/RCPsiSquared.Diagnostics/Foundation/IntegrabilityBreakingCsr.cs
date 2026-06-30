using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The F89 Door-C CSR sweep: does breaking the (SE,DE) Liouvillian block's free-fermion
/// additivity (XXZ anisotropy Δ) drive the fixed-q complex spacing ratio (CSR) from Poisson toward
/// Ginibre? The galoischaos witness is the Δ=0 baseline (a clean null: the Galois-S_n H_B-mixed half
/// reads Poisson, NOT GinUE); this sweep turns Δ on and re-reads the CSR with the methodologically
/// correct recipe (review round 2):
/// <list type="bullet">
///   <item>POOL the per-spectrum z-values over the q-grid — never concatenate raw eigenvalues across q
///         (that superimposes independent point processes and fakes Poisson, erasing a Ginibre signal);</item>
///   <item>BOOTSTRAP a 95% CI on ⟨|z|⟩, so a 0.08 Poisson-vs-Ginibre gap is callable;</item>
///   <item>compare against FINITE-SIZE-MATCHED Poisson/GinUE references (pooled over many draws at the
///         measurement's ~50-point per-spectrum size), not the asymptotic 0.658/0.738 that carry the
///         wrong edge bias.</item>
/// </list>
/// The (SE,DE) block is built by <see cref="XxzCoherenceBlock.BuildFull"/> (U(1)-closed under Δ; the
/// Δ·ZZ term is a diagonal IMAGINARY frequency shift, so the AT real-part split is unchanged under Δ).
/// At Δ=0 this reproduces the galoischaos witness baseline (the regression anchor).</summary>
public static class IntegrabilityBreakingCsr
{
    public enum Half { HbMixed, AtLocked, Full }

    private const double ImTol = 1e-6;
    private const double RateTol = 1e-6;

    /// <summary>⟨|z|⟩, ⟨cos arg z⟩, the pooled z count, and a 95% bootstrap CI on ⟨|z|⟩.</summary>
    public readonly record struct CsrReading(int ZCount, double MeanAbs, double MeanCos, double CiLo, double CiHi);

    /// <summary>Upper-half-plane eigenvalues of the (SE,DE) block at (q, Δ), filtered to the chosen half.
    /// AT-locked = Re ∈ {−2, −6} (the absorption-theorem rungs, free-fermion Bloch frequencies);
    /// H_B-mixed = the spread residue (chain Galois S_n). The split is by real part, which Δ leaves
    /// untouched (Δ·ZZ shifts only the imaginary frequency).</summary>
    private static List<Complex> HalfEigs(int n, double q, double delta, Half half, Domain domain)
        => Filter(Matrix<Complex>.Build.DenseOfArray(XxzCoherenceBlock.BuildFull(n, new Complex(q, 0), delta))
            .Evd().EigenValues, half, domain);

    /// <summary>Filter a spectrum to the chosen half (AT-locked Re∈{−2,−6} / H_B-mixed residue / full)
    /// inside the CSR domain (UpperHalf or OffReal).</summary>
    private static List<Complex> Filter(IEnumerable<Complex> vals, Half half, Domain domain)
    {
        var res = new List<Complex>();
        foreach (var lam in vals)
        {
            bool inDomain = domain == Domain.UpperHalf ? lam.Imaginary > ImTol : Math.Abs(lam.Imaginary) > ImTol;
            if (!inDomain) continue;
            bool locked = Math.Abs(lam.Real + 2) < RateTol || Math.Abs(lam.Real + 6) < RateTol;
            bool take = half == Half.Full || (half == Half.HbMixed ? !locked : locked);
            if (take) res.Add(lam);
        }
        return res;
    }

    /// <summary>The CSR fundamental domain. <see cref="UpperHalf"/> (Im &gt; tol) is correct ONLY when the
    /// spectrum is conjugation-symmetric (λ → λ*) — true at Δ=0 (free fermion, 100% conjugate matches) but
    /// NOT at Δ≠0 (the diagonal Δ·ZZ imaginary shift breaks it to 0%). <see cref="OffReal"/> (|Im| &gt;
    /// tol) is the correct domain with no conjugation symmetry: at Δ≠0 there are no real eigenvalues, so
    /// off-real is the full bulk. Using UpperHalf at Δ≠0 keeps an arbitrary non-fundamental ≈60-of-147
    /// subset and biases the CSR.</summary>
    public enum Domain { UpperHalf, OffReal }

    /// <summary>The full (SE,DE) block spectrum at (q, Δ) — all N·C(N,2) eigenvalues, no half-plane
    /// filter. For the symmetry-class diagnostic (is the spectrum conjugation-symmetric, validating the
    /// upper-half-plane CSR restriction? does it carry a reflection about the AT midpoint?).</summary>
    public static Complex[] FullSpectrum(int n, double q, double delta)
        => Matrix<Complex>.Build.DenseOfArray(XxzCoherenceBlock.BuildFull(n, new Complex(q, 0), delta))
            .Evd().EigenValues.ToArray();

    /// <summary>The pooled per-spectrum z-values of the chosen half over the q-grid.</summary>
    private static List<Complex> PooledZ(int n, double delta, double[] qs, Half half, Domain domain)
    {
        var pool = new List<Complex>();
        foreach (var q in qs) pool.AddRange(ComplexSpacingRatio.ZValues(HalfEigs(n, q, delta, half, domain)));
        return pool;
    }

    /// <summary>The pooled-z CSR of the chosen half at anisotropy Δ over the q-grid. Pass the CSR domain
    /// valid for this Δ: UpperHalf at Δ=0 (conjugation-symmetric), OffReal at Δ≠0 (no symmetry).</summary>
    public static CsrReading Sweep(int n, double delta, double[] qs, Half half,
        Domain domain = Domain.UpperHalf, int bootSeed = 1234)
        => Reduce(PooledZ(n, delta, qs, half, domain), bootSeed);

    /// <summary>Stage 2: the random-field disorder-ensemble pooled CSR. For each of <paramref name="realizations"/>
    /// realizations draw a per-site field w_k ~ U[−w, w], build the (SE,DE) block at (q, Δ) + field, and pool
    /// the chosen-half OffReal z-values across realizations. The random field breaks conjugation symmetry, so
    /// OffReal is the valid domain; pooling z's across realizations is both the correct ensemble and the
    /// large-sample source. Δ=0 is free fermion + disorder (1D Anderson, expected Poisson); Δ≠0 is interacting
    /// + disorder (the genuine non-integrability / MBL-ergodic test).</summary>
    public static CsrReading DisorderSweep(int n, double q, double delta, double w, int realizations, Half half, int seed)
    {
        var rng = new Random(seed);
        var pool = new List<Complex>();
        for (int r = 0; r < realizations; r++)
        {
            var field = new double[n];
            for (int k = 0; k < n; k++) field[k] = (2 * rng.NextDouble() - 1) * w;        // U[−w, w]
            var block = XxzCoherenceBlock.BuildFullWithField(n, new Complex(q, 0), delta, field);
            var vals = Matrix<Complex>.Build.DenseOfArray(block).Evd().EigenValues;
            pool.AddRange(ComplexSpacingRatio.ZValues(Filter(vals, half, Domain.OffReal)));
        }
        return Reduce(pool, seed + 7919);
    }

    /// <summary>Per-q ⟨|z|⟩ of the chosen half (the stationarity check: confirm it is flat across q
    /// before trusting the pool; near an EP/discriminant locus the local statistics shift).</summary>
    public static double[] PerQMeanAbs(int n, double delta, double[] qs, Half half,
        Domain domain = Domain.UpperHalf)
        => qs.Select(q =>
        {
            var z = ComplexSpacingRatio.ZValues(HalfEigs(n, q, delta, half, domain));
            return z.Count == 0 ? double.NaN : z.Average(c => c.Magnitude);
        }).ToArray();

    /// <summary>Finite-size-matched 2D-Poisson reference: pool the z's of <paramref name="draws"/> Poisson
    /// clouds, each of <paramref name="size"/> points (≈ the measurement's per-spectrum count), so the
    /// reference carries the SAME finite-size edge bias as the measurement.</summary>
    public static CsrReading PoissonReference(int size, int draws, int seed)
    {
        var pool = new List<Complex>();
        for (int d = 0; d < draws; d++) pool.AddRange(ComplexSpacingRatio.PoissonDiskZValues(size, seed + d));
        return Reduce(pool, seed + 9973);
    }

    /// <summary>Finite-size-matched GinUE reference: pool the z's of <paramref name="draws"/> GinUE spectra
    /// of <paramref name="size"/> eigenvalues each.</summary>
    public static CsrReading GinueReference(int size, int draws, int seed)
    {
        var pool = new List<Complex>();
        for (int d = 0; d < draws; d++) pool.AddRange(ComplexSpacingRatio.GinueZValues(size, seed + d));
        return Reduce(pool, seed + 9973);
    }

    /// <summary>⟨|z|⟩, ⟨cos θ⟩ + a 95% bootstrap CI on ⟨|z|⟩ from a pooled z-list. Shared with
    /// <see cref="FillingThresholdCsr"/> (the Door-C filling-threshold follow-up), the sibling CSR harness.</summary>
    internal static CsrReading Reduce(IReadOnlyList<Complex> zs, int bootSeed, int bootstraps = 400)
    {
        int nz = zs.Count;
        if (nz == 0) return new CsrReading(0, double.NaN, double.NaN, double.NaN, double.NaN);

        var abs = new double[nz];
        double sumCos = 0;
        for (int i = 0; i < nz; i++) { abs[i] = zs[i].Magnitude; sumCos += Math.Cos(zs[i].Phase); }
        double meanAbs = abs.Average();
        double meanCos = sumCos / nz;

        var r = new Random(bootSeed);
        var boot = new double[bootstraps];
        for (int b = 0; b < bootstraps; b++)
        {
            double s = 0;
            for (int i = 0; i < nz; i++) s += abs[r.Next(nz)];
            boot[b] = s / nz;
        }
        Array.Sort(boot);
        return new CsrReading(nz, meanAbs, meanCos, boot[(int)(0.025 * bootstraps)], boot[(int)(0.975 * bootstraps)]);
    }
}
