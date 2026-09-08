using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The F89 Door-C CSR sweep: does breaking the underlying XY Hamiltonian's free-fermion
/// additivity with XXZ anisotropy Δ drive the fixed-q Liouvillian complex spacing ratio (CSR) from Poisson toward
/// Ginibre? The galoischaos witness is the Δ=0 baseline (a clean null: the Galois-S_n H_B-mixed half
/// reads Poisson, NOT GinUE); this sweep turns Δ on and re-reads the CSR with the methodologically
/// correct recipe (review round 2):
/// <list type="bullet">
///   <item>POOL the per-spectrum z-values over the q-grid — never concatenate raw eigenvalues across q
///         (that superimposes independent point processes and fakes Poisson, erasing a Ginibre signal);</item>
///   <item>treat a fixed q-grid as a deterministic descriptive sweep with no sampling CI; for random
///         ensembles, bootstrap whole independently drawn spectra rather than correlated z-values;</item>
///   <item>compare against FINITE-SIZE-MATCHED Poisson/GinUE references (pooled over many draws at the
///         measurement's ~50-point per-spectrum size), not the asymptotic 0.658/0.738 that carry the
///         wrong edge bias.</item>
/// </list>
/// The (SE,DE) block is built by <see cref="XxzCoherenceBlock.BuildFull"/> (U(1)-closed under Δ; the
/// Δ·ZZ term adds -i*q*Delta*(zz(ket)-zz(bra)) to the matrix diagonal).
/// At real q this addition is imaginary; the unchanged dissipator does not fix eigenmode real parts.
/// At Δ=0 this reproduces the galoischaos witness baseline (the regression anchor).
/// Every CSR call retains finite-precision cluster representatives, one per 1e-9 rounded cluster. The resulting CSR and
/// representative count are tolerance-dependent and not an exact degeneracy census; sufficiently close
/// nondegenerate levels can be merged at this declared numerical resolution.</summary>
public static class IntegrabilityBreakingCsr
{
    public enum Half { HbMixed, AtLocked, Full }

    private const double ImTol = 1e-6;
    private const double RateTol = 1e-6;

    /// <summary>The numerical clustering boundary shared by every CSR consumer in this harness.</summary>
    public static string ClusteringScope => ComplexSpacingRatio.ClusteringScope;

    public enum UncertaintySemantics
    {
        NoneDeterministicGrid,
        InsufficientIndependentSpectra,
        SpectrumClusterBootstrap95
    }

    /// <summary>⟨|z|⟩ and ⟨cos arg z⟩ over the pooled z-values. A finite CI is present only for a
    /// whole-spectrum cluster bootstrap with at least two nonempty independently generated spectra.
    /// A deterministic q-grid has no sampling CI and reports zero independent spectra.</summary>
    public readonly record struct CsrReading(
        int ZCount,
        double MeanAbs,
        double MeanCos,
        double CiLo,
        double CiHi,
        int IndependentSpectrumCount,
        UncertaintySemantics Uncertainty);

    /// <summary>Upper-half-plane eigenvalues of the (SE,DE) block at (q, Δ), filtered to the chosen half.
    /// AT-locked = Re ∈ {−2, −6} (the absorption-theorem rungs, with frequencies inherited from the free-fermion XY Hamiltonian);
    /// H_B-mixed labels the remaining values. The filter selects recomputed eigenvalues by their current real parts;
    /// this is not an invariant AT/residual decomposition as Delta changes.</summary>
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
    /// spectrum is conjugation-symmetric (λ → λ*) — true at Δ=0 (100% conjugate matches) but
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

    /// <summary>The pooled per-spectrum z-values of the chosen half over the q-grid. Each spectrum first
    /// contributes one representative per finite-precision cluster under <see cref="ClusteringScope"/>.</summary>
    private static IReadOnlyList<Complex>[] ZGrid(int n, double delta, double[] qs, Half half, Domain domain)
        => qs.Select(q => (IReadOnlyList<Complex>)ComplexSpacingRatio.ZValues(
            HalfEigs(n, q, delta, half, domain))).ToArray();

    /// <summary>The pooled-z CSR of the chosen half at anisotropy Δ over a fixed q-grid. This is a
    /// deterministic descriptive sweep, not a random sample, so CiLo/CiHi are NaN. Pass the CSR domain
    /// valid for this Δ: UpperHalf at Δ=0 (conjugation-symmetric), OffReal at Δ≠0 (no symmetry).</summary>
    public static CsrReading Sweep(int n, double delta, double[] qs, Half half,
        Domain domain = Domain.UpperHalf)
        => ReduceDeterministicGrid(ZGrid(n, delta, qs, half, domain));

    /// <summary>Stage 2: the random-field disorder-ensemble pooled CSR. For each of <paramref name="realizations"/>
    /// realizations draw a per-site field w_k ~ U[−w, w], build the (SE,DE) block at (q, Δ) + field, and pool
    /// the chosen-half OffReal z-values across realizations. For w&gt;0 the 95% interval resamples these whole independent
    /// spectra; z-values within one spectrum are not treated as independent. At w=0 repeated copies are one
    /// deterministic spectrum and no sampling CI is reported. The random field breaks conjugation
    /// symmetry, so OffReal is the valid domain. At Δ=0 the random-field XY Hamiltonian is quadratic (1D Anderson, expected Poisson);
    /// this does not make the dephasing Liouvillian a quadratic free-fermion generator. At Δ≠0 the Hamiltonian is
    /// interacting and disordered (the genuine non-integrability / MBL-ergodic test).</summary>
    public static CsrReading DisorderSweep(int n, double q, double delta, double w, int realizations, Half half, int seed)
    {
        var rng = new Random(seed);
        var spectra = new List<IReadOnlyList<Complex>>();
        for (int r = 0; r < realizations; r++)
        {
            var field = new double[n];
            for (int k = 0; k < n; k++) field[k] = (2 * rng.NextDouble() - 1) * w;        // U[−w, w]
            var block = XxzCoherenceBlock.BuildFullWithField(n, new Complex(q, 0), delta, field);
            var vals = Matrix<Complex>.Build.DenseOfArray(block).Evd().EigenValues;
            spectra.Add(ComplexSpacingRatio.ZValues(Filter(vals, half, Domain.OffReal)));
        }
        return w == 0.0
            ? ReduceDeterministicGrid(spectra)
            : ReduceIndependentSpectra(spectra, seed + 7919);
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
        var spectra = new List<IReadOnlyList<Complex>>();
        for (int d = 0; d < draws; d++) spectra.Add(ComplexSpacingRatio.PoissonDiskZValues(size, seed + d));
        return ReduceIndependentSpectra(spectra, seed + 9973);
    }

    /// <summary>Finite-size-matched GinUE reference: pool the z's of <paramref name="draws"/> GinUE spectra
    /// of <paramref name="size"/> eigenvalues each.</summary>
    public static CsrReading GinueReference(int size, int draws, int seed)
    {
        var spectra = new List<IReadOnlyList<Complex>>();
        for (int d = 0; d < draws; d++) spectra.Add(ComplexSpacingRatio.GinueZValues(size, seed + d));
        return ReduceIndependentSpectra(spectra, seed + 9973);
    }

    internal static CsrReading ReduceDeterministicGrid(IReadOnlyList<IReadOnlyList<Complex>> spectra)
        => ReducePooled(spectra, independentSpectrumCount: 0,
            UncertaintySemantics.NoneDeterministicGrid, double.NaN, double.NaN);

    /// <summary>Pool z-values for the point estimate, but bootstrap the independent spectrum clusters.
    /// Empty filtered spectra contribute no statistic and are not counted as bootstrap units.</summary>
    internal static CsrReading ReduceIndependentSpectra(
        IReadOnlyList<IReadOnlyList<Complex>> spectra, int bootSeed, int bootstraps = 400)
    {
        var nonempty = spectra.Where(s => s.Count > 0).ToArray();
        if (nonempty.Length < 2)
            return ReducePooled(nonempty, nonempty.Length,
                UncertaintySemantics.InsufficientIndependentSpectra, double.NaN, double.NaN);

        var clusterSums = nonempty.Select(s => s.Sum(z => z.Magnitude)).ToArray();
        var clusterCounts = nonempty.Select(s => s.Count).ToArray();
        var rng = new Random(bootSeed);
        var boot = new double[bootstraps];
        for (int b = 0; b < bootstraps; b++)
        {
            double sum = 0.0;
            int count = 0;
            for (int i = 0; i < nonempty.Length; i++)
            {
                int pick = rng.Next(nonempty.Length);
                sum += clusterSums[pick];
                count += clusterCounts[pick];
            }
            boot[b] = sum / count;
        }
        Array.Sort(boot);
        return ReducePooled(nonempty, nonempty.Length,
            UncertaintySemantics.SpectrumClusterBootstrap95,
            boot[(int)(0.025 * bootstraps)], boot[(int)(0.975 * bootstraps)]);
    }

    private static CsrReading ReducePooled(
        IReadOnlyList<IReadOnlyList<Complex>> spectra,
        int independentSpectrumCount,
        UncertaintySemantics uncertainty,
        double ciLo,
        double ciHi)
    {
        int nz = spectra.Sum(s => s.Count);
        if (nz == 0)
            return new CsrReading(0, double.NaN, double.NaN, double.NaN, double.NaN,
                independentSpectrumCount, uncertainty);

        double sumAbs = 0.0;
        double sumCos = 0.0;
        foreach (var spectrum in spectra)
        {
            foreach (var z in spectrum)
            {
                sumAbs += z.Magnitude;
                sumCos += Math.Cos(z.Phase);
            }
        }
        return new CsrReading(nz, sumAbs / nz, sumCos / nz, ciLo, ciHi,
            independentSpectrumCount, uncertainty);
    }
}
