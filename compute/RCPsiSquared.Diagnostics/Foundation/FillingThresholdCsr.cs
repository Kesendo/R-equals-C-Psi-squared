using System;
using System.Collections.Generic;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The F89 Door-C filling comparison: finite executed CSR evidence compares dilute (SE,DE)=(1,2)
/// with dense blocks under the same interacting disorder at canonical Delta=1.
/// Nonzero Delta breaks the underlying Hamiltonian's free-fermion additivity, but uniform XXZ remains Bethe-integrable;
/// at Delta=0 the random-field XY Hamiltonian remains quadratic (Anderson/free fermions). This Hamiltonian label
/// does not classify the Z-dephasing Liouvillian as a quadratic free-fermion generator;
/// generic random field plus Delta!=0 is the interacting disordered nonintegrable test. This harness builds the GENERAL
/// (wKet,wBra) coherence block (<see cref="WeightCoherenceBlock.Build(int,int,int,Complex,double,double[])"/>) at
/// EXTENSIVE filling (wKet,wBra near N/2) and re-runs the same disordered CSR. If the DENSE block reaches GinUE
/// while the dilute one does not, that supports a finite-size filling dependence, not a universal
/// thermalization cause or a deduction of Hamiltonian integrability from spacings.
///
/// <para>GinUE is retained only as a comparison ensemble. Unequal weight (p,p+1) means the known F1 map leaves
/// the block, and <see cref="ConjugationMatchFraction"/> can show that one bare conjugation relation is absent.
/// Neither fact exhausts the unitary/antiunitary algebra after irreducible strong-symmetry reduction, so the
/// sector's full SRP class remains open. Methodology inherited from the Door-C harness: pool
/// per-spectrum z's (never raw eigenvalues), cluster-bootstrap whole independent disorder realizations,
/// treat a fixed q-grid as descriptive with no sampling CI, and read in the OffReal domain (|Im| &gt; tol) for a
/// like-for-like comparison with the complex GinUE cloud. Every spectrum contributes finite-precision cluster representatives,
/// one per 1e-9 rounded cluster; the resulting CSR and count are tolerance-dependent and not an exact degeneracy census.</para></summary>
public static class FillingThresholdCsr
{
    private const double ImTol = 1e-6;

    /// <summary>The numerical clustering boundary shared by every CSR consumer in this harness.</summary>
    public static string ClusteringScope => ComplexSpacingRatio.ClusteringScope;

    public static string ReferenceScope =>
        "GinUE comparison only; the full sector symmetry algebra and irreducible SRP class remain open. " +
        ClusteringScope;

    /// <summary>Off-real eigenvalues (|Im| &gt; tol) used only for the GinUE comparison.
    /// Selecting this subset does not determine the block's symmetry class.</summary>
    private static List<Complex> OffReal(IEnumerable<Complex> vals)
    {
        var res = new List<Complex>();
        foreach (var lam in vals)
            if (Math.Abs(lam.Imaginary) > ImTol) res.Add(lam);
        return res;
    }

    private static List<Complex> OffRealSpectrum(int n, int wKet, int wBra, double q, double delta, double[]? field)
        => OffReal(Matrix<Complex>.Build
            .DenseOfArray(WeightCoherenceBlock.Build(n, wKet, wBra, new Complex(q, 0), delta, field))
            .Evd().EigenValues);

    /// <summary>The disorder-ensemble pooled CSR of the (wKet,wBra) block at (q, Δ). For each of
    /// <paramref name="realizations"/> realizations draw a per-site field w_k ~ U[−w, w], build the block with that
    /// field, and pool the OffReal per-spectrum z's after finite-precision clustering; the 95% interval
    /// resamples whole spectra, not the correlated ratios inside them. w=0 is deterministic (the clean block, every
    /// realization identical), so repeated copies report no sampling CI. At Δ=0 the disordered XY Hamiltonian is quadratic/Anderson-like; this is not a
    /// quadratic-Liouvillian claim. At Δ≠0 the Hamiltonian is interacting and disordered.</summary>
    public static IntegrabilityBreakingCsr.CsrReading DisorderSweep(
        int n, int wKet, int wBra, double q, double delta, double w, int realizations, int seed)
    {
        var rng = new Random(seed);
        var spectra = new List<IReadOnlyList<Complex>>();
        for (int r = 0; r < realizations; r++)
        {
            double[]? field = null;
            if (w != 0.0)
            {
                field = new double[n];
                for (int k = 0; k < n; k++) field[k] = (2 * rng.NextDouble() - 1) * w;   // U[−w, w]
            }
            spectra.Add(ComplexSpacingRatio.ZValues(OffRealSpectrum(n, wKet, wBra, q, delta, field)));
        }
        return w == 0.0
            ? IntegrabilityBreakingCsr.ReduceDeterministicGrid(spectra)
            : IntegrabilityBreakingCsr.ReduceIndependentSpectra(spectra, seed + 7919);
    }

    /// <summary>The clean (disorder-free) pooled CSR of the (wKet,wBra) block at Δ, pooled over the q-grid. The
    /// integrable/no-disorder control: at Δ=0 the underlying XY Hamiltonian is free-fermion; at Δ≠0 the Hamiltonian
    /// is Bethe-integrable. Neither statement classifies the Z-dephasing Liouvillian as free-fermion.
    /// The q-grid is deterministic, so this method reports no sampling CI. OffReal throughout.</summary>
    public static IntegrabilityBreakingCsr.CsrReading CleanSweep(int n, int wKet, int wBra, double[] qs, double delta)
    {
        var spectra = new List<IReadOnlyList<Complex>>();
        foreach (var q in qs)
            spectra.Add(ComplexSpacingRatio.ZValues(OffRealSpectrum(n, wKet, wBra, q, delta, null)));
        return IntegrabilityBreakingCsr.ReduceDeterministicGrid(spectra);
    }

    /// <summary>The fraction of eigenvalues λ whose conjugate λ* is also in the spectrum (within tol). ≈ 1 ⟹
    /// conjugation-symmetric (Δ=0, no field), ≈ 0 means this particular spectral-conjugation
    /// match is absent. It does not exclude other antiunitary relations after irreducible-sector reduction and
    /// therefore does not assign a symmetry class.
    ///
    /// <para>The involution here is λ ↦ λ*, NOT the F1 palindrome λ ↦ −2σ − λ, so this is a different object
    /// from <c>F1SpectrumStatistics.MaxF1PairingDistance</c> and must not be replaced by it. What it does share
    /// with the F1 scorers is the matcher, and that is the part repaired: the match is now MULTISET matching WITH
    /// REMOVAL, each conjugate partner consumed once. Without removal a single self-conjugate (real) eigenvalue
    /// could answer for arbitrarily many λ, which inflates the fraction in one direction only. Both readers of
    /// this number test it against ≈ 0, so the old defect could never have manufactured their conclusion, only
    /// hidden a real symmetry. Do not read the repair as making the number exact in both directions:
    /// greedy-with-removal is not a maximum-cardinality matching and is order-dependent, so on a
    /// near-degenerate but genuinely conjugation-symmetric spectrum it can now come out BELOW 1. That
    /// is the ≈ 1 direction this summary opens with, and it is what the repair leaves open; both live
    /// consumers read against ≈ 0 and are unaffected.</para></summary>
    public static double ConjugationMatchFraction(int n, int wKet, int wBra, double q, double delta, double[]? field)
    {
        var spec = Matrix<Complex>.Build
            .DenseOfArray(WeightCoherenceBlock.Build(n, wKet, wBra, new Complex(q, 0), delta, field))
            .Evd().EigenValues;
        const double tol = 1e-6;
        var taken = new bool[spec.Count];
        int matched = 0;
        for (int i = 0; i < spec.Count; i++)
        {
            var c = Complex.Conjugate(spec[i]);
            int best = -1;
            double bestDist = double.MaxValue;
            for (int j = 0; j < spec.Count; j++)
            {
                if (taken[j]) continue;
                double d = (spec[j] - c).Magnitude;
                if (d < bestDist) { bestDist = d; best = j; }
            }
            if (best >= 0 && bestDist < tol) { taken[best] = true; matched++; }
        }
        return spec.Count == 0 ? double.NaN : (double)matched / spec.Count;
    }
}
