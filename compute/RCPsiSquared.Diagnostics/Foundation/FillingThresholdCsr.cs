using System;
using System.Collections.Generic;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The F89 Door-C DECISIVE follow-up: is fixed-q dissipative quantum chaos (GinUE complex spacing ratio)
/// a FILLING threshold rather than an integrability one? Door-C stages 1-2 (<see cref="IntegrabilityBreakingCsr"/>)
/// found the DILUTE (SE,DE)=(1,2) coherence block stays Poisson / non-GinUE under every integrability-breaking knob
/// (XXZ Δ, a random Z-field), because a 2-excitation sector cannot thermalize. This harness builds the GENERAL
/// (wKet,wBra) coherence block (<see cref="WeightCoherenceBlock.Build(int,int,int,Complex,double,double[])"/>) at
/// EXTENSIVE filling (wKet,wBra near N/2) and re-runs the same disordered CSR. If the DENSE block reaches GinUE
/// while the dilute one does not, Door-C's null is structural/kinematic (a FILLING threshold), not algebraic
/// (Galois/integrability).
///
/// <para>Class A is kept by using UNEQUAL weight (p, p+1): the F1 palindrome Π / conjugation maps the (p,p+1)
/// coherence block to the conjugate (p+1,p) block, NOT to itself, so there is no residual antiunitary — the GinUE
/// (class A) reference 0.738/−0.24 is the right target (methodology #5; confirmed empirically by
/// <see cref="ConjugationMatchFraction"/> ≈ 0 once a field/Δ breaks conjugation symmetry). Methodology inherited
/// from the Door-C harness: pool per-spectrum z's (never raw eigenvalues), bootstrap the CI (shared
/// <see cref="IntegrabilityBreakingCsr.Reduce"/>), and read in the OffReal domain (|Im| &gt; tol), valid once
/// conjugation symmetry is broken.</para></summary>
public static class FillingThresholdCsr
{
    private const double ImTol = 1e-6;

    /// <summary>Off-real eigenvalues (|Im| &gt; tol) of a coherence block — the valid CSR domain once conjugation
    /// symmetry is broken (Δ≠0 or a random field), apples-to-apples with the GinUE reference (a full complex cloud).</summary>
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
    /// field, and pool the OffReal per-spectrum z's; bootstrap a 95% CI. w=0 is deterministic (the clean block, every
    /// realization identical). Δ=0 is free-fermion + disorder (Anderson-like); Δ≠0 is interacting + disorder.</summary>
    public static IntegrabilityBreakingCsr.CsrReading DisorderSweep(
        int n, int wKet, int wBra, double q, double delta, double w, int realizations, int seed)
    {
        var rng = new Random(seed);
        var pool = new List<Complex>();
        for (int r = 0; r < realizations; r++)
        {
            double[]? field = null;
            if (w != 0.0)
            {
                field = new double[n];
                for (int k = 0; k < n; k++) field[k] = (2 * rng.NextDouble() - 1) * w;   // U[−w, w]
            }
            pool.AddRange(ComplexSpacingRatio.ZValues(OffRealSpectrum(n, wKet, wBra, q, delta, field)));
        }
        return IntegrabilityBreakingCsr.Reduce(pool, seed + 7919);
    }

    /// <summary>The clean (disorder-free) pooled CSR of the (wKet,wBra) block at Δ, pooled over the q-grid. The
    /// integrable/no-disorder control: at Δ=0 the block is free-fermion; at Δ≠0 the Hamiltonian is Bethe-integrable
    /// but the Liouvillian is not free-fermion (the Δ·ZZ breaks the additivity). OffReal throughout.</summary>
    public static IntegrabilityBreakingCsr.CsrReading CleanSweep(int n, int wKet, int wBra, double[] qs, double delta)
    {
        var pool = new List<Complex>();
        foreach (var q in qs)
            pool.AddRange(ComplexSpacingRatio.ZValues(OffRealSpectrum(n, wKet, wBra, q, delta, null)));
        return IntegrabilityBreakingCsr.Reduce(pool, 13337);
    }

    /// <summary>The fraction of eigenvalues λ whose conjugate λ* is also in the spectrum (within tol). ≈ 1 ⟹
    /// conjugation-symmetric (free-fermion Δ=0, no field), ≈ 0 ⟹ broken (a field or Δ≠0). The class-A guard: a
    /// near-zero fraction confirms no residual antiunitary, so the GinUE (class A) reference is the right target.</summary>
    public static double ConjugationMatchFraction(int n, int wKet, int wBra, double q, double delta, double[]? field)
    {
        var spec = Matrix<Complex>.Build
            .DenseOfArray(WeightCoherenceBlock.Build(n, wKet, wBra, new Complex(q, 0), delta, field))
            .Evd().EigenValues;
        const double tol = 1e-6;
        int matched = 0;
        for (int i = 0; i < spec.Count; i++)
        {
            var c = Complex.Conjugate(spec[i]);
            for (int j = 0; j < spec.Count; j++)
                if ((spec[j] - c).Magnitude < tol) { matched++; break; }
        }
        return spec.Count == 0 ? double.NaN : (double)matched / spec.Count;
    }
}
