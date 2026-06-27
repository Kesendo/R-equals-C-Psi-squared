using System;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The F89 Door-C CSR sweep (docs/superpowers/plans/2026-06-27-f89-door-c-csr-integrability-sweep.md):
/// does breaking the (SE,DE) Liouvillian block's free-fermion additivity (XXZ anisotropy Δ) drive the
/// fixed-q complex spacing ratio from Poisson toward Ginibre? Methodology per review round 2: pool the
/// per-spectrum z-values over the q-grid (never raw eigenvalues), bootstrap the CI, and compare against
/// finite-size-MATCHED Poisson/GinUE references (not the asymptotic 0.658/0.738).</summary>
public class IntegrabilityBreakingCsrTests
{
    private readonly ITestOutputHelper _out;
    public IntegrabilityBreakingCsrTests(ITestOutputHelper output) => _out = output;

    private static double[] Qs(int count)
    {
        var q = new double[count];
        for (int i = 0; i < count; i++) q[i] = 0.3 + i * (3.7 / (count - 1));   // q ∈ [0.3, 4.0], matches the witness
        return q;
    }

    /// <summary>Regression anchor: at Δ=0 the H_B-mixed pooled CSR reproduces the galoischaos witness
    /// baseline — Poisson-like / sub-Poisson (⟨|z|⟩ &lt; 0.72), NOT GinUE (0.738). Validates the pooled-z
    /// harness against the trusted Δ=0 witness before any Δ-sweep claim.</summary>
    [Fact]
    public void DeltaZero_HbMixed_ReproducesPoissonBaseline()
    {
        var r = IntegrabilityBreakingCsr.Sweep(n: 7, delta: 0.0, qs: Qs(20), half: IntegrabilityBreakingCsr.Half.HbMixed);
        Assert.True(r.ZCount > 200, $"need a real pool, got {r.ZCount}");
        Assert.InRange(r.MeanAbs, 0.40, 0.72);     // Poisson-like, NOT GinUE 0.738
    }

    /// <summary>Finite-size-matched references at the measurement's per-spectrum size must still SEPARATE
    /// the two classes (Poisson &lt; GinUE), so a Δ-sweep verdict reads against the right yardsticks rather
    /// than the asymptotic values that carry the wrong edge bias for ~50-point spectra.</summary>
    [Fact]
    public void FiniteSizeReferences_StillSeparate()
    {
        var pois = IntegrabilityBreakingCsr.PoissonReference(size: 45, draws: 60, seed: 1);
        var gin = IntegrabilityBreakingCsr.GinueReference(size: 45, draws: 60, seed: 2);
        Assert.True(gin.MeanAbs > pois.MeanAbs + 0.02,
            $"finite-size refs must separate: Poisson {pois.MeanAbs:F3} vs GinUE {gin.MeanAbs:F3}");
    }

    /// <summary>Reconnaissance: the headline sweep — pooled H_B-mixed ⟨|z|⟩ (with 95% bootstrap CI) vs Δ,
    /// against finite-size-matched Poisson/GinUE references. Does breaking free-fermion additivity drive
    /// the fixed-q CSR toward Ginibre, or do you need to break the Hamiltonian's integrability too?</summary>
    [Fact]
    public void Reconnaissance_HbMixedCsrVsDelta()
    {
        var qs = Qs(24);
        double[] deltas = { 0.0, 0.25, 0.5, 1.0, 2.0 };
        foreach (int n in new[] { 6, 7 })
        {
            _out.WriteLine($"N={n}  (domain valid per Δ: UpperHalf at Δ=0 [conj-symmetric], " +
                           "OffReal at Δ>0 [conj-symmetry broken]; refs size-matched per Δ)");
            _out.WriteLine("   Δ    | dom   | ⟨|z|⟩  [95% CI]        | ⟨cosθ⟩ | Npool | refs P/G");
            _out.WriteLine("  ------|-------|-----------------------|--------|-------|----------");
            foreach (double d in deltas)
            {
                var dom = d == 0.0 ? IntegrabilityBreakingCsr.Domain.UpperHalf : IntegrabilityBreakingCsr.Domain.OffReal;
                var r = IntegrabilityBreakingCsr.Sweep(n, d, qs, IntegrabilityBreakingCsr.Half.HbMixed, dom);
                int perSpec = Math.Max(10, r.ZCount / qs.Length);
                var pRef = IntegrabilityBreakingCsr.PoissonReference(perSpec, draws: 60, seed: 100 + n);
                var gRef = IntegrabilityBreakingCsr.GinueReference(perSpec, draws: 60, seed: 200 + n);
                _out.WriteLine($"  {d,5:F2} | {(dom == IntegrabilityBreakingCsr.Domain.UpperHalf ? "upper" : "offRl")} | " +
                               $"{r.MeanAbs:F3} [{r.CiLo:F3}, {r.CiHi:F3}] | {r.MeanCos,+6:F3} | {r.ZCount,5} | " +
                               $"{pRef.MeanAbs:F3}/{gRef.MeanAbs:F3}");
            }
            _out.WriteLine("");
        }
    }

    [Theory]
    [InlineData(5)]   // 5·C(5,2) = 50
    [InlineData(7)]   // 7·C(7,2) = 147
    public void FullSpectrum_HasNTimesCN2Eigenvalues(int n)
    {
        int expected = n * (n * (n - 1) / 2);
        Assert.Equal(expected, IntegrabilityBreakingCsr.FullSpectrum(n, q: 2.0, delta: 0.0).Length);
    }

    /// <summary>Symmetry-class check (review round 2 #5): the CSR upper-half-plane restriction is only
    /// valid if the (SE,DE) spectrum is conjugation-symmetric (λ → λ*). The reviewer's heuristic says Π
    /// maps (SE,DE) → the conjugate (DE,SE) block, so (SE,DE) alone may NOT be self-conjugate — which
    /// would put us in class A (GinUE reference) and make the upper-half restriction questionable. We
    /// measure it: the conjugation-match fraction, the real-axis fraction, and a reflection about the AT
    /// midpoint Re=−4 (λ → −8−λ and λ → −8−λ*), at the Δ=0 baseline and the Δ=0.5 peak.</summary>
    [Fact]
    public void Reconnaissance_SpectrumSymmetry_N7()
    {
        const double tol = 1e-6;
        foreach (var (q, d) in new[] { (2.0, 0.0), (2.0, 0.5), (0.5, 0.5) })
        {
            var spec = IntegrabilityBreakingCsr.FullSpectrum(7, q, d);
            int nReal = spec.Count(z => Math.Abs(z.Imaginary) < tol);
            int nUpper = spec.Count(z => z.Imaginary > tol);

            bool InSpec(Complex w) => spec.Any(z => (z - w).Magnitude < tol);
            double conjFrac = spec.Count(z => InSpec(Complex.Conjugate(z))) / (double)spec.Length;
            double reflFrac = spec.Count(z => InSpec(new Complex(-8, 0) - z)) / (double)spec.Length;
            double reflConjFrac = spec.Count(z => InSpec(new Complex(-8, 0) - Complex.Conjugate(z))) / (double)spec.Length;

            _out.WriteLine($"q={q}, Δ={d}: total={spec.Length}, real={nReal}, upper={nUpper} | " +
                           $"conj-sym(λ→λ*)={conjFrac:P0}, refl(−8−λ)={reflFrac:P0}, refl-conj(−8−λ*)={reflConjFrac:P0}");
        }
    }

    /// <summary>Stationarity insurance (review round 2 #4, physics-first #6): the pooled ⟨|z|⟩ is only
    /// meaningful if ⟨|z|⟩(q) is roughly flat across the q-grid. A large spread means the pool averages a
    /// q-varying population (e.g. near an EP/discriminant locus) and the pooled mean is suspect. We print
    /// the per-q ⟨|z|⟩ spread at the Δ=0 baseline, the Δ=0.5 peak, and the Δ=2 sub-Poisson tail.</summary>
    [Fact]
    public void Reconnaissance_Stationarity_PerQ_N7()
    {
        var qs = Qs(24);
        foreach (double d in new[] { 0.0, 0.5, 2.0 })
        {
            var dom = d == 0.0 ? IntegrabilityBreakingCsr.Domain.UpperHalf : IntegrabilityBreakingCsr.Domain.OffReal;
            var perQ = IntegrabilityBreakingCsr.PerQMeanAbs(7, d, qs, IntegrabilityBreakingCsr.Half.HbMixed, dom);
            var valid = perQ.Where(x => !double.IsNaN(x)).ToArray();
            double mean = valid.Average();
            double sd = Math.Sqrt(valid.Select(x => (x - mean) * (x - mean)).Average());
            _out.WriteLine($"Δ={d,4:F2}: per-q ⟨|z|⟩ mean={mean:F3} sd={sd:F3} " +
                           $"min={valid.Min():F3} max={valid.Max():F3} (n_q valid={valid.Length}/{qs.Length})");
            _out.WriteLine("   " + string.Join(" ", perQ.Select(x => double.IsNaN(x) ? " NaN" : $"{x:F2}")));
            _out.WriteLine("");
        }
    }
}
