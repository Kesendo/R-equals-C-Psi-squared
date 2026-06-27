using System;
using System.Linq;
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
            // size the finite-size references at this N's per-spectrum H_B point count.
            var probe = IntegrabilityBreakingCsr.Sweep(n, 0.0, qs, IntegrabilityBreakingCsr.Half.HbMixed);
            int perSpec = Math.Max(10, probe.ZCount / qs.Length);
            var pRef = IntegrabilityBreakingCsr.PoissonReference(perSpec, draws: 80, seed: 100 + n);
            var gRef = IntegrabilityBreakingCsr.GinueReference(perSpec, draws: 80, seed: 200 + n);

            _out.WriteLine($"N={n}  (per-spectrum H_B ≈ {perSpec} pts;  finite-size refs: " +
                           $"Poisson ⟨|z|⟩={pRef.MeanAbs:F3}, GinUE ⟨|z|⟩={gRef.MeanAbs:F3})");
            _out.WriteLine("   Δ    | ⟨|z|⟩  [95% CI]          | ⟨cosθ⟩  | Npool");
            _out.WriteLine("  ------|-------------------------|---------|------");
            foreach (double d in deltas)
            {
                var r = IntegrabilityBreakingCsr.Sweep(n, d, qs, IntegrabilityBreakingCsr.Half.HbMixed);
                _out.WriteLine($"  {d,5:F2} | {r.MeanAbs:F3} [{r.CiLo:F3}, {r.CiHi:F3}]      | " +
                               $"{r.MeanCos,+7:F3} | {r.ZCount}");
            }
            _out.WriteLine("");
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
            var perQ = IntegrabilityBreakingCsr.PerQMeanAbs(7, d, qs, IntegrabilityBreakingCsr.Half.HbMixed);
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
