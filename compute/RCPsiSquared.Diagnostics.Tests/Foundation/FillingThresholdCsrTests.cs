using System;
using System.Linq;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The F89 Door-C DECISIVE follow-up: is fixed-q dissipative chaos (GinUE) a FILLING threshold, not an
/// integrability one? Door-C stages 1-2 showed the DILUTE (SE,DE)=(1,2) coherence block stays Poisson / non-GinUE
/// under every integrability-breaking knob (Δ, disorder), because a 2-excitation sector cannot thermalize. This
/// harness builds the general (wKet,wBra) coherence block at extensive filling (wKet,wBra near N/2) and re-runs the
/// disordered CSR: if the DENSE block DOES reach GinUE under the SAME knobs while the dilute one does not, the
/// Door-C null is structural/kinematic (filling), not algebraic/integrability. Unequal weight (p,p+1) keeps the
/// non-Hermitian symmetry class A (Π maps (p,p+1)→the conjugate (p+1,p) block, not a self-symmetry), licensing the
/// GinUE 0.738/−0.24 target. Methodology inherited from IntegrabilityBreakingCsr (pool per-spectrum z's, bootstrap
/// CI, finite-size-matched references, OffReal domain once conjugation symmetry is broken).</summary>
public class FillingThresholdCsrTests
{
    private readonly ITestOutputHelper _out;
    public FillingThresholdCsrTests(ITestOutputHelper output) => _out = output;

    /// <summary>Zero-field anchor: w=0 ⟹ every realization is the identical (q,Δ) block, so the pool scales
    /// linearly with the realization count and ⟨|z|⟩ is unchanged. Catches RNG/pooling bugs.</summary>
    [Fact]
    public void DisorderSweep_ZeroField_IsDeterministicAcrossRealizations()
    {
        var r1 = FillingThresholdCsr.DisorderSweep(6, 3, 4, q: 1.0, delta: 1.0, w: 0.0, realizations: 1, seed: 1);
        var r3 = FillingThresholdCsr.DisorderSweep(6, 3, 4, q: 1.0, delta: 1.0, w: 0.0, realizations: 3, seed: 1);
        Assert.True(r1.ZCount > 100, $"a dense block must give a real per-spectrum pool, got {r1.ZCount}");
        Assert.Equal(r1.ZCount * 3, r3.ZCount);            // identical copies pool linearly
        Assert.Equal(r1.MeanAbs, r3.MeanAbs, 9);           // no disorder ⟹ identical statistic
    }

    /// <summary>The dilute control, reproduced through the GENERAL builder: the (1,2)=(SE,DE) block under
    /// interacting disorder (Δ=1, intermediate W) must stay NON-GinUE — no strong angular repulsion (⟨cosθ⟩ not
    /// near GinUE's −0.24). This is the Door-C Stage-2 result re-derived via WeightCoherenceBlock, the apples-to-
    /// apples baseline the dense block is compared against.</summary>
    [Fact]
    public void DiluteBlock_UnderInteractingDisorder_StaysNonGinUE()
    {
        var r = FillingThresholdCsr.DisorderSweep(7, 1, 2, q: 2.0, delta: 1.0, w: 1.0, realizations: 120, seed: 1000);
        Assert.True(r.ZCount > 200, $"need a real pool, got {r.ZCount}");
        Assert.True(r.MeanCos > -0.08, $"dilute block should NOT show GinUE angular repulsion, got ⟨cosθ⟩={r.MeanCos:F3}");
    }

    /// <summary>Symmetry-class guard (methodology #5): a per-site random field breaks conjugation symmetry, so the
    /// (p,p+1) block spectrum is NOT conjugation-symmetric (match fraction ≈ 0) ⟹ OffReal is the valid domain and
    /// class A (the GinUE reference) is licensed. Without this the CSR target would be AI+/AII+, not 0.738/−0.24.</summary>
    [Fact]
    public void RandomField_BreaksConjugationSymmetry_LicensingClassA()
    {
        var rng = new Random(7);
        var field = Enumerable.Range(0, 6).Select(_ => 2 * rng.NextDouble() - 1).ToArray();
        double frac = FillingThresholdCsr.ConjugationMatchFraction(6, 3, 4, q: 1.0, delta: 1.0, field: field);
        Assert.True(frac < 0.1, $"random field must break conjugation symmetry (class A), but match fraction was {frac:P0}");
    }

    /// <summary>The headline reconnaissance: the FILLING LADDER. At fixed N, walk (wKet,wBra) from the dilute
    /// (1,2) up to near-half-filling, at Δ∈{0,1} and a W-grid, and print ⟨|z|⟩/⟨cosθ⟩ against finite-size-matched
    /// Poisson/GinUE references. The question: does ⟨cosθ⟩ go NEGATIVE (GinUE repulsion) as filling rises — and
    /// does it need disorder, interactions, or both?</summary>
    [Fact]
    [Trait("Category", "SLOW_FILLCSR")]
    public void Reconnaissance_FillingLadder_N6()
    {
        int n = 6; double q = 1.0; int R = 40;
        (int wk, int wb)[] ladder = { (1, 2), (2, 3), (3, 4) };
        double[] ws = { 0.0, 0.5, 1.0, 2.0 };

        foreach (double delta in new[] { 0.0, 1.0 })
        {
            _out.WriteLine($"=== N={n}, q={q}, Δ={delta}, R={R} (disorder), refs finite-size-matched per block ===");
            _out.WriteLine("  (wk,wb) dim |  W   | ⟨|z|⟩ [95% CI]        | ⟨cosθ⟩ | Npool | ref P/G ⟨|z|⟩ | ref G ⟨cos⟩");
            foreach (var (wk, wb) in ladder)
            {
                int dim = Binom(n, wk) * Binom(n, wb);
                int perSpec = Math.Max(10, dim / 2);
                var pRef = IntegrabilityBreakingCsr.PoissonReference(perSpec, draws: 40, seed: 50);
                var gRef = IntegrabilityBreakingCsr.GinueReference(perSpec, draws: 40, seed: 60);
                foreach (double w in ws)
                {
                    int real = w == 0.0 ? 1 : R;
                    var r = FillingThresholdCsr.DisorderSweep(n, wk, wb, q, delta, w, real, seed: 2000);
                    _out.WriteLine($"  ({wk},{wb}) {dim,4} | {w,4:F1} | {r.MeanAbs:F3} [{r.CiLo:F3},{r.CiHi:F3}] | " +
                                   $"{r.MeanCos,+6:F3} | {r.ZCount,5} | {pRef.MeanAbs:F3}/{gRef.MeanAbs:F3}    | {gRef.MeanCos,+6:F3}");
                }
            }
            _out.WriteLine("");
        }
    }

    /// <summary>Does the dense block CONVERGE to GinUE as N grows? Walk the near-half-filling block at N=6,7
    /// (Δ=1, ergodic window W=0.75), with the Poisson/GinUE references size-matched to the ACTUAL per-spectrum
    /// z-count (methodology #3 — the off-real count is ≈ dim, not dim/2). If ⟨cosθ⟩ tracks the GinUE reference up
    /// as N grows while the dilute (1,2) stays flat at 0, the filling-threshold reading lands.</summary>
    [Fact]
    [Trait("Category", "SLOW_FILLCSR")]
    public void Reconnaissance_DenseSizeScaling()
    {
        double delta = 1.0, q = 1.0, W = 0.75;
        (int n, int wk, int wb, int R)[] cases =
        {
            (6, 1, 2, 60),   // the dilute control at N=6
            (6, 3, 4, 30),   // dense N=6
            (7, 1, 2, 60),   // the dilute control at N=7
            (7, 3, 4, 12),   // dense N=7 (1225-dim; R kept low — each spectrum already gives ~1200 z's)
        };
        _out.WriteLine($"Δ={delta}, q={q}, W={W}; references SIZE-MATCHED to the measured per-spectrum z-count");
        _out.WriteLine("  N (wk,wb) dim  | ⟨|z|⟩ [95% CI]        | ⟨cosθ⟩ | perSpec | Poisson z/cos | GinUE z/cos");
        foreach (var (n, wk, wb, R) in cases)
        {
            int dim = Binom(n, wk) * Binom(n, wb);
            var probe = FillingThresholdCsr.DisorderSweep(n, wk, wb, q, delta, W, 1, seed: 555);
            int perSpec = Math.Max(10, probe.ZCount);
            var pRef = IntegrabilityBreakingCsr.PoissonReference(perSpec, draws: 30, seed: 70);
            var gRef = IntegrabilityBreakingCsr.GinueReference(perSpec, draws: 30, seed: 80);
            var r = FillingThresholdCsr.DisorderSweep(n, wk, wb, q, delta, W, R, seed: 3000);
            _out.WriteLine($"  {n} ({wk},{wb}) {dim,4}  | {r.MeanAbs:F3} [{r.CiLo:F3},{r.CiHi:F3}] | " +
                           $"{r.MeanCos,+6:F3} | {perSpec,6} | {pRef.MeanAbs:F3}/{pRef.MeanCos,+5:F3}  | {gRef.MeanAbs:F3}/{gRef.MeanCos,+6:F3}");
        }
    }

    /// <summary>Timing/feasibility probe for the N=8 dense block (4,5)=3920-dim: one disordered spectrum, to gauge
    /// MathNet managed-EVD cost before any multi-realization run. Prints the per-spectrum CSR (one 3920-point cloud
    /// already gives a tight single-spectrum statistic).</summary>
    [Fact]
    [Trait("Category", "SLOW_FILLCSR")]
    public void Reconnaissance_DenseN8_SingleProbe()
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var r = FillingThresholdCsr.DisorderSweep(8, 4, 5, q: 1.0, delta: 1.0, w: 0.75, realizations: 1, seed: 99);
        sw.Stop();
        _out.WriteLine($"N=8 (4,5) dim=3920: one EVD+CSR in {sw.Elapsed.TotalSeconds:F1}s — " +
                       $"⟨|z|⟩={r.MeanAbs:F3} ⟨cosθ⟩={r.MeanCos:+0.000;-0.000} over {r.ZCount} z's");
    }

    /// <summary>The N=8 capstone of the size trend: dense (4,5)=3920 with a real CI (R realizations) and a
    /// size-matched GinUE reference, beside the cheap dilute (1,2) control. ~50s/EVD, so this is a minutes-long run.</summary>
    [Fact]
    [Trait("Category", "SLOW_FILLCSR")]
    public void Reconnaissance_DenseN8_Full()
    {
        double q = 1.0, delta = 1.0, W = 0.75;
        var dilute = FillingThresholdCsr.DisorderSweep(8, 1, 2, q, delta, W, 60, seed: 4001);
        var dense = FillingThresholdCsr.DisorderSweep(8, 4, 5, q, delta, W, 4, seed: 4002);
        var pRef = IntegrabilityBreakingCsr.PoissonReference(Math.Max(10, dense.ZCount / 4), draws: 20, seed: 71);
        var gRef = IntegrabilityBreakingCsr.GinueReference(Math.Max(10, dense.ZCount / 4), draws: 4, seed: 81);
        _out.WriteLine($"N=8, Δ={delta}, q={q}, W={W}");
        _out.WriteLine($"  dilute (1,2) 224 : ⟨|z|⟩={dilute.MeanAbs:F3} [{dilute.CiLo:F3},{dilute.CiHi:F3}] ⟨cosθ⟩={dilute.MeanCos:+0.000;-0.000} (N={dilute.ZCount})");
        _out.WriteLine($"  dense  (4,5) 3920: ⟨|z|⟩={dense.MeanAbs:F3} [{dense.CiLo:F3},{dense.CiHi:F3}] ⟨cosθ⟩={dense.MeanCos:+0.000;-0.000} (N={dense.ZCount})");
        _out.WriteLine($"  refs(size {dense.ZCount / 4}): Poisson ⟨|z|⟩={pRef.MeanAbs:F3} ⟨cos⟩={pRef.MeanCos:+0.000;-0.000} | GinUE ⟨|z|⟩={gRef.MeanAbs:F3} ⟨cos⟩={gRef.MeanCos:+0.000;-0.000}");
    }

    private static int Binom(int n, int k)
    {
        long r = 1;
        for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
        return (int)r;
    }
}
