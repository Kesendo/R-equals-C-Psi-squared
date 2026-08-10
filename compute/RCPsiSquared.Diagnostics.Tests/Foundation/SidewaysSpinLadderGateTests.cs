using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The C# gate of the <c>sideways_spin_ladder</c> open arc (cockpit rule 5: the evidence lives
/// as a witness/gate, not only as a Python script): asserts on the shared chain walk
/// <see cref="SidewaysSpinLadderChain"/> (one construction for gate AND live witness), the same checks and
/// order as <c>simulations/eta_ladder_chain.py</c> at odd N and <c>simulations/eta_ladder_chain_n4.py</c>
/// at N=4 (<see cref="Gate_N4"/>, walking all four complex-λ loci); the pass criteria are identical except
/// the image-eigenvector gate, where the fixed 1e-12 of the Python gates is replaced by the error-model
/// law below (the case-2 formulation the house tolerance rule asks for).
///
/// <para>CONTROL — Φ and S⁺ both intertwine L on every rung of the two chains p+q̃ = N∓1, residual
/// compared to 0.0 EXACTLY (the identity cancels pairwise; the Python gate holds the same comparison; a
/// nonzero here is a finding about the construction, not a tolerance case). SHAPE — the interior blocks
/// carrying the cross-fold partner −λ_A − 2N are exactly the two chain interiors. COUNT — fold + band
/// sectors = 4N−8, the F125 orbit size at odd N. LADDER — the transport norms are the Clebsch-Gordan
/// coefficients √(ℓ(ℓ+1) − m(m+1)), ℓ = (N−3)/2, at 1e-6 with no free parameter; the chain terminates at
/// the eigensolver floor into the boundary block (highest weight S⁺|ℓ,ℓ⟩ = 0, a RATIO gate) with the
/// survivors negative control. Plus the σ_min rejection: F125's pinned σ_min values are NOT the CG
/// coefficients and must not be quoted as a confirmation.</para></summary>
public class SidewaysSpinLadderGateTests
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private readonly ITestOutputHelper _out;
    public SidewaysSpinLadderGateTests(ITestOutputHelper o) => _out = o;

    private static string F(double v) => v.ToString("E1", Inv);

    [Fact(DisplayName = "the chain-walk seeds are recorded census entries (RealDefectiveSeeds, R-parity +1)")]
    [Trait("Category", "SIDEWAYS")]
    public void Seeds_AreRecordedCensusEntries()
    {
        foreach (int n in new[] { 5, 7 })
        {
            var (qStar, lamA) = SidewaysSpinLadderChain.Seed(n);
            Assert.Contains(RealDefectiveSeeds.All,
                s => s.N == n && s.QStar == qStar && s.LambdaA == lamA && s.RParity == +1);
        }
        Assert.Throws<ArgumentOutOfRangeException>(() => SidewaysSpinLadderChain.Seed(6));
    }

    [Fact(DisplayName = "sideways ladder gate, N=5: CONTROL exact, SHAPE, COUNT 4N−8, CG norms √2,√2")]
    [Trait("Category", "SIDEWAYS")]
    public void Gate_N5() => RunGate(5);

    [Fact(DisplayName = "sideways ladder gate, N=7: CONTROL exact, SHAPE, COUNT 4N−8, CG norms 2,√6,√6,2")]
    [Trait("Category", "SLOW_SIDEWAYS")]
    public void Gate_N7() => RunGate(7);

    private void RunGate(int n)
    {
        var run = SidewaysSpinLadderChain.Run(n);
        _out.WriteLine($"N = {n}, γ = 1 uniform, q* = {run.QStar.ToString(Inv)}, J_b = {(2 * run.QStar).ToString(Inv)}");
        _out.WriteLine($"  λ_A = {run.LambdaA.ToString(Inv)}, cross-fold partner = {run.Fold.ToString("0.####", Inv)}, centre −N");

        _out.WriteLine($"  CONTROL worst residual (both ladders, all rungs): {F(run.WorstControlResidual)}");
        Assert.True(run.WorstControlResidual == 0.0,
            $"CONTROL: a ladder failed to intertwine exactly, worst residual {run.WorstControlResidual:E3}");

        _out.WriteLine($"  SHAPE: {run.FoldSet.Count} fold sectors, tolerance {F(run.Tolerance)}");
        Assert.True(run.FoldSet.SequenceEqual(run.PredictedSet),
            $"SHAPE: fold sectors [{string.Join(", ", run.FoldSet)}] != chain interiors [{string.Join(", ", run.PredictedSet)}]");

        _out.WriteLine($"  COUNT: {run.FoldSet.Count} + {run.BandSet.Count} = {run.FoldSet.Count + run.BandSet.Count} (4N−8 = {4 * n - 8})");
        Assert.True(run.FoldSet.Count + run.BandSet.Count == 4 * n - 8,
            $"COUNT: fold + band = {run.FoldSet.Count} + {run.BandSet.Count} != 4N−8 = {4 * n - 8}");

        foreach (var chain in run.Chains)
        {
            foreach (var rung in chain.Rungs)
            {
                // Case-2 tolerance (no exact route: v comes from zgeev), so the gate is a LAW, not a bare
                // number: residual / (eps·√targetDim) must stay O(10) across N. Measured range ≈9–13 at
                // N=5 (2e-14 on target dims 100/50) and ≈8–30 at N=7 (4.6e-14..8.0e-14 on target dims
                // 735/1225/147; the small 147-dim rung into the boundary-adjacent block scales largest,
                // its residual being set by the 147² source's conditioning, not by its target); gate < 100.
                double scaled = rung.EigResidual / (Math.Pow(2, -52) * Math.Sqrt(rung.TargetDim));
                _out.WriteLine($"    ({rung.P},{rung.Q})→({rung.P + 1},{rung.Q - 1}) image eigenvector residual "
                               + $"{F(rung.EigResidual)} = {scaled.ToString("0.0", Inv)}·eps·√dim, ‖S⁺v‖ = {rung.Norm.ToString("0.000000", Inv)}");
                Assert.True(scaled < 100.0,
                    $"LADDER p+q={chain.Total}: image at ({rung.P},{rung.Q}) is not an eigenvector, residual "
                    + $"{rung.EigResidual:E3} = {scaled:0.0}·eps·√dim (law: O(10), measured ≈8–30 at N=5,7)");
            }
            var norms = chain.Rungs.Select(r => r.Norm).ToList();
            _out.WriteLine($"  LADDER p+q={chain.Total}: norms [{string.Join(", ", norms.Select(x => x.ToString("0.000000", Inv)))}]"
                           + $" vs CG [{string.Join(", ", run.CgPredicted.Select(x => x.ToString("0.000000", Inv)))}]");
            // 1e-6 is the ported Python gate's own criterion (eta_ladder_chain.py); the measured
            // agreement is ~1e-13, eigensolver-limited, so the ported threshold is seven decades slack.
            Assert.True(norms.Count == run.CgPredicted.Count
                        && norms.Zip(run.CgPredicted).All(t => Math.Abs(t.First - t.Second) < 1e-6),
                $"LADDER p+q={chain.Total}: norms [{string.Join(", ", norms)}] are not the CG coefficients for ℓ = {(n - 3) / 2.0}");

            _out.WriteLine($"  LADDER p+q={chain.Total}: terminal ‖S⁺v‖ = {F(chain.Terminal)}, ratio to interior = "
                           + $"{F(chain.TerminalRatio)}; survivors {chain.Survivors} of {chain.TerminalSourceDim} "
                           + $"(target dim {chain.TerminalTargetDim})");
            Assert.True(chain.Survivors >= chain.TerminalTargetDim,
                $"LADDER p+q={chain.Total}: only {chain.Survivors} survivors < target dim {chain.TerminalTargetDim} (negative control failed)");
            Assert.True(chain.TerminalRatio < 1e-12,
                $"LADDER p+q={chain.Total}: fold vector not among the dying, ratio {chain.TerminalRatio:E3} (highest weight S⁺|ℓ,ℓ⟩ = 0)");
        }
    }

    // ------------------------------------------------------------------------- the N=4 walk

    [Fact(DisplayName = "the N=4 walk inputs are the four real-q defective loci (ReferenceDefectiveLoci)")]
    [Trait("Category", "SLOW_SIDEWAYS")]
    public void SeedsN4_AreTheRealDefectiveLoci()
    {
        // The N=4 inputs are LOCI, not census seeds (their λ are complex; RealDefectiveSeeds is the
        // real-λ census object and lists odd N only), so the source of record is the scanned,
        // conjugation-closed defective set. 1e-4 is ReferenceDefectiveLoci's own merge tolerance.
        // SLOW_SIDEWAYS: For(4) runs the fine-cell EP scan (~10 s), while Gate_N4 itself is ~80 ms.
        var loci = ReferenceDefectiveLoci.For(4);
        Assert.Equal(4, SidewaysSpinLadderChain.SeedsN4.Count);
        foreach (var qStar in SidewaysSpinLadderChain.SeedsN4)
            Assert.Contains(loci, q => Math.Abs(q.Imaginary) < 1e-6 && Math.Abs(q.Real - qStar) < 1e-4);
        Assert.Equal(4, loci.Count(q => Math.Abs(q.Imaginary) < 1e-6));   // and none is missing
    }

    [Fact(DisplayName = "sideways ladder gate, N=4: band+fold conjugate pair share the 4-orbit, four ℓ=1/2 doublets")]
    [Trait("Category", "SIDEWAYS")]
    public void Gate_N4()
    {
        foreach (var qStar in SidewaysSpinLadderChain.SeedsN4)
        {
            var run = SidewaysSpinLadderChain.RunN4(qStar);
            _out.WriteLine($"N = 4, γ = 1 uniform, q* = {qStar.ToString(Inv)}, J_b = {(2 * qStar).ToString(Inv)}");
            _out.WriteLine($"  λ* = {run.LambdaStar.Real.ToString("0.000000", Inv)} "
                           + $"{(run.LambdaStar.Imaginary >= 0 ? "+" : "−")} {Math.Abs(run.LambdaStar.Imaginary).ToString("0.000000", Inv)}i "
                           + $"(pair split {F(run.PairSplit)}, tol {F(run.Tolerance)})");

            _out.WriteLine($"  CONTROL worst residual (both ladders, all rungs): {F(run.WorstControlResidual)}");
            Assert.True(run.WorstControlResidual == 0.0,
                $"CONTROL: a ladder failed to intertwine exactly, worst residual {run.WorstControlResidual:E3}");

            // The self-fold pin Re λ* = −4, read honestly (the RunN4 docstring carries the derivation
            // and the sweep numbers): it holds machine-exact (9e-16..1.3e-14) wherever the closest pair
            // is σ-stable, which is NOT generic at N=4 — over ~27% of q ∈ [0.05, 3] a non-stable pair
            // is closest, its σ-image ties it exactly, and the pin fails by O(1). It certifies
            // σ-stability of the closest pair, never defectiveness; the loci labels are inherited from
            // ReferenceDefectiveLoci. fold(λ*) = conj(λ*) is printed, not separately gated:
            // |fold − conj λ*| = 2·pin exactly, the pin restated.
            _out.WriteLine($"  self-fold pin |Re λ* + 4| = {F(run.SelfFoldPinResidual)}, "
                           + $"|fold − conj λ*| = {F(run.FoldVsConjDistance)}");
            Assert.True(run.SelfFoldPinResidual < 1e-6,
                $"self-fold pin: Re λ* = {run.LambdaStar.Real} is not −4 (residual {run.SelfFoldPinResidual:E3})");

            // BOTH VALUES in every orbit block, and COUNT as the union: fold set = band set = the
            // confined 4-orbit, |fold ∪ band| = 4N−12 = 4. The odd-N sum |fold| + |band| = 8 = 4N−8
            // still holds numerically at N=4, but it double-counts the shared sectors; the union is
            // the count that carries the even-N accounting. Read at real strength: block (1,2)'s own
            // carries-both check is automatic (λ* is its closest-pair mean, d = split/2 = tol/10 by
            // construction), and the fold half is a corollary of the λ half in EVERY block (spectra
            // are conjugation-closed, the pin makes fold = conj λ*, so d(fold) tracks d(λ*) to
            // ~1e-11) — the measured content is the λ-sharing across (2,1)/(2,3)/(3,2) and the union
            // gate below, the most q-selective assertion (pass-windows: ±1.3e-4 around the first two
            // loci, ±3e-5 around the twin 0.857458, −0.076/+0.101 around 1.738181; the pin adds
            // ±1.5e-3 at 0.460212 and one-sided ~2e-3..5e-3 lower bounds at the twins; the RunN4
            // docstring carries the calibration). foldSet == bandSet is kept as a consistency check, not as evidence: it
            // follows from the same corollary.
            foreach (var (block, dLam, dFold) in run.OrbitCarries)
            {
                _out.WriteLine($"  block {block}: d(λ*) = {F(dLam)}, d(fold) = {F(dFold)}");
                Assert.True(dLam < run.Tolerance && dFold < run.Tolerance,
                    $"block {block} does not carry both values (d(λ*) = {dLam:E3}, d(fold) = {dFold:E3}, tol {run.Tolerance:E3})");
            }
            var union = run.FoldSet.Union(run.BandSet).OrderBy(k => k).ToList();
            _out.WriteLine($"  COUNT: |fold ∪ band| = {union.Count} (4N−12 = 4); fold set [{string.Join(", ", run.FoldSet)}], "
                           + $"band set [{string.Join(", ", run.BandSet)}]");
            Assert.True(union.SequenceEqual(run.Orbit),
                $"COUNT: fold ∪ band = [{string.Join(", ", union)}] != the confined orbit [{string.Join(", ", run.Orbit)}]");
            Assert.True(run.FoldSet.SequenceEqual(run.BandSet),
                $"fold set [{string.Join(", ", run.FoldSet)}] != band set [{string.Join(", ", run.BandSet)}] (they share all sectors at N=4)");

            // The four ℓ = 1/2 doublets: two chains × two conjugate values, CG norm 1 from the closed
            // form, for both split-pair members AND their mix (the twins are near-parallel here, so
            // the mix agreement is a smoke check, not the N=9 robustness; the defect is printed).
            // Structure caveat, sharpened from the σ_min fence: S⁺(4,1,2) has singular spectrum
            // {1×20, 2×4} and S⁺(4,2,1) is rank 4 at √3, so these norm/terminal/survivors readings
            // hold for most eigenvectors and at any generic q ≳ 0.58 (below, the closest pair
            // transports at norm 2 — the 0.460212 locus is a ~1e-3 norm-1 island there, the one locus
            // where this reading is itself selective); the run-specific content is the λ-sharing and
            // the union gate above. Survivors gated == dim(target), the N=4 Python sibling's own
            // criterion (the odd-N gate's ≥ is ITS sibling's). Kept as port fidelity.
            double cg = SidewaysSpinLadderChain.CgCoefficients(4).Single();   // ℓ = 1/2: exactly 1
            Assert.Equal(4, run.Doublets.Count);
            foreach (var d in run.Doublets)
            {
                string name = $"{d.B0}→{d.B1} {(d.AtFoldValue ? "fold" : "λ*  ")}";
                _out.WriteLine($"    {name}: ‖S⁺m1‖ = {d.NormM1.ToString("0.000000000", Inv)}, "
                               + $"‖S⁺m2‖ = {d.NormM2.ToString("0.000000000", Inv)}, ‖S⁺mix‖ = {d.NormMix.ToString("0.000000000", Inv)}, "
                               + $"1−|overlap| = {F(d.OverlapDefect)}");
                foreach (var (tag, nrm) in new[] { ("m1", d.NormM1), ("m2", d.NormM2), ("mix", d.NormMix) })
                    Assert.True(Math.Abs(nrm - cg) < 1e-6,
                        $"{name} {tag}: ‖S⁺v‖ = {nrm:0.000000000} is not the CG norm {cg} (ℓ = 1/2)");

                // Case-2 tolerance as in the odd-N gate: residual / (eps·√targetDim), gated O(10).
                // At N=4 every target dim is 24, so nothing varies across decades here; the LAW
                // reading lives in the odd-N gate where dims span 50..1225, this is its threshold.
                double scaled = d.EigResidual / (Math.Pow(2, -52) * Math.Sqrt(d.RungTargetDim));
                _out.WriteLine($"    {name}: image eigenvector residual {F(d.EigResidual)} = "
                               + $"{scaled.ToString("0.0", Inv)}·eps·√dim; terminal ratio {F(d.TerminalRatio)}; "
                               + $"survivors {d.Survivors} (target dim {d.TerminalTargetDim})");
                Assert.True(scaled < 100.0,
                    $"{name}: image is not an eigenvector, residual {d.EigResidual:E3} = {scaled:0.0}·eps·√dim");
                Assert.True(d.TerminalRatio < 1e-12,
                    $"{name}: transported vector not among the dying at {d.B1}→{d.B2}, ratio {d.TerminalRatio:E3}");
                Assert.True(d.Survivors == d.TerminalTargetDim,
                    $"{name}: survivors {d.Survivors} != dim(target) = {d.TerminalTargetDim}");
            }
        }
    }

    // ---------------------------------------------------------------------- the σ_min rejection

    // A review round proposed that F125's pinned σ_min(W) values (1 at N=4, √2 at N=5) ARE the
    // η-side CG coefficients, a free two-point confirmation of the multiplet reading. Measured
    // before adopting, and it is false: σ_min is the smallest singular value of the WHOLE rung
    // map, a block carries several multiplets, so it reads the weakest direction present; it
    // coincides where the block is small. Gated because the open arc quotes the disagreement in
    // order to REJECT that proposal, and a rejection needs a route back too.
    [Fact(DisplayName = "σ_min of an η rung is NOT the CG coefficient (N=7 middle rungs: √2 vs √6)")]
    [Trait("Category", "SIDEWAYS")]
    public void SigmaMin_IsNotTheCgCoefficient()
    {
        var agree = new List<(int N, int P, double Sv, double Cg)>();
        var differ = new List<(int N, int P, double Sv, double Cg)>();
        foreach (int n in new[] { 4, 5, 7 })
        {
            double spinL = (n - 3) / 2.0;
            double mw = -spinL;
            for (int p = 1; p <= n - 3; p++)         // guard p+2 ≤ n−1, matching the Python loop
            {
                var sv = SpectatorIntertwiner.BuildW(n, p, p + 1).Svd(false).S;
                double sMin = sv[sv.Count - 1].Real;
                double cgVal = Math.Sqrt(spinL * (spinL + 1) - mw * (mw + 1));
                bool same = Math.Abs(sMin - cgVal) < 1e-9;
                (same ? agree : differ).Add((n, p, sMin, cgVal));
                _out.WriteLine($"  N={n} ({p},{p + 1})→({p + 1},{p + 2}): σ_min = {sMin.ToString("0.000000", Inv)}  "
                               + $"CG = {cgVal.ToString("0.000000", Inv)}  {(same ? "equal" : "DIFFERENT")}");
                mw += 1.0;
            }
        }
        Assert.True(differ.Count > 0, "σ_min unexpectedly equals the CG coefficient on every rung");
        // the last rung starts at p = N−3, not N−2: its target block is (N−2, N−1)
        Assert.True(agree.All(t => t.N is 4 or 5 || t.P == 1 || t.P == t.N - 3),
            "an interior large-N rung agrees: the σ_min ≠ CG rejection would be unsound");
    }
}
