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
/// order as <c>simulations/eta_ladder_chain.py</c>; the pass criteria are identical except the
/// image-eigenvector gate, where the fixed 1e-12 of the Python gate is replaced by the error-model law
/// below (the case-2 formulation the house tolerance rule asks for).
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
