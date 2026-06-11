using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live F116 witness: the assertions hold against what the witness computes at inspect time
/// (the router re-run and the root-found metallic ratio), not against the closed forms it then compares to.
/// The golden c = 1 ratio is φ, the silver c = 2 ratio is 1 + √2, the ceiling reaches zero on the canonical
/// Z-middle cases, and the witness surfaces the expected child nodes.</summary>
public class GoldenRouterWitnessTests
{
    private const double Phi = 1.6180339887498949;     // (1 + √5) / 2
    private const double Silver = 2.4142135623730951;  // 1 + √2

    private static List<IInspectable> Children(GoldenRouterWitness w) =>
        ((IInspectable)w).Children.ToList();

    [Fact]
    public void LiveGoldenRatio_c1_IsPhi_To1e12()
    {
        // The live ratio comes from where the router's window-summed anticommutator vanishes, NOT the
        // closed form; its agreement with φ is the genuine check.
        double live = KBodyPalindromeRouting.LiveMetallicRatio(1.0);
        Assert.Equal(Phi, live, 12);
        Assert.Equal(Phi, KBodyPalindromeRouting.MetallicMean(1.0), 12);
    }

    [Fact]
    public void LiveSilverRatio_c2_IsOnePlusRoot2_To1e12()
    {
        double live = KBodyPalindromeRouting.LiveMetallicRatio(2.0);
        Assert.Equal(Silver, live, 12);
        Assert.Equal(Silver, KBodyPalindromeRouting.MetallicMean(2.0), 12);
    }

    [Theory]
    [InlineData(0.0, 1.0)]                      // r(0) = 1, the 45° frame
    [InlineData(1.0, 1.6180339887498949)]      // golden φ
    [InlineData(2.0, 2.4142135623730951)]      // silver 1+√2
    [InlineData(3.0, 3.3027756377319946)]      // bronze (3+√13)/2
    public void LiveMetallicRatio_MatchesClosedForm_OnEveryStation(double c, double expected)
    {
        double live = KBodyPalindromeRouting.LiveMetallicRatio(c);
        Assert.Equal(expected, live, 9);
        Assert.Equal(expected, KBodyPalindromeRouting.MetallicMean(c), 12);
    }

    [Fact]
    public void WindowSummedResidual_VanishesAtTheMetallicMean_AndIsNonzeroOffIt()
    {
        // At r = r(c) the residual is essentially zero; a frame ratio off the locus does not palindromize.
        double r = KBodyPalindromeRouting.MetallicMean(1.0);
        Assert.True(KBodyPalindromeRouting.WindowSummedAnticommutatorNorm(1.0, r) < 1e-9);
        Assert.True(KBodyPalindromeRouting.WindowSummedAnticommutatorNorm(1.0, r + 0.3) > 1e-2);
    }

    [Fact]
    public void Ceiling_ReachesZero_OnTheCanonicalZMiddleCases_WithTheRouter()
    {
        // The 2 → 0 step live: the per-term lens declines, the full certifier certifies via the
        // window-summed golden router. The ceiling closes.
        var w = new GoldenRouterWitness();
        var close = w.CeilingClose();
        Assert.Equal(2, close.Count);
        foreach (var (_, perTermRoutes, certified, strategy) in close)
        {
            Assert.False(perTermRoutes);                 // per-term lens cannot see the cross-template cancellation
            Assert.True(certified);                      // the certifier closes it
            Assert.Equal(PalindromeSoftCertifier.SoftStrategy.RoutingWindowSummed, strategy);
        }
    }

    [Fact]
    public void Readings_AllMatch_ForTheDefaultWeights()
    {
        var w = new GoldenRouterWitness();
        var readings = w.Readings();
        Assert.Equal(4, readings.Count);                 // default {0, 1, 2, 3}
        Assert.All(readings, r => Assert.True(r.Matches, $"c={r.C}: live {r.LiveRatio} vs closed {r.ClosedForm}"));
    }

    [Fact]
    public void Witness_SurfacesArcAndMetallicChildren()
    {
        var labels = Children(new GoldenRouterWitness()).Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("ceiling arc 6 → 4 → 2 → 0"));
        Assert.Contains(labels, l => l.Contains("c = 1"));
        Assert.Contains(labels, l => l.Contains("c = 2"));
        Assert.Contains(labels, l => l.Contains("residual"));
        // The arc node carries the two live 2 → 0 rows plus the 6 → 4 and 4 → 2 steps.
        var arc = Children(new GoldenRouterWitness()).Single(c => c.DisplayName.Contains("ceiling arc"));
        Assert.Equal(4, arc.Children.Count());
    }

    [Fact]
    public void Witness_ChildCount_IsArcPlusOnePerWeightPlusResidual()
    {
        var weights = new[] { 1.0, 2.0 };
        var children = Children(new GoldenRouterWitness(weights));
        // 1 arc node + 2 per-weight nodes + 1 residual curve.
        Assert.Equal(1 + weights.Length + 1, children.Count);
    }

    [Fact]
    public void Witness_Summary_ReportsCeilingClosedAndMatchCount()
    {
        var summary = new GoldenRouterWitness().Summary;
        Assert.Contains("closes the ceiling", summary);
        Assert.Contains("RoutingWindowSummed", summary);
        Assert.Contains("4/4 MATCH", summary);
    }

    [Fact]
    public void Witness_RendersToJson()
    {
        var json = InspectionJsonExporter.ToJson(new GoldenRouterWitness());
        Assert.Contains("metallic", json);
        Assert.Contains("MATCH", json);
    }
}
