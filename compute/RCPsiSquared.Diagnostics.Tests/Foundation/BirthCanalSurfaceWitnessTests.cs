using System;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class BirthCanalSurfaceWitnessTests
{
    private static readonly BirthCanalSurfaceWitness W = new(n: 5);

    [Fact]
    public void ReadPoint_SterileAnchors_HaveZeroDeviation()
    {
        Assert.True(W.ReadPoint(1.0, 1.0).Deviation < 1e-4);    // uniform
        Assert.True(W.ReadPoint(0.25, 3.0).Deviation < 1e-4);   // peaked-V
        Assert.False(W.ReadPoint(1.0, 1.0).IsCanal);
        Assert.False(W.ReadPoint(0.25, 3.0).IsCanal);
    }

    [Fact]
    public void ReadPoint_CanalAnchor_DriftsAndMatchesThePinnedRates()
    {
        var p = W.ReadPoint(0.25, 1.5);                          // flat-bulk-edge
        Assert.True(p.IsCanal);
        Assert.True(p.Deviation > 1e-3);
        Assert.Equal(1.2483, p.Low!.SlowestRate, 3);            // rate(Q=1.5)
        Assert.Equal(4.0 / 3.0, p.High!.SlowestRate, 2);        // rate(Q=1000) -> 4/3
        Assert.True(p.DriftMax > 5e-2);                          // max per-site light drift ~7.8e-2
    }

    [Fact]
    public void ReadPoint_SterilePeakedV_LightAvoidsTheStrongCenter()
    {
        var p = W.ReadPoint(0.25, 3.0);                          // light should be ~[1/4,1/4,0,1/4,1/4]
        Assert.True(p.DriftMax < 1e-6);                          // genuine freeze (robust)
        Assert.True(p.Low!.PerSiteLight[2] < 1e-6);             // zero light on the strong center
    }

    [Fact]
    public void ReadPoint_DeviationEqualsPostEpFlowFieldMembership()
    {
        // M7: the grid value IS the codebase's BirthCanalDeviation, not a reimplementation.
        var profile = SymmetricGammaSlice.Profile(5, 0.25, 1.5);
        var field = new PostEpFlowField(5, new[] { 1.5, 1000.0 }, new[] { 0.0, 1.0 }, profile);
        Assert.Equal(field.BirthCanalDeviation, W.ReadPoint(0.25, 1.5).Deviation, 9);
    }

    [Fact]
    public void ReadPoint_AbsorptionResidualIsMachineZero_L3()
    {
        // L3 / M4: rate == 2*sum gamma_l*light_l at the low probe, both sterile and canal.
        Assert.True(Math.Abs(W.ReadPoint(1.0, 1.0).Low!.SlowestRate - W.ReadPoint(1.0, 1.0).Low!.AbsorptionRate) < 1e-6);
        Assert.True(Math.Abs(W.ReadPoint(0.25, 1.5).Low!.SlowestRate - W.ReadPoint(0.25, 1.5).Low!.AbsorptionRate) < 1e-6);
    }

    [Fact]
    public void Grid_HasANonEmptyBoundaryCurve_AtN5()
    {
        var (xs, ys) = W.BoundaryCurve();
        Assert.NotEmpty(xs);                       // at least one sterile->canal crossing in the box
        Assert.Equal(xs.Count, ys.Count);
        foreach (double x in xs) Assert.InRange(x, 0.2, 1.0);   // crossings live in the edge range
    }

    [Fact]
    public void TheSurfaceNode_RendersHeatmapAndBoundary()
    {
        var surface = System.Linq.Enumerable.First(W.Children);   // "the surface" is the first node
        Assert.Equal("the surface", surface.DisplayName);
        var kids = System.Linq.Enumerable.ToList(surface.Children);
        Assert.Contains(kids, k => k.Payload is RCPsiSquared.Core.Inspection.InspectablePayload.MatrixView);
        Assert.Contains(kids, k => k.Payload is RCPsiSquared.Core.Inspection.InspectablePayload.Curve);
    }
}
