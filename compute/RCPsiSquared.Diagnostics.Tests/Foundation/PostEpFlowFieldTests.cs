using System;
using System.Linq;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class PostEpFlowFieldTests
{
    private static double[] Linspace(double a, double b, int n)
    {
        var g = new double[n];
        for (int i = 0; i < n; i++) g[i] = a + (b - a) * i / (n - 1);
        return g;
    }

    [Theory]
    [InlineData(3, 0.333333333333)]
    [InlineData(4, 0.25)]
    [InlineData(5, 0.2)]
    [InlineData(6, 0.166666666667)]
    public void Target_IsOneOverN(int n, double expected)
    {
        var field = new PostEpFlowField(n, new[] { 1.0 }, Linspace(0, 6, 10));
        Assert.Equal(expected, field.Target, 9);
    }

    [Fact]
    public void SingleExcitation_IsConserved_AtEveryTau()
    {
        var taus = Linspace(0, 6, 30);
        var field = new PostEpFlowField(4, new[] { 1.0 }, taus);
        var q = field.Flows.Single();
        for (int t = 0; t < taus.Length; t++)
        {
            double sum = q.Sites.Sum(s => s.Occupation[t]);
            Assert.Equal(1.0, sum, 9);
        }
    }

    [Fact]
    public void InitialState_HasExcitationOnSiteZero()
    {
        var field = new PostEpFlowField(4, new[] { 2.5 }, Linspace(0, 6, 30));
        var q = field.Flows.Single();
        Assert.Equal(1.0, q.Sites[0].Occupation[0], 9);
        for (int s = 1; s < 4; s++)
            Assert.Equal(0.0, q.Sites[s].Occupation[0], 9);
    }

    [Fact]
    public void Transport_ExcitationHopsOffSiteZero()
    {
        // Above the EP (Q=2.5) the excitation sloshes: site 0 must lose population and a
        // neighbour must gain it. This pins the dynamics, not just the conserved sector.
        var field = new PostEpFlowField(4, new[] { 2.5 }, Linspace(0, 6, 60));
        var q = field.Flows.Single();
        double site0End = q.Sites[0].Occupation[^1];
        double site1Max = Enumerable.Range(0, q.Sites[1].Occupation.Count)
            .Select(i => q.Sites[1].Occupation[i]).Max();
        Assert.True(site0End < 0.9, $"site 0 should lose population, end={site0End:F3}");
        Assert.True(site1Max > 0.1, $"site 1 should receive population, max={site1Max:F3}");
    }
}
