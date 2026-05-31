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
}
