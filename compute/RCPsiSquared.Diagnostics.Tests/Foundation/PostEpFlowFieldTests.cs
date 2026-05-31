using System;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using System.Numerics;
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

    [Fact]
    public void LongTime_RelaxesToOneOverN()
    {
        // Slowest non-kernel rate ~ O(1) in τ; τ up to 50 is fully converged.
        var field = new PostEpFlowField(4, new[] { 2.5 }, Linspace(0, 50, 60));
        var q = field.Flows.Single();
        foreach (var s in q.Sites)
            Assert.Equal(0.25, s.Occupation[^1], 6);
    }

    [Fact]
    public void Kernel_FutureIsAlreadyPresent_AtTimeZero()
    {
        // THE_FLOW: the λ=0 component of ρ(0) already equals the uniform 1/N target, bit-exact.
        const int n = 4;
        var field = new PostEpFlowField(n, new[] { 2.5 }, Linspace(0, 6, 10));
        var L = field.DimensionlessLiouvillian(2.5);
        var rho0 = field.InitialStateVec();

        var evd = L.Evd();
        var R = evd.EigenVectors;
        var lambda = evd.EigenValues;
        var c0 = R.Solve(rho0);

        // Keep only the kernel modes (|λ| < 1e-7), reconstruct their contribution to vec(ρ).
        var masked = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>.Build.Dense(c0.Count);
        for (int i = 0; i < c0.Count; i++)
            if (lambda[i].Magnitude < 1e-7) masked[i] = c0[i];
        var kernelProjection = R * masked;

        // Target: uniform single-excitation ρ = (1/N) Σ_{popcount-1 |s⟩⟨s|}, row-major vec.
        int d = 1 << n;
        var target = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>.Build.Dense(d * d);
        for (int a = 0; a < d; a++)
            if (System.Numerics.BitOperations.PopCount((uint)a) == 1)
                target[a * d + a] = new System.Numerics.Complex(1.0 / n, 0.0);

        double diff = (kernelProjection - target).L2Norm();
        Assert.True(diff < 1e-9, $"kernel projection differs from 1/N target by {diff:E2}");
    }

    [Fact]
    public void Tree_HasQNodesWithSiteCurveLeaves()
    {
        var field = new PostEpFlowField(4, new[] { 0.5, 2.5 }, Linspace(0, 6, 20));
        var qNodes = field.Children.ToList();
        Assert.Equal(2, qNodes.Count);

        var firstQ = qNodes[0];
        var siteLeaves = firstQ.Children.ToList();
        Assert.Equal(4, siteLeaves.Count);
        Assert.IsType<RCPsiSquared.Core.Inspection.InspectablePayload.Curve>(siteLeaves[0].Payload);

        var curve = (RCPsiSquared.Core.Inspection.InspectablePayload.Curve)siteLeaves[0].Payload;
        Assert.Equal(20, curve.X.Count);
        Assert.Equal(20, curve.Y.Count);
    }

    [Fact]
    public void MatchesPythonPrototype_AtN4_Q2p5()
    {
        var taus = Linspace(0, 6, 60);
        var field = new PostEpFlowField(4, new[] { 2.5 }, taus);
        var q = field.Flows.Single();
        // Anchor recorded from `python simulations/post_ep_dynamics_4d.py 4` (Q=2.50 row, site 0 end).
        const double pythonSite0End = 0.25;
        Assert.Equal(pythonSite0End, q.Sites[0].Occupation[^1], 2);
    }

    [Fact]
    public void NormalizeToTotal_SumsToTargetAndKeepsRatios()
    {
        var norm = PostEpFlowField.NormalizeToTotal(new[] { 0.2, 0.6, 2.4, 0.6, 0.2 }, 5.0);
        Assert.Equal(5.0, norm.Sum(), 9);
        Assert.Equal(12.0, norm[2] / norm[0], 9);
    }

    [Fact]
    public void Constructor_RejectsWrongLengthProfile()
    {
        Assert.Throws<ArgumentException>(() =>
            new PostEpFlowField(4, new[] { 1.0 }, Linspace(0, 6, 10), gammaProfile: new[] { 1.0, 1.0, 1.0 }));
    }

    [Fact]
    public void Constructor_RejectsNonPositiveProfileEntry()
    {
        Assert.Throws<ArgumentException>(() =>
            new PostEpFlowField(4, new[] { 1.0 }, Linspace(0, 6, 10), gammaProfile: new[] { 1.0, 0.0, 1.0, 1.0 }));
    }
}
