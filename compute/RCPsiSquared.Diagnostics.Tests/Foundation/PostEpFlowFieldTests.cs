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

    [Fact]
    public void UniformProfile_MatchesScalarDefault()
    {
        var taus = Linspace(0, 6, 30);
        var dflt = new PostEpFlowField(4, new[] { 2.5 }, taus);
        var explicitUniform = new PostEpFlowField(4, new[] { 2.5 }, taus, gammaProfile: new[] { 1.0, 1.0, 1.0, 1.0 });
        var a = dflt.Flows.Single().Sites;
        var b = explicitUniform.Flows.Single().Sites;
        for (int s = 0; s < 4; s++)
            for (int t = 0; t < taus.Length; t++)
                Assert.Equal(a[s].Occupation[t], b[s].Occupation[t], 12);
    }

    [Fact]
    public void NonUniformProfile_StillRelaxesToOneOverN()
    {
        var profile = PostEpFlowField.NormalizeToTotal(new[] { 0.2, 0.6, 2.4, 0.6, 0.2 }, 5.0);
        var field = new PostEpFlowField(5, new[] { 2.5 }, Linspace(0, 50, 60), gammaProfile: profile);
        var q = field.Flows.Single();
        foreach (var s in q.Sites)
            Assert.Equal(0.2, s.Occupation[^1], 6);
    }

    [Fact]
    public void NonUniformProfile_ConservesExcitation()
    {
        var profile = PostEpFlowField.NormalizeToTotal(new[] { 0.2, 0.6, 2.4, 0.6, 0.2 }, 5.0);
        var taus = Linspace(0, 6, 30);
        var field = new PostEpFlowField(5, new[] { 1.5 }, taus, gammaProfile: profile);
        var q = field.Flows.Single();
        for (int t = 0; t < taus.Length; t++)
            Assert.Equal(1.0, q.Sites.Sum(s => s.Occupation[t]), 9);
    }

    [Fact]
    public void Profile_ChangesSpectralGap_AtFixedTotal()
    {
        var taus = Linspace(0, 6, 20);
        var vshape = PostEpFlowField.NormalizeToTotal(new[] { 0.2, 0.6, 2.4, 0.6, 0.2 }, 5.0);
        var uniform = new PostEpFlowField(5, new[] { 20.0 }, taus); // uniform, Sum gamma = 5
        var shaped = new PostEpFlowField(5, new[] { 20.0 }, taus, gammaProfile: vshape); // Sum gamma = 5
        double rateU = uniform.Flows.Single().SlowestRate;
        double rateV = shaped.Flows.Single().SlowestRate;
        Assert.True(rateU > 0, $"uniform slowest rate should be positive, got {rateU}");
        Assert.True(rateV > 0, $"shaped slowest rate should be positive, got {rateV}");
        Assert.True(Math.Abs(rateU - rateV) > 1e-6, $"shape should change the gap: uniform={rateU:F6}, vshape={rateV:F6}");
    }

    [Fact]
    public void Tree_QNodeSummary_IncludesSlowestRate()
    {
        var field = new PostEpFlowField(4, new[] { 2.5 }, Linspace(0, 6, 20));
        var qNode = field.Children.First();
        Assert.Contains("rate", qNode.Summary);
    }

    [Fact]
    public void Uniform_IsInSterileZone_ClosedFormRateIsTwo()
    {
        var field = new PostEpFlowField(5, new[] { 20.0 }, Linspace(0, 6, 10));
        Assert.True(field.IsInSterileZone);
        Assert.False(field.IsInBirthCanal);
        Assert.Equal(0.0, field.BirthCanalDeviation, 4);
        Assert.Equal(2.0, field.ClosedFormRate, 4);
    }

    [Fact]
    public void PeakedV_IsSterile_ClosedFormRateIsOne()
    {
        var field = new PostEpFlowField(5, new[] { 20.0 }, Linspace(0, 6, 10),
            gammaProfile: new[] { 0.25, 0.75, 3.0, 0.75, 0.25 });
        Assert.True(field.IsInSterileZone);
        Assert.Equal(1.0, field.ClosedFormRate, 4);
    }

    [Fact]
    public void FlatBulkEdge_IsInBirthCanal_ClosedFormRateThrows()
    {
        var field = new PostEpFlowField(5, new[] { 20.0 }, Linspace(0, 6, 10),
            gammaProfile: new[] { 0.25, 1.5, 1.5, 1.5, 0.25 });
        Assert.True(field.IsInBirthCanal);
        Assert.False(field.IsInSterileZone);
        Assert.True(field.BirthCanalDeviation > 0.05, $"expected positive Q-drift, got {field.BirthCanalDeviation}");
        Assert.Throws<InvalidOperationException>(() => field.ClosedFormRate);
    }

    [Fact]
    public void MaxSaturationCeiling_IsQuarter()
    {
        // C_block ≤ 1/4 for any state (BlockCoherenceContent, Theorem 2, the bilinear apex ¼).
        Assert.Equal(0.25, PostEpFlowField.MaxSaturationCeiling, 12);
    }

    [Fact]
    public void SingleExcitationFlow_StaysOutOfBirthChannel()
    {
        // The post-EP flow is a single excitation (definite number, number-conserving, even rung):
        // it carries NO between-block coherence, so it never enters the birth channel. Measured live.
        var field = new PostEpFlowField(4, new[] { 2.5 }, Linspace(0, 6, 30));
        Assert.True(field.PeakBetweenBlockSaturation < 1e-9,
            $"single-excitation flow should carry no between-block coherence, got {field.PeakBetweenBlockSaturation:E2}");
        Assert.False(field.FlowEntersBirthChannel);
    }

    [Fact]
    public void BirthChannel_IsSaturableToQuarter_ByDickeSuperposition()
    {
        // The odd-rung channel the flow avoids is real and reaches the ceiling: the Dicke
        // superposition (|D_n⟩+|D_{n+1}⟩)/√2 saturates C_block = MaxSaturationCeiling = 1/4.
        const int N = 4, n = 1;
        var v = (DickeState(N, n) + DickeState(N, n + 1)).Divide(new System.Numerics.Complex(Math.Sqrt(2.0), 0));
        var rho = v.OuterProduct(v.Conjugate());
        double cBlock = RCPsiSquared.Core.F86.BlockCoherenceContent.Compute(rho, n);
        Assert.Equal(PostEpFlowField.MaxSaturationCeiling, cBlock, 9);
    }

    [Fact]
    public void RingTopology_ChangesTheFlowDynamics()
    {
        // The wrap bond closes the loop: the excitation can go both ways around, so the ring
        // flow differs from the open chain (benzene C4 ring vs butadiene C4 chain).
        var taus = Linspace(0, 3, 30);
        var chain = new PostEpFlowField(4, new[] { 2.5 }, taus);
        var ring = new PostEpFlowField(4, new[] { 2.5 }, taus, topology: FlowTopology.Ring);
        var cOcc = chain.Flows.Single().Sites[2].Occupation;
        var rOcc = ring.Flows.Single().Sites[2].Occupation;
        bool differ = false;
        for (int t = 0; t < taus.Length; t++)
            if (Math.Abs(cOcc[t] - rOcc[t]) > 1e-3) differ = true;
        Assert.True(differ, "ring topology should change the flow dynamics vs the open chain");
    }

    [Fact]
    public void ReadAssembly_PutsThePiecesTogether_OnTheBirthRail()
    {
        // The slowest non-kernel mode (the birth-canal object) read through the whole assembly:
        // depth = light = rate (Absorption), the parity rung, the saturation ceiling, in one place.
        var field = new PostEpFlowField(5, new[] { 20.0 }, Linspace(0, 6, 10));
        var a = field.ReadAssembly(20.0);
        Assert.Equal(1.0, a.SlowestDepth, 4);          // one light quantum
        Assert.Equal(1, a.Rung);                        // odd rung
        Assert.True(a.OnBirthRail);
        Assert.Equal(a.SlowestRate, a.AbsorptionRate, 6);   // rate = 2·Σ γ_k·light_k (Absorption)
        Assert.Equal(2.0, a.SlowestRate, 4);            // = 2γ uniform
        Assert.Equal(0.25, a.MaxSaturationCeiling, 12);
        // the slowest rate is degenerate; the carrier is the gauge-invariant average over that
        // subspace. For the symmetric uniform chain it must be left-right symmetric, sum to the
        // depth (one quantum), and be roughly equipartitioned (no gauge-dependent single-mode spike).
        Assert.True(a.Degeneracy > 1, $"uniform slowest rate should be degenerate, got {a.Degeneracy}");
        Assert.Equal(a.PerSiteLight[0], a.PerSiteLight[4], 3);
        Assert.Equal(a.PerSiteLight[1], a.PerSiteLight[3], 3);
        Assert.Equal(1.0, a.PerSiteLight.Sum(), 4);
        foreach (double lk in a.PerSiteLight)
            Assert.InRange(lk, 0.10, 0.30);
    }

    [Fact]
    public void ReadAssembly_AbsorptionCrossCheck_HoldsForNonUniformProfile()
    {
        // The per-site light (carrier vector) reproduces the rate for a shaped profile too:
        // peaked-V puts the one quantum on the low-γ edges, so the rate drops to 1γ (still odd rung).
        var field = new PostEpFlowField(5, new[] { 20.0 }, Linspace(0, 6, 10),
            gammaProfile: new[] { 0.25, 0.75, 3.0, 0.75, 0.25 });
        var a = field.ReadAssembly(20.0);
        Assert.Equal(a.SlowestRate, a.AbsorptionRate, 6);
        Assert.Equal(1.0, a.SlowestDepth, 4);           // still one quantum on the odd rail
        Assert.Equal(1, a.Rung);
        Assert.Equal(1.0, a.SlowestRate, 4);            // γ-weighted share, not 2γ
    }

    private static MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex> DickeState(int N, int n)
    {
        int d = 1 << N;
        var v = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>.Build.Dense(d);
        for (int x = 0; x < d; x++)
            if (System.Numerics.BitOperations.PopCount((uint)x) == n) v[x] = System.Numerics.Complex.One;
        return v.Divide(new System.Numerics.Complex(v.L2Norm(), 0));
    }
}
