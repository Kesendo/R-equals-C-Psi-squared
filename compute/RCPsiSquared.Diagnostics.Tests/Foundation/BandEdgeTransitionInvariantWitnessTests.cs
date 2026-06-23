using System;
using System.Linq;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Topo = RCPsiSquared.Diagnostics.Foundation.BandEdgeTransitionInvariantWitness.Topo;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Gate-first (the C# twin of simulations/handshake_M_checksum.py and
/// simulations/handshake_F124_adversarial.py): for the open chain's band-edge carrier, the full
/// bond-transition matrix M[b,k]=⟨ψ_k|V_b|ψ_1⟩ (all N modes) has ‖M‖_F² + λ_min(MMᵀ) = z = 2 exactly,
/// with ‖M‖_F²=2−E and λ_min=E. The non-trivial half λ_min=E is the Dirichlet-edge coupling; the genuine
/// minimum, the frame identities, and the carrier/object/topology breakages are all gated here.</summary>
public class BandEdgeTransitionInvariantWitnessTests
{
    private static readonly int[] ChainSweep = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    /// <summary>The headline identity: ‖M‖_F² + λ_min(MMᵀ) = z = 2 exactly, for every N.</summary>
    [Fact]
    public void Identity_SumIsTwo_AllN()
    {
        foreach (int n in ChainSweep)
            Assert.Equal(2.0, BandEdgeTransitionInvariantWitness.Analyse(Topo.Chain, n).Sum, 9);
    }

    /// <summary>The real content: λ_min = E = (4/(N+1))sin²(π/(N+1)) (the Dirichlet-edge coupling).</summary>
    [Fact]
    public void LambdaMin_EqualsEndpointClosedForm_AllN()
    {
        foreach (int n in ChainSweep)
        {
            var r = BandEdgeTransitionInvariantWitness.Analyse(Topo.Chain, n);
            Assert.Equal(BandEdgeTransitionInvariantWitness.EndpointClosedForm(n), r.LamMin, 9);
        }
    }

    /// <summary>The trace half: ‖M‖_F² = 2 − E (the carrier's degree-weighted-norm deficit from z).</summary>
    [Fact]
    public void Frobenius_EqualsTwoMinusEndpoint_AllN()
    {
        foreach (int n in ChainSweep)
        {
            var r = BandEdgeTransitionInvariantWitness.Analyse(Topo.Chain, n);
            Assert.Equal(2.0 - BandEdgeTransitionInvariantWitness.EndpointClosedForm(n), r.Fro2, 9);
        }
    }

    /// <summary>The staggered bond wave (−1)^b is the GENUINE minimum at the band edge (an eigenvector
    /// whose Rayleigh quotient IS λ_min), and the band-edge Gram off-diagonals are all positive (Perron).</summary>
    [Fact]
    public void Staggered_IsGenuineMinimum_AtBandEdge_AllN()
    {
        foreach (int n in ChainSweep)
        {
            var r = BandEdgeTransitionInvariantWitness.Analyse(Topo.Chain, n);
            Assert.True(r.StaggeredIsGenuineMinimum, $"N={n}: staggered should be the genuine minimum");
            Assert.True(r.OffDiagAllPositive, $"N={n}: band-edge Gram off-diagonals should all be positive");
        }
    }

    /// <summary>The conserved discrete-energy envelope Q_a = c₀² along the whole chain (Part 2's
    /// bulk-cancellation heart, the reason λ_min = 2Q₀ = 2c₀² = E).</summary>
    [Fact]
    public void Envelope_IsConstantAndEqualsC0Squared_AllN()
    {
        foreach (int n in ChainSweep)
        {
            var (constant, eqC0) = BandEdgeTransitionInvariantWitness.EnvelopeCheck(n);
            Assert.True(constant, $"N={n}: envelope Q should be constant along the chain");
            Assert.True(eqC0, $"N={n}: envelope Q should equal c₀²");
        }
    }

    /// <summary>The carrier-selecting step (the genuine-minimum hazard): for an INTERIOR carrier the staggered
    /// mode is still an eigenvector but NOT the least, and the sum drops strictly below 2. The band edge is
    /// load-bearing for the "=2".</summary>
    [Theory]
    [InlineData(7)]
    [InlineData(9)]
    [InlineData(11)]
    public void InteriorCarrier_StaggeredEigenvectorButNotMinimum_SumBelowTwo(int n)
    {
        var r = BandEdgeTransitionInvariantWitness.Analyse(Topo.Chain, n, carrierRank: 1);
        Assert.True(r.StagResid < 1e-8, "staggered should remain an eigenvector for any carrier");
        Assert.True(r.StagRayleigh > r.LamMin + 1e-8, "staggered should NOT be the minimum for an interior carrier");
        Assert.True(r.Sum < 2.0 - 1e-8, "the sum should drop below 2 for an interior carrier");
        // tie the refutation to the Perron mechanism: a noded interior carrier loses the all-positive
        // Gram off-diagonals that make the staggered mode the minimum for the band edge.
        Assert.False(r.OffDiagAllPositive, "an interior (noded) carrier should lose the Perron positivity precondition");
    }

    /// <summary>The frame reading: λ_min = σ_min²(M) (the lower frame bound = the Eckart-Young distance² to
    /// rank-collapse) and the K-partner column ⟨ψ_N|V_b|ψ_1⟩ ≡ 0 (the exact kernel).</summary>
    [Fact]
    public void FrameReading_LambdaMinEqualsSigmaMinSquared_AndKPartnerColumnNull_AllN()
    {
        foreach (int n in ChainSweep)
        {
            var r = BandEdgeTransitionInvariantWitness.Analyse(Topo.Chain, n);
            Assert.Equal(r.LamMin, r.SigmaMinSq, 9);
            Assert.True(r.KPartnerColNorm < 1e-9, $"N={n}: the K-partner column should be null");
        }
    }

    /// <summary>The conditioner: κ = λ_max/λ_min grows monotonically with N (shorter chains better-conditioned;
    /// the long-chain limit E→0 goes singular).</summary>
    [Fact]
    public void ConditionNumber_GrowsWithN()
    {
        var kappas = new[] { 4, 5, 6, 7, 8 }
            .Select(n => BandEdgeTransitionInvariantWitness.Analyse(Topo.Chain, n).Kappa).ToArray();
        for (int i = 1; i < kappas.Length; i++)
            Assert.True(kappas[i] > kappas[i - 1], $"κ should grow: κ[{i}]={kappas[i]} ≤ κ[{i - 1}]={kappas[i - 1]}");
    }

    /// <summary>The object guard: the decoder's location dictionary k=2..N (strength channel dropped) has
    /// λ_min = 0: the K-partner null column makes M_loc M_locᵀ rank-deficient. The strength column lifts the
    /// floor from 0 to E.</summary>
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void LocationDictionary_LambdaMinIsZero(int n)
    {
        Assert.True(BandEdgeTransitionInvariantWitness.LocationDictionaryLamMin(n) < 1e-9,
            $"N={n}: the location dictionary k=2..N should have λ_min = 0");
    }

    /// <summary>Topology scope: the even ring holds DEGENERATELY (E=0, no boundary): sum=2 like the chain but
    /// with λ_min=0, so it is not a second instance of the Dirichlet-edge mechanism. The proof discloses this in
    /// its Scope section; this gate keeps the topology table honest.</summary>
    [Theory]
    [InlineData(4)]
    [InlineData(6)]
    public void EvenRing_HoldsDegenerately_SumIsTwoAndLambdaMinZero(int n)
    {
        var r = BandEdgeTransitionInvariantWitness.Analyse(Topo.Ring, n);
        Assert.Equal(2.0, r.Sum, 9);
        Assert.True(r.LamMin < 1e-9, $"N={n}: the even ring should have λ_min = 0 (no boundary, E=0)");
    }

    /// <summary>Topology scope: the odd ring frustrates the staggering (λ_min > 0, sum > 2): the "2" is the
    /// coordination number, not a universal constant.</summary>
    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void OddRing_FrustratesStaggering_SumAboveTwo(int n)
    {
        Assert.True(BandEdgeTransitionInvariantWitness.Analyse(Topo.Ring, n).Sum > 2.0 + 1e-6,
            $"N={n}: the odd ring should frustrate the staggering, sum > 2");
    }

    /// <summary>Topology scope: the star breaks the trace half (‖M‖_F² = N/2): the hub's high degree changes
    /// the degree-weighted norm.</summary>
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    public void Star_BreaksTraceHalf_FrobeniusIsNOverTwo(int n)
    {
        Assert.Equal(n / 2.0, BandEdgeTransitionInvariantWitness.Analyse(Topo.Star, n).Fro2, 9);
    }
}
