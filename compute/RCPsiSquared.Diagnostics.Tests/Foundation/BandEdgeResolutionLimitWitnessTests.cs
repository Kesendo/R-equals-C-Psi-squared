using System;
using System.Linq;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Gate-first (the C# twin of simulations/f124_inverse_problem_gate.py Stage 0): F124's bond-to-mode
/// map M, read as a bond-recovery inverse problem, has a defect-localization resolution limit. σ_min²=E is the
/// reconstruction floor; κ=λ_max/λ_min ~ N²; the contrast σ_max/σ_min=√κ ~ N; the worst-conditioned bond
/// direction is the staggered q=π mode; the floor E·(N+1)³ → 4π².</summary>
public class BandEdgeResolutionLimitWitnessTests
{
    private static readonly int[] Ladder = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
    private static readonly int[] Range = { 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    /// <summary>The worst-conditioned bond direction (the least-observable input) is the staggered q=π mode, the
    /// diffraction-limit detail, at every N.</summary>
    [Fact]
    public void WorstConditionedDirection_IsStaggered_AllN()
    {
        foreach (int n in Ladder)
            Assert.True(BandEdgeResolutionLimitWitness.WorstDirectionIsStaggered(n),
                $"N={n}: the worst-conditioned bond direction should be the staggered q=π mode");
    }

    /// <summary>The worst-case reconstruction floor σ_min²(M) equals F124's lower frame bound E.</summary>
    [Fact]
    public void ReconstructionFloor_SigmaMinSquared_EqualsEndpointE_AllN()
    {
        foreach (int n in Ladder)
        {
            double sminSq = Math.Pow(BandEdgeResolutionLimitWitness.SigmaMin(n), 2);
            Assert.Equal(BandEdgeTransitionInvariantWitness.EndpointClosedForm(n), sminSq, 9);
        }
    }

    /// <summary>The condition number κ = λ_max/λ_min grows ~ N² (the noise amplification of the inverse).</summary>
    [Fact]
    public void ConditionNumber_GrowsAsNSquared()
    {
        double slope = BandEdgeResolutionLimitWitness.LogLogSlope(Range, BandEdgeResolutionLimitWitness.Kappa);
        Assert.InRange(slope, 1.8, 2.2);
    }

    /// <summary>The contrast / resolution ratio σ_max/σ_min equals √κ exactly and grows ~ N (a staggered defect
    /// is √κ ~ N times harder to localize than a band-edge one).</summary>
    [Fact]
    public void ContrastRatio_IsSqrtKappa_AndGrowsAsN()
    {
        foreach (int n in Range)
            Assert.Equal(Math.Sqrt(BandEdgeResolutionLimitWitness.Kappa(n)),
                BandEdgeResolutionLimitWitness.ContrastRatio(n), 9);
        double slope = BandEdgeResolutionLimitWitness.LogLogSlope(Range, BandEdgeResolutionLimitWitness.ContrastRatio);
        Assert.InRange(slope, 0.8, 1.2);
    }

    /// <summary>The reconstruction floor vanishes as (N+1)^(−3): E·(N+1)³ → 4π² = 39.478 (the right-variable
    /// law; the naive small-N exponent in N is a pre-asymptotic artifact).</summary>
    [Fact]
    public void ReconstructionFloor_AsymptoticConstant_ApproachesFourPiSquared()
    {
        // E·(N+1)³ approaches 4π²=39.478 from below; at N=60 it is ~39.44 (still climbing). Tolerance, not rounding.
        Assert.True(System.Math.Abs(BandEdgeResolutionLimitWitness.FloorConstant(60) - 4 * Math.PI * Math.PI) < 0.1,
            $"E·(N+1)³ at N=60 = {BandEdgeResolutionLimitWitness.FloorConstant(60):0.####}, should approach 4π² = {4 * Math.PI * Math.PI:0.####}");
        // right-variable log-log slope over N=4..60 is ~ -3 (in N+1), NOT the naive ~ -2.5 in N
        double slope = BandEdgeResolutionLimitWitness.LogLogSlope(
            Enumerable.Range(4, 57), BandEdgeTransitionInvariantWitness.EndpointClosedForm, versusNPlus1: true);
        Assert.InRange(slope, -3.1, -2.9);
    }
}
