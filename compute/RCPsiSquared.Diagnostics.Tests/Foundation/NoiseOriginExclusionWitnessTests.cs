using System;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The light live witness for INCOMPLETENESS_PROOF's noise-origin 5-candidate elimination.
/// The recomputed core is Candidate 5, the dimension algebra: d²−2d = 0 has roots exactly {0, 2}
/// (two-sided — 0 and 2 ARE roots AND d=1, d≥3 are NOT), so nothing other than a qubit or nothing
/// can carry the mirror condition; the other candidates are typed/surfaced, not recomputed.</summary>
public class NoiseOriginExclusionWitnessTests
{
    [Fact]
    public void Constructor_RejectsBadArgs()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new NoiseOriginExclusionWitness(maxD: 1));
    }

    [Fact]
    public void MinMemoryPolynomial_RootsAreZeroAndTwo_OthersExcluded()   // two-sided
    {
        Assert.Equal(0, NoiseOriginExclusionWitness.MinMemoryPolynomial(0));   // nothing
        Assert.Equal(0, NoiseOriginExclusionWitness.MinMemoryPolynomial(2));   // qubit
        Assert.Equal(-1, NoiseOriginExclusionWitness.MinMemoryPolynomial(1));  // d=1 excluded
        Assert.Equal(3, NoiseOriginExclusionWitness.MinMemoryPolynomial(3));   // d=3 excluded
        Assert.Equal(8, NoiseOriginExclusionWitness.MinMemoryPolynomial(4));   // d=4 excluded
    }

    [Fact]
    public void AllowedDimensions_AreExactlyZeroAndTwo()
    {
        var w = new NoiseOriginExclusionWitness(maxD: 5);
        Assert.Equal(new[] { 0, 2 }, w.AllowedDimensions);
    }

    [Fact]
    public void FiveCandidates_Candidate5IsLive_Candidate2IsDeferred()
    {
        var w = new NoiseOriginExclusionWitness();
        Assert.Equal(5, w.Candidates.Count);
        Assert.True(w.Candidates.Single(c => c.Index == 5).Recomputed,
            "Candidate 5 (the dimension algebra) is recomputed live");
        Assert.False(w.Candidates.Single(c => c.Index == 2).Recomputed,
            "Candidate 2 (single-qubit decay) is the deferred heavy compute in this light witness");
    }

    [Fact]
    public void Witness_SurfacesDimensionAlgebraAndConclusion()
    {
        var labels = ((IInspectable)new NoiseOriginExclusionWitness()).Children.Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("dimension algebra") && l.Contains("Candidate 5"));
        Assert.Contains(labels, l => l.Contains("conclusion") && l.Contains("OUTSIDE"));
    }
}
