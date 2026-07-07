using System;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Guard for the live eigenvector-holonomy witness: the (1,2)-block defective seed's eigenvector
/// frame, encircling the EP, is the mod-4 memory loop i⁴=1 (M₁ eigenvalues ±i, M₂=−I, M₄=I), recomputed
/// at inspect time. Fast (N=5, ~1 s) — run by class-name filter, never bare (Diagnostics.Tests carries a
/// 44 GB SLOW_F104 test). Filter: <c>--filter "FullyQualifiedName~SeedHolonomyWitnessTests"</c>.</summary>
public class SeedHolonomyWitnessTests
{
    [Fact]
    public void N5_LiveHolonomy_IsTheMod4MemoryLoop_i4Equals1()
    {
        var w = new SeedHolonomyWitness(n: 5);
        Assert.Equal(50, w.Dim);

        // M₁ eigenvalues ±i: trace ~0, det ~1
        var ev = w.M1Eigenvalues;
        Assert.True((ev[0] + ev[1]).Magnitude < 1e-2, $"tr M₁ = {ev[0] + ev[1]} not ~0");
        Assert.True((ev[0] * ev[1] - Complex.One).Magnitude < 1e-2, $"det M₁ = {ev[0] * ev[1]} not ~1");

        // M₂ = −I, M₄ = I — the frame is single-valued only after four loops
        Assert.True(w.M2DistanceFromNegI < 1e-2, $"‖M₂+I‖ = {w.M2DistanceFromNegI:e2} not ~0");
        Assert.True(w.M4DistanceFromI < 1e-2, $"‖M₄−I‖ = {w.M4DistanceFromI:e2} not ~0");

        // the 2D span is preserved on every loop at N=5 (the full block tracks cleanly here)
        Assert.True(w.MaxSpanResidual < 1e-3, $"max span residual {w.MaxSpanResidual:e2} — tracking not clean");
    }

    [Fact]
    public void Guard_RejectsUnsupportedN()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new SeedHolonomyWitness(2));   // below range
        Assert.Throws<ArgumentOutOfRangeException>(() => new SeedHolonomyWitness(11));  // above the live-build guard
        Assert.Throws<ArgumentOutOfRangeException>(() => new SeedHolonomyWitness(6));   // no reference seed
        Assert.Throws<ArgumentOutOfRangeException>(() => new SeedHolonomyWitness(7));   // unverified N, not admitted
    }

    [Fact]
    public void WitnessAndRecomputedChildrenAreLive_ReadingNodeIsNot()
    {
        var w = new SeedHolonomyWitness(n: 5);
        Assert.Equal(NodeProvenance.Live, w.Provenance);
        foreach (var c in w.Children)
        {
            if (c.DisplayName.Contains("reading"))                    // interpretive prose, not a recomputed number
                Assert.NotEqual(NodeProvenance.Live, c.Provenance);
            else                                                      // M₁/M₂/M₄/span residual = live-recomputed
                Assert.Equal(NodeProvenance.Live, c.Provenance);
        }
    }
}
