using System;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class EnvelopeTheoremWitnessTests
{
    [Fact]
    public void Global_EnvelopeNonIncreasing_TheoremHoldsLive()
    {
        var w = new EnvelopeTheoremWitness(n: 3);
        Assert.Equal(0, w.GlobalBell.RiseCount);
        Assert.True(w.GlobalBell.IsNonIncreasing);
    }

    [Fact]
    public void LocalBell_EnvelopeRises_AboveTheGenuinenessBar_TheFreedom()
    {
        var w = new EnvelopeTheoremWitness(n: 3);
        Assert.True(w.LocalBell.RiseCount >= 1);
        Assert.True(w.LocalBell.MaxRiseMagnitude > EnvelopeTheoremWitness.GenuinenessBar);
    }

    [Fact]
    public void SingleExcitation_IsSubBarAndVanishesUnderRefinement_TheArtifactControl()
    {
        var w = new EnvelopeTheoremWitness(n: 3);
        Assert.True(w.SingleCoarse.MaxRiseMagnitude < EnvelopeTheoremWitness.GenuinenessBar);
        Assert.True(w.SingleFine.MaxRiseMagnitude < EnvelopeTheoremWitness.GenuinenessBar);
        Assert.Equal(0, w.SingleFine.RiseCount);
    }

    [Fact]
    public void BondingMode_IsSilent_TheHEigenstateControl()
    {
        var w = new EnvelopeTheoremWitness(n: 3);
        Assert.Equal(0, w.LocalBonding.RiseCount);
    }

    [Fact]
    public void N2_Rejected_WithTheIdentityCaseMessage()
    {
        var ex = Assert.Throws<ArgumentOutOfRangeException>(() => new EnvelopeTheoremWitness(n: 2));
        Assert.Contains("identity", ex.Message, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void Witness_SummaryNamesTheClaimAndTheBoundary_ChildrenSurfaceTheThreeStories()
    {
        var w = new EnvelopeTheoremWitness(n: 3);
        Assert.Contains("CpsiEnvelopeTheoremClaim", w.Summary);
        Assert.Contains("freedom", w.Summary);
        var labels = ((IInspectable)w).Children.Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("theorem"));
        Assert.Contains(labels, l => l.Contains("freedom"));
        Assert.Contains(labels, l => l.Contains("control"));
    }
}
