using System;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class EnvelopeTheoremWitnessTests
{
    [Fact]
    public void Witness_AtN4_IsHonest_GlobalRefutes_ControlDoesNotFakeVanishing()
    {
        // The envelope_n4_rise finding: at N=4 the full-state envelope GENUINELY rises (refinement-stable).
        // The witness must NOT lie: the Summary must say the live evolution REFUTES the theorem at this N,
        // never "CONFIRMS"; and the artifact control must not hardcode "vanish under refinement" when
        // SingleExcitation actually persists at N=4.
        var w = new EnvelopeTheoremWitness(n: 4);
        Assert.True(w.GlobalBell.RiseCount > 0, "precondition: the N=4 global envelope rises");
        Assert.DoesNotContain("CONFIRMS", w.Summary);
        Assert.Contains("REFUTES", w.Summary);

        var control = ((IInspectable)w).Children.Single(c => c.DisplayName.Contains("control"));
        if (w.SingleFine.RiseCount > 0)
            Assert.DoesNotContain("they vanish under refinement", control.Summary);
    }

    [Fact]
    public void Witness_AtN3_StillConfirms_TheoremHoldsAndControlVanishes()
    {
        var w = new EnvelopeTheoremWitness(n: 3);
        Assert.Equal(0, w.GlobalBell.RiseCount);
        Assert.Contains("CONFIRMS", w.Summary);
        Assert.DoesNotContain("REFUTES", w.Summary);
    }

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
