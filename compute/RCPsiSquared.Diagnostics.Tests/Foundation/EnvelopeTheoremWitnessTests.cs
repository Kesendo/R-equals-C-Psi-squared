using System;
using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class EnvelopeTheoremWitnessTests
{
    [Fact]
    public void GlobalRiseCount_StaticDetector_MatchesTheLiveWitness_NoDrift()
    {
        // The boundary sweep reuses a static detector lifted out of the witness; it must reproduce the
        // witness's own live GlobalBell reading at the witness regime (J=5, γ=0.01, tMax=25, 1600 pts)
        // exactly — the named N=3 window resolves 0, the N=4 window resolves rises. If the paths drift, the sweep stops measuring the
        // same thing the witness shows.
        foreach (int n in new[] { 3, 4 })
        {
            var w = new EnvelopeTheoremWitness(n);
            int viaStatic = EnvelopeTheoremWitness.GlobalRiseCount(n, 5.0, 0.01, 25.0, 1600);
            Assert.Equal(w.GlobalBell.RiseCount, viaStatic);
        }
    }

    [Fact]
    public void Witness_AtN4_ReportsFiniteRiseWithoutCallingItATheoremVerdict()
    {
        // The envelope_n4_rise row: at this N=4 grid the full-state envelope resolves rises above the bar.
        // That is a finite reading, not a theorem verdict or a causal refinement claim.
        var w = new EnvelopeTheoremWitness(n: 4);
        Assert.True(w.GlobalBell.RiseCount > 0, "precondition: the N=4 global envelope rises");
        Assert.Contains("finite", w.Summary);
        Assert.DoesNotContain("CONFIRMS", w.Summary);
        Assert.DoesNotContain("REFUTES", w.Summary);
        Assert.Contains("do not draw a theorem boundary", w.Summary);

        var control = ((IInspectable)w).Children.Single(c => c.DisplayName.Contains("control"));
        Assert.DoesNotContain("vanish under refinement", control.Summary);
        Assert.DoesNotContain("silent", control.Summary, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("without beating", control.Summary, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void Witness_AtN3_ReportsNoResolvedRiseWithoutAnAbsenceClaim()
    {
        var w = new EnvelopeTheoremWitness(n: 3);
        Assert.Equal(0, w.GlobalBell.RiseCount);
        Assert.Contains("no global predecessor rise", w.Summary);
        Assert.Contains("neither proves absence", w.Summary);
        Assert.DoesNotContain("REFUTES", w.Summary);
        Assert.DoesNotContain("CONFIRMS", w.Summary);
    }

    [Fact]
    public void Global_FiniteReadingHasNoResolvedRiseInNamedN3Window()
    {
        var w = new EnvelopeTheoremWitness(n: 3);
        Assert.Equal(0, w.GlobalBell.RiseCount);
    }

    [Fact]
    public void LocalBell_EnvelopeRises_AboveTheReportingBar_OnNamedGrid()
    {
        var w = new EnvelopeTheoremWitness(n: 3);
        Assert.True(w.LocalBell.RiseCount >= 1);
        Assert.True(w.LocalBell.MaxRiseMagnitude > EnvelopeTheoremWitness.RiseReportingBar);
    }

    [Fact]
    public void SingleExcitation_ReportsOnlyItsTwoNamedFiniteGridReadings()
    {
        var w = new EnvelopeTheoremWitness(n: 3);
        Assert.True(w.SingleCoarse.MaxRiseMagnitude < EnvelopeTheoremWitness.RiseReportingBar);
        Assert.True(w.SingleFine.MaxRiseMagnitude < EnvelopeTheoremWitness.RiseReportingBar);
        Assert.Equal(0, w.SingleFine.RiseCount);
    }

    [Fact]
    public void BondingMode_InitialStateIsHEigenstate_AndNamedGridHasNoAboveBarRise()
    {
        var w = new EnvelopeTheoremWitness(n: 3);
        Assert.Equal(0, w.LocalBonding.RiseCount);

        var s = new Symphony(n: 3, j: 5.0, gamma: 0.01,
            hType: HamiltonianType.XY, initialState: InitialStateKind.BondingMode,
            tMax: 25.0, tPoints: 1600);
        var rho0 = s.BuildInitialState(s.Dim);
        var h = new ChainSystem(s.N, s.J, s.Gamma, s.HType, s.Topology).BuildHamiltonian();
        Assert.InRange((h * rho0 - rho0 * h).FrobeniusNorm(), 0.0, 1e-12);
        var envelope = QuarterEnvelope.Of(s.States.Select(s.LocalCpsi).ToArray(),
            s.TimeGrid.ToArray(), riseTol: EnvelopeTheoremWitness.RiseReportingBar);
        Assert.Equal(0, envelope.RiseCount);
    }

    [Fact]
    public void N2_Rejected_WithTheIdentityCaseMessage()
    {
        var ex = Assert.Throws<ArgumentOutOfRangeException>(() => new EnvelopeTheoremWitness(n: 2));
        Assert.Contains("identity", ex.Message, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void Witness_SummaryNamesTheHistoricalClaimAndChildrenSurfaceTheThreeReadings()
    {
        var w = new EnvelopeTheoremWitness(n: 3);
        Assert.Contains("CpsiEnvelopeTheoremClaim", w.Summary);
        Assert.Contains("named grid reports", w.Summary);
        var labels = ((IInspectable)w).Children.Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("finite global"));
        Assert.Contains(labels, l => l.Contains("local carrier-pair"));
        Assert.Contains(labels, l => l.Contains("control"));
    }
}
