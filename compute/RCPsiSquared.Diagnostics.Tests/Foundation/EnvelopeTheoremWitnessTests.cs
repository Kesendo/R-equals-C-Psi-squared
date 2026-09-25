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

        // The control clause is chosen by a computed condition; test the branch that condition picks.
        // At N=4 the localized state's carrier pair also rises above the bar on both grids.
        var control = ((IInspectable)w).Children.Single(c => c.DisplayName.Contains("control"));
        bool singleBelowBarOnBoth = w.SingleCoarse.RiseCount == 0 && w.SingleFine.RiseCount == 0;
        Assert.False(singleBelowBarOnBoth, "precondition: at N=4 SingleExcitation rises above the bar");
        Assert.Contains("SingleExcitation = an above-bar finite row", control.Summary);
        Assert.DoesNotContain("sub-bar grid readings", control.Summary);
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

        // The N=3 control takes the other branch: the localized state stays below the bar on both grids.
        var control = ((IInspectable)w).Children.Single(c => c.DisplayName.Contains("control"));
        Assert.Contains("SingleExcitation = sub-bar grid readings", control.Summary);
    }

    [Fact]
    public void Global_NamedN3Window_MaximaFallStrictly()
    {
        // Raw, not at the reporting bar: every global maximum is below its predecessor (the closest pair
        // falls by 1.7e-3, twelve decades above rounding), so the largest positive predecessor delta is
        // exactly zero. A finite window, not an absence theorem.
        var w = new EnvelopeTheoremWitness(n: 3);
        Assert.Equal(0, w.GlobalBell.RiseCount);
        Assert.True(w.GlobalBell.IsNonIncreasing);
        Assert.True(w.GlobalBell.MaxRiseMagnitude == 0.0, $"raw max rise {w.GlobalBell.MaxRiseMagnitude:R}");
    }

    [Fact]
    public void LocalBell_EnvelopeRises_AboveTheReportingBar_OnNamedGrid()
    {
        // Five carrier-pair rises at 1600 points, the smallest 6.0e-3 and the largest 1.2e-2; they
        // persist under refinement (SymphonyTests' nested-grid control).
        var w = new EnvelopeTheoremWitness(n: 3);
        Assert.Equal(5, w.LocalBell.RiseCount);
        Assert.True(w.LocalBell.MaxRiseMagnitude > EnvelopeTheoremWitness.RiseReportingBar);
    }

    [Fact]
    public void SingleExcitation_SubBarRiseAt400Points_IsGoneAt1600()
    {
        // The refinement control on the witness's own readings. At 400 points the three-point apex
        // estimator leaves one sub-bar rise (5.5e-4); at 1600 points every maximum falls, raw.
        var w = new EnvelopeTheoremWitness(n: 3);
        Assert.True(w.SingleCoarse.MaxRiseMagnitude > 0.0, "expected the coarse-grid estimator rise");
        Assert.True(w.SingleCoarse.MaxRiseMagnitude < EnvelopeTheoremWitness.RiseReportingBar);
        Assert.True(w.SingleFine.IsNonIncreasing);
        Assert.True(w.SingleFine.MaxRiseMagnitude == 0.0, $"raw max rise {w.SingleFine.MaxRiseMagnitude:R}");
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
    public void Witness_SummaryNamesTheClaimAndChildrenSurfaceTheThreeReadings()
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
