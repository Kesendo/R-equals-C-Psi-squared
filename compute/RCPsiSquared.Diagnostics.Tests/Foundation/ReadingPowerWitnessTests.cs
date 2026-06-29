using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class ReadingPowerWitnessTests
{
    static ReadingPowerWitness W() => new ReadingPowerWitness(n: 4);

    [Fact]
    public void NoEpPeak_InAnyBasis_FiMonotoneInQ()
    {
        var w = W();
        foreach (var basis in new[] { ReadoutBasis.Z, ReadoutBasis.X, ReadoutBasis.Y })
            Assert.True(w.IsMonotoneInQ(basis),
                $"FI(Q) must be monotone increasing in {basis} (no EP peak in this readout)");
    }

    [Fact]
    public void ResolutionLaw_ZBasis_IsLinearInQ()
    {
        var w = W();
        Assert.True(w.ZLinearityRelativeResidual() < 0.25,
            $"residual {w.ZLinearityRelativeResidual()}");
    }

    [Fact]
    public void BasisOrdering_CoherenceReadoutFallsFasterThanPopulation()
    {
        var w = W();
        Assert.True(w.SpanRatio(ReadoutBasis.X) > 10 * w.SpanRatio(ReadoutBasis.Z));
    }

    [Fact]
    public void Children_CarryTheThreeVerdictNodes_AndCurves()
    {
        var labels = ((IInspectable)W()).Children.Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("resolution law"));
        Assert.Contains(labels, l => l.Contains("no EP peak"));
        Assert.Contains(labels, l => l.Contains("basis ordering"));
    }

    [Fact]
    public void Children_CarryTheDecodeDemo_ThreePlantedDefects_AllExactBond()
    {
        var demo = ((IInspectable)W()).Children.Single(c => c.DisplayName.Contains("decoder (read a planted defect)"));
        var plantCases = ((IInspectable)demo).Children.Where(c => c.DisplayName.StartsWith("planted:")).ToList();
        Assert.Equal(3, plantCases.Count);
        foreach (var c in plantCases)
            Assert.Contains("match", c.Summary);   // each planted case decodes to its true bond
    }

    [Fact]
    public void DecodeDemo_CarriesTheDeLossNode_AlphaAmbiguous_DeviationResolvesWithSign()
    {
        // The de-loss before/after (spec §4.3): at the N=5 mirror pair the α path is AMBIGUOUS and the
        // signed deviation path RESOLVES it with the correct (weakened) sign. The node reports the actual
        // computed numbers live.
        var demo = ((IInspectable)W()).Children.Single(c => c.DisplayName.Contains("decoder (read a planted defect)"));
        var deloss = ((IInspectable)demo).Children.Single(c => c.DisplayName.Contains("de-loss"));
        Assert.Contains("AMBIGUOUS", deloss.Summary);          // the α path flags it
        Assert.Contains("RESOLVED", deloss.Summary);           // the deviation path resolves it
        Assert.Contains("weakened, sign read", deloss.Summary); // the sign is recovered
    }

    [Fact]
    public void Sweep_ComputedOnce()
    {
        var w = W();
        _ = ((IInspectable)w).Summary;
        _ = ((IInspectable)w).Children.Select(c => c.Summary).ToList();
        _ = w.IsMonotoneInQ(ReadoutBasis.Z);
        Assert.Equal(1, w.SweepCount);
    }

    [Fact]
    public void Guard_RejectsNOutsideThreeToFive()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new ReadingPowerWitness(2));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ReadingPowerWitness(6));
    }
}
