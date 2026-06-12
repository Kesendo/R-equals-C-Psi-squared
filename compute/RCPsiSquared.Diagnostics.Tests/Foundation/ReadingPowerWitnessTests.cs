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
