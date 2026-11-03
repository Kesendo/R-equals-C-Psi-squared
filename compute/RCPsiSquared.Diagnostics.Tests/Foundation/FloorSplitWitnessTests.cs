using System;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>F86's K split by the free-fermion floor: in the c = 2 block the floor family's span is
/// invariant under every single bond and under L and L† alike, so the per-site readout amplitude
/// splits linearly into a floor part and a rest, and the quadratic K into floor², rest² and the
/// interference between them, at each bond's canonical Q_peak.</summary>
public class FloorSplitWitnessTests
{
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void TheThreeTerms_AddUpToTheProductionScannersEnvelope(int n)
    {
        var w = new FloorSplitWitness(n);
        Assert.Equal(n - 1, w.Bonds.Count);
        foreach (var b in w.Bonds)
        {
            Assert.True(b.FamilyClosedUnderThisBond, $"N={n} bond {b.Bond}: family not closed under the bond alone");
            double sum = b.KFloor + b.KRest + b.KCross;
            // two eigen routes (this witness and ResonanceScan) on the same Q: agreement at rounding scale
            Assert.True(Math.Abs(Math.Abs(sum) - b.ScanEnvelope) <= 1e-9 * b.ScanEnvelope,
                $"N={n} bond {b.Bond}: |floor²+rest²+cross| = {Math.Abs(sum)} vs scanner {b.ScanEnvelope}");
            Assert.True(Math.Abs(sum - b.K) <= 1e-12 * Math.Abs(b.K));
        }
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void TheFloorTerm_OnTwoRoutes_BlockProjectionMeetsTheSingleParticleClosedForm(int n)
    {
        var w = new FloorSplitWitness(n);
        double scale = w.Bonds.Max(b => Math.Abs(b.K));
        foreach (var b in w.Bonds)
            Assert.True(Math.Abs(b.KFloor - b.KFloorClosedForm) <= 1e-10 * scale,
                $"N={n} bond {b.Bond}: block {b.KFloor} vs closed form {b.KFloorClosedForm}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    public void AtTheEndpointBond_TheFloorIsTheLargestTerm(int n)
    {
        var e = new FloorSplitWitness(n).Bonds[0];
        Assert.True(e.FloorShare > e.RestShare && e.FloorShare > Math.Abs(e.CrossShare),
            $"N={n}: floor {e.FloorShare}, rest {e.RestShare}, cross {e.CrossShare}");
    }

    [Fact]
    public void AtTheEndpointBond_TheFloorShareFalls_AndTheInterferenceGrows()
    {
        var floor = Enumerable.Range(5, 4).Select(n => new FloorSplitWitness(n).Bonds[0].FloorShare).ToArray();
        for (int i = 1; i < floor.Length; i++)
            Assert.True(floor[i] < floor[i - 1], $"floor share N={5 + i}: {floor[i]} not below N={4 + i}: {floor[i - 1]}");
        var cross = Enumerable.Range(5, 3).Select(n => Math.Abs(new FloorSplitWitness(n).Bonds[0].CrossShare)).ToArray();
        for (int i = 1; i < cross.Length; i++)
            Assert.True(cross[i] > cross[i - 1], $"|cross share| N={5 + i}: {cross[i]} not above N={4 + i}: {cross[i - 1]}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void AtOddN_TheBondsIntoTheMiddleSite_DoNotSeeTheFloor(int n)
    {
        var w = new FloorSplitWitness(n);
        int mid = (n - 1) / 2;
        double scale = w.Bonds.Max(b => Math.Abs(b.K));
        foreach (var b in w.Bonds)
        {
            bool touchesMiddle = b.Bond == mid - 1 || b.Bond == mid;
            if (touchesMiddle)
                Assert.True(Math.Abs(b.KFloorClosedForm) <= 1e-12 * scale,
                    $"N={n} bond {b.Bond}: floor term {b.KFloorClosedForm} at a bond into the middle site");
            else
                Assert.True(Math.Abs(b.KFloorClosedForm) > 1e-4 * scale,
                    $"N={n} bond {b.Bond}: floor term {b.KFloorClosedForm} unexpectedly silent");
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void AtOddN_TheBondsIntoTheMiddleSite_KeepOnlyTheFloorStateAgainstTheRestResponse(int n)
    {
        var w = new FloorSplitWitness(n);
        int mid = (n - 1) / 2;
        double scale = w.Bonds.Max(b => Math.Abs(b.K));
        foreach (var b in w.Bonds)
        {
            bool touchesMiddle = b.Bond == mid - 1 || b.Bond == mid;
            if (touchesMiddle)
                // a symmetry zero on the block route, with no exact route (ρ and ∂ρ come through the eigen-Duhamel
                // route of a non-normal L); the budget 1e-10·scale is borrowed from AllRoutesMeet's gate of the block
                // floor² against its closed form, not derived from an error model; observed 1e-17 to 1e-14 of scale
                Assert.True(Math.Abs(b.KCrossRestStateFloorResponse) <= 1e-10 * scale,
                    $"N={n} bond {b.Bond}: rest state against floor response {b.KCrossRestStateFloorResponse} at a bond into the middle site");
            else
                Assert.True(Math.Abs(b.KCrossRestStateFloorResponse) > 1e-4 * scale,
                    $"N={n} bond {b.Bond}: rest state against floor response {b.KCrossRestStateFloorResponse} unexpectedly silent");
            Assert.True(Math.Abs(b.KCrossFloorStateRestResponse) > 1e-4 * scale,
                $"N={n} bond {b.Bond}: floor state against rest response {b.KCrossFloorStateRestResponse} unexpectedly silent");
        }
    }

    [Fact]
    public void Constructor_RejectsOutOfRange()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new FloorSplitWitness(3));
        Assert.Throws<ArgumentOutOfRangeException>(() => new FloorSplitWitness(9));
        Assert.Throws<ArgumentOutOfRangeException>(() => new FloorSplitWitness(5, -0.05));
    }

    [Fact]
    public void Inspect_ShowsTheGateAndOneNodePerBond()
    {
        var w = new FloorSplitWitness(5);
        var kids = ((IInspectable)w).Children.ToList();
        Assert.Equal(1 + 4, kids.Count);
        Assert.True(w.AllRoutesMeet);
        Assert.StartsWith("routes meet:", kids[0].Summary);
    }
}
