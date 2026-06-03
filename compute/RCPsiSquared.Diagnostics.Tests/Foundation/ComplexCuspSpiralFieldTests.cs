using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class ComplexCuspSpiralFieldTests
{
    private static System.Collections.Generic.List<IInspectable> Children(ComplexCuspSpiralField f) =>
        ((IInspectable)f).Children.ToList();

    [Fact]
    public void Constructor_RejectsNonPositiveGamma() =>
        Assert.Throws<ArgumentOutOfRangeException>(() => new ComplexCuspSpiralField(gamma: 0.0));

    [Fact]
    public void Constructor_RejectsNegativeOmega() =>
        Assert.Throws<ArgumentOutOfRangeException>(() => new ComplexCuspSpiralField(omega: -0.1));

    [Fact]
    public void Constructor_RejectsTooFewOmegaPoints() =>
        Assert.Throws<ArgumentOutOfRangeException>(() => new ComplexCuspSpiralField(omegaPoints: 1));

    [Fact]
    public void Constructor_RejectsTMaxFactorNotPastCrossing() =>
        Assert.Throws<ArgumentOutOfRangeException>(() => new ComplexCuspSpiralField(tMaxFactor: 1.0));

    [Fact]
    public void Field_HasFiveChildren() => Assert.Equal(5, Children(new ComplexCuspSpiralField()).Count);

    [Fact]
    public void EveryChild_HasNameAndSummary()
    {
        foreach (var c in Children(new ComplexCuspSpiralField()))
        {
            Assert.False(string.IsNullOrWhiteSpace(c.DisplayName));
            Assert.False(string.IsNullOrWhiteSpace(c.Summary));
        }
    }

    [Fact]
    public void Field_Surfaces_Circle_Spiral_Winding_Hardware_Kinship()
    {
        var labels = Children(new ComplexCuspSpiralField()).Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("circle"));
        Assert.Contains(labels, l => l.Contains("spiral"));
        Assert.Contains(labels, l => l.Contains("winding"));
        Assert.Contains(labels, l => l.Contains("hardware"));
        Assert.Contains(labels, l => l.Contains("kinship"));
    }

    [Fact]
    public void CircleChild_TracesRadiusOneQuarter()
    {
        var circle = Children(new ComplexCuspSpiralField()).First(c => c.DisplayName.Contains("circle"));
        var curve = Assert.IsType<InspectablePayload.Curve>(circle.Payload);
        for (int i = 0; i < curve.X.Count; i++)
        {
            double r = Math.Sqrt(curve.X[i] * curve.X[i] + curve.Y[i] * curve.Y[i]);
            Assert.Equal(0.25, r, 9);
        }
    }

    [Fact]
    public void SpiralChild_NeverExceedsTheStartMagnitude()
    {
        var spiral = Children(new ComplexCuspSpiralField(gamma: 0.05, omega: 0.4))
            .First(c => c.DisplayName.Contains("spiral"));
        var curve = Assert.IsType<InspectablePayload.Curve>(spiral.Payload);
        for (int i = 0; i < curve.X.Count; i++)
        {
            double r = Math.Sqrt(curve.X[i] * curve.X[i] + curve.Y[i] * curve.Y[i]);
            Assert.True(r <= 1.0 / 3.0 + 1e-9, $"radius {r} exceeds the t=0 magnitude 1/3");
        }
    }

    [Fact]
    public void WindingChild_CrossingAngleMovesWithOmega()
    {
        // The third reading's point: the crossing time is flat across Ω (radial law Ω-independent),
        // the crossing angle is not. The curve is Ω → crossing angle°, so its endpoints differ.
        var winding = Children(new ComplexCuspSpiralField(gamma: 0.05, omega: 0.4))
            .First(c => c.DisplayName.Contains("winding"));
        var curve = Assert.IsType<InspectablePayload.Curve>(winding.Payload);
        Assert.NotEqual(curve.Y[0], curve.Y[^1], 3);
    }

    [Fact]
    public void Field_RendersToJson_CarriesTheCircleStory()
    {
        var json = InspectionJsonExporter.ToJson(new ComplexCuspSpiralField());
        Assert.Contains("circle", json);
        Assert.Contains("spiral", json);
    }
}
