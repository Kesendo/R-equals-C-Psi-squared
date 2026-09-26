using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class ApproachFamilyFieldTests
{
    private static System.Collections.Generic.List<IInspectable> Children(ApproachFamilyField f) =>
        ((IInspectable)f).Children.ToList();

    [Fact]
    public void Constructor_RejectsNonPositiveGamma() =>
        Assert.Throws<ArgumentOutOfRangeException>(() => new ApproachFamilyField(gamma: 0.0));

    [Fact]
    public void Constructor_RejectsBadSRange()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new ApproachFamilyField(sHi: 1.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ApproachFamilyField(sLo: -0.01));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ApproachFamilyField(sLo: 0.0, sHi: 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ApproachFamilyField(sLo: 0.9, sHi: 0.4));
    }

    [Fact]
    public void Field_AllowsAndRendersTheZeroEndpoint()
    {
        var field = new ApproachFamilyField(sLo: 0.0, sHi: 1.0, sPoints: 5);
        Assert.Contains("s∈[0, 1]", field.DisplayName, StringComparison.Ordinal);
        Assert.Contains("s=0 has w₀=w₁=0", field.Summary, StringComparison.Ordinal);

        var starts = Assert.IsType<InspectablePayload.Curve>(
            Children(field).First(c => c.DisplayName.Contains("starts")).Payload);
        Assert.Equal(0.0, starts.X[0], 12);
        Assert.Equal(0.0, starts.Y[0], 12);

        var shape = Assert.IsType<InspectablePayload.Curve>(
            Children(field).First(c => c.DisplayName.Contains("shape parameter")).Payload);
        Assert.Equal(0.0, shape.X[0], 12);
        Assert.Equal(0.0, shape.Y[0], 12);
    }

    [Fact]
    public void Constructor_RejectsTooFewPointsOrBadTMax()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new ApproachFamilyField(sPoints: 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ApproachFamilyField(tMaxFactor: 1.0));
    }

    [Fact]
    public void Field_HasFiveChildren() => Assert.Equal(5, Children(new ApproachFamilyField()).Count);

    [Fact]
    public void Field_Surfaces_Carrier_Starts_Threshold_Shape_Collapse()
    {
        var labels = Children(new ApproachFamilyField()).Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("shared slow exponential"));
        Assert.Contains(labels, l => l.Contains("starts"));
        Assert.Contains(labels, l => l.Contains("threshold"));
        Assert.Contains(labels, l => l.Contains("shape parameter"));
        Assert.Contains(labels, l => l.Contains("Bell+ specialization"));
    }

    [Fact]
    public void StartsChild_IsSOverThree()
    {
        var field = new ApproachFamilyField(sLo: 0.3, sHi: 1.0, sPoints: 8);
        var starts = Children(field)
            .First(c => c.DisplayName.Contains("starts"));
        var curve = Assert.IsType<InspectablePayload.Curve>(starts.Payload);
        for (int i = 0; i < curve.X.Count; i++)
            Assert.Equal(curve.X[i] / 3.0, curve.Y[i], 9);

        Assert.Contains("pure-state concurrence", field.Summary, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("one third", field.Summary, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("pure-state concurrence", starts.Summary, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("one third", starts.Summary, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("the entanglement itself", field.Summary, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("start height is the entanglement", starts.Summary, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void ThresholdChild_CrossesZeroAtThreeQuarters()
    {
        var field = new ApproachFamilyField(sLo: 0.3, sHi: 1.0, sPoints: 9);
        var thr = Children(field)
            .First(c => c.DisplayName.Contains("threshold"));
        var curve = Assert.IsType<InspectablePayload.Curve>(thr.Payload);
        Assert.True(curve.Y[0] < 0, "low s starts below ¼");
        Assert.True(curve.Y[^1] > 0, "high s starts above ¼");
        Assert.Contains("iff γ>0 and s>3/4", field.Summary, StringComparison.Ordinal);
        Assert.Contains("γ=0", field.Summary, StringComparison.Ordinal);
        Assert.Contains("constant", field.Summary, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("iff γ>0 and s>3/4", thr.Summary, StringComparison.Ordinal);
        Assert.Contains("touches", thr.Summary, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("does not cross", thr.Summary, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void ShapeChild_IsHarmonicFractionSSquaredOverTwo()
    {
        var shape = Children(new ApproachFamilyField(sLo: 0.3, sHi: 1.0, sPoints: 8))
            .First(c => c.DisplayName.Contains("shape parameter"));
        var curve = Assert.IsType<InspectablePayload.Curve>(shape.Payload);
        for (int i = 0; i < curve.X.Count; i++)
            Assert.Equal(0.5 * curve.X[i] * curve.X[i], curve.Y[i], 9);
    }

    [Fact]
    public void Field_RestrictsTheCarrierToNonzeroFamilyMembers()
    {
        var field = new ApproachFamilyField();
        Assert.Contains("every nonzero member", field.Summary, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("every member shares", field.Summary, StringComparison.OrdinalIgnoreCase);

        var carrier = Children(field).First(c => c.DisplayName.Contains("shared slow exponential"));
        Assert.Contains("every nonzero member", carrier.Summary, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void Field_RendersToJson()
    {
        var json = InspectionJsonExporter.ToJson(new ApproachFamilyField());
        Assert.Contains("carrier", json);
        Assert.Contains("threshold", json);
    }
}
