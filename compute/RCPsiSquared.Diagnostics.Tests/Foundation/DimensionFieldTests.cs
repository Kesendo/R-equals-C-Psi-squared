using System.Linq;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class DimensionFieldTests
{
    [Fact]
    public void DimensionField_Surfaces_MarksFlat_InBetweenRising_Polarity()
    {
        var field = new DimensionField(DimensionAxis.Crossover(N: 3, gamma: 0.5, thetaPoints: 13));
        var children = ((RCPsiSquared.Core.Inspection.IInspectable)field).Children.ToList();
        var labels = children.Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("marks"));
        Assert.Contains(labels, l => l.Contains("in-between"));
        Assert.Contains(labels, l => l.Contains("polarity"));
        Assert.Contains(labels, l => l.Contains("mirror"));
        // smoke: it renders to JSON without throwing and mentions the in-between
        var json = RCPsiSquared.Core.Inspection.InspectionJsonExporter.ToJson(field);
        Assert.Contains("in-between", json);
    }
}
