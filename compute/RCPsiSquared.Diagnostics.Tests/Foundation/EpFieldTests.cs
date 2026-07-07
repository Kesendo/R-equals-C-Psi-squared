using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class EpFieldTests
{
    private static System.Collections.Generic.List<IInspectable> Children(EpField f) =>
        ((IInspectable)f).Children.ToList();

    [Fact]
    public void Constructor_RejectsNonPositiveGamma() =>
        Assert.Throws<ArgumentOutOfRangeException>(() => new EpField(gamma: 0.0));

    [Fact]
    public void Constructor_RejectsNonPositiveGEff() =>
        Assert.Throws<ArgumentOutOfRangeException>(() => new EpField(gEff: 0.0));

    [Fact]
    public void Constructor_RejectsBadQRangeOrPoints()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new EpField(qLo: -0.1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new EpField(qLo: 3.0, qHi: 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new EpField(qPoints: 1));
    }

    [Fact]
    public void Field_HasFiveChildren() => Assert.Equal(5, Children(new EpField()).Count);

    [Fact]
    public void Field_Surfaces_Marks_Takt_Rotation_Defectiveness_Hardware()
    {
        var labels = Children(new EpField()).Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("marks"));
        Assert.Contains(labels, l => l.Contains("Takt") || l.Contains("coalesce"));
        Assert.Contains(labels, l => l.Contains("Rotation"));
        Assert.Contains(labels, l => l.Contains("defectiveness"));
        Assert.Contains(labels, l => l.Contains("hardware"));
    }

    [Fact]
    public void TaktChild_PinsAtFourGamma()
    {
        // The decay curve never exceeds 4γ₀ and reaches it (pins) for Q ≥ Q_EP.
        var f = new EpField(gamma: 1.0, gEff: 4.0 / 3.0, qLo: 0.3, qHi: 3.0, qPoints: 40);
        var takt = Children(f).First(c => c.DisplayName.Contains("Takt") || c.DisplayName.Contains("coalesce"));
        var curve = Assert.IsType<InspectablePayload.Curve>(takt.Payload);
        Assert.All(curve.Y, y => Assert.True(y <= 4.0 + 1e-9, $"decay {y} should not exceed 4γ₀"));
        Assert.True(curve.Y.Max() > 3.99, "decay should reach ~4γ₀ (pin) by the high-Q end");
    }

    [Fact]
    public void DefectivenessChild_PeaksAtOneAtTheToyEp()
    {
        // The toy 2×2 rate-channel overlap = min(x, 1/x) ≤ 1, reaching ~1 near its EP at Q_EP. This is a
        // property of the toy reduction (genuinely defective at its EP), NOT the physical chain block's
        // behaviour at Q_peak (non-normal there; its own real-axis defective EPs are F89's scattered seeds,
        // the retraction's "no real-axis defective EP" corrected 2026-07-07; LocalGlobalEpLink).
        var f = new EpField(gamma: 1.0, gEff: 4.0 / 3.0, qLo: 0.3, qHi: 3.0, qPoints: 60);
        var def = Children(f).First(c => c.DisplayName.Contains("defectiveness"));
        var curve = Assert.IsType<InspectablePayload.Curve>(def.Payload);
        Assert.All(curve.Y, y => Assert.True(y <= 1.0 + 1e-9));
        Assert.True(curve.Y.Max() > 0.97, "the toy overlap should approach 1 near the toy EP");
    }

    [Fact]
    public void RotationChild_IsZeroBelowEp_PositiveAbove()
    {
        var f = new EpField(gamma: 1.0, gEff: 4.0 / 3.0, qLo: 0.3, qHi: 3.0, qPoints: 60);
        var rot = Children(f).First(c => c.DisplayName.Contains("Rotation"));
        var curve = Assert.IsType<InspectablePayload.Curve>(rot.Payload);
        Assert.True(curve.Y[0] < 1e-6, "θ ≈ 0 well below the toy EP");
        Assert.True(curve.Y[^1] > 10.0, "θ has lifted off well above the toy EP");
    }

    [Fact]
    public void Field_RendersToJson()
    {
        var json = InspectionJsonExporter.ToJson(new EpField());
        Assert.Contains("Rotation", json);
        Assert.Contains("defectiveness", json);
    }
}
