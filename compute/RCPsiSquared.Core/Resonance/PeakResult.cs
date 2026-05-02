namespace RCPsiSquared.Core.Resonance;

public sealed record PeakResult(double QPeak, double KMax, double? HwhmLeft, double? HwhmRight)
{
    public double? HwhmLeftOverQPeak => HwhmLeft.HasValue ? HwhmLeft.Value / QPeak : null;
    public double? HwhmRightOverQPeak => HwhmRight.HasValue ? HwhmRight.Value / QPeak : null;
    public double? Asymmetry => (HwhmLeft.HasValue && HwhmRight.HasValue)
        ? HwhmRight.Value / HwhmLeft.Value : null;
}
