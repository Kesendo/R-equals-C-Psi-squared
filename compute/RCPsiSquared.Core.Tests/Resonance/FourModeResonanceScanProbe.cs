using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Resonance;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Resonance;

public class FourModeResonanceScanProbe(ITestOutputHelper output)
{
    [Fact(Skip = "diagnostic — extends Q-grid past 4 to discover actual 4-mode Q_peak; remove Skip to read the numbers.")]
    public void Probe_FourMode_QPeakValues_AcrossWideGrid()
    {
        var qGrid = ResonanceScan.LinearQGrid(0.2, 30.1, 300);

        foreach (int N in new[] { 5, 6, 7, 8 })
        {
            var block = new CoherenceBlock(N, 1, 0.05);
            var fullCurve = new ResonanceScan(block).ComputeKCurve(qGrid);
            var fourCurve = new FourModeResonanceScan(block).ComputeKCurve(qGrid);

            var fI = fullCurve.Peak(BondClass.Interior);
            var fE = fullCurve.Peak(BondClass.Endpoint);
            var f4I = fourCurve.Peak(BondClass.Interior);
            var f4E = fourCurve.Peak(BondClass.Endpoint);

            output.WriteLine($"N={N} Interior full   Q_peak={fI.QPeak,7:F3} HWHM/Q={fI.HwhmLeftOverQPeak ?? double.NaN,6:F3} K_max={fI.KMax,8:F4}");
            output.WriteLine($"N={N} Interior 4-mode Q_peak={f4I.QPeak,7:F3} HWHM/Q={f4I.HwhmLeftOverQPeak ?? double.NaN,6:F3} K_max={f4I.KMax,8:F4}");
            output.WriteLine($"N={N} Endpoint full   Q_peak={fE.QPeak,7:F3} HWHM/Q={fE.HwhmLeftOverQPeak ?? double.NaN,6:F3} K_max={fE.KMax,8:F4}");
            output.WriteLine($"N={N} Endpoint 4-mode Q_peak={f4E.QPeak,7:F3} HWHM/Q={f4E.HwhmLeftOverQPeak ?? double.NaN,6:F3} K_max={f4E.KMax,8:F4}");
        }
    }
}
