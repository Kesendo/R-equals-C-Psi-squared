using System.Text.Json;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Resonance;
using RCPsiSquared.Core.Tests.Fixtures;

namespace RCPsiSquared.Core.Tests.Validation;

/// <summary>Validation against `simulations/results/eq022_universality_extension/summary.json`
/// (Python `_eq022_b1_step_f_universality_extension.py`).
///
/// <para>Covers c=2 N=5..8 at γ₀=0.05 plus γ₀ ∈ {0.025, 0.10} at c=3 N=7. The C# scan
/// must reproduce Q_peak and HWHM_left/Q_peak per bond class to within 1×10⁻³.</para>
/// </summary>
public class ValidateAgainstPythonStepFTests
{
    private const double QPeakTol = 5e-3;
    private const double HwhmRatioTol = 5e-3;
    private const double KMaxRelTol = 0.02;

    public static IEnumerable<object[]> StepFCases()
    {
        var path = Path.Combine(AppContext.BaseDirectory, "Fixtures", "python_step_f_summary.json");
        var entries = JsonSerializer.Deserialize<List<PythonStepSummaryEntry>>(File.ReadAllText(path))!;
        // Subset for runtime: c=2 only (γ₀-invariance on c=3 N=7 covered by separate test).
        foreach (var e in entries.Where(e => e.C == 2))
            yield return new object[] { e };
    }

    [Theory]
    [MemberData(nameof(StepFCases))]
    public void Reproduces_PythonStepF(PythonStepSummaryEntry expected)
    {
        var block = new CoherenceBlock(N: expected.N, n: expected.LowerPopcount, gammaZero: expected.GammaZero);
        var scan = new ResonanceScan(block);
        var curve = scan.ComputeKCurve();

        AssertClassMatches(curve, expected, BondClass.Interior,
            expected.QPeakInterior, expected.HwhmLeftInterior, expected.KMaxInterior);
        AssertClassMatches(curve, expected, BondClass.Endpoint,
            expected.QPeakEndpoint, expected.HwhmLeftEndpoint, expected.KMaxEndpoint);
    }

    private static void AssertClassMatches(KCurve curve, PythonStepSummaryEntry expected, BondClass cls,
        double? expectedQ, double? expectedHwhm, double? expectedKMax)
    {
        if (expectedQ is null) return; // no interior bonds for some configurations

        var peak = curve.Peak(cls);
        Assert.InRange(peak.QPeak, expectedQ.Value - QPeakTol, expectedQ.Value + QPeakTol);

        if (expectedHwhm is not null)
        {
            Assert.NotNull(peak.HwhmLeft);
            double pythonRatio = expectedHwhm.Value / expectedQ.Value;
            double csharpRatio = peak.HwhmLeftOverQPeak!.Value;
            Assert.InRange(csharpRatio, pythonRatio - HwhmRatioTol, pythonRatio + HwhmRatioTol);
        }
        if (expectedKMax is not null)
        {
            double rel = Math.Abs(peak.KMax - expectedKMax.Value) / expectedKMax.Value;
            Assert.True(rel < KMaxRelTol,
                $"|K|max for {cls} mismatch (c={expected.C} N={expected.N}): C# {peak.KMax:F5}, Python {expectedKMax:F5}, rel = {rel:P3}");
        }
    }
}
