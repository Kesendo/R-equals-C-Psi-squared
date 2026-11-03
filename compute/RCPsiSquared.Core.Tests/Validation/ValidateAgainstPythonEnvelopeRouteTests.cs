using System.Text.Json;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Resonance;
using RCPsiSquared.Core.Tests.Fixtures;

namespace RCPsiSquared.Core.Tests.Validation;

/// <summary>The C# envelope scan against an independent Python route,
/// <c>Fixtures/python_envelope_route_summary.json</c> (written by
/// <c>simulations/f86_envelope_second_route.py</c>).
///
/// <para>The two routes share the observable and the Q grid (153 points on [0.2, 4]) and the
/// peak finder, and nothing else on the way to max_t: the Python route builds the c = 2 block by
/// the F89 eigendecomposition (J halved) and the c ≥ 3 blocks by a numpy transcription of the
/// block builder, uses the spatial-sum kernel through its rank-N factor, evaluates Duhamel in a
/// vectorised closed form and maximises in t by a dense scan on [0.3, 40] with a bounded Brent
/// step, where C# seeds a 45-point grid on [0.2, 2.4]·t_peak and refines by golden section.</para>
///
/// <para>Gate 1e-10, the same as the per-bond pins: the residual between the routes is the
/// eigensolvers' round-off of K (~10⁻¹³) magnified by the parabolic vertex on a grid of step
/// 0.025; the t-maximisers' tolerance enters K only at second order. Its size is the spread
/// between mirror-equal bonds within one route: up to 1.9·10⁻¹¹ at c = 2 through N = 10, up to
/// 1.0·10⁻¹⁰ per bond at c = 3. The class values compared here, means over bonds, agree to about
/// 10⁻¹¹.</para>
///
/// <para>Blocks run here: those of dimension ≤ 300 (c = 2, N = 5..8; c = 3, N = 5, 6). The larger
/// entries stay in the fixture as the record of the second route.</para></summary>
public class ValidateAgainstPythonEnvelopeRouteTests
{
    private const double Tol = 1e-10;
    private const int MaxBlockDimHere = 300;

    private static IEnumerable<EnvelopeRouteEntry> Entries()
    {
        var path = Path.Combine(AppContext.BaseDirectory, "Fixtures", "python_envelope_route_summary.json");
        return JsonSerializer.Deserialize<List<EnvelopeRouteEntry>>(File.ReadAllText(path))!;
    }

    private static long BlockDim(EnvelopeRouteEntry e) =>
        Binomial(e.N, e.LowerPopcount) * Binomial(e.N, e.LowerPopcount + 1);

    private static long Binomial(int n, int k)
    {
        long r = 1;
        for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
        return r;
    }

    public static IEnumerable<object[]> Cases() =>
        Entries().Where(e => BlockDim(e) <= MaxBlockDimHere).Select(e => new object[] { e });

    [Fact]
    public void Fixture_HoldsTheExpectedBlocks()
    {
        var labels = Entries().Select(e => e.ToString()).ToHashSet();
        foreach (var l in new[] { "c2_N5", "c2_N6", "c2_N7", "c2_N8", "c3_N5", "c3_N6", "c3_N7", "c3_N8", "c4_N7", "c4_N8" })
            Assert.Contains(l, labels);
        Assert.All(Entries(), e => Assert.Equal(e.C - 1, e.LowerPopcount));
    }

    [Theory]
    [MemberData(nameof(Cases))]
    public void ClassPeaks_MatchTheSecondRoute(EnvelopeRouteEntry expected)
    {
        var block = new CoherenceBlock(N: expected.N, n: expected.LowerPopcount, gammaZero: expected.GammaZero);
        var curve = new ResonanceScan(block).ComputeKCurve();

        AssertClass(curve.Peak(BondClass.Interior), expected.QPeakInterior, expected.HwhmRatioInterior, expected.KMaxInterior);
        AssertClass(curve.Peak(BondClass.Endpoint), expected.QPeakEndpoint, expected.HwhmRatioEndpoint, expected.KMaxEndpoint);
    }

    private static void AssertClass(PeakResult peak, double q, double ratio, double kMax)
    {
        Assert.Equal(q, peak.QPeak, Tol);
        Assert.NotNull(peak.HwhmLeftOverQPeak);
        Assert.Equal(ratio, peak.HwhmLeftOverQPeak!.Value, Tol);
        Assert.Equal(kMax, peak.KMax, Tol);
    }
}
