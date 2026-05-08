using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class C2DirectionAFalsificationProbeTests
{
    private readonly ITestOutputHelper _out;

    public C2DirectionAFalsificationProbeTests(ITestOutputHelper output)
    {
        _out = output;
    }

    [Fact]
    public void Probe_AcrossN5To8_ProducesDefinitiveVerdict()
    {
        var probe = C2DirectionAFalsificationProbe.Build();

        // Emit the per-N witness table so the verdict is human-inspectable in CI logs.
        _out.WriteLine($"BareFloor = {C2DirectionAFalsificationProbe.BareFloor:F6}");
        _out.WriteLine("");
        _out.WriteLine("  N | HWHM/Q* (E)| HWHM/Q* (I)| Δ (E)   | Δ (I)   | r (E)   | r (I)   | Δ_E−Δ_I | r_E−r_I | sign");
        _out.WriteLine("  --|------------|------------|---------|---------|---------|---------|---------|---------|----------");
        foreach (var w in probe.Witnesses)
        {
            _out.WriteLine(
                $"  {w.N} | {w.EndpointHwhmOverQpeak,10:F4} | {w.InteriorHwhmOverQpeak,10:F4} | " +
                $"{w.EndpointDelta,+7:F4} | {w.InteriorDelta,+7:F4} | " +
                $"{w.EndpointMagnitudeRatio,7:F4} | {w.InteriorMagnitudeRatio,7:F4} | " +
                $"{w.DeltaSignDifference,+7:F4} | {w.RatioSignDifference,+7:F4} | {w.SignRelation}");
        }
        _out.WriteLine("");
        _out.WriteLine($"Verdict: {probe.Verdict}");

        Assert.Equal(4, probe.Witnesses.Count);
        // We don't assert a specific verdict — the verdict IS the experiment result.
        // We do assert it is one of the three defined states (no implementation bug).
        Assert.True(
            probe.Verdict == DirectionAVerdict.StructurallyValid ||
            probe.Verdict == DirectionAVerdict.Falsified ||
            probe.Verdict == DirectionAVerdict.Ambiguous,
            $"Unexpected verdict value: {probe.Verdict}");
    }

    [Fact]
    public void Probe_IsTier2Verified()
    {
        var probe = C2DirectionAFalsificationProbe.Build(new[] { 5 });
        Assert.Equal(Tier.Tier2Verified, probe.Tier);
    }

    [Fact]
    public void Probe_AnchorReferences_F86Item1()
    {
        var probe = C2DirectionAFalsificationProbe.Build(new[] { 5 });
        Assert.Contains("PROOF_F86_QPEAK", probe.Anchor);
    }
}
