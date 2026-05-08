using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class C2DirectionCFalsificationProbeTests
{
    private readonly ITestOutputHelper _out;

    public C2DirectionCFalsificationProbeTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Probe_AcrossN5To8_ProducesDefinitiveVerdict()
    {
        var probe = C2DirectionCFalsificationProbe.Build();

        _out.WriteLine($"BareFloor = {C2DirectionCFalsificationProbe.BareFloor:F6}");
        _out.WriteLine($"MatchTolerance = {C2DirectionCFalsificationProbe.MatchTolerance:F4}, " +
                       $"FloorTolerance = {C2DirectionCFalsificationProbe.FloorTolerance:F4}");
        _out.WriteLine("");
        _out.WriteLine("  N | 4-mode E | 4-mode I | empir. E | empir. I | (4mode − empir) E | (4mode − empir) I | (4mode − floor) E | (4mode − floor) I");
        _out.WriteLine("  --|----------|----------|----------|----------|-------------------|-------------------|-------------------|-------------------");
        foreach (var w in probe.Witnesses)
        {
            _out.WriteLine(
                $"  {w.N} | {w.FourModeEndpointHwhmRatio,8:F4} | {w.FourModeInteriorHwhmRatio,8:F4} | " +
                $"{w.EmpiricalEndpointHwhmRatio,8:F4} | {w.EmpiricalInteriorHwhmRatio,8:F4} | " +
                $"{w.FourModeEmpiricalGapEndpoint,+17:F4} | {w.FourModeEmpiricalGapInterior,+17:F4} | " +
                $"{w.FourModeFloorGapEndpoint,+17:F4} | {w.FourModeFloorGapInterior,+17:F4}");
        }
        _out.WriteLine("");
        _out.WriteLine($"Verdict: {probe.Verdict}");

        Assert.Equal(4, probe.Witnesses.Count);
        Assert.True(
            probe.Verdict == DirectionCVerdict.FourModeMatchesEmpirical ||
            probe.Verdict == DirectionCVerdict.FourModeAtBareFloor ||
            probe.Verdict == DirectionCVerdict.FourModePartial);
    }

    [Fact]
    public void Probe_IsTier2Verified()
    {
        var probe = C2DirectionCFalsificationProbe.Build(new[] { 5 });
        Assert.Equal(Tier.Tier2Verified, probe.Tier);
    }

    [Fact]
    public void Probe_AnchorReferences_F86Item1()
    {
        var probe = C2DirectionCFalsificationProbe.Build(new[] { 5 });
        Assert.Contains("PROOF_F86_QPEAK", probe.Anchor);
    }
}
