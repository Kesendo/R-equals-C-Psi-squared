using RCPsiSquared.Core.Calibration;

namespace RCPsiSquared.Core.Tests.Calibration;

/// <summary>Pre-experiment audit workflow tests, pinned to 2026-04-25 Marrakesh.
/// Exercises the composition of <see cref="QubitRegime"/>,
/// <see cref="IbmCalibration.ChainScore"/>, and the CZ-coupling check into the
/// single <see cref="RegimeSummary"/> record.</summary>
public class RegimeSummaryTests
{
    private static Lazy<IReadOnlyList<QubitData>> Marrakesh20260425 => CalibrationFixtures.Marrakesh20260425;

    [Fact]
    public void FrameworkSnapshotsPath_IsRegimeMixedAndAddressable()
    {
        var qubits = Marrakesh20260425.Value;
        var s = RegimeSummary.For(qubits, new[] { 0, 1, 2 });
        Assert.Equal(RegimeVerdict.RegimeMixed, s.Verdict);
        Assert.False(s.IsRegimeUniform);
        Assert.True(s.IsAddressable);
        Assert.Equal(1, s.QuantumCount);
        Assert.Equal(0, s.BoundaryCount);
        Assert.Equal(2, s.ClassicalCount);
        Assert.Equal(597.27, s.Score, precision: 1);
    }

    [Fact]
    public void SoftBreakPath_IsUniformClassical()
    {
        var qubits = Marrakesh20260425.Value;
        var s = RegimeSummary.For(qubits, new[] { 48, 49, 50 });
        Assert.Equal(RegimeVerdict.UniformClassical, s.Verdict);
        Assert.True(s.IsRegimeUniform);
        Assert.True(s.IsAddressable);
        Assert.Equal(3, s.ClassicalCount);
        Assert.Equal(682.50, s.Score, precision: 1);
    }

    [Fact]
    public void Apr25Best5Chain_IsUniformClassical()
    {
        var qubits = Marrakesh20260425.Value;
        var s = RegimeSummary.For(qubits, new[] { 1, 2, 3, 4, 5 });
        Assert.Equal(RegimeVerdict.UniformClassical, s.Verdict);
        Assert.True(s.IsRegimeUniform);
        Assert.Equal(5, s.ClassicalCount);
    }

    [Fact]
    public void PerQubitRegime_CarriesRawCalibrationAndDerivedRegime()
    {
        var qubits = Marrakesh20260425.Value;
        var s = RegimeSummary.For(qubits, new[] { 0, 1, 2 });
        var q0 = s.Qubits[0];
        Assert.Equal(0, q0.Qubit);
        Assert.Equal(Regime.QuantumSide, q0.Regime);
        Assert.True(q0.RParam < QubitRegime.R_STAR);
        var q1 = s.Qubits[1];
        Assert.Equal(1, q1.Qubit);
        Assert.Equal(Regime.ClassicalSide, q1.Regime);
        Assert.True(q1.RParam > QubitRegime.R_STAR);
    }

    [Fact]
    public void NonCzCoupledPath_IsNotAddressable_VerdictTakesPrecedence()
    {
        var qubits = Marrakesh20260425.Value;
        // Q0 and Q50 are far apart on the heavy-hex; not CZ-coupled.
        // Note both are operational so the regime classification still runs,
        // but the verdict must still flag the topology problem first.
        var s = RegimeSummary.For(qubits, new[] { 0, 50 });
        Assert.False(s.AllCzCoupled);
        Assert.True(s.AllOperational);
        Assert.False(s.IsAddressable);
        Assert.Equal(RegimeVerdict.NotAddressable, s.Verdict);
    }

    [Fact]
    public void UnknownQubit_Throws()
    {
        var qubits = Marrakesh20260425.Value;
        Assert.Throws<ArgumentException>(() =>
            RegimeSummary.For(qubits, new[] { 0, 9999 }));
    }

    [Fact]
    public void EmptyPath_Throws()
    {
        var qubits = Marrakesh20260425.Value;
        Assert.Throws<ArgumentException>(() =>
            RegimeSummary.For(qubits, Array.Empty<int>()));
    }

    [Fact]
    public void Headline_IsHumanReadableOneLine()
    {
        var qubits = Marrakesh20260425.Value;
        var s = RegimeSummary.For(qubits, new[] { 0, 1, 2 });
        string h = s.ToHeadline();
        Assert.Contains("[0, 1, 2]", h);
        Assert.Contains("regime-mixed", h);
        Assert.Contains("addressable", h);
        Assert.Contains("score 597", h);
    }

    [Fact]
    public void SelectBestChain_LoweredThroughRegimeSummary_AllUniformClassical()
    {
        var qubits = Marrakesh20260425.Value;
        var calChain = IbmCalibration.SelectBestChain(qubits, length: 5);
        var s = RegimeSummary.For(qubits, calChain.QubitIds);
        Assert.True(s.IsAddressable);
        Assert.Equal(RegimeVerdict.UniformClassical, s.Verdict);
        Assert.Equal(calChain.Score, s.Score, precision: 1);
    }

    [Fact]
    public void RelatedConfirmations_FrameworkSnapshotsPath_FindsBothEntries()
    {
        var qubits = Marrakesh20260425.Value;
        var s = RegimeSummary.For(qubits, new[] { 0, 1, 2 });
        var hits = s.RelatedConfirmations().Select(c => c.Name).ToList();
        Assert.Contains("pi_protected_xiz_yzzy", hits);
        Assert.Contains("lebensader_skeleton_trace_decoupling", hits);
    }

    [Fact]
    public void RelatedConfirmations_MachineFilter_RestrictsToBackend()
    {
        var qubits = Marrakesh20260425.Value;
        var s = RegimeSummary.For(qubits, new[] { 48, 49, 50 });
        var marrakesh = s.RelatedConfirmations(machine: "ibm_marrakesh").ToList();
        Assert.NotEmpty(marrakesh);
        Assert.All(marrakesh, c => Assert.Equal("ibm_marrakesh", c.Machine));
        Assert.Empty(s.RelatedConfirmations(machine: "ibm_kingston"));
    }

    [Fact]
    public void RelatedConfirmations_UntestedPath_ReturnsEmpty()
    {
        var qubits = Marrakesh20260425.Value;
        // [4, 3, 2] is the Apr-25 best 3-chain, never run as documented path.
        var s = RegimeSummary.For(qubits, new[] { 4, 3, 2 });
        Assert.Empty(s.RelatedConfirmations());
    }

}
