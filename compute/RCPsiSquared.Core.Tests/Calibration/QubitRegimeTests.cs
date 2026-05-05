using RCPsiSquared.Core.Calibration;

namespace RCPsiSquared.Core.Tests.Calibration;

/// <summary>Tests for <see cref="QubitRegime"/>, the calibration-to-framework
/// bridge primitive. Empirical anchors are pinned to the 2026-04-25 Marrakesh
/// calibration via the same shared loader as <see cref="IbmCalibrationTests"/>;
/// the Marrakesh 91-day biography review (commit 787854f) confirms Q0 lives
/// permanently quantum-side (r mean 0.086 over 91 days, walk 0.022) and the
/// soft_break path [48, 49, 50] is uniform-classical.</summary>
public class QubitRegimeTests
{
    private static readonly Lazy<IReadOnlyList<QubitData>> Marrakesh20260425 = new(() =>
        IbmCalibration.Load(Path.Combine(FindRepoRoot(),
            "data", "ibm_calibration_snapshots",
            "ibm_marrakesh_calibrations_2026-04-25T11_28_00Z.csv")));

    [Fact]
    public void RStar_MatchesPolynomialFoldCatastropheConstant()
    {
        Assert.Equal(0.212755, QubitRegime.R_STAR, precision: 6);
    }

    [Fact]
    public void RParam_ComputesT2Over2T1()
    {
        Assert.Equal(0.5, QubitRegime.RParam(t1Us: 100, t2Us: 100), precision: 6);
        Assert.Equal(1.0, QubitRegime.RParam(t1Us: 100, t2Us: 200), precision: 6);
        Assert.Equal(0.1, QubitRegime.RParam(t1Us: 100, t2Us: 20), precision: 6);
    }

    [Fact]
    public void RParam_ZeroT1_ReturnsZero()
    {
        Assert.Equal(0.0, QubitRegime.RParam(t1Us: 0, t2Us: 100));
    }

    [Fact]
    public void Classify_BinaryDefault_NeverReturnsBoundary()
    {
        Assert.Equal(Regime.QuantumSide, QubitRegime.Classify(t1Us: 100, t2Us: 40));
        Assert.Equal(Regime.ClassicalSide, QubitRegime.Classify(t1Us: 100, t2Us: 50));
        Assert.Equal(Regime.ClassicalSide,
            QubitRegime.Classify(t1Us: 100, t2Us: 2.0 * 100 * QubitRegime.R_STAR));
    }

    [Fact]
    public void Classify_WithEpsilon_SurfacesBoundaryBand()
    {
        double tBoundary = 2.0 * 100 * QubitRegime.R_STAR;
        Assert.Equal(Regime.Boundary,
            QubitRegime.Classify(t1Us: 100, t2Us: tBoundary, epsilon: 0.01));
        Assert.Equal(Regime.QuantumSide,
            QubitRegime.Classify(t1Us: 100, t2Us: tBoundary - 5.0, epsilon: 0.01));
        Assert.Equal(Regime.ClassicalSide,
            QubitRegime.Classify(t1Us: 100, t2Us: tBoundary + 5.0, epsilon: 0.01));
    }

    [Fact]
    public void Q0_OnMarrakesh20260425_IsQuantumSide()
    {
        var qubits = Marrakesh20260425.Value;
        var q0 = qubits.Single(q => q.Qubit == 0);
        Assert.True(q0.IsQuantumSide,
            $"Q0 has T1={q0.T1Us:F1}, T2={q0.T2Us:F1}, r={q0.RParam:F4}; should be quantum-side");
        Assert.Equal(Regime.QuantumSide, q0.Regime);
    }

    [Fact]
    public void Q1_OnMarrakesh20260425_IsClassicalSide()
    {
        var qubits = Marrakesh20260425.Value;
        var q1 = qubits.Single(q => q.Qubit == 1);
        Assert.False(q1.IsQuantumSide,
            $"Q1 has T1={q1.T1Us:F1}, T2={q1.T2Us:F1}, r={q1.RParam:F4}; should be classical-side");
        Assert.Equal(Regime.ClassicalSide, q1.Regime);
    }

    [Fact]
    public void FrameworkSnapshotsPath_IsRegimeMixed()
    {
        var qubits = Marrakesh20260425.Value;
        var (q, b, c) = QubitRegime.PathComposition(qubits, new[] { 0, 1, 2 });
        Assert.Equal(1, q);
        Assert.Equal(0, b);
        Assert.Equal(2, c);
    }

    [Fact]
    public void SoftBreakPath_IsUniformClassical()
    {
        var qubits = Marrakesh20260425.Value;
        var (q, b, c) = QubitRegime.PathComposition(qubits, new[] { 48, 49, 50 });
        Assert.Equal(0, q);
        Assert.Equal(0, b);
        Assert.Equal(3, c);
    }

    [Fact]
    public void Apr25Best5Chain_IsUniformClassical()
    {
        var qubits = Marrakesh20260425.Value;
        var (q, b, c) = QubitRegime.PathComposition(qubits, new[] { 1, 2, 3, 4, 5 });
        Assert.Equal(0, q);
        Assert.Equal(0, b);
        Assert.Equal(5, c);
    }

    [Fact]
    public void PathComposition_RejectsUnknownQubit()
    {
        var qubits = Marrakesh20260425.Value;
        Assert.Throws<ArgumentException>(() =>
            QubitRegime.PathComposition(qubits, new[] { 0, 9999 }));
    }

    [Fact]
    public void QubitDataDerivedRegime_AgreesWithStaticClassifier()
    {
        var qubits = Marrakesh20260425.Value;
        foreach (var q in qubits)
        {
            Assert.Equal(QubitRegime.Classify(q.T1Us, q.T2Us), q.Regime);
            Assert.Equal(QubitRegime.RParam(q.T1Us, q.T2Us), q.RParam);
            Assert.Equal(QubitRegime.IsQuantumSide(q.T1Us, q.T2Us), q.IsQuantumSide);
        }
    }

    private static string FindRepoRoot()
    {
        var dir = new DirectoryInfo(AppContext.BaseDirectory);
        while (dir != null)
        {
            if (File.Exists(Path.Combine(dir.FullName, "MIRROR_THEORY.md"))
             && Directory.Exists(Path.Combine(dir.FullName, "compute")))
                return dir.FullName;
            dir = dir.Parent;
        }
        throw new InvalidOperationException(
            $"could not locate repository root starting from {AppContext.BaseDirectory}");
    }
}
