using RCPsiSquared.Core.Calibration;
using RCPsiSquared.Core.ChainSystems;

namespace RCPsiSquared.Core.Tests.Calibration;

/// <summary>Regression + parser tests for <see cref="IbmCalibration"/>.
///
/// <para>Empirical anchors are pinned to the 2026-04-25 Marrakesh CSV
/// (<c>data/ibm_calibration_snapshots/ibm_marrakesh_calibrations_2026-04-25T11_28_00Z.csv</c>).
/// They mirror the Python reference output of <c>simulations/ibm_calibration.py</c>:
/// 156 operational qubits, best 3-chain = [4, 3, 2] with score ≈ 867.07,
/// best 5-chain = [1, 2, 3, 4, 5] with score ≈ 1246.58, and the documented
/// path-quality ordering <c>ChainScore([48, 49, 50]) &gt; ChainScore([0, 1, 2])</c>
/// (~682.5 vs ~597.3) which underlies the 23× cleaner truly-baseline observed
/// on path [48, 49, 50] in the 2026-04-26 soft_break replication.</para>
/// </summary>
public class IbmCalibrationTests
{
    private static readonly Lazy<IReadOnlyList<QubitData>> Marrakesh20260425 = new(() =>
        IbmCalibration.Load(Path.Combine(FindRepoRoot(),
            "data", "ibm_calibration_snapshots",
            "ibm_marrakesh_calibrations_2026-04-25T11_28_00Z.csv")));

    [Fact]
    public void Load_Marrakesh20260425_Has156OperationalQubits()
    {
        var qubits = Marrakesh20260425.Value;
        Assert.Equal(156, qubits.Count);
        Assert.All(qubits, q => Assert.True(q.Operational));
    }

    [Fact]
    public void Path_48_49_50_OutscoresPath_0_1_2()
    {
        var qubits = Marrakesh20260425.Value;
        double good = IbmCalibration.ChainScore(qubits, new[] { 48, 49, 50 });
        double first = IbmCalibration.ChainScore(qubits, new[] { 0, 1, 2 });
        Assert.True(good > first,
            $"path [48,49,50] (score={good:F2}) should outscore [0,1,2] (score={first:F2}) " +
            "per the 2026-04-26 soft_break replication finding");
        Assert.Equal(682.50, good, precision: 1);
        Assert.Equal(597.27, first, precision: 1);
    }

    [Fact]
    public void BestChain_Length3_ReturnsContiguousCzCoupledTriple()
    {
        var qubits = Marrakesh20260425.Value;
        var (path, score) = IbmCalibration.BestChain(qubits, length: 3);
        Assert.Equal(3, path.Count);
        Assert.Equal(new[] { 4, 3, 2 }, path);
        Assert.Equal(867.07, score, precision: 1);
    }

    [Fact]
    public void BestChain_Length5_ReturnsContiguousCzCoupledQuintuple()
    {
        var qubits = Marrakesh20260425.Value;
        var (path, score) = IbmCalibration.BestChain(qubits, length: 5);
        Assert.Equal(5, path.Count);
        Assert.Equal(new[] { 1, 2, 3, 4, 5 }, path);
        Assert.Equal(1246.58, score, precision: 1);
    }

    [Fact]
    public void BestChain_PathBondsAreCzNeighbours()
    {
        var qubits = Marrakesh20260425.Value;
        var byId = qubits.ToDictionary(q => q.Qubit);
        var (path, _) = IbmCalibration.BestChain(qubits, length: 5);
        for (int i = 0; i < path.Count - 1; i++)
        {
            int a = path[i], b = path[i + 1];
            bool bonded = byId[a].CzNeighbours.ContainsKey(b)
                       || byId[b].CzNeighbours.ContainsKey(a);
            Assert.True(bonded, $"qubits {a} and {b} are not CZ-coupled in the calibration graph");
        }
    }

    [Fact]
    public void BestQubits_TopOne_IsQubit35_OnMarrakesh20260425()
    {
        var qubits = Marrakesh20260425.Value;
        var top = IbmCalibration.BestQubits(qubits, k: 1);
        Assert.Single(top);
        Assert.Equal(35, top[0].Qubit);
        Assert.Equal(453.70, IbmCalibration.ScoreQubit(top[0]), precision: 1);
    }

    [Fact]
    public void ScoreQubit_NonOperational_ReturnsNegativeInfinity()
    {
        var down = new QubitData(
            Qubit: 0, T1Us: 100, T2Us: 100, ReadoutError: 0.01,
            SxError: 0.001, PauliXError: 0.001, Operational: false,
            CzNeighbours: new Dictionary<int, double>(),
            RzzNeighbours: new Dictionary<int, double>());
        Assert.Equal(double.NegativeInfinity, IbmCalibration.ScoreQubit(down));
    }

    [Fact]
    public void ScoreQubit_T2CappedAt2T1_DampensAnomalousT2Reads()
    {
        var anomalous = new QubitData(
            Qubit: 0, T1Us: 100, T2Us: 1000, ReadoutError: 0.0,
            SxError: 0.0, PauliXError: 0.0, Operational: true,
            CzNeighbours: new Dictionary<int, double>(),
            RzzNeighbours: new Dictionary<int, double>());
        Assert.Equal(200.0, IbmCalibration.ScoreQubit(anomalous), precision: 6);
    }

    [Fact]
    public void Load_HandlesEmptyGateErrorFields_AsZero()
    {
        var path = WriteTempCsv(@"""Qubit"",""T1 (us)"",""T2 (us)"",""Readout assignment error"",""Prob meas0 prep1"",""Prob meas1 prep0"",""Readout length (ns)"",""ID error"",""Single-qubit gate length (ns)"",""RX error"",""Z-axis rotation (rz) error"",""√x (sx) error"",""Pauli-X error"",""CZ error"",""Gate length (ns)"",""RZZ error"",""MEASURE error"",""Operational""
""0"",""100.0"",""80.0"",""0.01"",""0.01"",""0.01"",""1000"","""",""36"","""","""","""","""",""1:0.005"",""1:68"",""1:0.006"",""0.01"",""Yes""");
        try
        {
            var qubits = IbmCalibration.Load(path);
            Assert.Single(qubits);
            Assert.Equal(0.0, qubits[0].SxError);
            Assert.Equal(0.0, qubits[0].PauliXError);
            Assert.True(qubits[0].Operational);
        }
        finally { File.Delete(path); }
    }

    [Fact]
    public void Load_ParsesNeighbourField_AsDirectedEdgeMap()
    {
        var path = WriteTempCsv(@"""Qubit"",""T1 (us)"",""T2 (us)"",""Readout assignment error"",""Prob meas0 prep1"",""Prob meas1 prep0"",""Readout length (ns)"",""ID error"",""Single-qubit gate length (ns)"",""RX error"",""Z-axis rotation (rz) error"",""√x (sx) error"",""Pauli-X error"",""CZ error"",""Gate length (ns)"",""RZZ error"",""MEASURE error"",""Operational""
""1"",""100"",""100"",""0.01"",""0.01"",""0.01"",""1000"",""0.001"",""36"",""0.001"",""0"",""0.001"",""0.001"",""2:0.000990685953658782;0:0.002133297688672098"",""2:68;0:68"",""2:0.001083;0:0.002170"",""0.0074"",""Yes""");
        try
        {
            var qubits = IbmCalibration.Load(path);
            Assert.Single(qubits);
            var nbrs = qubits[0].CzNeighbours;
            Assert.Equal(2, nbrs.Count);
            Assert.Equal(0.000990685953658782, nbrs[2], precision: 12);
            Assert.Equal(0.002133297688672098, nbrs[0], precision: 12);
        }
        finally { File.Delete(path); }
    }

    [Fact]
    public void BestChain_RejectsInvalidLength()
    {
        var qubits = Marrakesh20260425.Value;
        Assert.Throws<ArgumentOutOfRangeException>(() => IbmCalibration.BestChain(qubits, length: 0));
    }

    [Fact]
    public void SelectBestChain_BridgesToCalibrationChainRecord()
    {
        var qubits = Marrakesh20260425.Value;
        var chain = IbmCalibration.SelectBestChain(qubits, length: 5);
        Assert.Equal(5, chain.Length);
        Assert.Equal(new[] { 1, 2, 3, 4, 5 }, chain.QubitIds);
        Assert.Equal(1246.58, chain.Score, precision: 1);
        Assert.Equal(5, chain.Qubits.Count);
        Assert.Equal(1, chain.Qubits[0].Qubit);
        Assert.Equal(5, chain.Qubits[4].Qubit);
    }

    [Fact]
    public void CalibrationChain_ToChainSystem_BuildsLogicalChainOfMatchingSize()
    {
        var qubits = Marrakesh20260425.Value;
        var calChain = IbmCalibration.SelectBestChain(qubits, length: 5);
        var sys = calChain.ToChainSystem(J: 1.0, gammaZero: 0.05);
        Assert.Equal(5, sys.N);
        Assert.Equal(1.0, sys.J);
        Assert.Equal(0.05, sys.GammaZero);
        Assert.Equal(HamiltonianType.XY, sys.HType);
        Assert.Equal(TopologyKind.Chain, sys.Topology);
        Assert.Equal(4, sys.Bonds.Count);
        Assert.Equal(new Bond(0, 1, 1.0), sys.Bonds[0]);
        Assert.Equal(new Bond(3, 4, 1.0), sys.Bonds[3]);
    }

    [Fact]
    public void CalibrationChain_ToChainSystemWithT2Gamma_PerSiteGammaMatchesT2()
    {
        var qubits = Marrakesh20260425.Value;
        var calChain = IbmCalibration.SelectBestChain(qubits, length: 3);
        var (sys, gammaPerSite) = calChain.ToChainSystemWithT2Gamma(J: 1.0);
        Assert.Equal(3, sys.N);
        Assert.Equal(3, gammaPerSite.Length);
        for (int i = 0; i < 3; i++)
        {
            double expected = 1.0 / (2.0 * calChain.Qubits[i].T2Us);
            Assert.Equal(expected, gammaPerSite[i], precision: 12);
        }
        Assert.Equal(gammaPerSite.Average(), sys.GammaZero, precision: 12);
    }

    private static string WriteTempCsv(string contents)
    {
        string path = Path.Combine(Path.GetTempPath(), $"ibm_calib_test_{Guid.NewGuid():N}.csv");
        File.WriteAllText(path, contents);
        return path;
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
            $"could not locate repository root (MIRROR_THEORY.md + compute/) " +
            $"starting from {AppContext.BaseDirectory}");
    }
}
