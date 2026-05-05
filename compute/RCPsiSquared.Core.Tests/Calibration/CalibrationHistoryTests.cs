using RCPsiSquared.Core.Calibration;

namespace RCPsiSquared.Core.Tests.Calibration;

/// <summary>Tests for <see cref="CalibrationHistory.Load"/>, anchored to the
/// 91-day Marrakesh dataset committed at
/// <c>data/ibm_history/results/ibm_marrakesh_history.csv</c> (commit 787854f,
/// 14196 records, 156 qubits × 91 days, 2026-02-04 to 2026-05-05).</summary>
public class CalibrationHistoryTests
{
    private static readonly Lazy<IReadOnlyDictionary<int, QubitTimeline>> Marrakesh91d = new(() =>
        CalibrationHistory.Load(Path.Combine(FindRepoRoot(),
            "data", "ibm_history", "results", "ibm_marrakesh_history.csv")));

    [Fact]
    public void Load_Marrakesh91d_HasAll156Qubits()
    {
        var h = Marrakesh91d.Value;
        Assert.Equal(156, h.Count);
        for (int q = 0; q < 156; q++)
            Assert.True(h.ContainsKey(q), $"qubit {q} missing from history");
    }

    [Fact]
    public void Load_Marrakesh91d_HasNinetyOneDaysPerQubit()
    {
        var h = Marrakesh91d.Value;
        foreach (var (qid, timeline) in h)
            Assert.Equal(91, timeline.Days.Count);
    }

    [Fact]
    public void Load_Marrakesh91d_DaysAreSortedChronologically()
    {
        var h = Marrakesh91d.Value;
        var q0 = h[0];
        Assert.Equal("2026-02-04", q0.Days[0].Date);
        Assert.Equal("2026-05-05", q0.Days[^1].Date);
        for (int i = 1; i < q0.Days.Count; i++)
            Assert.True(string.CompareOrdinal(q0.Days[i].Date, q0.Days[i - 1].Date) >= 0,
                $"day {i} ({q0.Days[i].Date}) precedes day {i - 1} ({q0.Days[i - 1].Date})");
    }

    [Fact]
    public void CalibrationDay_DerivedRegime_AgreesWithQubitRegime()
    {
        var day = new CalibrationDay("2026-02-04", T1Us: 100, T2Us: 40);
        Assert.Equal(QubitRegime.RParam(100, 40), day.RParam, precision: 6);
        Assert.Equal(QubitRegime.Classify(100, 40), day.Regime);
    }

    [Fact]
    public void Load_RejectsMissingColumns()
    {
        string path = Path.Combine(Path.GetTempPath(), $"hist_test_{Guid.NewGuid():N}.csv");
        File.WriteAllText(path, "date,qubit\n2026-01-01,0\n");
        try
        {
            Assert.Throws<InvalidDataException>(() => CalibrationHistory.Load(path));
        }
        finally { File.Delete(path); }
    }

    [Fact]
    public void Load_SkipsRowsWithNonPositiveT1OrT2()
    {
        string path = Path.Combine(Path.GetTempPath(), $"hist_test_{Guid.NewGuid():N}.csv");
        File.WriteAllText(path,
            "date,qubit,T1_us,T2_us\n" +
            "2026-01-01,0,100.0,80.0\n" +
            "2026-01-02,0,0.0,80.0\n" +       // skipped (T1 <= 0)
            "2026-01-03,0,100.0,-1.0\n" +     // skipped (T2 <= 0)
            "2026-01-04,0,150.0,60.0\n");
        try
        {
            var h = CalibrationHistory.Load(path);
            Assert.Single(h);
            Assert.Equal(2, h[0].Days.Count);
            Assert.Equal("2026-01-01", h[0].Days[0].Date);
            Assert.Equal("2026-01-04", h[0].Days[1].Date);
        }
        finally { File.Delete(path); }
    }

    private static string FindRepoRoot()
    {
        var dir = new DirectoryInfo(AppContext.BaseDirectory);
        while (dir != null)
        {
            if (Directory.Exists(Path.Combine(dir.FullName, "ClaudeTasks"))
             && File.Exists(Path.Combine(dir.FullName, "MIRROR_THEORY.md")))
                return dir.FullName;
            dir = dir.Parent;
        }
        throw new InvalidOperationException(
            $"could not locate repository root starting from {AppContext.BaseDirectory}");
    }
}
