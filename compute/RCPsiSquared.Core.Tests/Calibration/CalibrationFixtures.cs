using RCPsiSquared.Core.Calibration;

namespace RCPsiSquared.Core.Tests.Calibration;

/// <summary>Shared test infrastructure for the Calibration test files:
/// repo-root finder, lazy data loaders pinned to in-repo CSVs, and synthetic
/// timeline factories. Single point of edit if either anchor moves or a third
/// archetype-shape needs covering.</summary>
internal static class CalibrationFixtures
{
    /// <summary>Walks up from the test binary directory until it finds the
    /// repo root (identified by MIRROR_THEORY.md plus the compute/ directory).</summary>
    public static string FindRepoRoot()
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

    /// <summary>Loaded once per test session: the 2026-04-25 Marrakesh
    /// calibration snapshot from <c>data/ibm_calibration_snapshots/</c>.</summary>
    public static readonly Lazy<IReadOnlyList<QubitData>> Marrakesh20260425 = new(() =>
        IbmCalibration.Load(Path.Combine(FindRepoRoot(),
            "data", "ibm_calibration_snapshots",
            "ibm_marrakesh_calibrations_2026-04-25T11_28_00Z.csv")));

    /// <summary>Loaded once per test session: the 91-day Marrakesh history
    /// from <c>data/ibm_history/results/ibm_marrakesh_history.csv</c>
    /// (2026-02-04 to 2026-05-05, 156 qubits × 91 days).</summary>
    public static readonly Lazy<IReadOnlyDictionary<int, QubitTimeline>> Marrakesh91d = new(() =>
        CalibrationHistory.Load(Path.Combine(FindRepoRoot(),
            "data", "ibm_history", "results", "ibm_marrakesh_history.csv")));

    /// <summary>Synthetic timeline of <paramref name="days"/> calibration
    /// entries with constant T1/T2; used to lock specific archetypes.</summary>
    public static QubitTimeline StableTimeline(int qid, int days, double t1Us, double t2Us)
    {
        var list = new List<CalibrationDay>(days);
        for (int i = 0; i < days; i++)
            list.Add(new CalibrationDay($"2026-01-{i + 1:D2}", t1Us, t2Us));
        return new QubitTimeline(qid, list);
    }

    /// <summary>Synthetic timeline alternating T2 between 20 and 80 (with
    /// fixed T1) day by day, producing walk rate ≈ 1.0 → twitch archetype.</summary>
    public static QubitTimeline AlternatingTimeline(int qid, int days, double t1Us = 100.0)
    {
        var list = new List<CalibrationDay>(days);
        for (int i = 0; i < days; i++)
        {
            double t2 = (i % 2 == 0) ? 20.0 : 80.0;
            list.Add(new CalibrationDay($"2026-01-{i + 1:D2}", t1Us, t2));
        }
        return new QubitTimeline(qid, list);
    }
}
