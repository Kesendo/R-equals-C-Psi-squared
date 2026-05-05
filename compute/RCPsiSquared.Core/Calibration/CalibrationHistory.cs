using System.Globalization;

namespace RCPsiSquared.Core.Calibration;

/// <summary>One day's calibration entry for a single qubit. Companion to
/// <see cref="QubitData"/> (which is the full snapshot record); this
/// stripped-down form holds only what's needed for multi-day lifecycle
/// analysis (date + T1/T2 + derived r and regime).</summary>
public sealed record CalibrationDay(string Date, double T1Us, double T2Us)
{
    public double RParam => QubitRegime.RParam(T1Us, T2Us);
    public Regime Regime => QubitRegime.Classify(T1Us, T2Us);
}

/// <summary>A single qubit's daily-calibration time series, in date order.
/// All derived statistics (mean, std, walk rate, archetype) are pure
/// functions of <see cref="Days"/>, computed on demand by
/// <see cref="QubitLifecycle"/>.</summary>
public sealed record QubitTimeline(int Qubit, IReadOnlyList<CalibrationDay> Days);

/// <summary>Loader for IBM-style daily-calibration history CSVs (the format
/// produced by <c>data/ibm_history/ibm_history_analysis.py --mode collect</c>).
/// Parses date + qubit + T1/T2 columns, ignores the derived columns the Python
/// pipeline pre-computes (we re-derive r and regime from <see cref="QubitRegime"/>
/// instead of trusting the file).
///
/// <para>Used by <see cref="LifecycleSummary"/> for path-quality auditing across
/// a multi-day window. Calibrated empirical anchors:
/// <c>data/ibm_history/results/ibm_marrakesh_history.csv</c> (91 days,
/// 156 qubits, 2026-02-04 to 2026-05-05) and
/// <c>data/ibm_history/ibm_torino_history.csv</c> (180 days, 133 qubits,
/// 2025-08-14 to 2026-02-10, retired Eagle r3 backend).</para></summary>
public static class CalibrationHistory
{
    /// <summary>Parse a daily-calibration CSV into per-qubit timelines, sorted
    /// by date. Skips rows with non-positive T1 or T2 (non-operational days).</summary>
    public static IReadOnlyDictionary<int, QubitTimeline> Load(string csvPath)
    {
        using var reader = new StreamReader(csvPath);
        string? header = reader.ReadLine() ?? throw new InvalidDataException("empty CSV");
        string[] headers = header.Split(',');
        int iDate = Array.IndexOf(headers, "date");
        int iQubit = Array.IndexOf(headers, "qubit");
        int iT1 = Array.IndexOf(headers, "T1_us");
        int iT2 = Array.IndexOf(headers, "T2_us");
        if (iDate < 0 || iQubit < 0 || iT1 < 0 || iT2 < 0)
            throw new InvalidDataException(
                $"missing required columns (date, qubit, T1_us, T2_us) in {csvPath}");

        var byQubit = new Dictionary<int, List<CalibrationDay>>();
        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            if (string.IsNullOrWhiteSpace(line)) continue;
            string[] f = line.Split(',');
            if (f.Length <= Math.Max(Math.Max(iDate, iQubit), Math.Max(iT1, iT2))) continue;
            if (!int.TryParse(f[iQubit], NumberStyles.Integer, CultureInfo.InvariantCulture,
                              out int qid)) continue;
            if (!double.TryParse(f[iT1], NumberStyles.Any, CultureInfo.InvariantCulture,
                                 out double t1)) continue;
            if (!double.TryParse(f[iT2], NumberStyles.Any, CultureInfo.InvariantCulture,
                                 out double t2)) continue;
            if (t1 <= 0 || t2 <= 0) continue;

            if (!byQubit.TryGetValue(qid, out var list))
            {
                list = new List<CalibrationDay>();
                byQubit[qid] = list;
            }
            list.Add(new CalibrationDay(f[iDate], t1, t2));
        }

        var result = new Dictionary<int, QubitTimeline>(byQubit.Count);
        foreach (var (qid, days) in byQubit)
        {
            days.Sort((a, b) => string.CompareOrdinal(a.Date, b.Date));
            result[qid] = new QubitTimeline(qid, days);
        }
        return result;
    }
}
