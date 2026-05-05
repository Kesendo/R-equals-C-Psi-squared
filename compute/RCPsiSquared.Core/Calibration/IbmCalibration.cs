using System.Globalization;

namespace RCPsiSquared.Core.Calibration;

/// <summary>IBM Heron r2 calibration CSV loader and qubit/chain quality scoring.
/// C# counterpart to <c>simulations/ibm_calibration.py</c>; both are validated
/// against the 2026-04-26 framework_snapshots / soft_break runs (path
/// [48, 49, 50] scores higher than [0, 1, 2], matching the empirical 23×
/// truly-baseline ratio).
///
/// <para>Score function: T2-dominated, with multiplicative penalties for
/// single-qubit-gate-error and readout error. Rationale: F88-Lens experiments
/// at t ~ 1/J ~ 1us are dephasing-limited, so T2 is the dominant per-qubit
/// noise driver of the truly Π²-odd-baseline.</para>
///
/// <para>Chain finder: DFS over the CZ-coupled graph; chain score = sum of
/// per-qubit scores plus −1000·CZ-error per bond (weaker errors raise the
/// score). Returns the global optimum over all paths of the requested length
/// without revisits.</para></summary>
public static class IbmCalibration
{
    /// <summary>Parse an IBM calibration CSV into <see cref="QubitData"/> records.
    /// Handles empty fields (non-operational qubits sometimes have missing gate-
    /// error entries) by defaulting to 0.</summary>
    public static IReadOnlyList<QubitData> Load(string csvPath)
    {
        using var reader = new StreamReader(csvPath);
        string? headerLine = reader.ReadLine() ?? throw new InvalidDataException("empty CSV");
        string[] headers = ParseCsvLine(headerLine);
        var col = new Dictionary<string, int>();
        for (int i = 0; i < headers.Length; i++) col[headers[i]] = i;

        int idx(string name) =>
            col.TryGetValue(name, out int i)
                ? i
                : throw new InvalidDataException($"missing column '{name}' in CSV");

        int iQ = idx("Qubit");
        int iT1 = idx("T1 (us)");
        int iT2 = idx("T2 (us)");
        int iReadout = idx("Readout assignment error");
        int iSx = idx("√x (sx) error");
        int iPx = idx("Pauli-X error");
        int iCz = idx("CZ error");
        int iRzz = idx("RZZ error");
        int iOp = idx("Operational");

        var qubits = new List<QubitData>();
        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            if (string.IsNullOrWhiteSpace(line)) continue;
            string[] f = ParseCsvLine(line);
            qubits.Add(new QubitData(
                Qubit: int.Parse(f[iQ], CultureInfo.InvariantCulture),
                T1Us: SafeFloat(f[iT1]),
                T2Us: SafeFloat(f[iT2]),
                ReadoutError: SafeFloat(f[iReadout]),
                SxError: SafeFloat(f[iSx]),
                PauliXError: SafeFloat(f[iPx]),
                Operational: f[iOp].Trim().Equals("Yes", StringComparison.OrdinalIgnoreCase),
                CzNeighbours: ParseNeighbourField(f[iCz]),
                RzzNeighbours: ParseNeighbourField(f[iRzz])));
        }
        return qubits;
    }

    /// <summary>Composite quality score; higher is better. Rewards high T1/T2,
    /// penalises gate and readout errors. Non-operational qubits return
    /// <see cref="double.NegativeInfinity"/>.</summary>
    public static double ScoreQubit(QubitData q)
    {
        if (!q.Operational) return double.NegativeInfinity;
        double coherence = Math.Min(q.T2Us, 2.0 * q.T1Us);
        double gateQuality = Math.Pow(1.0 - q.SxError, 4) * Math.Pow(1.0 - q.PauliXError, 4);
        double readoutQuality = 1.0 - q.ReadoutError;
        return coherence * gateQuality * readoutQuality;
    }

    /// <summary>Top-k qubits by <see cref="ScoreQubit"/>.</summary>
    public static IReadOnlyList<QubitData> BestQubits(IReadOnlyList<QubitData> qubits, int k) =>
        qubits.OrderByDescending(ScoreQubit).Take(k).ToList();

    /// <summary>Find the best contiguous CZ-coupled path of <paramref name="length"/>
    /// qubits via DFS. Score = qubit-sum + Σ −1000·CZ-error along bonds.
    /// Returns (path, score). Path is a list of qubit indices in chain order.</summary>
    public static (IReadOnlyList<int> Path, double Score) BestChain(
        IReadOnlyList<QubitData> qubits, int length)
    {
        if (length < 1) throw new ArgumentOutOfRangeException(nameof(length), "length must be ≥ 1");
        var byId = qubits.ToDictionary(q => q.Qubit);

        List<int> bestPath = new();
        double bestScore = double.NegativeInfinity;

        void Dfs(List<int> path, HashSet<int> visited, double partial)
        {
            if (path.Count == length)
            {
                if (partial > bestScore)
                {
                    bestScore = partial;
                    bestPath = new List<int>(path);
                }
                return;
            }
            var last = byId[path[^1]];
            foreach (var (nbrId, czErr) in last.CzNeighbours)
            {
                if (visited.Contains(nbrId) || !byId.TryGetValue(nbrId, out var nbr)) continue;
                if (!nbr.Operational) continue;
                double bondScore = -1000.0 * czErr;
                double qubitScore = ScoreQubit(nbr);
                visited.Add(nbrId);
                path.Add(nbrId);
                Dfs(path, visited, partial + qubitScore + bondScore);
                path.RemoveAt(path.Count - 1);
                visited.Remove(nbrId);
            }
        }

        foreach (var q in qubits.Where(x => x.Operational))
        {
            var path = new List<int> { q.Qubit };
            var visited = new HashSet<int> { q.Qubit };
            Dfs(path, visited, ScoreQubit(q));
        }
        return (bestPath, bestScore);
    }

    /// <summary>Score a specific chain (for comparison against
    /// <see cref="BestChain"/>). Same scoring rule.</summary>
    public static double ChainScore(IReadOnlyList<QubitData> qubits, IReadOnlyList<int> path)
    {
        var byId = qubits.ToDictionary(q => q.Qubit);
        double score = path.Sum(p => ScoreQubit(byId[p]));
        for (int i = 0; i < path.Count - 1; i++)
        {
            int a = path[i], b = path[i + 1];
            double czErr = byId[a].CzNeighbours.GetValueOrDefault(b,
                byId[b].CzNeighbours.GetValueOrDefault(a, 1.0));
            score += -1000.0 * czErr;
        }
        return score;
    }

    /// <summary>One-call bridge: load calibration CSV, find the best CZ-coupled chain
    /// of <paramref name="length"/>, and return a <see cref="CalibrationChain"/> ready
    /// to lower into a <see cref="ChainSystems.ChainSystem"/>. Equivalent to
    /// <c>SelectBestChain(Load(csvPath), length)</c>.</summary>
    public static CalibrationChain SelectBestChain(string csvPath, int length) =>
        SelectBestChain(Load(csvPath), length);

    /// <summary>Bridge an already-loaded calibration to a <see cref="CalibrationChain"/>.
    /// Returns the global optimum over all CZ-coupled paths of the requested length;
    /// throws if no operational chain of that length exists.</summary>
    public static CalibrationChain SelectBestChain(IReadOnlyList<QubitData> qubits, int length)
    {
        var (path, score) = BestChain(qubits, length);
        if (path.Count != length)
            throw new InvalidOperationException(
                $"no operational CZ-coupled chain of length {length} found in calibration");
        var byId = qubits.ToDictionary(q => q.Qubit);
        var ordered = path.Select(id => byId[id]).ToArray();
        return new CalibrationChain(QubitIds: path, Score: score, Qubits: ordered);
    }

    // ──────────── private helpers ────────────

    private static double SafeFloat(string value) =>
        double.TryParse(value, NumberStyles.Any, CultureInfo.InvariantCulture, out double v) ? v : 0.0;

    private static IReadOnlyDictionary<int, double> ParseNeighbourField(string field)
    {
        var result = new Dictionary<int, double>();
        if (string.IsNullOrEmpty(field)) return result;
        foreach (string entry in field.Split(';'))
        {
            int colon = entry.IndexOf(':');
            if (colon < 0) continue;
            if (int.TryParse(entry[..colon], NumberStyles.Integer, CultureInfo.InvariantCulture, out int nbr) &&
                double.TryParse(entry[(colon + 1)..], NumberStyles.Any, CultureInfo.InvariantCulture, out double err))
            {
                result[nbr] = err;
            }
        }
        return result;
    }

    private static string[] ParseCsvLine(string line)
    {
        var fields = new List<string>();
        var current = new System.Text.StringBuilder();
        bool inQuotes = false;
        foreach (char c in line)
        {
            if (c == '"') inQuotes = !inQuotes;
            else if (c == ',' && !inQuotes)
            {
                fields.Add(current.ToString());
                current.Clear();
            }
            else current.Append(c);
        }
        fields.Add(current.ToString());
        return fields.ToArray();
    }
}
