using RCPsiSquared.Core.Confirmations;

namespace RCPsiSquared.Core.Calibration;

/// <summary>Descriptive R* band composition of an addressable path.</summary>
public enum RegimeVerdict
{
    AllBelowRStar,
    AllAtOrAboveRStar,
    RStarBandMixed,
    NotAddressable,
}

/// <summary>Stable text labels for <see cref="RegimeVerdict"/>.</summary>
public static class RegimeVerdictLabels
{
    public static string Label(this RegimeVerdict value) => value switch
    {
        RegimeVerdict.AllBelowRStar => "all-below-r-star",
        RegimeVerdict.AllAtOrAboveRStar => "all-at-or-above-r-star",
        RegimeVerdict.RStarBandMixed => "r-star-band-mixed",
        RegimeVerdict.NotAddressable => "not-addressable",
        _ => throw new ArgumentOutOfRangeException(nameof(value)),
    };
}

/// <summary>One path row with measured T1/T2 and its derived R* band.</summary>
public sealed record PerQubitRegime(
    int Qubit,
    double T1Us,
    double T2Us,
    double RParam,
    Regime RStarBand,
    bool Operational);

/// <summary>Snapshot audit that keeps three readings separate: hardware
/// addressability, an empirical path score, and R* band composition.
///
/// <para>The reported roughly 13.5x contrast between earlier runs is a
/// confounded cross-backend, cross-path, and cross-date association. It is not
/// encoded as a causal rule or an action recommendation here. <see cref="Score"/>
/// is likewise the existing empirical <see cref="IbmCalibration.ChainScore"/>
/// heuristic, not a physical invariant.</para></summary>
public sealed record RegimeSummary(
    IReadOnlyList<int> Path,
    IReadOnlyList<PerQubitRegime> Qubits,
    int BelowRStarCount,
    int NearRStarCount,
    int AtOrAboveRStarCount,
    double Score,
    bool AllCzCoupled,
    bool AllOperational)
{
    /// <summary>True when all rows occupy one non-near R* band.</summary>
    public bool HasSingleNonBoundaryRStarBand =>
        BelowRStarCount == Path.Count || AtOrAboveRStarCount == Path.Count;

    /// <summary>True when each row is operational and every consecutive path
    /// pair has a recorded CZ edge.</summary>
    public bool IsAddressable => AllCzCoupled && AllOperational;

    /// <summary>Exact conjunction of the two named descriptive conditions.</summary>
    public bool IsAddressableAndSingleRStarBand =>
        IsAddressable && HasSingleNonBoundaryRStarBand;

    /// <summary>Addressability takes precedence; otherwise report the R* band
    /// composition without assigning a universal phase.</summary>
    public RegimeVerdict RStarBandVerdict =>
        !IsAddressable ? RegimeVerdict.NotAddressable :
        BelowRStarCount == Path.Count ? RegimeVerdict.AllBelowRStar :
        AtOrAboveRStarCount == Path.Count ? RegimeVerdict.AllAtOrAboveRStar :
        RegimeVerdict.RStarBandMixed;

    /// <summary>Build a snapshot summary. No measured calibration row is altered.</summary>
    public static RegimeSummary For(
        IReadOnlyList<QubitData> qubits,
        IReadOnlyList<int> path,
        double epsilon = 0.0)
    {
        if (path.Count == 0)
            throw new ArgumentException("path must not be empty", nameof(path));

        var byId = qubits.ToDictionary(q => q.Qubit);
        var perQubit = new List<PerQubitRegime>(path.Count);
        foreach (int qid in path)
        {
            if (!byId.TryGetValue(qid, out var data))
                throw new ArgumentException($"qubit {qid} not in calibration", nameof(path));

            perQubit.Add(new PerQubitRegime(
                data.Qubit,
                data.T1Us,
                data.T2Us,
                data.RParam,
                QubitRegime.Classify(data.T1Us, data.T2Us, epsilon),
                data.Operational));
        }

        int below = perQubit.Count(x => x.RStarBand == Regime.BelowRStar);
        int near = perQubit.Count(x => x.RStarBand == Regime.NearRStar);
        int atOrAbove = perQubit.Count(x => x.RStarBand == Regime.AtOrAboveRStar);

        bool allCz = true;
        for (int i = 0; i < path.Count - 1; i++)
        {
            int a = path[i];
            int b = path[i + 1];
            if (!byId[a].CzNeighbours.ContainsKey(b) && !byId[b].CzNeighbours.ContainsKey(a))
            {
                allCz = false;
                break;
            }
        }

        return new RegimeSummary(
            path,
            perQubit,
            below,
            near,
            atOrAbove,
            IbmCalibration.ChainScore(byId, path),
            allCz,
            perQubit.All(x => x.Operational));
    }

    /// <summary>One-line descriptive summary for logs.</summary>
    public string ToHeadline()
    {
        string addressability = IsAddressable ? "addressable" : "not-addressable";
        string path = "[" + string.Join(", ", Path) + "]";
        return $"path {path} | {RStarBandVerdict.Label()} "
             + $"({BelowRStarCount} below, {NearRStarCount} near, "
             + $"{AtOrAboveRStarCount} at-or-above) | empirical score {Score:F2} "
             + $"| {addressability}";
    }

    /// <summary>Existing confirmation rows whose path matches exactly. This is
    /// a read-only lookup and does not turn the summary into a run decision.</summary>
    public IEnumerable<Confirmation> RelatedConfirmations(string? machine = null) =>
        machine == null
            ? ConfirmationsRegistry.ByPath(Path)
            : ConfirmationsRegistry.ByMachineAndPath(machine, Path);
}
