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
/// <para>The observation behind the band composition, read through the F88b-Lens
/// Π²-odd-memory of the truly (XX + YY) category, every hardware run per path with its
/// runner script. soft_break: all-at-or-above-R* Marrakesh [48, 49, 50] 0.0013 (2026-04-26),
/// all-at-or-above-R* Marrakesh [4, 5, 6] 0.0017 (2026-04-30), all-below-R* Kingston
/// [43, 56, 63] 0.0022 (2026-05-05, job d7sqjpiudops73976960, Confirmation
/// <c>regime_uniformity_kingston_uniform_quantum</c>). framework_snapshots, Marrakesh,
/// 2026-04-26: [48, 49, 50] 0.0116, band-mixed [0, 1, 2] 0.0190 eleven minutes later and
/// 0.0297 nine hours after that. Producer: <c>simulations/f88b_lens_ibm_kingston_uniform_quantum.py</c>
/// (the lens) over those run files.</para>
///
/// <para>What the record carries. Within one script, minutes apart, mixed against uniform
/// is 1.64×; the uniform path [48, 49, 50] reads 8.6× apart between the two scripts on
/// one day, and the mixed path 1.57× apart between two runs of one script. The 22.1×
/// (0.0297 against 0.0013) and Kingston's 13.3× against 0.0297 each pair a soft_break run
/// with a framework_snapshots run, so they measure the scripts far more than the paths.
/// What is left for the band is a 1.64× inside run-to-run spread of the same size, and
/// nothing to carry a mechanism: R* is a purity-proxy threshold, not a physical phase
/// boundary. The nearest calibration (2026-04-25) singles out no culprit either: the mixed
/// path is worse on T2 (Q0 at 41 μs, which is what puts it below R*), not worse on CZ, and
/// not uniformly worse on readout. The band composition is therefore reported, not encoded
/// as a causal rule or an action recommendation. <see cref="Score"/>
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
