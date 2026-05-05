using RCPsiSquared.Core.Confirmations;

namespace RCPsiSquared.Core.Calibration;

/// <summary>Pre-submit audit label for a candidate qubit-path on hardware.
/// <see cref="NotAddressable"/> always wins (a path that can't run can't be
/// classified by regime); the three regime labels apply to addressable paths.</summary>
public enum RegimeVerdict
{
    UniformQuantum,
    UniformClassical,
    RegimeMixed,
    NotAddressable,
}

/// <summary><see cref="RegimeVerdict"/> ↔ kebab-case label conversion for logs
/// and headlines. Kept stable so callers can serialise to text, then back via
/// <see cref="Parse"/>.</summary>
public static class RegimeVerdictLabels
{
    public static string Label(this RegimeVerdict v) => v switch
    {
        RegimeVerdict.UniformQuantum => "uniform-quantum",
        RegimeVerdict.UniformClassical => "uniform-classical",
        RegimeVerdict.RegimeMixed => "regime-mixed",
        RegimeVerdict.NotAddressable => "not-addressable",
        _ => throw new ArgumentOutOfRangeException(nameof(v)),
    };
}

/// <summary>One row of a <see cref="RegimeSummary"/>: per-qubit calibration plus
/// the regime classification derived from <see cref="QubitRegime"/>.</summary>
public sealed record PerQubitRegime(
    int Qubit,
    double T1Us,
    double T2Us,
    double RParam,
    Regime Regime,
    bool Operational);

/// <summary>Pre-experiment audit of a candidate qubit-path on an IBM-style
/// backend: composes <see cref="IbmCalibration"/>, <see cref="QubitRegime"/>,
/// and the CZ-coupling graph into a single workflow object.
///
/// <para>Answers the recurring pre-submit question: "is this path addressable
/// on hardware (operational + CZ-coupled), and what regime mix does it sample?"
/// The 2026-04-26 framework_snapshots [0, 1, 2] vs soft_break [48, 49, 50]
/// comparison is the motivating case; the former is regime-mixed (Q0 quantum,
/// Q1+Q2 classical) and reads truly-baseline 0.030, the latter is
/// uniform-classical and reads 0.0013. <see cref="Verdict"/> labels exactly
/// this distinction so it can be flagged before any QPU minutes are spent.</para>
///
/// <para>Uses <see cref="QubitRegime.PathComposition"/> for the regime counts
/// and <see cref="IbmCalibration.ChainScore"/> for the path-quality scalar.</para></summary>
public sealed record RegimeSummary(
    IReadOnlyList<int> Path,
    IReadOnlyList<PerQubitRegime> Qubits,
    int QuantumCount,
    int BoundaryCount,
    int ClassicalCount,
    double Score,
    bool AllCzCoupled,
    bool AllOperational)
{
    /// <summary>True iff every qubit on the path lives on the same side of R*
    /// (all quantum or all classical). Boundary qubits do not count as either,
    /// so a path with a boundary qubit is not regime-uniform.</summary>
    public bool IsRegimeUniform =>
        QuantumCount == Path.Count || ClassicalCount == Path.Count;

    /// <summary>True iff the path can actually be submitted: every consecutive
    /// pair is CZ-coupled and every qubit is operational on this calibration.</summary>
    public bool IsAddressable => AllCzCoupled && AllOperational;

    /// <summary>Audit label for pre-submit gating.
    /// <see cref="RegimeVerdict.NotAddressable"/> overrides the regime labels
    /// (a non-addressable path can't run, so its regime composition is moot).
    /// For the addressable case, the regime-uniformity check picks
    /// <see cref="RegimeVerdict.UniformQuantum"/> /
    /// <see cref="RegimeVerdict.UniformClassical"/> /
    /// <see cref="RegimeVerdict.RegimeMixed"/>. Use <see cref="RegimeVerdictLabels.Label"/>
    /// for kebab-case strings in logs.</summary>
    public RegimeVerdict Verdict =>
        !IsAddressable ? RegimeVerdict.NotAddressable :
        QuantumCount == Path.Count ? RegimeVerdict.UniformQuantum :
        ClassicalCount == Path.Count ? RegimeVerdict.UniformClassical :
        RegimeVerdict.RegimeMixed;

    /// <summary>Recommended pre-submit gate for F88-Lens / F87-trichotomy
    /// experiments where the truly-baseline matters. True iff
    /// <see cref="IsAddressable"/> AND <see cref="IsRegimeUniform"/>.
    ///
    /// <para>Rule grounded in the 2026-05-05 Kingston [43, 56, 63] hardware
    /// run (Confirmation <c>regime_uniformity_kingston_uniform_quantum</c>):
    /// regime-uniform paths give clean truly-baselines (uniform-classical
    /// 0.0013, uniform-quantum 0.0022, ratio 1.69×) while regime-mixed
    /// paths add ~13.5× extra noise (Marrakesh [0, 1, 2] truly = 0.0297
    /// vs uniform-classical baseline). Mixed paths still produce signal
    /// but with reduced trichotomy contrast; uniform paths are the safer
    /// default for clean readings.</para></summary>
    public bool IsRecommendedForSubmit => IsAddressable && IsRegimeUniform;

    /// <summary>Build the summary for <paramref name="path"/> against an
    /// already-loaded calibration. Throws <see cref="ArgumentException"/> if
    /// any qubit on the path is missing from the calibration.</summary>
    public static RegimeSummary For(
        IReadOnlyList<QubitData> qubits,
        IReadOnlyList<int> path,
        double epsilon = 0.0)
    {
        if (path.Count == 0) throw new ArgumentException("path must not be empty", nameof(path));
        var byId = qubits.ToDictionary(q => q.Qubit);

        var perQubit = new List<PerQubitRegime>(path.Count);
        foreach (int qid in path)
        {
            if (!byId.TryGetValue(qid, out var d))
                throw new ArgumentException($"qubit {qid} not in calibration", nameof(path));
            perQubit.Add(new PerQubitRegime(
                Qubit: d.Qubit,
                T1Us: d.T1Us,
                T2Us: d.T2Us,
                RParam: d.RParam,
                Regime: QubitRegime.Classify(d.T1Us, d.T2Us, epsilon),
                Operational: d.Operational));
        }

        int q = perQubit.Count(x => x.Regime == Regime.QuantumSide);
        int b = perQubit.Count(x => x.Regime == Regime.Boundary);
        int c = perQubit.Count(x => x.Regime == Regime.ClassicalSide);

        bool allCz = true;
        for (int i = 0; i < path.Count - 1; i++)
        {
            int a = path[i], bb = path[i + 1];
            bool bonded = byId[a].CzNeighbours.ContainsKey(bb)
                       || byId[bb].CzNeighbours.ContainsKey(a);
            if (!bonded) { allCz = false; break; }
        }

        bool allOp = perQubit.All(x => x.Operational);
        double score = IbmCalibration.ChainScore(byId, path);

        return new RegimeSummary(
            Path: path,
            Qubits: perQubit,
            QuantumCount: q,
            BoundaryCount: b,
            ClassicalCount: c,
            Score: score,
            AllCzCoupled: allCz,
            AllOperational: allOp);
    }

    /// <summary>One-line human-readable summary for logs, e.g.
    /// <c>"path [0, 1, 2] | regime-mixed (1q, 0b, 2c) | score 597.27 | addressable"</c>.</summary>
    public string ToHeadline()
    {
        string addr = IsAddressable ? "addressable" : "not-addressable";
        string pathStr = "[" + string.Join(", ", Path) + "]";
        return $"path {pathStr} | {Verdict.Label()} ({QuantumCount}q, {BoundaryCount}b, {ClassicalCount}c) "
             + $"| score {Score:F2} | {addr}";
    }

    /// <summary>Confirmations from <see cref="ConfirmationsRegistry"/> whose
    /// <see cref="Confirmation.QubitPath"/> matches this summary's path exactly.
    /// Pass <paramref name="machine"/> to restrict to one backend (e.g.
    /// <c>"ibm_marrakesh"</c>).</summary>
    public IEnumerable<Confirmation> RelatedConfirmations(string? machine = null) =>
        machine == null
            ? ConfirmationsRegistry.ByPath(Path)
            : ConfirmationsRegistry.ByMachineAndPath(machine, Path);
}
