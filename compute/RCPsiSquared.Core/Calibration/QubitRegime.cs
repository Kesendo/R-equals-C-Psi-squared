namespace RCPsiSquared.Core.Calibration;

/// <summary>Regime label for a qubit relative to the CΨ = ¼ fold-catastrophe
/// boundary, expressed at the calibration-data level via r = T2 / (2·T1).</summary>
public enum Regime
{
    /// <summary>r &lt; R*: T2 dephasing dominates, CΨ_min crosses ¼ → bridge open
    /// in the framework's d=0 ↔ d=2 sense.</summary>
    QuantumSide,
    /// <summary>|r − R*| ≤ ε (only returned when an explicit ε is passed):
    /// qubit is at the fold-catastrophe boundary, neither regime is committed.</summary>
    Boundary,
    /// <summary>r &gt; R*: T1 relaxation dominates relative to T2, CΨ_min stays
    /// above ¼ → bridge stays closed in the framework's classical-side sense.</summary>
    ClassicalSide,
}

/// <summary>Per-qubit regime classifier: bridges raw hardware calibration
/// (T1, T2) to the framework's CΨ = ¼ fold-catastrophe boundary.
///
/// <para>The framework's polynomial trunk d² − 2d = 0 selects d=2 as the qubit
/// dimension; the same polynomial gives ¼ as the discriminant boundary of the
/// purity-coherence form CΨ. Solving "where does the time-minimum of CΨ touch ¼"
/// for a single transmon under T1 amplitude damping + T2 pure dephasing yields the
/// critical ratio <see cref="R_STAR"/> ≈ 0.213 in r = T2 / (2·T1). Below R*, CΨ_min
/// crosses ¼ (bridge open); above R*, CΨ_min stays &gt; ¼ (bridge closed).</para>
///
/// <para>r is what IBM publishes daily for every qubit on every backend. This
/// class exposes it as a typed primitive instead of letting every analysis
/// script re-derive 0.212755 by hand. See
/// <c>docs/BOTH_SIDES_VISIBLE.md</c> for the original 6-month Torino crossing
/// analysis; see <c>simulations/_qubit_biography.py</c> for the lifecycle
/// archetype derived from <see cref="Classify"/> applied to a daily history.</para>
/// </summary>
public static class QubitRegime
{
    /// <summary>Critical T2/(2·T1) ratio where CΨ_min(t) just touches ¼. The
    /// numerical value comes from minimising the single-transmon purity
    /// 1 − ½e^(−t/T1) + ½e^(−2t/T2) under T2 ≤ 2T1 and solving CΨ_min = ¼;
    /// it is hardware-independent (a property of the d² − 2d = 0 polynomial).</summary>
    public const double R_STAR = 0.212755;

    /// <summary>r = T2 / (2·T1), the single ratio that controls CΨ_min on a
    /// single transmon. Returns <see cref="double.PositiveInfinity"/> when
    /// T1 ≤ 0 so non-operational qubits classify as
    /// <see cref="Regime.ClassicalSide"/> (the safe default; invalid hardware
    /// must not look like a clean quantum-side reading). Callers needing strict
    /// validation should check operationality before classification.</summary>
    public static double RParam(double t1Us, double t2Us) =>
        t1Us > 0 ? t2Us / (2.0 * t1Us) : double.PositiveInfinity;

    /// <summary>Classify a qubit by its T1/T2. Default ε = 0 returns binary
    /// QuantumSide/ClassicalSide; pass ε &gt; 0 to surface the
    /// <see cref="Regime.Boundary"/> band where |r − R*| ≤ ε (useful for
    /// stability analysis where small drift can flip the classification).</summary>
    public static Regime Classify(double t1Us, double t2Us, double epsilon = 0.0)
    {
        double r = RParam(t1Us, t2Us);
        if (epsilon > 0 && Math.Abs(r - R_STAR) <= epsilon) return Regime.Boundary;
        return r < R_STAR ? Regime.QuantumSide : Regime.ClassicalSide;
    }

    /// <summary>Convenience: true iff this qubit's CΨ_min crosses ¼ (r &lt; R*).</summary>
    public static bool IsQuantumSide(double t1Us, double t2Us) =>
        RParam(t1Us, t2Us) < R_STAR;

    /// <summary>Regime composition of a path on a calibration snapshot:
    /// (Quantum, Boundary, Classical) counts. Useful for F88-Lens experimental
    /// design where regime-mixed paths read different physics from regime-uniform
    /// paths (cf. the 2026-04-26 framework_snapshots [0, 1, 2] mixed regime
    /// truly-baseline 0.030 vs soft_break [48, 49, 50] uniform-classical 0.0013).</summary>
    public static (int Quantum, int Boundary, int Classical) PathComposition(
        IReadOnlyList<QubitData> qubits, IReadOnlyList<int> path, double epsilon = 0.0)
    {
        var byId = qubits.ToDictionary(q => q.Qubit);
        int q = 0, b = 0, c = 0;
        foreach (int qid in path)
        {
            if (!byId.TryGetValue(qid, out var data))
                throw new ArgumentException($"qubit {qid} not in calibration", nameof(path));
            switch (Classify(data.T1Us, data.T2Us, epsilon))
            {
                case Regime.QuantumSide: q++; break;
                case Regime.Boundary: b++; break;
                case Regime.ClassicalSide: c++; break;
            }
        }
        return (q, b, c);
    }
}
