using System.Globalization;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live reading-power lab: how well a basis readout resolves a bond defect as a
/// function of the resonator quality factor Q = J/γ. It runs one sweep over Q (by varying γ at
/// fixed coupling J) and reports three things a stranger can read off the output:
///
/// <para>(1) the RESOLUTION LAW — on the population (Z) readout the classical Fisher information of
/// the defect grows linearly with Q: more coherent cycles per lifetime means more resolving power
/// per dose. A least-squares line through the origin, FI ≈ c·Q, is fit and the worst relative
/// residual reported (PASS under 25%).</para>
///
/// <para>(2) NO EXCEPTIONAL-POINT PEAK — the Fisher information is monotone increasing in Q in
/// every readout basis tested (Z, X, Y). The exceptional point Q = 1 is the WORST reading point,
/// not a resonance: nothing peaks there.</para>
///
/// <para>(3) BASIS ORDERING — coherence-basis readouts (X, Y) lose resolving power toward low Q far
/// faster than the population readout (Z). The span ratio FI(Q_max)/FI(Q_min) of the X readout
/// exceeds ten times that of Z; near Q = 1 only the population basis still reads.</para>
///
/// <para>The sweep is built once (lazy <c>??=</c>) and every basis is evaluated on the same shared
/// trajectories via <see cref="ReadoutFisher.FiMax"/> / <see cref="ReadoutFisher.DiscriminationMax"/>.
/// Honest scope: a fixed dose window K ∈ (0, 1], a forward δJ difference (O(δJ) bias ≈ 1%), N fixed
/// per instance, the law is per-DOSE. Design:
/// docs/superpowers/specs/2026-06-12-handshake-decoder-reading-grammar-design.md</para></summary>
public sealed class ReadingPowerWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    public int N { get; }

    /// <summary>The dephasing rates swept, giving Q = J/γ = 20 .. 1 at J = 1 (γ ascending ⟹ Q descending).</summary>
    static readonly double[] Gammas = { 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0 };
    const double J = 1.0, DeltaJ = 0.02, KMax = 1.0;
    const int Points = 80;

    /// <summary>The PASS threshold for the Z-readout linearity (max relative residual).</summary>
    public const double LinearityThreshold = 0.25;

    public int SweepCount { get; private set; }

    sealed record SweepRow(double Q, double FiZ, double FiX, double FiY, double DLocZ);
    List<SweepRow>? _rows;

    public ReadingPowerWitness(int n = 4)
    {
        if (n < 3 || n > 5) throw new ArgumentOutOfRangeException(nameof(n),
            "the reading-power sweep is specified for N in 3..5 (location needs an interior bond; cost grows fast above)");
        N = n;
    }

    List<SweepRow> Rows()
    {
        if (_rows != null) return _rows;
        SweepCount++;
        _rows = Gammas.Select(g =>
        {
            var ts = ReadoutFisher.KGrid(g, KMax, Points);
            var clean = ReadoutFisher.Trajectory(N, J, g, null, 0.0, ts);
            var d0    = ReadoutFisher.Trajectory(N, J, g, 0, DeltaJ, ts);
            var dFar  = ReadoutFisher.Trajectory(N, J, g, N - 2, DeltaJ, ts);
            return new SweepRow(
                Q: J / g,
                FiZ: ReadoutFisher.FiMax(clean, d0, DeltaJ, ReadoutBasis.Z),
                FiX: ReadoutFisher.FiMax(clean, d0, DeltaJ, ReadoutBasis.X),
                FiY: ReadoutFisher.FiMax(clean, d0, DeltaJ, ReadoutBasis.Y),
                DLocZ: ReadoutFisher.DiscriminationMax(clean, d0, dFar, ReadoutBasis.Z));
        }).ToList();
        return _rows;
    }

    static double Fi(SweepRow r, ReadoutBasis b) => b switch
    { ReadoutBasis.X => r.FiX, ReadoutBasis.Y => r.FiY, _ => r.FiZ };

    /// <summary>FI(Q) is monotone increasing in Q for the given readout basis (no EP peak).</summary>
    public bool IsMonotoneInQ(ReadoutBasis b)
    {
        var rows = Rows().OrderBy(r => r.Q).ToList();
        for (int i = 1; i < rows.Count; i++)
            if (Fi(rows[i], b) <= Fi(rows[i - 1], b)) return false;
        return true;
    }

    /// <summary>FI(Q_max) / FI(Q_min) for the given basis. Gammas ascend ⟹ rows[0] = Q_max (=20),
    /// rows[^1] = Q_min (=1).</summary>
    public double SpanRatio(ReadoutBasis b)
    { var r = Rows(); return Fi(r[0], b) / Fi(r[^1], b); }

    /// <summary>Fit FI_Z ≈ c·Q through the origin (least squares) and return the worst relative residual.</summary>
    public double ZLinearityRelativeResidual()
    {
        var rows = Rows();
        double c = rows.Sum(r => r.FiZ * r.Q) / rows.Sum(r => r.Q * r.Q);
        return rows.Max(r => Math.Abs(r.FiZ - c * r.Q) / Math.Max(r.FiZ, 1e-12));
    }

    /// <summary>The fitted slope c in FI_Z ≈ c·Q (origin-constrained least squares).</summary>
    public double ZLinearitySlope()
    {
        var rows = Rows();
        return rows.Sum(r => r.FiZ * r.Q) / rows.Sum(r => r.Q * r.Q);
    }

    public string DisplayName =>
        $"ReadingPowerWitness (FI of a bond defect vs Q = J/γ, N = {N}, readout bases Z/X/Y)";

    public string Summary
    {
        get
        {
            double c = ZLinearitySlope();
            double resid = ZLinearityRelativeResidual();
            bool noPeak = IsMonotoneInQ(ReadoutBasis.Z)
                       && IsMonotoneInQ(ReadoutBasis.X)
                       && IsMonotoneInQ(ReadoutBasis.Y);
            return $"resolving power grows with Q = J/γ: FI ≈ c·Q on the Z readout " +
                   $"(c = {c.ToString("0.####", Inv)}, residual {(resid * 100).ToString("0.#", Inv)} %" +
                   $"{(resid < LinearityThreshold ? "" : " ABOVE 25% threshold")}), " +
                   $"{(noPeak ? "no basis peaks at the exceptional point" : "WARNING: a basis peaks at the exceptional point")}, " +
                   $"and coherence readouts fade fastest toward low Q. " +
                   $"Design: docs/superpowers/specs/2026-06-12-handshake-decoder-reading-grammar-design.md";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var rows = Rows();

            // 1. The sweep table over Q, plus a Curve of FI_Z over Q.
            var qAxis = rows.Select(r => r.Q).ToArray();
            var fiZAxis = rows.Select(r => r.FiZ).ToArray();
            var tableLines = rows.Select(r =>
                $"Q = {r.Q.ToString("0.##", Inv),5}: " +
                $"FI_Z = {r.FiZ.ToString("0.####e+00", Inv)}, " +
                $"FI_X = {r.FiX.ToString("0.####e+00", Inv)}, " +
                $"FI_Y = {r.FiY.ToString("0.####e+00", Inv)}, " +
                $"D_location = {r.DLocZ.ToString("0.####e+00", Inv)}");
            yield return new InspectableNode(
                displayName: "the sweep (FI of the defect readout over Q = J/γ)",
                summary: $"{rows.Count} rows, Q from {rows.Max(r => r.Q).ToString("0.##", Inv)} down to " +
                         $"{rows.Min(r => r.Q).ToString("0.##", Inv)}. " + string.Join("  |  ", tableLines),
                payload: new InspectablePayload.Curve("FI_Z vs Q", qAxis, fiZAxis, "Q = J/γ", "Fisher information (Z readout)"));

            // 2. The resolution law FI = c·Q on the Z readout.
            double c = ZLinearitySlope();
            double resid = ZLinearityRelativeResidual();
            bool lawPasses = resid < LinearityThreshold;
            yield return new InspectableNode(
                displayName: "the resolution law (FI = c·Q on the Z readout)",
                summary: $"Fisher information of the defect readout grows linearly with Q = J/γ " +
                         $"(the resonator quality factor): more coherent cycles per lifetime, more resolving " +
                         $"power per dose. Fitted slope c = {c.ToString("0.######", Inv)}, " +
                         $"worst relative residual = {(resid * 100).ToString("0.##", Inv)} % ⟹ " +
                         $"{(lawPasses ? "PASS" : "FAIL")} at the 25% bar.",
                payload: new InspectablePayload.Real("max relative residual", resid));

            // 3. No exceptional-point peak: monotone in all three bases.
            bool monoZ = IsMonotoneInQ(ReadoutBasis.Z);
            bool monoX = IsMonotoneInQ(ReadoutBasis.X);
            bool monoY = IsMonotoneInQ(ReadoutBasis.Y);
            bool noPeak = monoZ && monoX && monoY;
            yield return new InspectableNode(
                displayName: "no EP peak (FI monotone in Q, every basis)",
                summary: $"no readout basis peaks at the exceptional point Q = 1; it is the worst reading " +
                         $"point in every basis tested. Monotone increasing in Q: " +
                         $"Z = {monoZ}, X = {monoX}, Y = {monoY} ⟹ {(noPeak ? "PASS" : "FAIL")}.");

            // 4. Basis ordering: coherence readouts fall faster than population.
            double spanZ = SpanRatio(ReadoutBasis.Z);
            double spanX = SpanRatio(ReadoutBasis.X);
            double spanY = SpanRatio(ReadoutBasis.Y);
            bool ordering = spanX > 10 * spanZ;
            yield return new InspectableNode(
                displayName: "basis ordering (coherence readouts fade faster than population)",
                summary: $"coherence-basis readouts (X/Y) lose resolving power toward low Q faster than the " +
                         $"population basis (Z); near Q = 1 only the population basis still reads. " +
                         $"Span ratio FI(Q_max)/FI(Q_min): Z = {spanZ.ToString("0.###", Inv)}, " +
                         $"X = {spanX.ToString("0.###", Inv)}, Y = {spanY.ToString("0.###", Inv)} ⟹ " +
                         $"X span / Z span = {(spanX / spanZ).ToString("0.#", Inv)}× " +
                         $"({(ordering ? "PASS" : "FAIL")} at the 10× bar).");

            // 5. Honest scope.
            yield return new InspectableNode(
                displayName: "honest scope",
                summary: $"fixed dose window K ∈ (0, 1] ({Points} grid points); forward difference at " +
                         $"δJ = {DeltaJ.ToString(Inv)} (O(δJ) bias ≈ 1%); N = {N} fixed per instance. " +
                         $"The law is per-DOSE: a fixed lab-time budget rates the two routes to high Q " +
                         $"(raise J, or lower γ) differently. Defect strength read at bond 0; location " +
                         $"discrimination is bond 0 vs bond {N - 2}.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
