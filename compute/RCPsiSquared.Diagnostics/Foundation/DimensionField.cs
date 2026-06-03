using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Object Manager's telescope onto a <see cref="DimensionAxis"/>: a single live
/// IInspectable node that, at every θ on the axis, reads the marks (the Liouvillian eigenvalue
/// multiset, the fixed contract) against the in-between (how the slow eigen-subspace rotates, the
/// content). The C# home of the dimension-field sweep. The in-between is read at two eyepieces: the
/// coarse <see cref="DimensionSweepResult.CumulativeRotation"/> (the single largest principal angle,
/// which saturates on a degenerate slow manifold once one direction rotates fully out), and the
/// resolving <see cref="DimensionSweepResult.PrincipalAngleSpectrum"/> (all k angles per θ). The
/// spectrum is surfaced directly as the fan: it separates the invariant core (angles that stay ≈ 0,
/// the part of the slow manifold the rotation fixes, on the crossover axis the {I, Z} shadow) from
/// the rotating directions (angles that grow, the {X, Y} lit in-between). The slow manifold carries
/// its own marks-and-in-between split, and the fan makes it visible.
///
/// <para>On the <see cref="DimensionAxis.Crossover"/> axis the two readings split cleanly: the
/// marks do not move (L(θ) is a similarity transform of L(0), so the sorted eigenvalues are
/// θ-invariant to the dense-Evd floor), while the slow subspace rotates with R_z(θ). The θ-grid
/// passes three framework gates: θ = 0° is pure XZ, θ = 45° is the T-gate symmetric crossover
/// where the polarity ladder α = sin²θ/2 reads ¼, and θ = 90° is the S-gate
/// (<see cref="RotationSuperoperator"/>'s NinetyDegreeMirror) where α reads ½.</para>
///
/// <para>A plain <see cref="IInspectable"/> (not a Claim): a live reading, not a typed-knowledge
/// assertion. The sweep computes lazily on first access and is cached. Built on demand from a
/// <see cref="DimensionAxis"/> plus the slow-mode count.</para></summary>
public sealed class DimensionField : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private readonly DimensionAxis _axis;
    private readonly int _slowCount;

    /// <summary>The axis this field reads; its θ-grid, N, γ-profile, and H(θ).</summary>
    public DimensionAxis Axis => _axis;

    /// <summary>How many of the slowest (smallest |Re λ|) Liouvillian modes span the in-between's
    /// slow subspace at each θ.</summary>
    public int SlowCount => _slowCount;

    public DimensionField(DimensionAxis axis, int slowCount = 16)
    {
        _axis = axis ?? throw new ArgumentNullException(nameof(axis));
        if (slowCount < 1) throw new ArgumentOutOfRangeException(nameof(slowCount), $"need at least one slow mode; got {slowCount}");
        _slowCount = slowCount;
    }

    private DimensionSweepResult? _sweep;
    /// <summary>The cached sweep: the marks (eigenvalues), the polarity ladder, and both in-between
    /// eyepieces (the coarse chordal <see cref="DimensionSweepResult.SubspaceRotation"/> and the
    /// resolving <see cref="DimensionSweepResult.CumulativeRotation"/>).</summary>
    public DimensionSweepResult Sweep => _sweep ??= DimensionSweep.Compute(_axis, _slowCount);

    /// <summary>The θ-grid in degrees, the x-axis of every curve.</summary>
    private double[] ThetaDegrees()
    {
        var theta = _axis.Theta;
        var deg = new double[theta.Count];
        for (int p = 0; p < theta.Count; p++)
            deg[p] = theta[p] * 180.0 / Math.PI;
        return deg;
    }

    /// <summary>Per-θ marks drift: driftPerTheta[p] = max_k |λ_k(θ_p) − λ_k(θ_0)| over the
    /// index-aligned sorted eigenvalues. The contract curve; sits at ≈0 on the crossover axis.</summary>
    private double[] DriftPerTheta()
    {
        var evals = Sweep.Eigenvalues;
        int points = evals.Count;
        var reference = evals[0];
        var drift = new double[points];
        for (int p = 0; p < points; p++)
        {
            var sorted = evals[p];
            double maxD = 0.0;
            int len = Math.Min(sorted.Length, reference.Length);
            for (int i = 0; i < len; i++)
            {
                double d = (sorted[i] - reference[i]).Magnitude;
                if (d > maxD) maxD = d;
            }
            drift[p] = maxD;
        }
        return drift;
    }

    /// <summary>The coarse in-between eyepiece in degrees: the largest principal angle between the
    /// slow subspace at θ₀ and at θ[p]. Saturates on a degenerate manifold; the fan resolves it.</summary>
    private double[] CumulativeRotationDegrees()
    {
        var cum = Sweep.CumulativeRotation;
        var deg = new double[cum.Count];
        for (int p = 0; p < cum.Count; p++)
            deg[p] = cum[p] * 180.0 / Math.PI;
        return deg;
    }

    /// <summary>The principal-angle fan: a (k × points) real matrix, entry [i, p] = the i-th
    /// principal angle (degrees, ascending) of the slow subspace at θ[p] from θ₀. Rows are the angle
    /// index (the invariant core at the top, near 0; the rotating directions toward the bottom),
    /// columns are θ. Drawn as a magnitude heatmap it shows the fan open: the top rows stay blank
    /// (the core the rotation fixes), the lower rows fill in as θ grows (the in-between that turns).</summary>
    private ComplexMatrix SpectrumFanDegrees()
    {
        var spectrum = Sweep.PrincipalAngleSpectrum;
        int points = spectrum.Count;
        int k = spectrum[0].Length; // reference dimension: one angle per slow-basis column at θ₀
        return ComplexMatrix.Build.Dense(k, points, (i, p) =>
        {
            double[] angles = spectrum[p];
            double deg = i < angles.Length ? angles[i] * 180.0 / Math.PI : 0.0;
            return new Complex(deg, 0.0);
        });
    }

    /// <summary>The split read at the final θ: how many of the k principal angles stay in the
    /// invariant core (below <paramref name="tolDeg"/>) versus rotate out, and how far the rotating
    /// ones turned. The core is the part of the slow manifold the rotation fixes; the rest is the
    /// content. tolDeg = 1° absorbs the dense-Evd floor on the core without catching any genuinely
    /// rotating direction (which turns by tens of degrees).
    ///
    /// <para>The core count is lens-dependent in <see cref="SlowCount"/>. On the N = 3 crossover it
    /// reads true (4) in the window slowCount ≈ 8 to 16: below it the slow subspace is too small to
    /// hold the whole core, above it the subspace fills enough of the operator space that the θ₀ and θ
    /// subspaces re-include each other's rotated images and the apparent core inflates (18 at
    /// slowCount 24, 27 at 32). slowCount = 16, the default, sits at the top of the true window.</para></summary>
    private (int Core, int Rotating, double MaxRotDeg) CoreSplitAtFinalTheta(double tolDeg = 1.0)
    {
        double[] finalAngles = Sweep.PrincipalAngleSpectrum[^1];
        int core = 0, rotating = 0;
        double maxRotDeg = 0.0;
        foreach (double a in finalAngles)
        {
            double deg = a * 180.0 / Math.PI;
            if (deg < tolDeg) core++;
            else { rotating++; if (deg > maxRotDeg) maxRotDeg = deg; }
        }
        return (core, rotating, maxRotDeg);
    }

    public string DisplayName =>
        $"DimensionField (axis={_axis.Name}, N={_axis.N}, {_axis.Theta.Count} samples)";

    public string Summary
    {
        get
        {
            double drift = Sweep.MaxEigenvalueDriftAcrossTheta;
            var (core, rotating, maxRotDeg) = CoreSplitAtFinalTheta();
            double aLo = Sweep.Polarity[0];
            double aHi = Sweep.Polarity[^1];
            return $"α: {aLo.ToString("0.###", Inv)} → {aHi.ToString("0.###", Inv)}; " +
                   $"marks flat at {drift.ToString("E1", Inv)}; in-between: {rotating} of {core + rotating} slow " +
                   $"directions rotate to {maxRotDeg.ToString("0.#", Inv)}°, {core} stay fixed (the core). " +
                   "Gates on the θ-grid: 0°→XZ, 45°→T-gate→α=1/4, 90°→S-gate→α=1/2.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var thetaDeg = ThetaDegrees();

            // 1. The marks (the contract): per-θ eigenvalue drift from θ₀, ≈0 on the crossover axis.
            yield return new InspectableNode(
                displayName: "marks (the contract)",
                summary: $"max |Δλ| across θ = {Sweep.MaxEigenvalueDriftAcrossTheta.ToString("E2", Inv)} (the marks do not move)",
                payload: new InspectablePayload.Curve(
                    "eigenvalue drift from θ₀", thetaDeg, DriftPerTheta(), "θ°", "max |Δλ|"));

            // 2. The in-between (the content), coarse eyepiece: the largest principal angle from θ₀.
            // Saturates on the degenerate slow manifold; the fan (next child) resolves what it hides.
            var cumDeg = CumulativeRotationDegrees();
            double maxRotDeg = cumDeg.Length == 0 ? 0.0 : cumDeg.Max();
            yield return new InspectableNode(
                displayName: "in-between (the content)",
                summary: $"slow subspace rotates to {maxRotDeg.ToString("0.#", Inv)}° (largest principal angle from θ₀; saturates on degenerate manifolds, the fan below resolves the split)",
                payload: new InspectablePayload.Curve(
                    "slow-subspace rotation from θ₀", thetaDeg, cumDeg, "θ°", "principal angle°"));

            // 3. The fan (resolving eyepiece): the full principal-angle spectrum, k angles per θ. The
            // heatmap shows the invariant core (top rows, ≈0, blank) split from the rotating
            // directions (lower rows, filling as θ grows), the structure the single largest angle hides.
            var (core, rotating, fanMaxDeg) = CoreSplitAtFinalTheta();
            yield return new InspectableNode(
                displayName: "the fan (principal-angle spectrum)",
                summary: $"{core + rotating} angles at the final θ: {core} stay ≈0 (the invariant core), " +
                         $"{rotating} rotate to {fanMaxDeg.ToString("0.#", Inv)}° (the in-between)",
                payload: new InspectablePayload.MatrixView(
                    "principal-angle fan (rows: angle index, cols: θ; degrees)", SpectrumFanDegrees()));

            // 4. The polarity ladder α = sin²θ/2: ¼ at the T-gate (45°), ½ at the S-gate (90°).
            yield return new InspectableNode(
                displayName: "polarity on the ladder",
                summary: $"α = sin²θ/2: {Sweep.Polarity[0].ToString("0.###", Inv)} → {Sweep.Polarity[^1].ToString("0.###", Inv)} (¼ at 45°, ½ at 90°)",
                payload: new InspectablePayload.Curve(
                    "α = sin²θ/2", thetaDeg, Sweep.Polarity, "θ°", "α"));

            // 5. The mirror here at θ = 45° (the T-gate): Ad_{R_z(π/4)} on one qubit, the √-of-90°.
            yield return new InspectableNode(
                displayName: "the mirror here (θ=45°, the T-gate)",
                summary: "Ad_{R_z(π/4)}: the continuous mirror at the symmetric crossover, the √ of the 90° S-gate",
                payload: new InspectablePayload.MatrixView(
                    "Ad_{R_z(π/4)}", RotationSuperoperator.AdRzVec(Math.PI / 4, 1, new[] { 0 })));
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
