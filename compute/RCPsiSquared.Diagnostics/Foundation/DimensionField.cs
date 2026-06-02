using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Object Manager's telescope onto a <see cref="DimensionAxis"/>: a single live
/// IInspectable node that, at every θ on the axis, reads the marks (the Liouvillian eigenvalue
/// multiset, the fixed contract) against the in-between (how the slow eigen-subspace rotates, the
/// content). The C# home of the dimension-field sweep; the in-between is read as the principal-angle
/// rotation of the slow subspace (<see cref="DimensionSweepResult.CumulativeRotation"/>). Caveat: the
/// largest-angle reading saturates on a degenerate slow manifold (one direction rotates fully out at
/// the first step), so it is a coarse eyepiece; the full principal-angle spectrum via
/// <see cref="DimensionSweepResult.SlowBasis"/> is the finer reading.
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

    /// <summary>The resolving in-between eyepiece in degrees: the largest principal angle between
    /// the slow subspace at θ₀ and at θ[p].</summary>
    private double[] CumulativeRotationDegrees()
    {
        var cum = Sweep.CumulativeRotation;
        var deg = new double[cum.Count];
        for (int p = 0; p < cum.Count; p++)
            deg[p] = cum[p] * 180.0 / Math.PI;
        return deg;
    }

    public string DisplayName =>
        $"DimensionField (axis={_axis.Name}, N={_axis.N}, {_axis.Theta.Count} samples)";

    public string Summary
    {
        get
        {
            double drift = Sweep.MaxEigenvalueDriftAcrossTheta;
            double maxRotDeg = 0.0;
            foreach (double a in Sweep.CumulativeRotation)
            {
                double deg = a * 180.0 / Math.PI;
                if (deg > maxRotDeg) maxRotDeg = deg;
            }
            double aLo = Sweep.Polarity[0];
            double aHi = Sweep.Polarity[^1];
            return $"α: {aLo.ToString("0.###", Inv)} → {aHi.ToString("0.###", Inv)}; " +
                   $"marks flat at {drift.ToString("E1", Inv)}; in-between rotates to {maxRotDeg.ToString("0.#", Inv)}°. " +
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

            // 2. The in-between (the content): the principal-angle rotation from θ₀ (largest angle;
            // saturates on a degenerate slow manifold, the full spectrum via SlowBasis is finer).
            var cumDeg = CumulativeRotationDegrees();
            double maxRotDeg = cumDeg.Length == 0 ? 0.0 : cumDeg.Max();
            yield return new InspectableNode(
                displayName: "in-between (the content)",
                summary: $"slow subspace rotates to {maxRotDeg.ToString("0.#", Inv)}° (largest principal angle from θ₀; saturates on degenerate manifolds, see SlowBasis)",
                payload: new InspectablePayload.Curve(
                    "slow-subspace rotation from θ₀", thetaDeg, cumDeg, "θ°", "principal angle°"));

            // 3. The polarity ladder α = sin²θ/2: ¼ at the T-gate (45°), ½ at the S-gate (90°).
            yield return new InspectableNode(
                displayName: "polarity on the ladder",
                summary: $"α = sin²θ/2: {Sweep.Polarity[0].ToString("0.###", Inv)} → {Sweep.Polarity[^1].ToString("0.###", Inv)} (¼ at 45°, ½ at 90°)",
                payload: new InspectablePayload.Curve(
                    "α = sin²θ/2", thetaDeg, Sweep.Polarity, "θ°", "α"));

            // 4. The mirror here at θ = 45° (the T-gate): Ad_{R_z(π/4)} on one qubit, the √-of-90°.
            yield return new InspectableNode(
                displayName: "the mirror here (θ=45°, the T-gate)",
                summary: "Ad_{R_z(π/4)}: the continuous mirror at the symmetric crossover, the √ of the 90° S-gate",
                payload: new InspectablePayload.MatrixView(
                    "Ad_{R_z(π/4)}", RotationSuperoperator.AdRzVec(Math.PI / 4, 1, new[] { 0 })));
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
