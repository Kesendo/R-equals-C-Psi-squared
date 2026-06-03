using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Ptf;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Object Manager's telescope onto the PTF J-defect axis, the contrast to
/// <see cref="DimensionField"/>'s crossover. There the in-between is a pure rotation and the marks are
/// frozen by an exact similarity. Here the perturbation (a single detuned bond, δJ) is Π-invariant but
/// NOT a similarity, so the telescope reads a different shape:
///
/// <list type="bullet">
///   <item><b>The contract holds.</b> The palindrome residual ‖Π L Π⁻¹ + L + 2Σγ‖ stays ~10⁻¹⁵ across
///   δJ: the spectrum is kept mirror-symmetric (the mirror is the contract, not the eigenvalues).</item>
///   <item><b>The spectrum moves.</b> Unlike the crossover's frozen marks, the eigenvalues genuinely
///   drift with δJ. It is a real perturbation, not a similarity.</item>
///   <item><b>The in-between is mixing, not rotation.</b> The first-order matrix elements
///   ⟨W_s|V_L|M_{s'}⟩ have a protected kernel (the steady states, shift ≈ 0 by U(1)) but a live
///   off-diagonal (the eigenvectors mix, the PTF α_i). A pure rotation would leave the off-diagonal
///   near zero.</item>
/// </list>
///
/// <para>A plain <see cref="IInspectable"/> (a live reading, not a Claim). N ≤ 6 (dense Liouvillian +
/// eigenvectors; the mixing also needs the eigenvector inverse). The α_i painter-rate phenomenology
/// (the Python PTF workflow) is the deferred twin of this operator-algebra reading.</para></summary>
public sealed class JDefectField : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private readonly DimensionAxis _axis;
    private readonly int _defectBond;
    private readonly int _slowCount;

    /// <summary>The J-defect axis this field reads (δJ grid, N, γ-profile, H(δJ)).</summary>
    public DimensionAxis Axis => _axis;

    /// <summary>The detuned bond index b (the bond on sites b, b+1); canonical PTF defect is 0.</summary>
    public int DefectBond => _defectBond;

    /// <summary>How many of the slowest Liouvillian modes the mixing and fan readings span.</summary>
    public int SlowCount => _slowCount;

    public JDefectField(int N, double gamma, int defectBond,
        double deltaJMax = 0.1, int points = 25, int slowCount = 16)
    {
        if (slowCount < 1) throw new ArgumentOutOfRangeException(nameof(slowCount), $"need at least one slow mode; got {slowCount}");
        _axis = DimensionAxis.JDefect(N, gamma, defectBond, deltaJMax, points);
        _defectBond = defectBond;
        _slowCount = slowCount;
    }

    private DimensionSweepResult? _sweep;
    /// <summary>The cached sweep over δJ: eigenvalues (marks), the slow basis, and the principal-angle
    /// fan of the slow subspace.</summary>
    public DimensionSweepResult Sweep => _sweep ??= DimensionSweep.Compute(_axis, _slowCount);

    private SlowModeMixing.Reading? _mixing;
    /// <summary>The cached first-order mixing reading ⟨W_s|V_L|M_{s'}⟩ at δJ = 0, V_L the defect
    /// bond's ∂L/∂J.</summary>
    public SlowModeMixing.Reading Mixing => _mixing ??= ComputeMixing();

    private SlowModeMixing.Reading ComputeMixing()
    {
        var lA = PauliDephasingDissipator.BuildZ(_axis.Hamiltonian(0.0), _axis.GammaPerSite);
        var vL = BondPerturbation.Build(_axis.N, _defectBond, _defectBond + 1, BondPerturbation.Kind.XY);
        return SlowModeMixing.Compute(lA, vL, _slowCount);
    }

    private double[]? _residualCurve;
    /// <summary>The palindrome residual ‖Π L Π⁻¹ + L + 2Σγ‖_F at each δJ; flat at ~10⁻¹⁵ (the
    /// contract). Cached, since it rebuilds a Liouvillian and a Pauli-basis conjugation per point.</summary>
    private double[] ResidualCurve => _residualCurve ??= ComputeResidualCurve();

    private double[] ComputeResidualCurve()
    {
        double sigmaGamma = 0.0;
        foreach (double g in _axis.GammaPerSite) sigmaGamma += g;
        var theta = _axis.Theta;
        var res = new double[theta.Count];
        for (int p = 0; p < theta.Count; p++)
        {
            var L = PauliDephasingDissipator.BuildZ(_axis.Hamiltonian(theta[p]), _axis.GammaPerSite);
            res[p] = PalindromeResidual.Build(L, _axis.N, sigmaGamma).FrobeniusNorm();
        }
        return res;
    }

    private double[]? _movementCurve;
    /// <summary>The cached reorder-robust spectrum-movement curve (the directed Hausdorff per δJ).</summary>
    private double[] MovementCurve => _movementCurve ??= SpectrumMovementCurve();

    /// <summary>Per-δJ spectrum movement, reorder-robust: the directed Hausdorff distance
    /// max_i min_j |λ_i(0) − λ_j(δJ)|, how far the farthest original eigenvalue's nearest neighbour
    /// has moved. Nonzero and rising (the spectrum genuinely moves), unlike the crossover's dead-flat
    /// similarity. An index-aligned drift would be inflated by sort-order swaps here, because, unlike
    /// the crossover, the eigenvalues actually move and cross; the nearest-neighbour distance does
    /// not care about ordering.</summary>
    private double[] SpectrumMovementCurve()
    {
        var evals = Sweep.Eigenvalues;
        var reference = evals[0];
        var move = new double[evals.Count];
        for (int p = 0; p < evals.Count; p++)
        {
            var cur = evals[p];
            double maxNearest = 0.0;
            foreach (var lam0 in reference)
            {
                double nearest = double.PositiveInfinity;
                foreach (var lam in cur)
                {
                    double d = (lam0 - lam).Magnitude;
                    if (d < nearest) nearest = d;
                }
                if (nearest > maxNearest) maxNearest = nearest;
            }
            move[p] = maxNearest;
        }
        return move;
    }

    /// <summary>The principal-angle fan as a (k × points) real matrix, entry [i, p] = the i-th
    /// principal angle (degrees) of the slow subspace at δJ[p] from δJ = 0. The slow subspace moving
    /// as the eigenvectors mix.</summary>
    private ComplexMatrix SpectrumFanDegrees()
    {
        var spectrum = Sweep.PrincipalAngleSpectrum;
        int points = spectrum.Count;
        int k = spectrum[0].Length;
        return ComplexMatrix.Build.Dense(k, points, (i, p) =>
        {
            double[] angles = spectrum[p];
            double deg = i < angles.Length ? angles[i] * 180.0 / Math.PI : 0.0;
            return new Complex(deg, 0.0);
        });
    }

    /// <summary>The mixing matrix as real magnitudes |⟨W_s|V_L|M_{s'}⟩| for the heatmap.</summary>
    private ComplexMatrix MixingMagnitude()
    {
        var m = Mixing.Mixing;
        return ComplexMatrix.Build.Dense(m.RowCount, m.ColumnCount,
            (i, j) => new Complex(m[i, j].Magnitude, 0.0));
    }

    /// <summary>The kernel protection: how many slow modes are kernel (steady states, |Re λ| &lt; 10⁻⁶)
    /// and the largest first-order shift among them (≈ 0 by U(1) conservation).</summary>
    private (int Count, double MaxShift) KernelProtection()
    {
        var shifts = Mixing.DiagonalShiftMagnitudes;
        var evals = Mixing.SlowEigenvalues;
        int count = 0;
        double maxShift = 0.0;
        for (int s = 0; s < shifts.Count; s++)
            if (Math.Abs(evals[s].Real) < 1e-6)
            {
                count++;
                if (shifts[s] > maxShift) maxShift = shifts[s];
            }
        return (count, maxShift);
    }

    public string DisplayName =>
        $"JDefectField (N={_axis.N}, defect bond {_defectBond}=({_defectBond},{_defectBond + 1}), {_axis.Theta.Count} δJ samples)";

    public string Summary
    {
        get
        {
            double maxRes = ResidualCurve.Max();
            double movement = MovementCurve.Max();
            var (kernelCount, maxKernelShift) = KernelProtection();
            return $"palindrome holds across δJ (residual ≤ {maxRes.ToString("E1", Inv)}, the contract); " +
                   $"spectrum moves (to {movement.ToString("0.###", Inv)}, no similarity); " +
                   $"kernel ({kernelCount} steady states) protected to {maxKernelShift.ToString("E1", Inv)}, " +
                   $"eigenvectors mix (off-diagonal {Mixing.OffDiagonalMass.ToString("0.00", Inv)}).";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var deltaJ = _axis.Theta.ToArray();

            // 1. The contract: the palindrome residual stays ~1e-15 across δJ (the mirror is kept).
            var res = ResidualCurve;
            yield return new InspectableNode(
                displayName: "the contract (the palindrome holds)",
                summary: $"palindrome residual ≤ {res.Max().ToString("E1", Inv)} across δJ (Π-invariant: the mirror is kept, even as the spectrum moves)",
                payload: new InspectablePayload.Curve(
                    "palindrome residual", deltaJ, res, "δJ", "‖Π L Π⁻¹ + L + 2Σγ‖"));

            // 2. The spectrum moves: the nearest-eigenvalue movement rises with δJ, unlike the
            // crossover's frozen marks (reorder-robust, so sort swaps do not inflate it).
            var movement = MovementCurve;
            yield return new InspectableNode(
                displayName: "the spectrum moves (no similarity)",
                summary: $"eigenvalues move to {movement.Max().ToString("0.###", Inv)} across δJ (a real perturbation, not the crossover's frozen marks)",
                payload: new InspectablePayload.Curve(
                    "spectrum movement from δJ=0", deltaJ, movement, "δJ", "max nearest |Δλ|"));

            // 3. The mixing: the kernel is protected (U(1)) but the off-diagonal is alive (the in-between
            // is eigenvector mixing, not a rotation). Two readings: the per-mode shift profile, the matrix.
            var (kernelCount, maxKernelShift) = KernelProtection();
            var shifts = Mixing.DiagonalShiftMagnitudes;
            var shiftVec = ComplexVector.Build.Dense(shifts.Count, i => new Complex(shifts[i], 0.0));
            var shiftLabels = Mixing.SlowEigenvalues
                .Select(e => $"Re{e.Real.ToString("0.0", Inv)}").ToList();
            yield return new InspectableNode(
                displayName: "the mixing (more than a rotation)",
                summary: $"kernel ({kernelCount} steady states) protected to {maxKernelShift.ToString("E1", Inv)} (U(1)); " +
                         $"the coherences shift; off-diagonal mixing {Mixing.OffDiagonalMass.ToString("0.00", Inv)} alive (the eigenvectors mix, not a rotation)",
                children: new IInspectable[]
                {
                    new InspectableNode(
                        displayName: "the shift profile (kernel protected, coherences move)",
                        summary: $"first-order |δλ_s| per slow mode: {kernelCount} kernel ≈ 0, the rest shift",
                        payload: new InspectablePayload.Vector("first-order |δλ_s|", shiftVec, shiftLabels)),
                    new InspectableNode(
                        displayName: "the mixing matrix (eigenvectors mix)",
                        summary: "|⟨W_s|V_L|M_s'⟩| on the slow modes: dark diagonal kernel, live off-diagonal",
                        payload: new InspectablePayload.MatrixView("mixing |⟨W|V_L|M⟩|", MixingMagnitude())),
                });

            // 4. The slow subspace moves: the principal-angle fan across δJ (the in-between moving).
            yield return new InspectableNode(
                displayName: "the slow subspace moves (the fan)",
                summary: "principal-angle spectrum of the slow subspace across δJ (the in-between in motion)",
                payload: new InspectablePayload.MatrixView(
                    "principal-angle fan (rows: angle index, cols: δJ; degrees)", SpectrumFanDegrees()));
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
