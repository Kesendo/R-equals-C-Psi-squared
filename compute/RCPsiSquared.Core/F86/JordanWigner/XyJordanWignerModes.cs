using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>F86 Item 1' Direction (b'') JW track foundation: open-boundary-condition
/// sine-mode basis ψ_k(j) = √(2/(N+1))·sin(π·k·(j+1)/(N+1)) for k = 1..N (1-indexed),
/// with single-particle dispersion ε_k = 2·J·cos(π·k/(N+1)).
///
/// <para><b>Tier outcome: Tier1Derived.</b> Textbook XY Jordan-Wigner identity:
/// the XY chain Hamiltonian H_XY = J·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}) under JW
/// transformation yields free fermions whose single-particle eigenstates are these sine
/// modes with the cosine dispersion. The witnesses:</para>
/// <list type="bullet">
///   <item><b>Row-orthonormality of <see cref="SineModeMatrix"/>.</b> The N×N matrix with
///   <c>M[k-1, j] = ψ_k(j)</c> satisfies <c>M · Mᵀ = I</c> to machine precision.</item>
///   <item><b>Dispersion match against direct hopping-matrix EVD.</b> Diagonalising the
///   tridiagonal hopping matrix <c>h[i, i±1] = J</c> yields the same {ε_k} set.</item>
///   <item><b>Chiral pairing ε_{N+1−k} = −ε_k.</b> Verified at <see cref="Build"/> via
///   <see cref="ChiralPairingResidual"/> ≤ <see cref="ChiralPairingTolerance"/>.</item>
/// </list>
///
/// <para><b>Chirality named 2026-06-10 (fresh-eyes pass):</b> the dispersion satisfies
/// ε_{N+1−k} = −ε_k, verbatim the AZ-class-BDI spectrum inversion of
/// <see cref="RCPsiSquared.Core.Symmetry.ChiralKClaim"/> (K = diag((−1)^ℓ), K·H·K = −H for the
/// tridiagonal hopping matrix this primitive diagonalises), wired as this claim's typed parent
/// in <c>RegisterF86JordanWignerLight</c>; downstream, the cluster differences
/// δ = ε_k − ε_{k₁} − ε_{k₂} of <see cref="JwDispersionStructure"/> flip sign under the global
/// k-reflection, so the cluster tables driving the bond-selection effect are chirally ± paired;
/// the eigenvector-side sibling of the same sublattice chirality is
/// <c>ChiralMirrorTrajectoryClaim</c> (Diagnostics/Ptf, registered next to ChiralKClaim).</para>
///
/// <para>This primitive is the foundation for the JW free-fermion track in Item 1'
/// Direction (b'') — the single surviving 2026-05-06 letter direction (the (a''), (c''),
/// (e'') directions of the same taxonomy are retired). After the (c'') falsification
/// (FourModePartial verdict in <see cref="Item1Derivation.C2DirectionCFalsificationProbe"/>:
/// 4-mode reduction reproduces Interior HWHM exactly but fails Endpoint by 0.28-0.50), the
/// JW track became the structurally honest path for deriving the closed-form Endpoint
/// HWHM/Q_peak constant; the 2026-05-11 F90 bridge identity then numerically closed
/// Direction (b'') by exhibiting the F89 path-(N−1) equivalence per bond, so the remaining
/// JW-track work targets the analytical AT-locked F_a/F_b derivation rather than the
/// numerical structure itself. T2 (<c>C2BlockJwDecomposition</c>, per-bond (k1, k2) bilinear
/// coefficients) and T3 (<c>C2BondKModeProfile</c>, per-bond k-mode profile) build on this
/// primitive.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track
/// foundation; the proof's "Item 1 (status update 2026-05-05)" section frames the open
/// closed-form HWHM/Q_peak target this primitive feeds) + textbook XY-JW identity.</para>
/// </summary>
public sealed class XyJordanWignerModes : Claim
{
    /// <summary>Machine-precision ceiling for the chiral pairing witness |ε_k + ε_{N+1−k}|.
    /// The identity is algebraically exact (cos(π − x) = −cos(x)); the residual is
    /// transcendental-evaluation FP drift, observed ≲ 7e-16 for N ≤ 12.</summary>
    public const double ChiralPairingTolerance = 1e-12;

    public int N { get; }
    public double J { get; }

    /// <summary>Single-particle dispersion <c>ε_k = 2·J·cos(π·k/(N+1))</c> for k = 1..N.
    /// Stored 0-indexed: <c>Dispersion[k-1] = ε_k</c>.</summary>
    public IReadOnlyList<double> Dispersion { get; }

    /// <summary>max_k |ε_k + ε_{N+1−k}|: the construction witness for the chiral pairing
    /// ε_{N+1−k} = −ε_k (the <see cref="RCPsiSquared.Core.Symmetry.ChiralKClaim"/> BDI
    /// spectrum inversion). Verified ≤ <see cref="ChiralPairingTolerance"/> at
    /// <see cref="Build"/>.</summary>
    public double ChiralPairingResidual { get; }

    /// <summary>N×N row-orthonormal matrix with <c>SineModeMatrix[k-1, j] = ψ_k(j)</c>.
    /// Row k-1 is the k-th sine mode evaluated on sites j = 0..N−1.</summary>
    public Matrix<double> SineModeMatrix { get; }

    public static XyJordanWignerModes Build(int N, double J = 1.0)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 1; got {N}.");

        var dispersion = new double[N];
        var matrix = Matrix<double>.Build.Dense(N, N);
        double normalisation = Math.Sqrt(2.0 / (N + 1));
        for (int k = 1; k <= N; k++)
        {
            dispersion[k - 1] = 2.0 * J * Math.Cos(Math.PI * k / (N + 1));
            for (int j = 0; j < N; j++)
                matrix[k - 1, j] = normalisation * Math.Sin(Math.PI * k * (j + 1) / (N + 1));
        }

        // Chiral pairing witness (ChiralKClaim spectrum inversion): ε_{N+1−k} = −ε_k.
        // Dispersion[N − k] is ε_{N+1−k} in the 0-indexed storage. The residual scales
        // linearly with J, so the ceiling does too (J is O(1) throughout the repo).
        double chiralResidual = 0.0;
        for (int k = 1; k <= N; k++)
            chiralResidual = Math.Max(chiralResidual, Math.Abs(dispersion[k - 1] + dispersion[N - k]));
        if (chiralResidual > ChiralPairingTolerance * Math.Max(1.0, Math.Abs(J)))
            throw new InvalidOperationException(
                $"Chiral pairing witness failed at N={N}, J={J:G4}: max_k |ε_k + ε_(N+1−k)| = " +
                $"{chiralResidual:E3} exceeds {ChiralPairingTolerance:E0}·max(1, |J|).");

        return new XyJordanWignerModes(N, J, dispersion, matrix, chiralResidual);
    }

    private XyJordanWignerModes(int n, double j, IReadOnlyList<double> dispersion, Matrix<double> matrix,
        double chiralPairingResidual)
        : base("XY Jordan-Wigner OBC sine-mode basis + dispersion",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1' Direction (b'') (JW track) + textbook XY-JW")
    {
        N = n;
        J = j;
        Dispersion = dispersion;
        SineModeMatrix = matrix;
        ChiralPairingResidual = chiralPairingResidual;
    }

    /// <summary>ψ_k(j) for k ∈ [1, N], j ∈ [0, N−1]. Throws on out-of-range.</summary>
    public double SineMode(int k, int j)
    {
        if (k < 1 || k > N) throw new ArgumentOutOfRangeException(nameof(k), $"k must be in [1, {N}]; got {k}.");
        if (j < 0 || j >= N) throw new ArgumentOutOfRangeException(nameof(j), $"j must be in [0, {N - 1}]; got {j}.");
        return SineModeMatrix[k - 1, j];
    }

    public override string DisplayName => $"XY Jordan-Wigner modes (N = {N}, J = {J:G4})";

    public override string Summary =>
        $"OBC sine-mode basis: {N} modes, ε_min = {Dispersion.Min():F4}, ε_max = {Dispersion.Max():F4} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: N.ToString());
            yield return InspectableNode.RealScalar("J", J, "G4");
            yield return InspectableNode.RealScalar("ε_min", Dispersion.Min(), "F4");
            yield return InspectableNode.RealScalar("ε_max", Dispersion.Max(), "F4");
            yield return new InspectableNode("dispersion (ε_k for k = 1..N)",
                summary: string.Join(", ", Dispersion.Select(e => e.ToString("F4"))));
            yield return new InspectableNode("chirality (ChiralKClaim spectrum inversion)",
                summary: "the JW dispersion is ChiralKClaim's BDI spectrum inversion: " +
                         "ε_{N+1−k} = −ε_k (K = diag((−1)^ℓ), K·H·K = −H on the tridiagonal " +
                         "hopping matrix); the cluster differences δ = ε_k − ε_{k₁} − ε_{k₂} " +
                         "flip sign under global k-reflection, so the JwDispersionStructure " +
                         "cluster table is chirally ± paired; " +
                         $"max_k |ε_k + ε_{{N+1−k}}| = {ChiralPairingResidual:E2}");
        }
    }
}
