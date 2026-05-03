using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F71;

/// <summary>F71 chain-mirror operator R as a typed claim. Wraps
/// <see cref="ChainMirror.Build"/> with an inspectable representation: matrix dimension,
/// involution check (R² = I), and references to the symmetric / antisymmetric projectors.
/// </summary>
public sealed class F71MirrorOperator : Claim
{
    public int N { get; }
    public ComplexMatrix R { get; }

    private readonly Lazy<double> _involutionResidual;

    /// <summary>Frobenius norm of (R² − I). Should be ~0 for the F71 mirror; bit-exact since
    /// R is a permutation matrix.</summary>
    public double InvolutionResidual => _involutionResidual.Value;

    public F71MirrorOperator(int N)
        : base("F71 chain-mirror operator R",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md")
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 2; got {N}");
        this.N = N;
        R = ChainMirror.Build(N);
        _involutionResidual = new Lazy<double>(() =>
            (R * R - ComplexMatrix.Build.SparseIdentity(R.RowCount)).FrobeniusNorm());
    }

    public override string DisplayName => $"R (N={N}, dim {R.RowCount})";

    public override string Summary =>
        $"chain-mirror operator on {1 << N}-dim Hilbert space; R² = I (residual {InvolutionResidual:E2})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("N", N);
            yield return InspectableNode.RealScalar("dim", R.RowCount);
            yield return InspectableNode.RealScalar("‖R² − I‖_F (involution residual)", InvolutionResidual, "E3");
            yield return new InspectableNode("symmetric projector",
                summary: "P_sym = (I + R)/2, F71 +1 eigenspace");
            yield return new InspectableNode("antisymmetric projector",
                summary: "P_asym = (I − R)/2, F71 −1 eigenspace");
        }
    }
}
