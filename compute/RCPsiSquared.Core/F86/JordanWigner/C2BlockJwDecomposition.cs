using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>F86 Item 1' Direction (b'') JW track, T2 of 4: per-bond bilinear-fermion
/// coefficients in the OBC sine-mode basis.
///
/// <para>For each bond b ∈ [0, N−2] the c=2 hopping term H_b = J·(X_b X_{b+1} + Y_b Y_{b+1})
/// becomes, under Jordan-Wigner + sine-mode change of basis (T1 <see cref="XyJordanWignerModes"/>),
/// a sum of bilinear fermion operators η_{k1}^† η_{k2} weighted by</para>
///
/// <para><c>C_b(k1, k2) = ψ_{k1}(b)·ψ_{k2}(b+1) + ψ_{k1}(b+1)·ψ_{k2}(b)</c></para>
///
/// <para>which is symmetric in (k1, k2) by the hopping-symmetric construction. Indices k1,
/// k2 run over 1..N at the math level; in the matrix representation
/// <see cref="BondJwCoefficients.Coefficients"/> they are stored 0-indexed, so the matrix
/// entry <c>[k1−1, k2−1]</c> holds <c>C_b(k1, k2)</c>.</para>
///
/// <para><b>Class-level Tier: <c>Tier2Verified</c>.</b> The single-particle decomposition
/// is exact algebra (textbook XY-JW); the F73 sum-rule check at construction is the
/// verification:</para>
/// <list type="bullet">
///   <item><see cref="F73SumRuleResidual"/> = <c>‖Σ_b C_b − C_total‖_F</c> with
///   <c>C_total[k1, k2] = (ε_{k1}/J)·δ_{k1, k2}</c> (textbook XY identity:
///   <c>Σ_b [ψ_{k1}(b)ψ_{k2}(b+1) + ψ_{k1}(b+1)ψ_{k2}(b)] = 2 cos(π k1 / (N+1)) δ_{k1,k2}
///   = (ε_{k1}/J) δ_{k1,k2}</c>) below <see cref="F73Tolerance"/> = 1e-10 across N = 5..8.</item>
///   <item>Mapping to the full bilinear-fermion superoperator structure on the (n=1, n+1=2)
///   coherence block requires a full XY-block-L mapping that is documented but not promoted
///   to <c>Tier1Derived</c> in this primitive. The Tier 1 promotion path is the open work
///   noted alongside the JW track in <c>F86OpenQuestions</c>.</item>
/// </list>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track).</para>
///
/// <para>F90 status (2026-05-11): the F86 c=2 ↔ F89 path-(N−1) bridge identity
/// achieves numerical Tier-1 for Direction (b'') via per-bond Hellmann-Feynman
/// (bit-exact 20/22 bonds at N=5..8). The JW-track primitives in this file remain
/// active as the alternative analytical route toward the closed-form HWHM_left/Q_peak
/// constants; the per-bond numerical answer itself is no longer the open piece.
/// See <c>docs/proofs/PROOF_F90_F86C2_BRIDGE.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F90F86C2BridgeIdentity.cs</c>.</para>
/// </summary>
public sealed class C2BlockJwDecomposition : Claim
{
    /// <summary>F73 sum-rule residual tolerance: <c>1e-10</c>. The sum
    /// <c>Σ_b C_b</c> matches the diagonal dispersion identity <c>(ε_{k1}/J)·δ_{k1, k2}</c>
    /// to this precision under the textbook XY-JW identity.</summary>
    public const double F73Tolerance = 1e-10;

    public CoherenceBlock Block { get; }

    /// <summary>Composition: T1 sine-mode basis + dispersion. Built once at <see cref="Build"/>
    /// time and reused for every bond.</summary>
    public XyJordanWignerModes Modes { get; }

    /// <summary>One <see cref="BondJwCoefficients"/> per bond, in bond-index order, tagged
    /// with <see cref="BondClass"/>. Each entry exposes the N×N coefficient matrix
    /// <c>C_b(k1, k2)</c>.</summary>
    public IReadOnlyList<BondJwCoefficients> Bonds { get; }

    /// <summary><c>‖Σ_b C_b − C_total‖_F</c> where
    /// <c>C_total[k1, k2] = (Modes.Dispersion[k1−1]/J) · δ_{k1, k2}</c>. The textbook XY-JW
    /// identity guarantees this is zero analytically; the residual measures floating-point
    /// drift. Below <see cref="F73Tolerance"/> = 1e-10 the decomposition is verified.</summary>
    public double F73SumRuleResidual { get; }

    /// <summary>Public factory: validates c=2, builds the T1 sine-mode basis once at the
    /// block's effective hopping J (default J = 1.0), computes per-bond N×N coefficient
    /// matrices, and verifies the F73 sum-rule residual against the textbook dispersion
    /// identity. Static-factory pattern keeps future Tier 1 promotion type-safe.</summary>
    public static C2BlockJwDecomposition Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2BlockJwDecomposition applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        // CoherenceBlock does not expose a hopping J; the c=2 ChainSystem convention
        // used elsewhere in F86 is J = 1.0 (textbook XY-JW unit normalisation).
        const double J = 1.0;
        int N = block.N;
        int numBonds = block.NumBonds;

        var modes = XyJordanWignerModes.Build(N, J);
        var bonds = new BondJwCoefficients[numBonds];
        for (int b = 0; b < numBonds; b++)
        {
            var coeffs = Matrix<double>.Build.Dense(N, N);
            for (int k1 = 0; k1 < N; k1++)
                for (int k2 = 0; k2 < N; k2++)
                {
                    double psik1_b   = modes.SineMode(k1 + 1, b);
                    double psik1_bp1 = modes.SineMode(k1 + 1, b + 1);
                    double psik2_b   = modes.SineMode(k2 + 1, b);
                    double psik2_bp1 = modes.SineMode(k2 + 1, b + 1);
                    coeffs[k1, k2] = psik1_b * psik2_bp1 + psik1_bp1 * psik2_b;
                }
            var bondClass = BondClassExtensions.OfBond(b, numBonds);
            bonds[b] = new BondJwCoefficients(Bond: b, BondClass: bondClass, Coefficients: coeffs);
        }

        var sum = Matrix<double>.Build.Dense(N, N);
        foreach (var bond in bonds) sum += bond.Coefficients;

        var expected = Matrix<double>.Build.Dense(N, N);
        for (int k = 0; k < N; k++) expected[k, k] = modes.Dispersion[k] / J;

        double residual = (sum - expected).FrobeniusNorm();

        return new C2BlockJwDecomposition(block, modes, bonds, residual);
    }

    private C2BlockJwDecomposition(
        CoherenceBlock block,
        XyJordanWignerModes modes,
        IReadOnlyList<BondJwCoefficients> bonds,
        double f73SumRuleResidual)
        : base("c=2 block JW decomposition: per-bond bilinear coefficients C_b(k1, k2)",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1' Direction (b'') (JW track)")
    {
        Block = block;
        Modes = modes;
        Bonds = bonds;
        F73SumRuleResidual = f73SumRuleResidual;
    }

    public override string DisplayName =>
        $"c=2 block JW decomposition (N={Block.N}, bonds={Block.NumBonds})";

    public override string Summary =>
        $"per-bond C_b(k1, k2) = ψ_{{k1}}(b)·ψ_{{k2}}(b+1) + ψ_{{k1}}(b+1)·ψ_{{k2}}(b); " +
        $"F73 sum-rule residual = {F73SumRuleResidual:G3} (tol {F73Tolerance:G3}) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("NumBonds", summary: Block.NumBonds.ToString());
            yield return Modes;
            yield return InspectableNode.RealScalar("F73 sum-rule residual", F73SumRuleResidual, "G3");
            yield return InspectableNode.RealScalar("F73 tolerance", F73Tolerance, "G3");
            yield return InspectableNode.Group("Bonds", Bonds, Bonds.Count);
        }
    }
}

/// <summary>Per-bond JW bilinear coefficients C_b(k1, k2) for the c=2 block, plus the bond's
/// <see cref="BondClass"/>. Visible in the inspection tree under
/// <see cref="C2BlockJwDecomposition.Bonds"/>.
///
/// <para>The N×N matrix <see cref="Coefficients"/> stores
/// <c>Coefficients[k1−1, k2−1] = ψ_{k1}(b)·ψ_{k2}(b+1) + ψ_{k1}(b+1)·ψ_{k2}(b)</c>
/// for k1, k2 = 1..N (0-indexed in the matrix). Symmetric in (k1, k2) by hopping-symmetric
/// construction.</para>
/// </summary>
public sealed record BondJwCoefficients(
    int Bond,
    BondClass BondClass,
    Matrix<double> Coefficients
) : IInspectable
{
    public string DisplayName => $"bond {Bond} JW coefficients ({BondClass})";

    public string Summary =>
        $"N×N C_b(k1, k2), ‖C_b‖_F = {Coefficients.FrobeniusNorm():F4}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("bond", summary: Bond.ToString());
            yield return new InspectableNode("bond class", summary: BondClass.ToString());
            yield return InspectableNode.RealScalar("‖C_b‖_F", Coefficients.FrobeniusNorm(), "F4");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
