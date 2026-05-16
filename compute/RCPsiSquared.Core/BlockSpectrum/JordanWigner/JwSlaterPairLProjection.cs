using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.BlockSpectrum.JordanWigner;

/// <summary>Full L = L_H + L_D Liouvillian block of the chain XY + Z-dephasing dynamics
/// transformed into the JW Slater-pair basis of <see cref="JwSlaterPairBasis"/>. The
/// Hamiltonian part is diagonal (as the basis already witnesses); the dissipator part is
/// the new structural content of this primitive — its off-diagonal support per row is
/// bounded by the Slater-pair reach of a single Z_l ⊗ Z_l swap:
///
/// <code>
///   nnz_off(row) ≤ (1 + p_r·(N − p_r)) · (1 + p_c·(N − p_c)) − 1
/// </code>
///
/// <para>This bound is structural: Z_l = 1 − 2·Σ_{k,k'} ψ_k(l)·ψ_{k'}(l)·η_k†·η_{k'} acts
/// on a Slater determinant |L⟩ by either leaving it unchanged or producing a single
/// orbital-swap L → (L \ {k'}) ∪ {k} with k' ∈ L, k ∉ L — at most 1 + |L|·(N − |L|) target
/// Slater determinants. Z_l ⊗ Z_l on |L⟩⟨K| acts on both sides independently; the support
/// is the Cartesian product, minus the source itself.</para>
///
/// <para><b>Tier outcome: <see cref="Tier.Tier1Derived"/>.</b> Pure basis change of the
/// dense L block; the sparsity bound is algebraic from the JW expansion of Z_l. Witnesses:</para>
/// <list type="bullet">
///   <item><see cref="TheoreticalNnzBoundPerRow"/> = (1 + p_r(N−p_r))·(1 + p_c(N−p_c)) − 1.</item>
///   <item><see cref="MaxOffDiagonalNonzeroPerRow"/>, <see cref="MeanOffDiagonalNonzeroPerRow"/>:
///         empirical counts of |L_JW[α, β]| > <see cref="NonzeroThreshold"/> for α ≠ β.</item>
///   <item><see cref="RespectsTheoreticalBound"/> := empirical max ≤ theoretical bound.</item>
/// </list>
///
/// <para>Stepping-stone for the N=10 sparse-eig path (Phase 2 of the N=10 push plan): the
/// sparsity bound proven here justifies a future sparse-construction primitive that builds
/// L_JW directly element-by-element without the dense U^†·L·U detour, then feeds a sparse
/// Arnoldi/Lanczos eig at sector dimensions beyond <see cref="JwSlaterPairBasis.MaxSectorDimForDenseWitness"/>.</para>
///
/// <para>Anchor: <see cref="JwSlaterPairBasis"/> (basis transformation) +
/// <see cref="PerBlockLiouvillianBuilder"/> (computational-basis block) +
/// textbook XY-JW Slater-determinant action of Z_l.</para>
/// </summary>
public sealed class JwSlaterPairLProjection : Claim
{
    /// <summary>Magnitude threshold below which an off-diagonal entry is treated as zero
    /// for sparsity counting. Chosen to swallow FP round-off in the dense U^†·L·U product
    /// (typical residual ~1e-14) while still treating any real algebraic non-zero as such.</summary>
    public const double NonzeroThreshold = 1e-12;

    public JwSlaterPairBasis Basis { get; }
    public IReadOnlyList<double> GammaPerSite { get; }

    /// <summary>L in the original computational basis on the sector defined by
    /// <see cref="JwSlaterPairBasis.FlatIndices"/>. Kept for unitary-invariance witnesses
    /// at small N; do not depend on this property beyond the dense witness regime.</summary>
    public ComplexMatrix LBlockComputational { get; }

    /// <summary>L in the Slater-pair basis: <c>L_JW = U^† · L_block · U</c>. Diagonal
    /// holds the Hamiltonian eigenvalues −i·(Σε(L) − Σε(K)) plus the JW-basis dissipator
    /// diagonal; off-diagonal is the dissipator's Slater-swap content.</summary>
    public ComplexMatrix LJw { get; }

    public int MaxOffDiagonalNonzeroPerRow { get; }
    public double MeanOffDiagonalNonzeroPerRow { get; }
    public int TheoreticalNnzBoundPerRow { get; }
    public bool RespectsTheoreticalBound => MaxOffDiagonalNonzeroPerRow <= TheoreticalNnzBoundPerRow;

    public static JwSlaterPairLProjection Build(JwSlaterPairBasis basis, IReadOnlyList<double> gammaPerSite)
    {
        if (basis is null) throw new ArgumentNullException(nameof(basis));
        if (gammaPerSite is null) throw new ArgumentNullException(nameof(gammaPerSite));
        if (gammaPerSite.Count != basis.N)
            throw new ArgumentException(
                $"gammaPerSite length {gammaPerSite.Count} != basis.N {basis.N}", nameof(gammaPerSite));

        var H = PauliHamiltonian.XYChain(basis.N, basis.J).ToMatrix();
        var lBlock = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, basis.FlatIndices);
        var lJw = basis.Uinv * lBlock * basis.U;

        int dim = lJw.RowCount;
        int maxNnz = 0;
        long totalNnz = 0;
        for (int alpha = 0; alpha < dim; alpha++)
        {
            int rowNnz = 0;
            for (int beta = 0; beta < dim; beta++)
            {
                if (beta == alpha) continue;
                Complex z = lJw[alpha, beta];
                if (z.Real * z.Real + z.Imaginary * z.Imaginary > NonzeroThreshold * NonzeroThreshold)
                    rowNnz++;
            }
            if (rowNnz > maxNnz) maxNnz = rowNnz;
            totalNnz += rowNnz;
        }
        double meanNnz = (double)totalNnz / dim;

        int bound = (1 + basis.PRow * (basis.N - basis.PRow))
                  * (1 + basis.PCol * (basis.N - basis.PCol))
                  - 1;

        return new JwSlaterPairLProjection(basis, gammaPerSite.ToArray(), lBlock, lJw, maxNnz, meanNnz, bound);
    }

    private JwSlaterPairLProjection(
        JwSlaterPairBasis basis,
        double[] gammaPerSite,
        ComplexMatrix lBlock,
        ComplexMatrix lJw,
        int maxOffDiagonalNonzeroPerRow,
        double meanOffDiagonalNonzeroPerRow,
        int theoreticalNnzBoundPerRow)
        : base($"L_JW = U^† · L · U on JW Slater-pair sector (p_c={basis.PCol}, p_r={basis.PRow}, N={basis.N}): " +
               $"dim {lJw.RowCount}, max-off-diagonal-nnz/row {maxOffDiagonalNonzeroPerRow} ≤ " +
               $"theoretical {theoreticalNnzBoundPerRow} = (1 + p_r·(N−p_r))·(1 + p_c·(N−p_c)) − 1.",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/BlockSpectrum/JordanWigner/JwSlaterPairBasis.cs (basis transformation U) + " +
               "compute/RCPsiSquared.Core/BlockSpectrum/PerBlockLiouvillianBuilder.cs (computational-basis L block); " +
               "JW expansion Z_l = 1 − 2·Σ ψ_k(l)·ψ_{k'}(l)·η_k†·η_{k'} bounds the per-side Slater-pair " +
               "reach by 1 + p·(N − p), and Z_l ⊗ Z_l on |L⟩⟨K| factors as the Cartesian product.")
    {
        Basis = basis;
        GammaPerSite = gammaPerSite;
        LBlockComputational = lBlock;
        LJw = lJw;
        MaxOffDiagonalNonzeroPerRow = maxOffDiagonalNonzeroPerRow;
        MeanOffDiagonalNonzeroPerRow = meanOffDiagonalNonzeroPerRow;
        TheoreticalNnzBoundPerRow = theoreticalNnzBoundPerRow;
    }

    public override string DisplayName =>
        $"L_JW Slater-pair projection (p_c={Basis.PCol}, p_r={Basis.PRow}, N={Basis.N}); " +
        $"dim {LJw.RowCount}";

    public override string Summary =>
        $"max-nnz/row={MaxOffDiagonalNonzeroPerRow}/{TheoreticalNnzBoundPerRow} (bound), " +
        $"mean={MeanOffDiagonalNonzeroPerRow:F1}, density={(double)MaxOffDiagonalNonzeroPerRow / LJw.RowCount:P1} " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return Basis;
            yield return InspectableNode.RealScalar("N", Basis.N);
            yield return InspectableNode.RealScalar("p_c", Basis.PCol);
            yield return InspectableNode.RealScalar("p_r", Basis.PRow);
            yield return InspectableNode.RealScalar("sector dimension", LJw.RowCount);
            yield return InspectableNode.RealScalar("max off-diagonal nnz/row", MaxOffDiagonalNonzeroPerRow);
            yield return InspectableNode.RealScalar("mean off-diagonal nnz/row", MeanOffDiagonalNonzeroPerRow, "F2");
            yield return InspectableNode.RealScalar("theoretical nnz bound/row", TheoreticalNnzBoundPerRow);
            yield return InspectableNode.RealScalar("respects theoretical bound", RespectsTheoreticalBound ? 1 : 0);
        }
    }
}
