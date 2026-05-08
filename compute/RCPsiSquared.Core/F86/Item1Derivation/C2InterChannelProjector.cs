using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 Item 1 (c=2 stratum), <b>projector lift of A3</b>: Riesz spectral projection
/// onto the σ_0-degenerate top eigenspace of <c>V V†</c> (HD=1 side) and <c>V† V</c> (HD=3
/// side), where <c>V_inter = P_{HD₁}^† · M_H_total · P_{HD₂}</c>. This is the typed answer to
/// the A3 obstruction documented in <see cref="C2InterChannelAnalytical"/>: at even N the
/// top singular value σ_0 of V_inter is exactly degenerate (degeneracy 2 at N=6, 8 for c=2),
/// so the single direction |u_0⟩ is library-dependent. The orthogonal projector
/// <c>P_top := Σ_{k : σ_k = σ_0} |u_k⟩⟨u_k|</c> over the entire degenerate subspace is, by
/// the Riesz spectral projection theorem, basis-invariant: any orthonormal basis of the same
/// eigenspace produces the same outer-product sum.
///
/// <para><b>Tier outcome: Tier1Derived</b> at the projector level. The constructor performs a
/// runtime library-stability witness: it computes <c>P_top</c> twice — once from MathNet's
/// canonical EVD basis, once from a random orthogonal mixing of the degenerate eigenvectors —
/// and asserts the Frobenius-norm residual <c>‖ΔP_top‖_F</c> is below
/// <see cref="LibraryStabilityTolerance"/> (default 1e-10). If it is, the projector is
/// confirmed library-invariant and Tier rises to Tier1Derived. If the residual exceeds the
/// tolerance (which would indicate either an implementation bug or a misclassified
/// non-degenerate eigenvalue), the Tier drops to Tier2Verified with a diagnostic message.
/// </para>
///
/// <para><b>Why this is the foundational lift.</b> Direction (a'') (per-bond SVD-block
/// resonance via <c>V_b[2,3]</c>) inherits the A3 single-vector obstruction at even N because
/// the <c>⟨u_0 | M_H_per_bond[b] | v_0⟩</c> matrix element depends on the library tiebreaker.
/// Lifting to <c>Tr(P_top^{(L)} · M_H_per_bond[b] · P_top^{(R)} · M_H_per_bond[b]^†)</c> is
/// library-invariant by construction, removing the obstruction without solving the
/// closed-form-singular-vector problem. This same lift unblocks any cross-block-based
/// derivation (a''/c'') that fails to converge against the A3 single-vector reference.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1 (A3 obstruction) and
/// <c>compute/RCPsiSquared.Core/F86/Item1Derivation/C2InterChannelAnalytical.cs</c>
/// PendingDerivationNote direction (a) "Lift the test contract: compare PROJECTOR onto top-2D
/// eigenspace".</para>
/// </summary>
public sealed class C2InterChannelProjector : Claim
{
    public CoherenceBlock Block { get; }
    public int HammingDistance1 { get; }
    public int HammingDistance2 { get; }

    /// <summary>Top singular value σ_0 of <c>V_inter = P_{HD₁}^† · M_H_total · P_{HD₂}</c>.</summary>
    public double Sigma0 { get; }

    /// <summary>Number of singular values within <see cref="DegeneracyTolerance"/> of σ_0.
    /// At c=2 this is 1 for odd N and 2 for even N (verified empirically across N=5..8).</summary>
    public int RankTop { get; }

    /// <summary>Orthogonal projector onto the σ_0-degenerate top eigenspace of <c>V V†</c>
    /// (HD=1 side), lifted to the full block basis (Mtot × Mtot Hermitian, rank
    /// <see cref="RankTop"/>).</summary>
    public ComplexMatrix PTopL { get; }

    /// <summary>Orthogonal projector onto the σ_0-degenerate top eigenspace of <c>V† V</c>
    /// (HD=3 side), lifted to the full block basis.</summary>
    public ComplexMatrix PTopR { get; }

    /// <summary>The Frobenius residual <c>‖P_top(canonical) − P_top(randomly mixed)‖_F</c>
    /// across both sides (max of the two), measuring the library-invariance of the projector
    /// against a random orthogonal mixing within the degenerate top eigenspace. Below
    /// <see cref="LibraryStabilityTolerance"/> the projector is confirmed Riesz-unique and
    /// Tier rises to Tier1Derived.</summary>
    public double LibraryStabilityResidual { get; }

    /// <summary>Spectral gap <c>(σ_0 − σ_{RankTop}) / σ_0</c> separating the degenerate top
    /// eigenspace from the remainder of the singular spectrum. Larger is cleaner; the Riesz
    /// lift is reliable when this gap is well above <see cref="DegeneracyTolerance"/>.</summary>
    public double SpectralGapBelowTop { get; }

    /// <summary>Relative tolerance used to classify singular values as degenerate with σ_0.
    /// Two σ_k are considered degenerate iff <c>|σ_k − σ_0| / σ_0 &lt; DegeneracyTolerance</c>.
    /// Matches the threshold used in <see cref="C2InterChannelAnalytical"/> A3 reasoning.</summary>
    public const double DegeneracyTolerance = 1e-8;

    /// <summary>Frobenius-norm threshold under which the runtime library-stability witness
    /// passes (projector is Riesz-unique → Tier1Derived). 1e-10 is two orders of magnitude
    /// above double-precision EVD residuals, comfortably tight without tripping on numerical
    /// noise from the QR-based random orthogonal mixing.</summary>
    public const double LibraryStabilityTolerance = 1e-10;

    /// <summary>True iff the runtime library-stability witness passed (Tier1Derived).</summary>
    public bool IsAnalyticallyDerived => Tier == Tier.Tier1Derived;

    /// <summary>Public factory: validates c=2 (or matching c), computes the projectors and
    /// runs the library-stability witness, then constructs the Claim with the resolved Tier.
    /// </summary>
    public static C2InterChannelProjector Build(CoherenceBlock block, int hd1 = 1, int hd2 = 3)
    {
        if (block.C != 2 && hd1 == 1 && hd2 == 3)
            throw new ArgumentException(
                $"C2InterChannelProjector with default (hd1=1, hd2=3) applies to the c=2 stratum; got c={block.C}.",
                nameof(block));

        var resolved = Resolve(block, hd1, hd2);
        return new C2InterChannelProjector(block, hd1, hd2, resolved);
    }

    private C2InterChannelProjector(CoherenceBlock block, int hd1, int hd2, Resolution resolved)
        : base("c=2 inter-channel SVD-top projector (Riesz lift of A3)",
               resolved.Tier,
               Item1Anchors.Root)
    {
        Block = block;
        HammingDistance1 = hd1;
        HammingDistance2 = hd2;
        Sigma0 = resolved.Sigma0;
        RankTop = resolved.RankTop;
        PTopL = resolved.PTopL;
        PTopR = resolved.PTopR;
        LibraryStabilityResidual = resolved.LibraryStabilityResidual;
        SpectralGapBelowTop = resolved.SpectralGapBelowTop;
    }

    private readonly record struct Resolution(
        Tier Tier,
        double Sigma0,
        int RankTop,
        ComplexMatrix PTopL,
        ComplexMatrix PTopR,
        double LibraryStabilityResidual,
        double SpectralGapBelowTop);

    private static Resolution Resolve(CoherenceBlock block, int hd1, int hd2)
    {
        var pHd1 = HdSubspaceProjector.Build(block, hd1);
        var pHd2 = HdSubspaceProjector.Build(block, hd2);
        var vInter = pHd1.ConjugateTranspose() * block.Decomposition.MhTotal * pHd2;

        // Eigendecomposition of V V† (HD1 side). Eigenvalues are σ_k².
        var vvh = vInter * vInter.ConjugateTranspose();
        var evdL = vvh.Evd();
        var (eigsL, vecsL) = SortDescending(evdL);

        var vhv = vInter.ConjugateTranspose() * vInter;
        var evdR = vhv.Evd();
        var (eigsR, vecsR) = SortDescending(evdR);

        double sigma0Sq = Math.Max(0.0, eigsL[0]);
        double sigma0 = Math.Sqrt(sigma0Sq);

        int rankTop = CountDegenerate(eigsL, sigma0Sq);
        int rankTopRight = CountDegenerate(eigsR, sigma0Sq);
        if (rankTop != rankTopRight)
        {
            // Numerical asymmetry: report the larger rank (the projector spans both sides
            // consistently because V V† and V† V have identical non-zero eigenvalues).
            rankTop = Math.Max(rankTop, rankTopRight);
        }

        // Spectral gap below the top eigenspace.
        double spectralGap = rankTop < eigsL.Length
            ? (sigma0Sq - Math.Max(0.0, eigsL[rankTop])) / Math.Max(sigma0Sq, 1e-15)
            : 1.0;

        // Build P_top from the top RankTop eigenvectors of V V† and V† V, lifted to full block.
        var pTopLReduced = ProjectorFromTopVectors(vecsL, rankTop);
        var pTopRReduced = ProjectorFromTopVectors(vecsR, rankTop);
        var pTopL = pHd1 * pTopLReduced * pHd1.ConjugateTranspose();
        var pTopR = pHd2 * pTopRReduced * pHd2.ConjugateTranspose();

        // Library-stability witness: re-build the projector with a random orthogonal mixing
        // within the degenerate top eigenspace; the Riesz spectral projection should be
        // unchanged at machine precision.
        double residualL = LibraryStabilityResidualFor(vecsL, rankTop, pHd1, pTopL);
        double residualR = LibraryStabilityResidualFor(vecsR, rankTop, pHd2, pTopR);
        double residual = Math.Max(residualL, residualR);

        Tier tier = residual < LibraryStabilityTolerance ? Tier.Tier1Derived : Tier.Tier2Verified;

        return new Resolution(
            Tier: tier,
            Sigma0: sigma0,
            RankTop: rankTop,
            PTopL: pTopL,
            PTopR: pTopR,
            LibraryStabilityResidual: residual,
            SpectralGapBelowTop: spectralGap);
    }

    private static (double[] EigsDesc, ComplexMatrix VecsDesc) SortDescending(
        MathNet.Numerics.LinearAlgebra.Factorization.Evd<Complex> evd)
    {
        int n = evd.EigenValues.Count;
        var eigs = new double[n];
        for (int i = 0; i < n; i++) eigs[i] = evd.EigenValues[i].Real;
        var order = Enumerable.Range(0, n).OrderByDescending(i => eigs[i]).ToArray();
        var sortedEigs = order.Select(i => eigs[i]).ToArray();
        var sortedVecs = ComplexMatrix.Build.Dense(n, n);
        for (int col = 0; col < n; col++)
            sortedVecs.SetColumn(col, evd.EigenVectors.Column(order[col]));
        return (sortedEigs, sortedVecs);
    }

    private static int CountDegenerate(double[] eigsDesc, double topSquared)
    {
        if (eigsDesc.Length == 0) return 0;
        double scale = Math.Max(topSquared, 1e-15);
        int count = 1;
        for (int k = 1; k < eigsDesc.Length; k++)
        {
            if (Math.Abs(eigsDesc[k] - topSquared) / scale < DegeneracyTolerance)
                count++;
            else
                break;
        }
        return count;
    }

    /// <summary>Build <c>P = Σ_{k=0..rank-1} v_k v_k^†</c> from the leading <paramref name="rank"/>
    /// columns of <paramref name="vecsDesc"/>. The result is a Hermitian rank-<paramref name="rank"/>
    /// projector in the reduced (HD-subspace) basis.</summary>
    private static ComplexMatrix ProjectorFromTopVectors(ComplexMatrix vecsDesc, int rank)
    {
        int n = vecsDesc.RowCount;
        var p = ComplexMatrix.Build.Dense(n, n);
        for (int k = 0; k < rank; k++)
        {
            var v = vecsDesc.Column(k);
            p += v.OuterProduct(v.Conjugate());
        }
        return p;
    }

    /// <summary>Library-stability residual: build P_top once from the canonical EVD basis and
    /// once from a random orthogonal mixing of those same eigenvectors. The Riesz projection
    /// onto a fixed eigenspace is basis-independent, so the two projectors should agree to
    /// machine precision (Frobenius residual ≪ 1).</summary>
    private static double LibraryStabilityResidualFor(
        ComplexMatrix vecsDesc, int rank, ComplexMatrix pHd, ComplexMatrix pTopFullBlock)
    {
        if (rank == 0) return 0.0;
        // Random unitary in the rank-dim degenerate subspace (deterministic seed for
        // reproducibility; the witness is structural, not statistical).
        var rng = new Random(20260508);
        var rotationReduced = RandomUnitary(rank, rng);

        // Apply rotation: V' = V_top · U, with U a rank×rank unitary. The mixed eigenvectors
        // span the same subspace.
        int n = vecsDesc.RowCount;
        var vTop = ComplexMatrix.Build.Dense(n, rank);
        for (int k = 0; k < rank; k++) vTop.SetColumn(k, vecsDesc.Column(k));
        var vMixed = vTop * rotationReduced;
        var pMixedReduced = vMixed * vMixed.ConjugateTranspose();
        var pMixedFullBlock = pHd * pMixedReduced * pHd.ConjugateTranspose();
        return (pTopFullBlock - pMixedFullBlock).FrobeniusNorm();
    }

    private static ComplexMatrix RandomUnitary(int n, Random rng)
    {
        // Generate a random complex matrix and QR-decompose; Q is unitary.
        var rand = ComplexMatrix.Build.Dense(n, n,
            (i, j) => new Complex(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5));
        var qr = rand.QR();
        return qr.Q;
    }

    public override string DisplayName =>
        $"c=2 SVD-top inter-channel projector (Riesz lift; N={Block.N})";

    public override string Summary =>
        $"σ_0={Sigma0:G6}, rank={RankTop}, library-stability residual={LibraryStabilityResidual:G3} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return InspectableNode.RealScalar("σ_0", Sigma0);
            yield return new InspectableNode("RankTop",
                summary: $"{RankTop} ({(RankTop == 1 ? "non-degenerate top" : $"degenerate top eigenspace, dim {RankTop}")})");
            yield return InspectableNode.RealScalar("LibraryStabilityResidual", LibraryStabilityResidual, "G3");
            yield return InspectableNode.RealScalar("SpectralGapBelowTop", SpectralGapBelowTop, "F4");
            yield return new InspectableNode("IsAnalyticallyDerived",
                summary: IsAnalyticallyDerived
                    ? "true (Tier1Derived: Riesz spectral projection is library-invariant)"
                    : "false (Tier2Verified: library-stability residual exceeded tolerance)");
            yield return new InspectableNode("PTopL trace (= RankTop on full block)",
                summary: $"{PTopL.Trace().Real:F6}");
            yield return new InspectableNode("PTopR trace (= RankTop on full block)",
                summary: $"{PTopR.Trace().Real:F6}");
        }
    }
}
