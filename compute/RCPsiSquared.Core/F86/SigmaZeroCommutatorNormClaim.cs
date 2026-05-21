using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 Tier-1-derived: the F86e inter-channel σ_0 IS a commutator operator norm.
///
/// <para><b>Identity (c=2 stratum):</b> for the c=2 coherence block (n=1, popcount-1 ⊗
/// popcount-2) of an N-qubit XY chain under Z-dephasing,</para>
/// <code>σ_0 = ‖[Π_HD1, M_H]‖</code>
/// <para>where σ_0 is the inter-channel SVD-top singular value (the F86e quantity, see
/// <see cref="SigmaZeroChromaticityScaling"/> and <see cref="InterChannelSvd"/>), Π_HD1 is
/// the Hamming-distance-1 subspace projector, M_H is the block Hamiltonian super-operator
/// (<c>block.Decomposition.MhTotal</c>), and ‖·‖ is the operator (top-singular-value) norm.</para>
///
/// <para><b>Why Tier1Derived.</b> On the c=2 block the Hamming distance HD(p, q) takes only
/// values {1, 3} (chromaticity c = 2 ⇒ HD ∈ {1, …, 2c−1} = {1, 3}). Therefore the two
/// HD-subspace projectors are complementary:</para>
/// <code>Π_HD1 + Π_HD3 = I.</code>
/// <para>The F86e inter-channel coupling is V_inter = Π_HD1 · M_H · Π_HD3, so substituting
/// Π_HD3 = I − Π_HD1 gives V_inter = Π_HD1 · M_H · (I − Π_HD1). For any orthogonal
/// projector P and Hermitian M, the standard lemma ‖P·M·(1−P)‖ = ‖[P, M]‖ holds: the
/// commutator [P, M] = P·M·(1−P) − (1−P)·M·P is anti-Hermitian, block-anti-diagonal in the
/// P / (1−P) split, and its two off-diagonal blocks P·M·(1−P) and (1−P)·M·P = (P·M·(1−P))†
/// are adjoints, so [P, M] and P·M·(1−P) share the same set of singular values; hence the
/// top singular values (the operator norms) coincide. Thus σ_0 = ‖V_inter‖ = ‖[Π_HD1, M_H]‖
/// exactly. Verified bit-exact (Python, residual ~1e-15) for N = 4..12; the C# witness here
/// reproduces it through N = 8 with <see cref="Residual"/> &lt; 1e-10.</para>
///
/// <para>The identity is c=2-specific: at c ≥ 3 the HD spectrum is {1, 3, …, 2c−1} with
/// more than two values, so Π_HD1 + Π_HD3 ≠ I and the substitution Π_HD3 = I − Π_HD1 fails.
/// The constructor therefore rejects any block with <c>block.C != 2</c>.</para>
///
/// <para><b>Bloch-basis Hadamard form.</b> In the F89 (SE, DE) Bloch / OBC-sine operator
/// basis M_H is diagonal (multiplication by Δ = E_k − E_{k₁} − E_{k₂}), so the commutator
/// is the Hadamard (Schur) product [Π_HD1, M_H] = Π̃_HD1 ⊙ ΔDiff with ΔDiff[a, b] = Δ_b − Δ_a
/// (verified bit-exact). σ_0 is thus a Schur-multiplier norm. The Δ-ordered commutator is
/// neither Toeplitz (diagonal CV ≈ 0.37, non-vanishing) nor Hankel (anti-diagonal CV grows
/// with N), so the σ_0(N→∞) asymptote is neither a Fourier-symbol supremum nor a Nehari
/// symbol distance; it is a genuine Schur-multiplier-norm constant. This characterises the
/// non-elementarity of the F86e asymptote rather than removing it.</para>
///
/// <para>Reference σ_0 values (γ-independent, V_inter uses only M_H, not the dissipator D):
/// σ_0(N=5) = 2.7650951722, σ_0(N=7) = 2.8284271247 (= 2√2, the
/// <see cref="SigmaZeroChromaticityScaling.SigmaZeroSweetSpotIdentity_C2N7"/> sweet-spot
/// crossing), σ_0(N=8) = 2.8393460323.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1 (c=2): this claim sharpens
/// the F86e σ_0 quantity from "top singular value of a coupling block" to "operator norm of
/// a commutator", a basis-free algebraic characterisation.</para>
/// </summary>
public sealed class SigmaZeroCommutatorNormClaim : Claim
{
    /// <summary>The c=2 coherence block this claim is computed on.</summary>
    public CoherenceBlock Block { get; }

    /// <summary>The F86e inter-channel SVD-top singular value σ_0 = ‖Π_HD1·M_H·Π_HD3‖,
    /// from <see cref="InterChannelSvd.Build"/> with HD₁ = 1, HD₂ = 3.</summary>
    public double Sigma0 { get; }

    /// <summary>The operator norm ‖[Π_HD1, M_H]‖ of the commutator of the (square) HD=1
    /// projector with the block Hamiltonian super-operator. Equal to <see cref="Sigma0"/>
    /// by the c=2 complementarity identity.</summary>
    public double CommutatorNorm { get; }

    /// <summary>|σ_0 − ‖[Π_HD1, M_H]‖|. The identity predicts 0; the live witness lands at
    /// ~1e-15 (machine precision).</summary>
    public double Residual { get; }

    /// <summary>Idempotency residual ‖Π·Π − Π‖ of the square projector built as P·P† from
    /// the rectangular indicator matrix returned by <see cref="HdSubspaceProjector.Build"/>.
    /// Sanity check that the square projector is genuinely a projector; predicts 0.</summary>
    public double ProjectorIdempotencyResidual { get; }

    /// <summary>Public factory: validates the c=2 stratum, computes σ_0 and the commutator
    /// norm, then constructs the instance. Mirrors
    /// <see cref="Item1Derivation.C2InterChannelAnalytical.Build"/>: a factory rather than a
    /// public constructor so the c=2 guard runs before any field is set.</summary>
    /// <exception cref="ArgumentException">Thrown when <paramref name="block"/> is not the
    /// c=2 stratum (the identity needs Π_HD1 + Π_HD3 = I, true only at c = 2).</exception>
    public static SigmaZeroCommutatorNormClaim Build(CoherenceBlock block)
    {
        if (block is null) throw new ArgumentNullException(nameof(block));
        if (block.C != 2)
            throw new ArgumentException(
                $"SigmaZeroCommutatorNormClaim applies only to the c=2 stratum; got c={block.C} " +
                $"(N={block.N}, n={block.LowerPopcount}). The identity σ_0 = ‖[Π_HD1, M_H]‖ relies on " +
                "Π_HD1 + Π_HD3 = I, which holds only when HD ∈ {1, 3}, i.e. only at chromaticity 2.",
                nameof(block));

        double sigma0 = InterChannelSvd.Build(block, hd1: 1, hd2: 3).Sigma0;

        // HdSubspaceProjector.Build returns a RECTANGULAR Mtot × n_HD1 indicator matrix
        // (one nonzero per column), NOT the square projector. The square projector onto the
        // HD=1 subspace is Π_HD1 = P · P†; because P has orthonormal columns, Π is a genuine
        // orthogonal projector (Π = Π† = Π²).
        var pRect = HdSubspaceProjector.Build(block, hammingDistance: 1);
        var piSquare = pRect * pRect.ConjugateTranspose();   // Mtot × Mtot
        var mh = block.Decomposition.MhTotal;

        var comm = piSquare * mh - mh * piSquare;            // [Π_HD1, M_H]
        double commutatorNorm = comm.L2Norm();               // top singular value = operator norm

        // Sanity check: Π must be idempotent. Computed and exposed, not asserted, so the
        // value is visible in the inspection tree; the test class asserts it is ~0.
        var idemResidual = (piSquare * piSquare - piSquare).L2Norm();

        return new SigmaZeroCommutatorNormClaim(block, sigma0, commutatorNorm, idemResidual);
    }

    private SigmaZeroCommutatorNormClaim(
        CoherenceBlock block, double sigma0, double commutatorNorm, double projectorIdempotencyResidual)
        : base("σ_0 = ‖[Π_HD1, M_H]‖: F86e inter-channel singular value as a commutator norm",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1 (c=2): σ_0 = ‖Π_HD1·M_H·(I−Π_HD1)‖ = ‖[Π_HD1, M_H]‖ " +
               "via Π_HD1 + Π_HD3 = I on the c=2 block (HD ∈ {1, 3})")
    {
        Block = block;
        Sigma0 = sigma0;
        CommutatorNorm = commutatorNorm;
        Residual = Math.Abs(sigma0 - commutatorNorm);
        ProjectorIdempotencyResidual = projectorIdempotencyResidual;
    }

    public override string DisplayName =>
        $"σ_0 = ‖[Π_HD1, M_H]‖ (c=2, N={Block.N})";

    public override string Summary =>
        $"σ_0={Sigma0:G8}, ‖[Π_HD1, M_H]‖={CommutatorNorm:G8}, residual={Residual:G3} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("c", summary: Block.C.ToString());
            yield return InspectableNode.RealScalar("σ_0 (InterChannelSvd, HD 1↔3)", Sigma0, "F10");
            yield return InspectableNode.RealScalar("‖[Π_HD1, M_H]‖ (commutator operator norm)", CommutatorNorm, "F10");
            yield return InspectableNode.RealScalar("residual |σ_0 − ‖[Π_HD1, M_H]‖|", Residual, "G3");
            yield return InspectableNode.RealScalar("projector idempotency residual ‖Π·Π − Π‖",
                ProjectorIdempotencyResidual, "G3");
        }
    }
}
