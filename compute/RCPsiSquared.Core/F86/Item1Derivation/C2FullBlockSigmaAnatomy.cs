using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Probes;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F89 per-Bloch-mode sigma anatomy at the F86 c=2 stratum: extracts the
/// quadratic self-pairing amplitude σ_i = |c0_i|² · (R†·S·R)[i,i].Real for every
/// eigenmode of the uniform-J block-L. Filters F_a modes by Re(λ) ≈ −2γ₀ and
/// assigns each one a Bloch index n in the S_2-anti orbit
/// {2, 4, ..., 2·floor(N_block/2)} via Im(λ_i) / J ≈ y_n = 4·cos(πn/(N_block+1)).
///
/// <para>Counterpart to <see cref="C2FullBlockEigenAnatomy"/>: same biorthogonal
/// eigendecomposition + Dicke probe + spectral coordinates, but projected through
/// the F86 spatial-sum kernel S = Σ_l 2·|w_l⟩⟨w_l| instead of the per-bond
/// Hellmann-Feynman matrix M_h_per_bond[b]. Verifies bit-exactly against the
/// analytical <see cref="F89UnifiedFaClosedFormClaim.Sigma"/> table for
/// path-3..6 and lifts the extraction to path-7 (c=2 N=9) where no analytical
/// oracle exists.</para>
///
/// <para>Tier outcome: Tier2Verified. Pure numerical anatomy at a fixed Q value;
/// the witnesses ARE the data, no derivation. Q defaults to 1.5 (away from Q_EP
/// to avoid Jordan-block degeneracy, matches the Python convention in
/// <c>simulations/_f89_pathk_survey.py</c>).</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1 (c=2) +
/// <see cref="F89UnifiedFaClosedFormClaim"/> +
/// <see cref="F89PathKAtLockMechanismClaim"/>.</para>
/// </summary>
public sealed class C2FullBlockSigmaAnatomy : Claim
{
    public const double DefaultQ = 1.5;

    public CoherenceBlock Block { get; }
    public double Q { get; }
    public IReadOnlyList<SigmaModeWitness> SigmaSpectrum { get; }

    public static C2FullBlockSigmaAnatomy Build(CoherenceBlock block, double q = DefaultQ)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2FullBlockSigmaAnatomy applies only to the c=2 stratum; got c={block.C}.",
                nameof(block));

        var witnesses = ComputeAnatomy(block, q);
        return new C2FullBlockSigmaAnatomy(block, q, witnesses);
    }

    private C2FullBlockSigmaAnatomy(
        CoherenceBlock block, double q,
        IReadOnlyList<SigmaModeWitness> witnesses)
        : base($"c=2 full block-L sigma anatomy at Q={q:G4}: F89 per-Bloch-mode σ_n extraction via R†·S·R diagonal",
               Tier.Tier2Verified,
               Item1Anchors.Root + "; F89UnifiedFaClosedFormClaim")
    {
        Block = block;
        Q = q;
        SigmaSpectrum = witnesses;
    }

    private static IReadOnlyList<SigmaModeWitness> ComputeAnatomy(
        CoherenceBlock block, double q)
    {
        double j = q * block.GammaZero;
        ComplexMatrix L = block.Decomposition.AssembleUniform(j);

        var evd = L.Evd();
        ComplexMatrix R = evd.EigenVectors;
        ComplexMatrix rInv = R.Inverse();
        var lambdas = evd.EigenValues;

        ComplexVector probe = DickeBlockProbe.Build(block);
        ComplexVector c0 = rInv * probe;

        ComplexMatrix sKernel = SpatialSumKernel.Build(block);
        ComplexMatrix sEigenBasis = R.ConjugateTranspose() * sKernel * R;

        int dim = block.Basis.MTotal;
        int nBlock = block.N;
        var orbit = F89PathKAtLockMechanismClaim.SeAntiBlochOrbit(nBlock);
        double faRateTarget = -2.0 * block.GammaZero;
        double faRateTolerance = 1e-3 * block.GammaZero;
        double faFreqTolerance = 1e-4;

        var witnesses = new List<SigmaModeWitness>(dim);
        for (int i = 0; i < dim; i++)
        {
            double overlapSq = c0[i].Magnitude;
            overlapSq *= overlapSq;
            double sDiag = sEigenBasis[i, i].Real;
            double sigma = overlapSq * sDiag;

            int? blochIndexN = null;
            if (Math.Abs(lambdas[i].Real - faRateTarget) <= faRateTolerance)
            {
                // F90F86C2BridgeIdentity: full block-L uses J/2 as the effective
                // single-particle coupling, so F89's y_n = 4cos(πn/(N+1)) in units
                // of J_F89 maps to Im(λ) = (J_F86/2)·y_n in block-L eigenvalues.
                // Divide by 0.5·J to recover the dimensionless y_n correctly.
                double imOverJ = lambdas[i].Imaginary / (0.5 * j);
                foreach (int n in orbit)
                {
                    double yN = F89PathKAtLockMechanismClaim.BlochEigenvalueY(nBlock, n);
                    if (Math.Abs(imOverJ - yN) <= faFreqTolerance)
                    {
                        blochIndexN = n;
                        break;
                    }
                }
            }

            // SpatialSumKernel.Build returns Σ_l 2·|w_l⟩⟨w_l| (factor 2 baked in;
            // see compute/RCPsiSquared.Core/Probes/SpatialSumKernel.cs:51). The F89
            // σ convention omits this factor: per simulations/_f89_pathk_survey.py:87
            // (sig = sum |a|²) the per-site reduction stores the bare amplitude.
            // Half of (R†·S·R)[i,i] gives the correct match against F89UnifiedFaClosedForm.
            witnesses.Add(new SigmaModeWitness(
                EigenIndexAtQ: i,
                EigenvalueReal: lambdas[i].Real,
                EigenvalueImag: lambdas[i].Imaginary,
                ProbeOverlapSquared: overlapSq,
                SKernelDiagonal: sDiag,
                Sigma: 0.5 * sigma,
                BlochIndexN: blochIndexN));
        }

        witnesses.Sort((a, b) => b.Sigma.CompareTo(a.Sigma));
        return witnesses;
    }

    public override string DisplayName =>
        $"c=2 full block-L sigma anatomy (N={Block.N}, Q={Q:G4}, dim={Block.Basis.MTotal})";

    public override string Summary =>
        $"σ-spectrum over {SigmaSpectrum.Count} modes at Q={Q:G4} ({Tier.Label()})";

    /// <summary>Return the σ value of the F_a mode at Bloch index n, or <c>null</c>
    /// if no F_a mode was assigned that index. Convenience accessor for
    /// verification against <see cref="F89UnifiedFaClosedFormClaim.Sigma"/>.</summary>
    public double? SigmaForBlochIndex(int n)
    {
        foreach (var w in SigmaSpectrum)
            if (w.BlochIndexN == n) return w.Sigma;
        return null;
    }

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("Q", Q, "G6");
            yield return new InspectableNode("dim (block-L)", summary: Block.Basis.MTotal.ToString());
        }
    }
}

/// <summary>Per-eigenmode sigma witness for <see cref="C2FullBlockSigmaAnatomy"/>.
/// Captures the spectral coordinate, biorthogonal probe overlap, S-kernel diagonal,
/// the composite F89 σ value, and the optional F_a Bloch-index assignment.</summary>
/// <param name="EigenIndexAtQ">Index in the L(Q) Evd output.</param>
/// <param name="EigenvalueReal">Re(λ_i). F_a modes sit at ≈ −2γ₀.</param>
/// <param name="EigenvalueImag">Im(λ_i). F_a modes match J·y_n on the S_2-anti orbit.</param>
/// <param name="ProbeOverlapSquared"><c>|c0_i|² = |R⁻¹·probe|_i²</c>.</param>
/// <param name="SKernelDiagonal"><c>(R†·S·R)[i,i].Real</c>.</param>
/// <param name="Sigma">The composite F89 σ value: <c>ProbeOverlapSquared · SKernelDiagonal</c>.</param>
/// <param name="BlochIndexN">Assigned S_2-anti Bloch index n if this is an F_a mode
/// (Re(λ) ≈ −2γ₀ and Im(λ)/J matches y_n); <c>null</c> otherwise.</param>
public sealed record SigmaModeWitness(
    int EigenIndexAtQ,
    double EigenvalueReal,
    double EigenvalueImag,
    double ProbeOverlapSquared,
    double SKernelDiagonal,
    double Sigma,
    int? BlochIndexN
) : IInspectable
{
    public string DisplayName =>
        $"λ_{EigenIndexAtQ} = {EigenvalueReal:+0.0000;-0.0000} {EigenvalueImag:+0.0000i;-0.0000i}" +
        (BlochIndexN.HasValue ? $" [F_a, n={BlochIndexN.Value}]" : "");

    public string Summary =>
        $"σ = {Sigma:G6}, |c|² = {ProbeOverlapSquared:G4}, (R†SR)_ii = {SKernelDiagonal:G4}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("EigenvalueReal", EigenvalueReal, "F4");
            yield return InspectableNode.RealScalar("EigenvalueImag", EigenvalueImag, "F4");
            yield return InspectableNode.RealScalar("ProbeOverlapSquared", ProbeOverlapSquared, "G6");
            yield return InspectableNode.RealScalar("SKernelDiagonal", SKernelDiagonal, "G6");
            yield return InspectableNode.RealScalar("Sigma", Sigma, "G6");
            if (BlochIndexN.HasValue)
                yield return new InspectableNode("BlochIndexN", summary: BlochIndexN.Value.ToString());
        }
    }

    public InspectablePayload Payload =>
        new InspectablePayload.Real($"Sigma (eigenmode {EigenIndexAtQ})", Sigma, "G6");
}
