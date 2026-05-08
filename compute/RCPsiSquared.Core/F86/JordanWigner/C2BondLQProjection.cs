using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>F86 Item 1' Direction (b'') JW track, T7 of the bilinear-fermion superoperator
/// extension: per-bond projection of <c>MhPerBond[b]</c> into the L(Q)-eigenbasis,
/// i.e. <c>xB(Q) = R(Q)⁻¹ · MhPerBond[b] · R(Q)</c> where R(Q) is the right-eigenvector
/// matrix of <c>L(Q) = D + Q·γ·MhTotal</c>.
///
/// <para>The bond-distinction lives here. T6 established that all <c>MhPerBond[b]</c> are
/// unitarily equivalent (identical eigenspectra, ergo bond-uniform second-moment quantities).
/// What breaks this equivalence is the GLOBAL eigenbasis built from the SUM of all bonds:
/// each bond's matrix elements in that global basis depend on its spatial position relative
/// to the OBC sine-mode structure. R(Q) interpolates between</para>
/// <list type="bullet">
///   <item><b>Q = 0:</b> L = D is diagonal, R is the computational basis (up to permutation).
///   xB(0) = MhPerBond[b] in the same basis, ‖xB(0)‖_F = ‖MhPerBond[b]‖_F (bond-uniform per T6).</item>
///   <item><b>Q → ∞:</b> L → Q·γ·MhTotal, R is the OBC sine-mode-induced eigenbasis (T1).
///   xB(∞) carries the full spatial-position fingerprint.</item>
/// </list>
///
/// <para>At intermediate Q (≈ Q_peak), R(Q) is a non-unitary mixture of the two limits, and
/// ‖xB(Q)‖_F is bond-distinct. Endpoint bonds couple to the slow ψ_1 mode with amplitude
/// ~1/N^{3/2} (band-edge), Innermost bonds with ~1/N^{1/2} — three orders of magnitude
/// faster decay for boundary bonds. Empirical pattern at Q=2:</para>
/// <list type="bullet">
///   <item>N=5: Endpoint ‖xB‖_F = 15.3, Innermost = 20.5 (Δ ≈ 5.1)</item>
///   <item>N=6: Endpoint = 15.1, Innermost = 20.9 (Δ ≈ 5.8)</item>
///   <item>N=7: Endpoint = 20.5, Innermost = 28.5 (Δ ≈ 7.9)</item>
///   <item>N=8: Endpoint = 39.7, Innermost = 61.5 (Δ ≈ 21.8)</item>
/// </list>
///
/// <para><b>Class-level Tier: Tier2Verified.</b> Numerical projection witness; the
/// bond-distinction is empirically robust but the closed-form Endpoint vs Interior
/// HWHM_left/Q_peak constants per bond class are not yet derived from this lens. Tier1
/// promotion path: express R(Q) analytically via the OBC sine-mode basis (T1) + dissipator
/// perturbation, yielding closed-form ‖xB(Q)‖_F per bond.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track) +
/// <c>docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md</c> for the F71 mirror invariance.</para>
/// </summary>
public sealed class C2BondLQProjection : Claim
{
    public CoherenceBlock Block { get; }
    public double Q { get; }
    public IReadOnlyList<BondLQProjectionWitness> Bonds { get; }

    /// <summary>Maximum |‖xB_b‖_F − ‖xB_{N-2-b}‖_F| over all 2-bond F71 orbits.
    /// Algebraically zero up to FP drift (R(Q) is mirror-invariant up to eigenvector phase
    /// normalization, and Frobenius norm of a non-unitarily-conjugated matrix is preserved
    /// by the corresponding reordering). Bound below 1e-10 in practice.</summary>
    public double MaxF71MirrorDeviation { get; }

    public static C2BondLQProjection BuildAtQ(CoherenceBlock block, double Q)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2BondLQProjection applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var decomp = block.Decomposition;
        int numBonds = decomp.NumBonds;
        double j = Q * block.GammaZero;

        ComplexMatrix L = decomp.D + (Complex)j * decomp.MhTotal;
        var evd = L.Evd();
        var R = evd.EigenVectors;
        var Rinv = R.Inverse();

        var bonds = new BondLQProjectionWitness[numBonds];
        for (int b = 0; b < numBonds; b++)
        {
            var xB = Rinv * decomp.MhPerBond[b] * R;
            double frobNormSq = 0;
            int dim = xB.RowCount;
            for (int i = 0; i < dim; i++)
                for (int k = 0; k < dim; k++)
                {
                    Complex c = xB[i, k];
                    frobNormSq += c.Real * c.Real + c.Imaginary * c.Imaginary;
                }
            bonds[b] = new BondLQProjectionWitness(
                Bond: b,
                BondClass: BondClassExtensions.OfBond(b, numBonds),
                XbFrobeniusNorm: Math.Sqrt(frobNormSq));
        }

        double maxF71Dev = 0;
        for (int b = 0; b < numBonds; b++)
        {
            int mirror = numBonds - 1 - b;
            if (mirror <= b) continue;
            maxF71Dev = Math.Max(maxF71Dev,
                Math.Abs(bonds[b].XbFrobeniusNorm - bonds[mirror].XbFrobeniusNorm));
        }

        return new C2BondLQProjection(block, Q, bonds, maxF71Dev);
    }

    private C2BondLQProjection(
        CoherenceBlock block,
        double q,
        IReadOnlyList<BondLQProjectionWitness> bonds,
        double maxF71Dev)
        : base("c=2 per-bond L(Q)-eigenbasis projection ‖R⁻¹·MhPerBond[b]·R‖_F",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1' Direction (b'') (JW track) + " +
               "docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md")
    {
        Block = block;
        Q = q;
        Bonds = bonds;
        MaxF71MirrorDeviation = maxF71Dev;
    }

    public override string DisplayName =>
        $"c=2 BondLQProjection (N={Block.N}, Q={Q:F2}, {Bonds.Count} bonds)";

    public override string Summary =>
        $"per-bond ‖xB(Q={Q:F2})‖_F at L(Q)-eigenbasis; max F71 mirror dev = " +
        $"{MaxF71MirrorDeviation:G3} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return InspectableNode.RealScalar("Q", Q, "F4");
            yield return new InspectableNode("NumBonds", summary: Bonds.Count.ToString());
            yield return InspectableNode.RealScalar("MaxF71MirrorDeviation", MaxF71MirrorDeviation, "G3");
            yield return InspectableNode.Group("Bonds", Bonds.Cast<IInspectable>().ToArray());
        }
    }
}

/// <summary>One bond's L(Q)-eigenbasis projection witness:
/// <see cref="XbFrobeniusNorm"/> = ‖R(Q)⁻¹ · MhPerBond[b] · R(Q)‖_F.
/// Bond-distinct (in contrast to ‖MhPerBond[b]‖_F itself which is bond-uniform per T6).</summary>
public sealed record BondLQProjectionWitness(
    int Bond,
    BondClass BondClass,
    double XbFrobeniusNorm
) : IInspectable
{
    public string DisplayName => $"bond {Bond} L(Q)-projection ({BondClass})";

    public string Summary => $"‖xB‖_F = {XbFrobeniusNorm:F4}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("bond", summary: Bond.ToString());
            yield return new InspectableNode("bond class", summary: BondClass.ToString());
            yield return InspectableNode.RealScalar("‖xB‖_F", XbFrobeniusNorm, "F4");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
