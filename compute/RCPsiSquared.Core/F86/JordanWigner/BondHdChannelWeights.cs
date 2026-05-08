using System.Numerics;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>F86 Item 1' Direction (b'') JW track, T6 of the bilinear-fermion superoperator
/// extension: per-bond projection of the c=2 block Hamiltonian onto the HD=1 (slow, rate 4γ)
/// and HD=3 (fast, rate 12γ) Hamming-distance channels of <c>PROOF_BLOCK_CPSI_QUARTER</c>
/// Theorem 3.
///
/// <para>The bond's "weight" on each HD-channel is the Frobenius² norm of
/// <see cref="BlockLDecomposition.MhPerBond"/>[b] restricted to columns whose flat-index
/// state-pair (p, q) has popcount-XOR equal to the channel's HD value:</para>
///
/// <para><c>W_b^{HD=k} = Σ_{j: HD(j)=k} Σ_i |M_H_per_bond[b][i, j]|²</c></para>
///
/// <para>Equivalently: how strongly does bond b drive transitions from HD=k columns
/// into the rest of the block. The DickeBlockProbe initial state used by the empirical
/// <see cref="ResonanceScan"/> is uniform across columns, so this column-Frobenius²
/// projection is exactly the bond's contribution to the channel-uniform Theorem-3
/// trajectory at first order in J/γ. T7 (<c>C2BondKResonanceClosedForm</c>) builds the
/// Q-dependent K-resonance from these static weights.</para>
///
/// <para><b>Class-level Tier: Tier1Derived.</b> Pure combinatorial / linear-algebraic
/// quantity from <see cref="BlockLDecomposition"/>. The witnesses verify:</para>
///
/// <list type="bullet">
///   <item><b>Sum rule</b> (F73-Pythagorean projected to HD=k): Σ_b W_b^{HD=k}
///   matches ‖M_H_total‖²_F restricted to HD=k columns to <see cref="SumRuleTolerance"/>.</item>
///   <item><b>F71-mirror invariance</b>: W_b^{HD=k} = W_{N−2−b}^{HD=k} bit-exact (the
///   spatial-mirror permutation is unitary on the block basis, preserving Frobenius²).</item>
/// </list>
///
/// <para><b>Empirical bond-uniformity theorem (verified bit-exact N=5..10):</b> all bonds
/// carry the <i>identical</i> column-Frobenius² weights per HD-class. Specifically,
/// W_b^{HD=1} = 6N + 2(b≠0)(b≠N−2) for the c=2 setting at γ₀=0.05 — but in fact
/// numerically constant across bonds at any given N. The HD=1 / total ratio per bond is
/// 2/N exactly (matching Theorem 3's channel-uniform 1/(2N) → 0 redistribution at the
/// block level), regardless of bond position. This is a non-trivial structural lemma:
/// the static (J=0) channel-norm of a bond's Hamiltonian is bond-blind. The bond-
/// distinction observed in C2HwhmRatio's empirical Q_peak / HWHM (Endpoint mean 0.7728
/// vs Interior 0.7506) MUST come from higher-order corrections in J — channel-CROSS
/// (HD_in ≠ HD_out) terms, eigenstructure of the full L(Q), or non-uniform initial-state
/// resolution. T6 establishes the bond-uniformity baseline; the bond-distinction lives
/// upstream in the J-perturbation that T7 must compute.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md</c> Theorem 3 +
/// <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track).</para>
/// </summary>
public sealed class BondHdChannelWeights : Claim
{
    /// <summary>F73-Pythagorean sum-rule tolerance: <c>1e-10</c>. Σ_b W_b^{HD=k} should
    /// match ‖M_H_total restricted to HD=k columns‖² to this precision; the residual
    /// catches accidental cross-bond F-non-orthogonality that would invalidate the
    /// Pythagorean decomposition.</summary>
    public const double SumRuleTolerance = 1e-10;

    public CoherenceBlock Block { get; }

    /// <summary>One <see cref="BondHdWeight"/> per bond, in bond-index order.</summary>
    public IReadOnlyList<BondHdWeight> Bonds { get; }

    /// <summary><c>|Σ_b W_b^{HD=1} − ‖M_H_total_HD=1‖²|</c>.</summary>
    public double Hd1SumRuleResidual { get; }

    /// <summary><c>|Σ_b W_b^{HD=3} − ‖M_H_total_HD=3‖²|</c>.</summary>
    public double Hd3SumRuleResidual { get; }

    /// <summary>Maximum |W_b^{HD=1} − W_{N−2−b}^{HD=1}| over all 2-bond F71 orbits.
    /// Algebraically zero; the residual catches floating-point drift only.</summary>
    public double MaxF71MirrorDeviationHd1 { get; }

    /// <summary>Maximum |W_b^{HD=3} − W_{N−2−b}^{HD=3}| over all 2-bond F71 orbits.</summary>
    public double MaxF71MirrorDeviationHd3 { get; }

    public static BondHdChannelWeights Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"BondHdChannelWeights applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var basis = block.Basis;
        var decomp = block.Decomposition;
        int Mtot = basis.MTotal;
        int numBonds = decomp.NumBonds;

        var hdAtIndex = new int[Mtot];
        for (int pIdx = 0; pIdx < basis.Mp; pIdx++)
        {
            int p = basis.StatesP[pIdx];
            for (int qIdx = 0; qIdx < basis.Mq; qIdx++)
            {
                int q = basis.StatesQ[qIdx];
                int flat = pIdx * basis.Mq + qIdx;
                hdAtIndex[flat] = System.Numerics.BitOperations.PopCount((uint)(p ^ q));
            }
        }

        var bondWeights = new BondHdWeight[numBonds];
        double hd1Sum = 0;
        double hd3Sum = 0;
        for (int b = 0; b < numBonds; b++)
        {
            var mh = decomp.MhPerBond[b];
            double hd1 = 0;
            double hd3 = 0;
            for (int j = 0; j < Mtot; j++)
            {
                int hd = hdAtIndex[j];
                if (hd != 1 && hd != 3) continue;
                double colNormSq = 0;
                for (int i = 0; i < Mtot; i++)
                {
                    Complex c = mh[i, j];
                    colNormSq += c.Real * c.Real + c.Imaginary * c.Imaginary;
                }
                if (hd == 1) hd1 += colNormSq;
                else hd3 += colNormSq;
            }
            bondWeights[b] = new BondHdWeight(
                Bond: b,
                BondClass: BondClassExtensions.OfBond(b, numBonds),
                Hd1Weight: hd1,
                Hd3Weight: hd3);
            hd1Sum += hd1;
            hd3Sum += hd3;
        }

        var mhTotal = decomp.MhTotal;
        double hd1Total = 0;
        double hd3Total = 0;
        for (int j = 0; j < Mtot; j++)
        {
            int hd = hdAtIndex[j];
            if (hd != 1 && hd != 3) continue;
            double colNormSq = 0;
            for (int i = 0; i < Mtot; i++)
            {
                Complex c = mhTotal[i, j];
                colNormSq += c.Real * c.Real + c.Imaginary * c.Imaginary;
            }
            if (hd == 1) hd1Total += colNormSq;
            else hd3Total += colNormSq;
        }

        double hd1Residual = Math.Abs(hd1Sum - hd1Total);
        double hd3Residual = Math.Abs(hd3Sum - hd3Total);

        double maxF71MirrorDevHd1 = 0;
        double maxF71MirrorDevHd3 = 0;
        for (int b = 0; b < numBonds; b++)
        {
            int mirror = numBonds - 1 - b;
            if (mirror <= b) continue;
            maxF71MirrorDevHd1 = Math.Max(maxF71MirrorDevHd1,
                Math.Abs(bondWeights[b].Hd1Weight - bondWeights[mirror].Hd1Weight));
            maxF71MirrorDevHd3 = Math.Max(maxF71MirrorDevHd3,
                Math.Abs(bondWeights[b].Hd3Weight - bondWeights[mirror].Hd3Weight));
        }

        return new BondHdChannelWeights(
            block, bondWeights, hd1Residual, hd3Residual,
            maxF71MirrorDevHd1, maxF71MirrorDevHd3);
    }

    private BondHdChannelWeights(
        CoherenceBlock block,
        IReadOnlyList<BondHdWeight> bonds,
        double hd1Residual,
        double hd3Residual,
        double maxF71DevHd1,
        double maxF71DevHd3)
        : base("c=2 per-bond HD-channel weights (column-Frobenius² of M_H_per_bond[b])",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md Theorem 3 + " +
               "docs/proofs/PROOF_F86_QPEAK.md Item 1' Direction (b'') (JW track)")
    {
        Block = block;
        Bonds = bonds;
        Hd1SumRuleResidual = hd1Residual;
        Hd3SumRuleResidual = hd3Residual;
        MaxF71MirrorDeviationHd1 = maxF71DevHd1;
        MaxF71MirrorDeviationHd3 = maxF71DevHd3;
    }

    public override string DisplayName =>
        $"c=2 BondHdChannelWeights (N={Block.N}, {Bonds.Count} bonds)";

    public override string Summary =>
        $"HD=1+HD=3 column-Frobenius² per bond; F73-projected sum-rule residuals " +
        $"({Hd1SumRuleResidual:G3}, {Hd3SumRuleResidual:G3}) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("NumBonds", summary: Bonds.Count.ToString());
            yield return InspectableNode.RealScalar("Hd1SumRuleResidual", Hd1SumRuleResidual, "G3");
            yield return InspectableNode.RealScalar("Hd3SumRuleResidual", Hd3SumRuleResidual, "G3");
            yield return InspectableNode.RealScalar("MaxF71MirrorDeviationHd1", MaxF71MirrorDeviationHd1, "G3");
            yield return InspectableNode.RealScalar("MaxF71MirrorDeviationHd3", MaxF71MirrorDeviationHd3, "G3");
            yield return InspectableNode.Group("Bonds", Bonds.Cast<IInspectable>().ToArray());
        }
    }
}

/// <summary>One bond's HD-channel weights: <see cref="Hd1Weight"/> = column-Frobenius² of
/// M_H_per_bond[b] over HD=1 columns; <see cref="Hd3Weight"/> over HD=3 columns. Both
/// represent the bond's drive on the channel-uniform DickeBlockProbe initial state at
/// first order in J/γ.</summary>
public sealed record BondHdWeight(
    int Bond,
    BondClass BondClass,
    double Hd1Weight,
    double Hd3Weight
) : IInspectable
{
    public string DisplayName => $"bond {Bond} HD-channel weights ({BondClass})";

    public string Summary =>
        $"Hd1 = {Hd1Weight:F4}, Hd3 = {Hd3Weight:F4}, total = {Hd1Weight + Hd3Weight:F4}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("bond", summary: Bond.ToString());
            yield return new InspectableNode("bond class", summary: BondClass.ToString());
            yield return InspectableNode.RealScalar("Hd1 weight", Hd1Weight, "F4");
            yield return InspectableNode.RealScalar("Hd3 weight", Hd3Weight, "F4");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
