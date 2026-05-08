using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Probes;
using RCPsiSquared.Core.Resonance;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 Item 1 (c=2 stratum), <b>Direction (b'') pair-level anatomy</b>: full Duhamel
/// pair decomposition of K_b(Q, t) at fixed (Q, t). Companion to
/// <see cref="C2FullBlockEigenAnatomy"/> (which only sees the diagonal Duhamel slice). The
/// K-resonance is structurally bilinear in eigenmodes:
/// <code>
///     K_b(Q, t) = 2·Re Σ_{i,j} c̄_i · S_ij · X_b[i,j] · I_ij(t)
/// </code>
/// where <c>c = R⁻¹·probe</c>, <c>S = R†·S_kernel·R</c>, <c>X_b = R⁻¹·M_h_per_bond[b]·R</c>,
/// and the Duhamel kernel
/// <code>
///     I_ij(t) = (e^{λ_j t} − e^{λ_i t}) / (λ_j − λ_i)   for i ≠ j
///     I_ii(t) = t · e^{λ_i t}                            (degenerate limit)
/// </code>
///
/// <para>The diagonal anatomy in <see cref="C2FullBlockEigenAnatomy"/> only counts the
/// <c>i = j</c> contributions (and even then, summed without the S_ij factor); cross-pair
/// contributions can dominate when two slow modes interfere coherently. This primitive
/// computes the full per-pair contribution magnitude and reports K_90_pair / K_99_pair as
/// the decisive read for whether HWHM is governed by a small N-stable set of mode pairs.
/// </para>
///
/// <para>Default <c>(Q, t) = (Q_peak, t_peak)</c>. <c>Q_peak ≈ 2.197 · Q_EP</c> per the
/// universal post-EP location <see cref="C2HwhmRatio.BareDoubledPtfXPeak"/>; this is where
/// HWHM is actually read, not at Q_EP itself. <c>t_peak = 1/(4γ₀)</c> per Statement 1.</para>
///
/// <para><b>Tier outcome: Tier2Verified.</b> Pair-level numerical anatomy at one (Q, t).
/// The witnesses ARE the data. Tier1 promotion would require an analytical identification
/// of the dominant pair set as a function of (c, N).</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Statement 1 +
/// <see cref="C2HwhmRatio.PendingDerivationNote"/> Direction (b'').</para>
/// </summary>
public sealed class C2FullBlockPairAnatomy : Claim
{
    public CoherenceBlock Block { get; }

    /// <summary>The Q value at which the pair anatomy was taken (default Q_peak).</summary>
    public double Q { get; }

    /// <summary>The t value at which the pair anatomy was taken (default t_peak = 1/(4γ₀)).</summary>
    public double T { get; }

    /// <summary>Number of dim² eigenvalue pairs (i, j) that together capture 90 % of
    /// <c>Σ_{ij,b} |contrib_{ij}(b)|</c>. Small + N-stable K_90_pair means HWHM lives in a
    /// low-rank pair truncation; growing K_90_pair means the K-mass spreads coherently
    /// across many pairs.</summary>
    public int K90Pair { get; }

    /// <summary>K_99 analogue at the 99 % cumulative threshold.</summary>
    public int K99Pair { get; }

    /// <summary>Total <c>Σ_{ij,b} |contrib_{ij}(b)|</c>; normalisation reference.</summary>
    public double TotalPairContribution { get; }

    /// <summary>Top <c>TopPairsToShow</c> eigenvalue-pair witnesses by |contribution|,
    /// sorted descending. Truncated for inspect-tree readability; the K_90/K_99 metrics
    /// summarise the long tail.</summary>
    public IReadOnlyList<EigenPairAnatomyWitness> TopPairs { get; }

    public const int TopPairsToShow = 20;

    public static C2FullBlockPairAnatomy Build(
        CoherenceBlock block, double? q = null, double? t = null)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2FullBlockPairAnatomy applies only to the c=2 stratum; got c={block.C}.",
                nameof(block));

        double qEp = new QEpLaw(InterChannelSvd.Build(block, 1, 3).Sigma0).Value;
        double qValue = q ?? C2HwhmRatio.BareDoubledPtfXPeak * qEp;
        double tValue = t ?? EpAlgebra.TPeak(block.GammaZero);

        var (k90, k99, total, topPairs) = ComputePairAnatomy(block, qValue, tValue);
        return new C2FullBlockPairAnatomy(block, qValue, tValue, k90, k99, total, topPairs);
    }

    private C2FullBlockPairAnatomy(
        CoherenceBlock block, double q, double t,
        int k90, int k99, double total,
        IReadOnlyList<EigenPairAnatomyWitness> topPairs)
        : base($"c=2 full block-L pair anatomy at (Q={q:G4}, t={t:G4}) (Direction (b'') pair view)",
               Tier.Tier2Verified,
               Item1Anchors.Root)
    {
        Block = block;
        Q = q;
        T = t;
        K90Pair = k90;
        K99Pair = k99;
        TotalPairContribution = total;
        TopPairs = topPairs;
    }

    private static (int K90, int K99, double Total, IReadOnlyList<EigenPairAnatomyWitness> TopPairs)
        ComputePairAnatomy(CoherenceBlock block, double q, double t)
    {
        double j = q * block.GammaZero;
        ComplexMatrix L = block.Decomposition.D + (Complex)j * block.Decomposition.MhTotal;
        var evd = L.Evd();
        ComplexMatrix R = evd.EigenVectors;
        ComplexMatrix rInv = R.Inverse();
        var lambdas = evd.EigenValues;
        int dim = block.Basis.MTotal;

        ComplexVector probe = DickeBlockProbe.Build(block);
        ComplexVector c0 = rInv * probe;

        ComplexMatrix sKernel = SpatialSumKernel.Build(block);
        // S_eig = R† · S_kernel · R
        ComplexMatrix sEig = R.ConjugateTranspose() * sKernel * R;

        int numBonds = block.NumBonds;
        var xBEig = new ComplexMatrix[numBonds];
        for (int b = 0; b < numBonds; b++)
        {
            xBEig[b] = rInv * block.Decomposition.MhPerBond[b] * R;
        }

        // Pre-compute exp(λ_i · t) for the Duhamel kernel.
        var expLam = new Complex[dim];
        for (int i = 0; i < dim; i++) expLam[i] = Complex.Exp(lambdas[i] * t);

        // I_ij(t) = (e^{λ_j t} − e^{λ_i t}) / (λ_j − λ_i), or t·e^{λ_i t} at degeneracy.
        // Stored row-major in a flat array of size dim·dim for cache efficiency in the
        // pair scan below. Memory at N=8 dim=224: 224² complex = 100K entries ≈ 1.6 MB. Fine.
        // Total summed |contribution| across all (i, j, b).
        double total = 0.0;

        // We collect ALL pairs (i, j) with their summed-over-bonds magnitude, then sort.
        // Memory: dim² doubles ≈ 50K @ N=5..N=8. Trivial.
        int pairCount = dim * dim;
        var pairContrib = new double[pairCount];
        var pairIndex = new (int I, int J)[pairCount];
        int p = 0;

        for (int i = 0; i < dim; i++)
        {
            for (int jj = 0; jj < dim; jj++)
            {
                Complex iIj;
                if (i == jj)
                {
                    iIj = t * expLam[i];
                }
                else
                {
                    Complex denom = lambdas[jj] - lambdas[i];
                    iIj = denom.Magnitude > 1e-12
                        ? (expLam[jj] - expLam[i]) / denom
                        : t * expLam[i];
                }

                Complex cBar = Complex.Conjugate(c0[i]);
                Complex sIj = sEig[i, jj];

                // Sum |2·Re[c̄_i · S_ij · X_b[i,j] · I_ij(t)]| over bonds.
                double summed = 0.0;
                for (int b = 0; b < numBonds; b++)
                {
                    Complex contrib = cBar * sIj * xBEig[b][i, jj] * iIj;
                    summed += Math.Abs(2.0 * contrib.Real);
                }

                pairContrib[p] = summed;
                pairIndex[p] = (i, jj);
                total += summed;
                p++;
            }
        }

        // Sort descending by |contribution|.
        var order = Enumerable.Range(0, pairCount).OrderByDescending(k => pairContrib[k]).ToArray();

        int k90 = CountToCumulativeFraction(order, pairContrib, 0.90, total);
        int k99 = CountToCumulativeFraction(order, pairContrib, 0.99, total);

        // Build top witnesses.
        int topCount = Math.Min(TopPairsToShow, pairCount);
        var top = new EigenPairAnatomyWitness[topCount];
        for (int rank = 0; rank < topCount; rank++)
        {
            int idx = order[rank];
            (int i, int jj) = pairIndex[idx];
            top[rank] = new EigenPairAnatomyWitness(
                IndexI: i,
                IndexJ: jj,
                EigenvalueRealI: lambdas[i].Real,
                EigenvalueRealJ: lambdas[jj].Real,
                ProbeOverlapI: c0[i].Magnitude * c0[i].Magnitude,
                AbsContribution: pairContrib[idx],
                FractionOfTotal: total > 0 ? pairContrib[idx] / total : 0);
        }

        return (k90, k99, total, top);
    }

    private static int CountToCumulativeFraction(
        int[] order, double[] values, double fraction, double total)
    {
        if (total <= 0.0) return 0;
        double target = fraction * total;
        double running = 0.0;
        for (int k = 0; k < order.Length; k++)
        {
            running += values[order[k]];
            if (running >= target) return k + 1;
        }
        return order.Length;
    }

    public override string DisplayName =>
        $"c=2 full block-L pair anatomy (N={Block.N}, Q={Q:G4}, t={T:G4})";

    public override string Summary =>
        $"K_90_pair = {K90Pair} / {Block.Basis.MTotal * Block.Basis.MTotal} pairs, " +
        $"K_99_pair = {K99Pair}, total Σ|contrib| = {TotalPairContribution:G4} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("dim (block-L)", summary: Block.Basis.MTotal.ToString());
            yield return new InspectableNode("dim² (total pairs)",
                summary: (Block.Basis.MTotal * Block.Basis.MTotal).ToString());
            yield return InspectableNode.RealScalar("Q", Q, "G6");
            yield return InspectableNode.RealScalar("t", T, "G6");
            yield return new InspectableNode("K_90_pair", summary: K90Pair.ToString());
            yield return new InspectableNode("K_99_pair", summary: K99Pair.ToString());
            yield return InspectableNode.RealScalar("TotalPairContribution", TotalPairContribution, "G4");
            yield return InspectableNode.Group($"top {TopPairs.Count} pairs by |contribution|",
                TopPairs.Cast<IInspectable>().ToArray());
        }
    }
}

/// <summary>Per-pair anatomy witness for <see cref="C2FullBlockPairAnatomy"/>.
/// Captures the Duhamel pair contribution magnitude summed over bonds, plus the spectral
/// coordinates of both modes.</summary>
public sealed record EigenPairAnatomyWitness(
    int IndexI,
    int IndexJ,
    double EigenvalueRealI,
    double EigenvalueRealJ,
    double ProbeOverlapI,
    double AbsContribution,
    double FractionOfTotal
) : IInspectable
{
    public string DisplayName =>
        $"({IndexI}, {IndexJ})  Re(λ) = ({EigenvalueRealI:+0.000;-0.000}, {EigenvalueRealJ:+0.000;-0.000})";

    public string Summary =>
        $"|contrib| = {AbsContribution:G4} ({FractionOfTotal * 100:F2}% of total), " +
        $"|c_i|² = {ProbeOverlapI:G4}, " +
        $"diagonal: {(IndexI == IndexJ ? "yes" : "no")}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("IndexI", summary: IndexI.ToString());
            yield return new InspectableNode("IndexJ", summary: IndexJ.ToString());
            yield return InspectableNode.RealScalar("EigenvalueRealI", EigenvalueRealI, "F4");
            yield return InspectableNode.RealScalar("EigenvalueRealJ", EigenvalueRealJ, "F4");
            yield return InspectableNode.RealScalar("ProbeOverlapI", ProbeOverlapI, "G4");
            yield return InspectableNode.RealScalar("AbsContribution", AbsContribution, "G4");
            yield return InspectableNode.RealScalar("FractionOfTotal", FractionOfTotal, "F4");
        }
    }

    public InspectablePayload Payload =>
        new InspectablePayload.Real($"|contrib| pair ({IndexI},{IndexJ})", AbsContribution, "G4");
}
