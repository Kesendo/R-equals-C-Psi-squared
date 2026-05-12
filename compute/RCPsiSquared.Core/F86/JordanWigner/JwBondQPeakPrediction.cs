using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>Per-bond Q_peak prediction via the JW-track closed-form architecture: bond's
/// top cluster-pair affinity (T13) selects a 2×2 sub-block; T12's cluster eigenstructure
/// + D's inter-cluster X̃ in W-eigenbasis yield Q_EP per the 2×2 closed form
/// <c>Q_EP^{(λ)} = 2·|X̃[λ, λ]|/(γ·|Δδ|)</c>; Q_peak ≈ 2.197·Q_EP per BareDoubledPtfXPeak.
///
/// <para>Architecture: bond-distinction lives in the SELECTION (which cluster-pair does
/// the bond couple to most strongly via MhPerBond[b]), not in the Q_EP formula itself
/// (which uses D-only X̃[λ,λ] and is bond-independent given the pair). Endpoint-bonds
/// select small-|Δδ| pairs → large Q_EP → large Q_peak; Innermost-bonds select
/// large-|Δδ| pairs → small Q_peak.</para>
///
/// <para>Tier2Verified runtime prediction; expected residual against empirical
/// <see cref="C2HwhmRatio"/> Q_peak is order 15-50% (varies per bond — multi-cluster
/// effects beyond top-pair contribute residual). Open Tier1-Promotion path: closed-form
/// for ‖M^{[b]}_{c1, c2}‖_F² selecting the dominant pair + multi-cluster combination
/// for residual closure.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track) +
/// <c>docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md</c>.</para>
///
/// <para>F90 status (2026-05-11): the F86 c=2 ↔ F89 path-(N−1) bridge identity
/// achieves numerical Tier-1 for Direction (b'') via per-bond Hellmann-Feynman
/// (bit-exact 20/22 bonds at N=5..8). The JW-track primitives in this file remain
/// active as the alternative analytical route toward the closed-form HWHM_left/Q_peak
/// constants; the per-bond numerical answer itself is no longer the open piece.
/// See <c>docs/proofs/PROOF_F90_F86C2_BRIDGE.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F90F86C2BridgeIdentity.cs</c>.</para>
/// </summary>
public sealed class JwBondQPeakPrediction : Claim
{
    public const double BareDoubledPtfXPeak = 2.196910;

    public CoherenceBlock Block { get; }
    public JwBondClusterPairAffinity Affinity { get; }
    public JwClusterDEigenstructure ClusterEigenstructure { get; }
    public C2HwhmRatio EmpiricalAnchor { get; }
    public IReadOnlyList<BondQPeakPrediction> Bonds { get; }

    /// <summary>Maximum |relative residual| = |(Q_peak_pred − Q_peak_empirical)/Q_peak_empirical|
    /// across bonds where empirical Q_peak is finite (skips orbit-escaped bonds).</summary>
    public double MaxRelativeResidual { get; }

    public static JwBondQPeakPrediction Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"JwBondQPeakPrediction applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var affinity = JwBondClusterPairAffinity.Build(block);
        var es = JwClusterDEigenstructure.Build(block);
        var empirical = C2HwhmRatio.Build(block);
        var jw = es.JwBasis;
        var disp = es.DispersionStructure;
        double γ = block.GammaZero;

        var DJw = jw.Uinv * block.Decomposition.D * jw.U;

        var tripleToAlpha = new Dictionary<(int, int, int), int>();
        for (int alpha = 0; alpha < jw.Triples.Count; alpha++)
        {
            var t = jw.Triples[alpha];
            tripleToAlpha[(t.K, t.K1, t.K2)] = alpha;
        }

        // Cache W eigenbases per cluster
        var clusterEigenBases = new Dictionary<int, (double[] λ, ComplexMatrix U)>();

        ComplexMatrix BuildWForCluster(int c)
        {
            int k = disp.Clusters[c].Triples.Count;
            int[] indices = disp.Clusters[c].Triples
                .Select(t => tripleToAlpha[(t.K, t.K1, t.K2)]).ToArray();
            var W = ComplexMatrix.Build.Dense(k, k);
            for (int i = 0; i < k; i++)
                for (int j = 0; j < k; j++)
                    W[i, j] = DJw[indices[i], indices[j]];
            return W;
        }

        (double[] λ, ComplexMatrix U) GetEigenBasis(int c)
        {
            if (!clusterEigenBases.TryGetValue(c, out var cached))
            {
                var W = BuildWForCluster(c);
                var evd = W.Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Hermitian);
                cached = (evd.EigenValues.Select(z => z.Real).ToArray(), evd.EigenVectors);
                clusterEigenBases[c] = cached;
            }
            return cached;
        }

        ComplexMatrix BuildXForPair(int c1, int c2)
        {
            int k1 = disp.Clusters[c1].Triples.Count;
            int k2 = disp.Clusters[c2].Triples.Count;
            int[] idx1 = disp.Clusters[c1].Triples.Select(t => tripleToAlpha[(t.K, t.K1, t.K2)]).ToArray();
            int[] idx2 = disp.Clusters[c2].Triples.Select(t => tripleToAlpha[(t.K, t.K1, t.K2)]).ToArray();
            var X = ComplexMatrix.Build.Dense(k1, k2);
            for (int i = 0; i < k1; i++)
                for (int j = 0; j < k2; j++)
                    X[i, j] = DJw[idx1[i], idx2[j]];
            return X;
        }

        // Per cluster-pair, compute the dominant 2x2 EP Q_EP via the GENERAL formula
        // Q_EP = √(4|x|² − (a−b)²) / (γ·|Δδ|) where x = X̃[i, j], a = λ_c1[i], b = λ_c2[j].
        // Reduces to T12's a=b form when λ_c1[i] = λ_c2[j]. Allows a≠b sub-block contributions.
        // Returns (max |x|, corresponding Q_EP). Cached per cluster-pair.
        var pairBestEpCache = new Dictionary<(int, int), (double xMag, double qEp)>();
        (double xMag, double qEp) GetPairBestEp(int c1, int c2)
        {
            var key = (c1, c2);
            if (pairBestEpCache.TryGetValue(key, out var cached)) return cached;
            var (λ1, U1) = GetEigenBasis(c1);
            var (λ2, U2) = GetEigenBasis(c2);
            var X = BuildXForPair(c1, c2);
            var Xtilde = U1.ConjugateTranspose() * X * U2;
            double absΔδ = Math.Abs(disp.Clusters[c1].Delta - disp.Clusters[c2].Delta);

            double bestX = 0;
            double bestQEp = double.NaN;
            for (int i = 0; i < λ1.Length; i++)
                for (int j = 0; j < λ2.Length; j++)
                {
                    double mag = Xtilde[i, j].Magnitude;
                    if (mag < 1e-12) continue;
                    double aMinusB = λ1[i] - λ2[j];
                    double disc = 4 * mag * mag - aMinusB * aMinusB;
                    if (disc <= 0) continue;  // No real-Q EP from this sub-block
                    if (absΔδ < 1e-8) continue;
                    double qEp = Math.Sqrt(disc) / (γ * absΔδ);
                    // Track entry with max |x| (dominant coupling); ties broken by smaller Q_EP
                    if (mag > bestX || (mag == bestX && (double.IsNaN(bestQEp) || qEp < bestQEp)))
                    {
                        bestX = mag;
                        bestQEp = qEp;
                    }
                }
            cached = (bestX, bestQEp);
            pairBestEpCache[key] = cached;
            return cached;
        }

        var bondPredictions = new BondQPeakPrediction[affinity.Bonds.Count];
        double maxRelResid = 0;
        for (int b = 0; b < affinity.Bonds.Count; b++)
        {
            var bondAff = affinity.Bonds[b];

            // Walk the bond's affinity ranking and pick the first pair where the GENERAL 2x2
            // EP formula yields a real-Q solution (i.e. some |x| with 4|x|² > (a-b)² exists).
            ClusterPairCouplingEntry? selectedPair = null;
            double selectedX = 0;
            double selectedQEp = double.NaN;
            foreach (var candidate in bondAff.RankedPairs)
            {
                var (xMag, qEpCandidate) = GetPairBestEp(candidate.Cluster1Index, candidate.Cluster2Index);
                if (!double.IsNaN(qEpCandidate))
                {
                    selectedPair = candidate;
                    selectedX = xMag;
                    selectedQEp = qEpCandidate;
                    break;
                }
            }

            if (selectedPair is null)
            {
                bondPredictions[b] = new BondQPeakPrediction(
                    Bond: b,
                    BondClass: bondAff.BondClass,
                    TopPairAbsDeltaδ: double.NaN,
                    XTildeMagnitude: double.NaN,
                    QEpPredicted: double.NaN,
                    QPeakPredicted: double.NaN,
                    QPeakEmpirical: empirical.Witnesses[b].QPeak,
                    RelativeResidual: double.NaN);
                continue;
            }

            double absΔδ = selectedPair.AbsoluteDeltaδ;
            double bestX = selectedX;
            double qEp = selectedQEp;
            double qPeakPred = double.IsNaN(qEp) ? double.NaN : BareDoubledPtfXPeak * qEp;
            double qPeakEmp = empirical.Witnesses[b].QPeak;
            double relResid = (double.IsNaN(qPeakPred) || qPeakEmp <= 0)
                ? double.NaN
                : Math.Abs(qPeakPred - qPeakEmp) / qPeakEmp;

            // Only count residual for non-escape bonds (escape bonds have Q_peak == grid upper bound)
            bool escaped = qPeakEmp >= 4.0 - 1e-6;  // canonical Q-grid upper bound is 4.0
            if (!double.IsNaN(relResid) && !escaped)
                maxRelResid = Math.Max(maxRelResid, relResid);

            bondPredictions[b] = new BondQPeakPrediction(
                Bond: b,
                BondClass: bondAff.BondClass,
                TopPairAbsDeltaδ: absΔδ,
                XTildeMagnitude: bestX,
                QEpPredicted: qEp,
                QPeakPredicted: qPeakPred,
                QPeakEmpirical: qPeakEmp,
                RelativeResidual: relResid);
        }

        return new JwBondQPeakPrediction(block, affinity, es, empirical, bondPredictions, maxRelResid);
    }

    private JwBondQPeakPrediction(
        CoherenceBlock block,
        JwBondClusterPairAffinity affinity,
        JwClusterDEigenstructure clusterEigenstructure,
        C2HwhmRatio empiricalAnchor,
        IReadOnlyList<BondQPeakPrediction> bonds,
        double maxRelResid)
        : base("c=2 per-bond Q_peak prediction via JW-track 2x2 closed form on top cluster-pair",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1' Direction (b'') (JW track) + " +
               "docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md")
    {
        Block = block;
        Affinity = affinity;
        ClusterEigenstructure = clusterEigenstructure;
        EmpiricalAnchor = empiricalAnchor;
        Bonds = bonds;
        MaxRelativeResidual = maxRelResid;
    }

    public override string DisplayName =>
        $"c=2 BondQPeakPrediction (N={Block.N}, max-rel-residual = {MaxRelativeResidual:P1})";

    public override string Summary =>
        $"per-bond Q_peak prediction; max |Δ|/Q_peak (non-escape) = {MaxRelativeResidual:P1} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("NumBonds", summary: Bonds.Count.ToString());
            yield return InspectableNode.RealScalar("MaxRelativeResidual", MaxRelativeResidual, "P1");
            yield return InspectableNode.Group("Bonds", Bonds.Cast<IInspectable>().ToArray());
        }
    }
}

/// <summary>One bond's Q_peak prediction: top-pair |Δδ|, |X̃[λ,λ]|, computed Q_EP / Q_peak,
/// empirical Q_peak from <see cref="C2HwhmRatio"/>, and relative residual.</summary>
public sealed record BondQPeakPrediction(
    int Bond,
    BondClass BondClass,
    double TopPairAbsDeltaδ,
    double XTildeMagnitude,
    double QEpPredicted,
    double QPeakPredicted,
    double QPeakEmpirical,
    double RelativeResidual
) : IInspectable
{
    public string DisplayName => $"bond {Bond} Q_peak prediction ({BondClass})";

    public string Summary =>
        $"Q_peak pred = {QPeakPredicted:F4}, emp = {QPeakEmpirical:F4}, rel-residual = {RelativeResidual:P1}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("bond", summary: Bond.ToString());
            yield return new InspectableNode("bond class", summary: BondClass.ToString());
            yield return InspectableNode.RealScalar("|Δδ| (top pair)", TopPairAbsDeltaδ, "F4");
            yield return InspectableNode.RealScalar("|X̃[λ, λ]|", XTildeMagnitude, "F4");
            yield return InspectableNode.RealScalar("Q_EP predicted", QEpPredicted, "F4");
            yield return InspectableNode.RealScalar("Q_peak predicted", QPeakPredicted, "F4");
            yield return InspectableNode.RealScalar("Q_peak empirical", QPeakEmpirical, "F4");
            if (!double.IsNaN(RelativeResidual))
                yield return InspectableNode.RealScalar("Relative residual", RelativeResidual, "P1");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
