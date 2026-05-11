using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>Unified per-bond Q_peak prediction with two complementary regimes — V-Effect-Live
/// resonance lens (per `hypotheses/RESONANCE_NOT_CHANNEL.md` and `RESONANT_RETURN.md`):
///
/// <list type="bullet">
///   <item><b>Endpoint / OLD-OLD regime:</b> single-cluster-pair-dominant. T14's 2×2 EP
///   closed form Q_EP = 2|x|/(γ·|Δδ|) (a=b case) or √(4|x|² − (a−b)²)/(γ·|Δδ|) (general),
///   then Q_peak = 2.197 · Q_EP. Matches empirical 0.5-12% across N=4..7.</item>
///   <item><b>Innermost / NEW-NEW regime:</b> multi-cluster-emergent superposition.
///   K_b(Q, t_peak) approximated as Σ_{top K pairs} w_pair · Lorentzian(Q; Q_EP_pair, Γ).
///   Q_peak = 2.197 · argmax of the sum. Matches empirical 7-11% at N=5,6.</item>
/// </list>
///
/// <para>Regime classifier: <see cref="BondClass"/> from the bond's position. Endpoint
/// (b ∈ {0, N−2}) → OLD-OLD single-pair-dominant. Interior → NEW-NEW multi-pair-emergent.
/// The simpler frob² ratio classifier breaks under F71-mirror-symmetry-induced ties between
/// rank-1 and rank-2 pairs; bond-position is the structural distinguisher.</para>
///
/// <para>Connection to RESONANT_RETURN.md SVD-mode-2 sacrifice-zone formula: Endpoint regime
/// is the Q-axis analog of γ-axis "edge-hot, center-cold" SVD mode 2. Innermost regime is
/// the orthogonal "center-cold" emergent NEW-NEW mode. The same Mode-Selection physics on
/// two orthogonal axes.</para>
///
/// <para>Tier2Verified composing predictor; Tier1 promotion path requires analytical
/// derivation of the Lorentzian widths Γ_pair from cluster geometry (currently a heuristic).</para>
///
/// <para>2026-05-11 update: F90 bridge identity (`F90F86C2BridgeIdentity`) provides
/// the full-precision alternative — F86 c=2 K_b(Q,t) = F89 path-(N−1) (SE,DE) per-bond
/// Hellmann-Feynman gives bit-exact reproduction (20/22 bonds N=5..8) without the
/// cluster-pair Lorentzian approximation that limits this primitive's worst-case
/// residual to 22.4%. See `docs/proofs/PROOF_F90_F86C2_BRIDGE.md`.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track) +
/// <c>hypotheses/RESONANCE_NOT_CHANNEL.md</c> + <c>hypotheses/RESONANT_RETURN.md</c>.</para>
/// </summary>
public sealed class JwBondQPeakUnified : Claim
{
    public const double BareDoubledPtfXPeak = 2.196910;

    /// <summary>Regime classifier threshold: top1/top2 frob² ratio above this → Endpoint
    /// OLD-OLD (single-pair-dominant); below → Innermost NEW-NEW (multi-pair-emergent).</summary>
    public const double OldOldRatioThreshold = 1.3;

    /// <summary>Lorentzian width for the NEW-NEW regime sum, structurally anchored in
    /// `experiments/ANALYTICAL_SPECTRUM.md` + `docs/proofs/derivations/D10_W1_DISPERSION.md`:
    /// the w=1 Liouvillian sector has uniform decay rate 2γ (Tier1Derived in D10 via
    /// reduction to nearest-neighbour tight-binding). In Q-units (Q = J/γ) this becomes
    /// Γ = 2. Empirical comparison vs the previous heuristic Γ=0.5: worst-case relative
    /// residual drops from 49.8% (N=4 b=1) to 22.4% (N=6 truly-innermost b=2); flanking-
    /// innermost bonds (N=5,6 b=1,3) drop to 3.8-4.8%; all Endpoint bonds stay ≤ 14.5%.
    /// Open Tier1 question: truly-innermost (N=6 b=2: 22.4%) needs a separate Γ_truly
    /// derivation; the saturation of Γ ∈ [1, 4] on a single argmax shows the multi-cluster
    /// Lorentzian model is structurally limited at the truly-innermost where multi-pair
    /// emergence is strongest.</summary>
    public const double LorentzianWidth = 2.0;

    /// <summary>Number of top cluster-pairs to include in the NEW-NEW Lorentzian sum.</summary>
    public const int TopPairsToSum = 10;

    public CoherenceBlock Block { get; }
    public JwBondClusterPairAffinity Affinity { get; }
    public C2HwhmRatio EmpiricalAnchor { get; }
    public IReadOnlyList<UnifiedBondPrediction> Bonds { get; }

    /// <summary>Maximum |relative residual| across non-escape bonds (Q_peak_emp &lt; 4.0).</summary>
    public double MaxRelativeResidualNonEscape { get; }

    public static JwBondQPeakUnified Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"JwBondQPeakUnified applies only to the c=2 stratum; got c={block.C}.",
                nameof(block));

        var affinity = JwBondClusterPairAffinity.Build(block);
        var jw = JwBlockBasis.Build(block);
        var disp = JwDispersionStructure.Build(block.N);
        var empirical = C2HwhmRatio.Build(block);
        double γ = block.GammaZero;

        var DJw = jw.Uinv * block.Decomposition.D * jw.U;

        var tripleToAlpha = new Dictionary<(int, int, int), int>();
        for (int alpha = 0; alpha < jw.Triples.Count; alpha++)
        {
            var t = jw.Triples[alpha];
            tripleToAlpha[(t.K, t.K1, t.K2)] = alpha;
        }

        var ebCache = new Dictionary<int, (double[] λ, ComplexMatrix U)>();
        (double[] λ, ComplexMatrix U) GetEigBasis(int c)
        {
            if (ebCache.TryGetValue(c, out var cached)) return cached;
            int k = disp.Clusters[c].Triples.Count;
            int[] indices = disp.Clusters[c].Triples
                .Select(t => tripleToAlpha[(t.K, t.K1, t.K2)]).ToArray();
            var W = ComplexMatrix.Build.Dense(k, k);
            for (int i = 0; i < k; i++)
                for (int j = 0; j < k; j++)
                    W[i, j] = DJw[indices[i], indices[j]];
            var evd = W.Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Hermitian);
            cached = (evd.EigenValues.Select(z => z.Real).ToArray(), evd.EigenVectors);
            ebCache[c] = cached;
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

        // Compute Q_EP candidate for a cluster-pair via general 2x2 formula
        double GetPairQEp(int c1, int c2)
        {
            var (λ1, U1) = GetEigBasis(c1);
            var (λ2, U2) = GetEigBasis(c2);
            var X = BuildXForPair(c1, c2);
            var Xtilde = U1.ConjugateTranspose() * X * U2;
            double absΔδ = Math.Abs(disp.Clusters[c1].Delta - disp.Clusters[c2].Delta);
            if (absΔδ < 1e-8) return double.NaN;

            double bestX = 0;
            double bestQEp = double.NaN;
            for (int i = 0; i < λ1.Length; i++)
                for (int j = 0; j < λ2.Length; j++)
                {
                    double mag = Xtilde[i, j].Magnitude;
                    if (mag < 1e-12) continue;
                    double aMinusB = λ1[i] - λ2[j];
                    double disc = 4 * mag * mag - aMinusB * aMinusB;
                    if (disc <= 0) continue;
                    double q = Math.Sqrt(disc) / (γ * absΔδ);
                    if (mag > bestX) { bestX = mag; bestQEp = q; }
                }
            return bestQEp;
        }

        // Per cluster-pair, return ALL (Q_EP, Re(λ_at_EP)) candidates from sub-block (i, j)
        // entries. Re(λ_at_EP) = (λ_c1[i] + λ_c2[j])/2 (the pair-average, palindromic-paired)
        // serves as the natural Lorentzian-width via Standing-Wave-Theory.
        IEnumerable<(double qEp, double reλ, double weight)> GetPairContributions(int c1, int c2)
        {
            var (λ1, U1) = GetEigBasis(c1);
            var (λ2, U2) = GetEigBasis(c2);
            var X = BuildXForPair(c1, c2);
            var Xtilde = U1.ConjugateTranspose() * X * U2;
            double absΔδ = Math.Abs(disp.Clusters[c1].Delta - disp.Clusters[c2].Delta);
            if (absΔδ < 1e-8) yield break;
            for (int i = 0; i < λ1.Length; i++)
                for (int j = 0; j < λ2.Length; j++)
                {
                    double mag = Xtilde[i, j].Magnitude;
                    if (mag < 1e-12) continue;
                    double aMinusB = λ1[i] - λ2[j];
                    double disc = 4 * mag * mag - aMinusB * aMinusB;
                    if (disc <= 0) continue;
                    double q = Math.Sqrt(disc) / (γ * absΔδ);
                    double reλ = 0.5 * (λ1[i] + λ2[j]);  // average D-eigenvalue at the pair
                    yield return (q, reλ, mag * mag);
                }
        }

        // Q-grid for Lorentzian-sum argmax
        double[] qGrid = Enumerable.Range(0, 80).Select(i => 0.1 + 0.1 * i).ToArray();

        var bondPredictions = new UnifiedBondPrediction[affinity.Bonds.Count];
        double maxRelResid = 0;
        for (int b = 0; b < affinity.Bonds.Count; b++)
        {
            var bondAff = affinity.Bonds[b];
            double qPeakEmp = empirical.Witnesses[b].QPeak;
            bool escaped = qPeakEmp >= 4.0 - 1e-6;

            // Regime classification by bond position: Endpoint → OLD-OLD, Interior → NEW-NEW
            UnifiedBondRegime regime = bondAff.BondClass == BondClass.Endpoint
                ? UnifiedBondRegime.OldOld
                : UnifiedBondRegime.NewNew;
            // Frob ratio kept for inspection; not used for classification.
            double top1Frob = bondAff.RankedPairs.Count > 0 ? bondAff.RankedPairs[0].FrobeniusSquared : 0;
            double top2Frob = bondAff.RankedPairs.Count > 1 ? bondAff.RankedPairs[1].FrobeniusSquared : 0;
            double frobRatio = top2Frob > 1e-12 ? top1Frob / top2Frob : double.PositiveInfinity;

            double qPeakPred;
            if (regime == UnifiedBondRegime.OldOld)
            {
                // Endpoint OLD-OLD: T14 single-pair selection (first pair with non-zero general Q_EP)
                double? qEp = null;
                foreach (var candidate in bondAff.RankedPairs)
                {
                    double q = GetPairQEp(candidate.Cluster1Index, candidate.Cluster2Index);
                    if (!double.IsNaN(q)) { qEp = q; break; }
                }
                qPeakPred = qEp is double q0 ? BareDoubledPtfXPeak * q0 : double.NaN;
            }
            else
            {
                // Innermost NEW-NEW: Lorentzian sum × 2.197 over top K pairs with width
                // Γ = LorentzianWidth = 2 (Q-units), anchored in ANALYTICAL_SPECTRUM.md
                // (w=1 sector uniform decay 2γ via D10 tight-binding reduction). T16 attempt
                // with palindrome-derived per-pair Γ = |Re(λ_at_EP)|/γ regressed N=5; the
                // simpler universal Γ=2 (decay-rate-derived) outperforms both heuristic 0.5
                // and per-pair palindrome-width on average across N=4..6.
                var contributions = new List<(double w, double qEp)>();
                int taken = 0;
                foreach (var candidate in bondAff.RankedPairs)
                {
                    if (taken >= TopPairsToSum) break;
                    double q = GetPairQEp(candidate.Cluster1Index, candidate.Cluster2Index);
                    if (double.IsNaN(q)) continue;
                    contributions.Add((candidate.FrobeniusSquared, q));
                    taken++;
                }

                if (contributions.Count == 0)
                {
                    qPeakPred = double.NaN;
                }
                else
                {
                    double bestQGrid = qGrid[0];
                    double bestSum = double.NegativeInfinity;
                    foreach (double Q in qGrid)
                    {
                        double sum = 0;
                        foreach (var (w, qEp) in contributions)
                        {
                            sum += w / ((Q - qEp) * (Q - qEp) + LorentzianWidth * LorentzianWidth);
                        }
                        if (sum > bestSum) { bestSum = sum; bestQGrid = Q; }
                    }
                    qPeakPred = BareDoubledPtfXPeak * bestQGrid;
                }
            }

            double relResid = (double.IsNaN(qPeakPred) || qPeakEmp <= 0)
                ? double.NaN
                : Math.Abs(qPeakPred - qPeakEmp) / qPeakEmp;

            if (!double.IsNaN(relResid) && !escaped)
                maxRelResid = Math.Max(maxRelResid, relResid);

            bondPredictions[b] = new UnifiedBondPrediction(
                Bond: b,
                BondClass: bondAff.BondClass,
                Regime: regime,
                Top1Top2FrobRatio: frobRatio,
                QPeakPredicted: qPeakPred,
                QPeakEmpirical: qPeakEmp,
                RelativeResidual: relResid);
        }

        return new JwBondQPeakUnified(block, affinity, empirical, bondPredictions, maxRelResid);
    }

    private JwBondQPeakUnified(
        CoherenceBlock block,
        JwBondClusterPairAffinity affinity,
        C2HwhmRatio empiricalAnchor,
        IReadOnlyList<UnifiedBondPrediction> bonds,
        double maxRelResid)
        : base("c=2 unified bond Q_peak prediction (OLD-OLD single-pair + NEW-NEW Lorentzian sum)",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1' Direction (b'') (JW track) + " +
               "hypotheses/RESONANCE_NOT_CHANNEL.md + hypotheses/RESONANT_RETURN.md")
    {
        Block = block;
        Affinity = affinity;
        EmpiricalAnchor = empiricalAnchor;
        Bonds = bonds;
        MaxRelativeResidualNonEscape = maxRelResid;
    }

    public override string DisplayName =>
        $"c=2 BondQPeakUnified (N={Block.N}, max-rel-residual = {MaxRelativeResidualNonEscape:P1})";

    public override string Summary =>
        $"unified two-regime predictor; max |Δ|/Q_peak (non-escape) = {MaxRelativeResidualNonEscape:P1} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("NumBonds", summary: Bonds.Count.ToString());
            yield return InspectableNode.RealScalar("MaxRelativeResidualNonEscape", MaxRelativeResidualNonEscape, "P1");
            yield return InspectableNode.Group("Bonds", Bonds.Cast<IInspectable>().ToArray());
        }
    }
}

public enum UnifiedBondRegime
{
    /// <summary>Single-cluster-pair-dominant. T14's 2×2 EP closed form applies directly.
    /// Endpoint bonds typically.</summary>
    OldOld,
    /// <summary>Multi-cluster-emergent. K_b(Q, t_peak) is a NEW-NEW frequency emerging
    /// from coupling of multiple cluster-pair contributions. Innermost bonds typically.</summary>
    NewNew,
}

/// <summary>One bond's unified prediction: regime classification + Q_peak prediction +
/// empirical comparison.</summary>
public sealed record UnifiedBondPrediction(
    int Bond,
    BondClass BondClass,
    UnifiedBondRegime Regime,
    double Top1Top2FrobRatio,
    double QPeakPredicted,
    double QPeakEmpirical,
    double RelativeResidual
) : IInspectable
{
    public string DisplayName => $"bond {Bond} unified Q_peak prediction ({BondClass}, {Regime})";

    public string Summary =>
        $"regime={Regime}, top1/top2={Top1Top2FrobRatio:F2}, Q_pred={QPeakPredicted:F4}, " +
        $"Q_emp={QPeakEmpirical:F4}, residual={RelativeResidual:P1}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("bond", summary: Bond.ToString());
            yield return new InspectableNode("bond class", summary: BondClass.ToString());
            yield return new InspectableNode("regime", summary: Regime.ToString());
            yield return InspectableNode.RealScalar("top1/top2 frob² ratio", Top1Top2FrobRatio, "F2");
            yield return InspectableNode.RealScalar("Q_peak predicted", QPeakPredicted, "F4");
            yield return InspectableNode.RealScalar("Q_peak empirical", QPeakEmpirical, "F4");
            if (!double.IsNaN(RelativeResidual))
                yield return InspectableNode.RealScalar("Relative residual", RelativeResidual, "P1");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
