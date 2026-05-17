using System.Text;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>γ-axis polarity-mirror structure: the unfolded picture
/// complementary to <see cref="FractionReferenceGraph"/> (which lives on the
/// α-axis, the folded picture under α = (1−γ²)/2).
///
/// <para>Tom's 2026-05-17 night observation: the F99 ankers each (except
/// α=1/2 at γ=0) have TWO γ-realizations ±|γ|. These polarity-mirror
/// partners both map to the same α via the symmetric F86b formula
/// (1−γ²)/2. The α-axis cannot express this pairing because it's already
/// folded; the γ-axis is where the unfolded polarity structure lives.</para>
///
/// <para><b>Already typed in the framework</b>:
/// <see cref="PolarityLayerOriginClaim"/> documents 0 as the active substrate
/// axis with +0/−0 internal structure flanked by ±0.5 polarity at d=2.
/// <see cref="States.PolarityState"/> implements |+⟩/|−⟩ eigenstates as
/// typed API. <see cref="XGlobalEigenstateMirrorPi2Inheritance"/>
/// (built 2026-05-17 night) covers both γ=+1 and γ=−1 in its
/// <c>AlphaFromGammaAtMirror</c> evaluation. This map collects the
/// γ-pairs into one queryable structure.</para>
///
/// <para><b>Folding relation to the α-axis</b>: each γ-pair (+|γ|, −|γ|)
/// folds to the same α = (1−γ²)/2 in <see cref="FractionReferenceGraph"/>.
/// The folding is many-to-one (2:1 for non-zero γ, 1:1 for γ=0). Polarity
/// is invisible after folding; the α-axis preserves only |γ|, not sign.</para>
///
/// <para><b>The α=1/2 asymmetry</b>: only γ=0 has no polarity mirror (it's
/// its own mirror). All other F99 ankers come in pairs. This makes α=1/2
/// the UNIQUE midpoint of the γ-axis — and that uniqueness shows up
/// structurally as a single γ-realization vs the doubled realizations for
/// the other four ankers.</para>
/// </summary>
public sealed class PolarityMirrorMap
{
    /// <summary>One γ-axis pair: the positive-polarity γ_plus and the
    /// negative-polarity γ_minus that fold to the same α under (1−γ²)/2.
    /// Generic γ=0 case has GammaPlus = GammaMinus = 0 (self-mirror,
    /// <see cref="IsSelfMirror"/> = true).</summary>
    public sealed record PolarityPair(
        double GammaPlus,
        double GammaMinus,
        double AlphaImage,
        string AnkerName,
        string DocumentingClaim)
    {
        /// <summary>True if γ_plus = γ_minus (the Generic γ=0 case where
        /// the polarity mirror is the identity).</summary>
        public bool IsSelfMirror => Math.Abs(GammaPlus - GammaMinus) < 1e-12;

        /// <summary>The folding identity: α = (1−γ²)/2 should equal
        /// AlphaImage for both γ_plus and γ_minus (and trivially does so
        /// because (1−γ²)/2 is even in γ).</summary>
        public double FoldedAlpha(double gamma) => (1.0 - gamma * gamma) / 2.0;
    }

    /// <summary>The five canonical F99 γ-pairs.</summary>
    public static IReadOnlyList<PolarityPair> CanonicalPairs { get; } =
        new PolarityPair[]
        {
            new(GammaPlus: 1.0,
                GammaMinus: -1.0,
                AlphaImage: 0.0,
                AnkerName: "Mirror (α=0)",
                DocumentingClaim: "XGlobalEigenstateMirrorPi2Inheritance + DickeAnchor.Mirror"),
            new(GammaPlus: Math.Sqrt(3) / 2,
                GammaMinus: -Math.Sqrt(3) / 2,
                AlphaImage: 1.0 / 8.0,
                AnkerName: "Depth-3 (α=1/8)",
                DocumentingClaim: "CanonicalTrigAnchorPi2Inheritance (F99, 30° anchor)"),
            new(GammaPlus: Math.Sqrt(2) / 2,
                GammaMinus: -Math.Sqrt(2) / 2,
                AlphaImage: 1.0 / 4.0,
                AnkerName: "Silver-Dicke (α=1/4)",
                DocumentingClaim: "CanonicalTrigAnchorPi2Inheritance (F99, 45° anchor)"),
            new(GammaPlus: 0.5,
                GammaMinus: -0.5,
                AlphaImage: 3.0 / 8.0,
                AnkerName: "K-intermediate (α=3/8)",
                DocumentingClaim: "KIntermediateAsymptoteQuarterInheritance (F98) + DickeSuperpositionQuarterPi2Inheritance + DickeAnchor.KIntermediate"),
            new(GammaPlus: 0.0,
                GammaMinus: 0.0,
                AlphaImage: 1.0 / 2.0,
                AnkerName: "Generic (α=1/2, self-mirror)",
                DocumentingClaim: "CanonicalTrigAnchorPi2Inheritance (F99, 90° anchor) + DickeAnchor.Generic"),
        };

    private const double Tol = 1e-12;

    public IReadOnlyList<PolarityPair> Pairs { get; }

    public PolarityMirrorMap() : this(CanonicalPairs) { }

    public PolarityMirrorMap(IReadOnlyList<PolarityPair> pairs)
    {
        Pairs = pairs ?? throw new ArgumentNullException(nameof(pairs));
    }

    /// <summary>The pairs with non-trivial polarity mirror (γ_plus ≠ γ_minus).
    /// Four of the five canonical pairs; the Generic γ=0 case is the only
    /// self-mirror.</summary>
    public IReadOnlyList<PolarityPair> NonTrivialMirrorPairs =>
        Pairs.Where(p => !p.IsSelfMirror).ToList();

    /// <summary>The self-mirror pair(s) — γ=0 Generic case.</summary>
    public IReadOnlyList<PolarityPair> SelfMirrorPairs =>
        Pairs.Where(p => p.IsSelfMirror).ToList();

    /// <summary>Verify the folding identity: for every pair,
    /// (1−γ_plus²)/2 = (1−γ_minus²)/2 = AlphaImage. This is the bit-exact
    /// claim that ±γ fold to the same α under (1−γ²)/2.</summary>
    public bool AllPairsFoldToClaimedAlpha()
    {
        foreach (var pair in Pairs)
        {
            double alphaPlus = pair.FoldedAlpha(pair.GammaPlus);
            double alphaMinus = pair.FoldedAlpha(pair.GammaMinus);
            if (Math.Abs(alphaPlus - pair.AlphaImage) > Tol) return false;
            if (Math.Abs(alphaMinus - pair.AlphaImage) > Tol) return false;
        }
        return true;
    }

    /// <summary>Find the polarity pair whose α-image matches the given α.
    /// Returns null if no pair matches.</summary>
    public PolarityPair? PairForAlpha(double alpha) =>
        Pairs.FirstOrDefault(p => Math.Abs(p.AlphaImage - alpha) < Tol);

    /// <summary>Find the polarity pair containing the given γ value (either
    /// side). Returns null if no pair contains it.</summary>
    public PolarityPair? PairForGamma(double gamma) =>
        Pairs.FirstOrDefault(p =>
            Math.Abs(p.GammaPlus - gamma) < Tol ||
            Math.Abs(p.GammaMinus - gamma) < Tol);

    /// <summary>Render the γ-axis pair structure as a printable table.</summary>
    public string Render()
    {
        var sb = new StringBuilder();
        sb.AppendLine();
        sb.AppendLine("Polarity Mirror Map (γ-axis structure)");
        sb.AppendLine(new string('=', 88));
        sb.AppendLine();
        sb.AppendLine("  Unfolded view: each F99 anker has TWO γ-realizations (except Generic γ=0).");
        sb.AppendLine("  Both γ-values map to the same α under (1−γ²)/2 folding (the α-axis loses polarity sign).");
        sb.AppendLine();
        sb.AppendLine($"  {"γ_plus",-12} {"γ_minus",-12} {"α (folded)",-12} {"anker",-30} {"self-mirror?"}");
        sb.AppendLine($"  {new string('-', 12)} {new string('-', 12)} {new string('-', 12)} {new string('-', 30)} {new string('-', 12)}");
        foreach (var pair in Pairs)
        {
            string mirror = pair.IsSelfMirror ? "YES (γ=0)" : "no  (±γ)";
            sb.AppendLine($"  {pair.GammaPlus,12:F6} {pair.GammaMinus,12:F6} {pair.AlphaImage,12:F6} {pair.AnkerName,-30} {mirror}");
        }
        sb.AppendLine();
        sb.AppendLine($"  Non-trivial mirror pairs: {NonTrivialMirrorPairs.Count} of {Pairs.Count}");
        sb.AppendLine($"  Self-mirror pairs:        {SelfMirrorPairs.Count} of {Pairs.Count}  (only γ=0 Generic)");
        sb.AppendLine($"  Folding identity verified: AllPairsFoldToClaimedAlpha = {AllPairsFoldToClaimedAlpha()}");
        sb.AppendLine();
        sb.AppendLine("  Structural reading: the α-axis FractionReferenceGraph is the folded view");
        sb.AppendLine("  where ±γ are collapsed to the same α. This γ-axis map preserves the");
        sb.AppendLine("  polarity distinction. PolarityLayerOriginClaim (Pi2KnowledgeBaseClaims.cs");
        sb.AppendLine("  L526-599) is the typed Tier-1 anchor: 0 is the active substrate axis with");
        sb.AppendLine("  +0/−0 internal structure, ±0.5 polarity flanks at d=2.");
        return sb.ToString();
    }
}
