using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 surface for the Liouvillian-free two-term Q-pair router
/// (<see cref="TwoTermPalindromeRouting"/>) as a single registry-queryable Claim. The router decides
/// a two-term bond bilinear's spectral fate (truly / soft / hard) AND routes the hidden palindrome
/// symmetry Q into a letter-based family (P1 / Uniform / Alternating / Continuous / None_) with no
/// eigensolve. This Claim is the non-diagonal counterpart to
/// <see cref="F87DiagonalCellBipartiteWitnessSet"/>: where that set reads the diagonal Klein-cell
/// (P1-family) special case through a chiral K, this Claim covers the full hidden-Q routing across
/// every Klein cell, including the off-diagonal Uniform / Alternating / Continuous escapes.
///
/// <para>The asserted fact is empirical: over ALL two-term bilinear pairs (the nine single bilinears
/// XX, YY, ZZ, XY, YX, XZ, ZX, YZ, ZY taken as unordered pairs, including self-pairs) the router's
/// fate verdict matches the spectral authority <see cref="PauliPairTrichotomy.Classify(ChainSystem,
/// IReadOnlyList{PauliTerm}, double, double, PauliLetter)"/> bit-exactly. The router is N-independent,
/// so the authority cross-check runs on a small internal chain (N=4, matching the proven
/// <c>TwoTermPalindromeRoutingTests</c>); the registry chain is kept only for the Claim's identity.</para>
///
/// <para>The comparison is built LAZILY: construction touches no eigensolve. The all-pairs
/// routing-vs-authority sweep (45 small N=4 eigendecompositions) runs only at first access of
/// <see cref="MatchCount"/> / <see cref="AllPairsMatchAuthority"/> / <see cref="Summary"/> /
/// <see cref="ExtraChildren"/>.</para>
///
/// <para>Tier: Tier2Empirical, like <see cref="F87StandardWitnessSet"/>. The router is a
/// routing-rule viewpoint empirically verified against the authority, not a derived theorem. Typed
/// parents: <see cref="F87TrichotomyClassification"/> (the spectral authority the router is verified
/// against, Tier1Derived), <see cref="F87DiagonalCellBipartiteWitnessSet"/> (the diagonal P1-family
/// special case this generalises, Tier1Candidate), and
/// <see cref="Core.Symmetry.CrossoverMirrorSqrtNinetyClaim"/> (the Continuous-family crossover
/// mirror's existing C# realisation, Tier1Derived). All three are at least as strong as
/// Tier2Empirical, so the strength-inheritance check (parent ≥ child) passes.</para>
///
/// <para>Anchor: <c>experiments/TWO_TERM_PALINDROME_KLEIN_ROUTING.md</c> +
/// <c>experiments/NON_HEISENBERG_PALINDROME.md</c> +
/// <c>simulations/framework/diagnostics/q_family_routing.py</c> (the Python reference) +
/// <see cref="TwoTermPalindromeRouting"/> + <c>TwoTermPalindromeRoutingTests</c>.</para></summary>
public sealed class TwoTermPalindromeRoutingClaim : Claim
{
    /// <summary>The nine canonical two-letter bilinears (the only letter pairs with both sites lit).
    /// The unordered pairs are the upper triangle (a, b) with b at or after a, including self-pairs.</summary>
    private static readonly string[] Bilinears = { "XX", "YY", "ZZ", "XY", "YX", "XZ", "ZX", "YZ", "ZY" };

    /// <summary>N used for the authority cross-check. The router is N-independent; N=4 matches the
    /// proven <c>TwoTermPalindromeRoutingTests.Routing_AgreesWithSpectralAuthority_OverAll36Combos_N4</c>
    /// and keeps the 45 eigendecompositions small (256×256 Liouvillian).</summary>
    private const int AuthorityN = 4;

    /// <summary>One routing-vs-authority comparison: a pair (a, b), its router verdict (fate + Q
    /// family + reason), and the spectral authority's fate. <see cref="Matches"/> iff the router fate
    /// equals the authority fate.</summary>
    public readonly record struct PairComparison(
        string Term1, string Term2, TrichotomyClass RoutedFate, QFamily Family, string Reason,
        TrichotomyClass AuthorityFate)
    {
        public bool Matches => RoutedFate == AuthorityFate;
    }

    /// <summary>The chain providing this Claim's identity (N for the display). The authority
    /// cross-check itself runs on an internal <see cref="AuthorityN"/>=4 chain.</summary>
    public ChainSystem Chain { get; }

    private readonly Lazy<IReadOnlyList<PairComparison>> _comparisons;

    public TwoTermPalindromeRoutingClaim(ChainSystem chain)
        : base("Q-pair routing: hidden-Q (P1/P4) two-term Klein classifier, verified vs the spectral authority",
               Tier.Tier2Empirical,
               "experiments/TWO_TERM_PALINDROME_KLEIN_ROUTING.md + " +
               "experiments/NON_HEISENBERG_PALINDROME.md + " +
               "simulations/framework/diagnostics/q_family_routing.py")
    {
        Chain = chain ?? throw new ArgumentNullException(nameof(chain));
        _comparisons = new Lazy<IReadOnlyList<PairComparison>>(BuildComparisons);
    }

    /// <summary>All <see cref="PairCount"/> routing-vs-authority comparisons (lazy: built at first
    /// access, runs the N=4 authority eigensolves once).</summary>
    public IReadOnlyList<PairComparison> Comparisons => _comparisons.Value;

    /// <summary>Total number of unordered two-term pairs covered (45: C(9,2) + 9 self-pairs).</summary>
    public int PairCount => Comparisons.Count;

    /// <summary>How many pairs the router classifies exactly as the spectral authority does.</summary>
    public int MatchCount => Comparisons.Count(c => c.Matches);

    /// <summary>True iff the router agrees with the spectral authority on every pair.</summary>
    public bool AllPairsMatchAuthority => MatchCount == PairCount;

    /// <summary>The family distribution over all pairs: counts of Truly (P1) / Uniform / Alternating /
    /// Continuous / None_ (the routed Q family for each pair).</summary>
    public IReadOnlyDictionary<QFamily, int> FamilyDistribution =>
        Comparisons.GroupBy(c => c.Family).ToDictionary(g => g.Key, g => g.Count());

    /// <summary>The canonical soft example XX+XZ, which routes to the Uniform Q-family.</summary>
    public PairComparison SoftExample => Find("XX", "XZ");

    /// <summary>A hard example XY+XZ, which routes to None_ and carries its Q7 Klein-cell reason.</summary>
    public PairComparison HardExample => Find("XY", "XZ");

    private PairComparison Find(string a, string b) =>
        Comparisons.First(c => c.Term1 == a && c.Term2 == b);

    public override string DisplayName =>
        $"Q-pair routing classifier (two-term, all pairs vs spectral authority, identity N={Chain.N})";

    public override string Summary
    {
        get
        {
            var d = FamilyDistribution;
            string dist = $"Truly={Count(d, QFamily.P1)}, Uniform={Count(d, QFamily.Uniform)}, " +
                          $"Alternating={Count(d, QFamily.Alternating)}, Continuous={Count(d, QFamily.Continuous)}, " +
                          $"None={Count(d, QFamily.None_)}";
            return $"hidden-Q two-term router: [{dist}]; {MatchCount}/{PairCount} match authority " +
                   $"(N={AuthorityN} cross-check) ({Tier.Label()})";
        }
    }

    private static int Count(IReadOnlyDictionary<QFamily, int> d, QFamily f) => d.TryGetValue(f, out var n) ? n : 0;

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("chain (identity)",
                summary: $"N={Chain.N}, {Chain.HType}, {Chain.Topology}; authority cross-check at N={AuthorityN}");

            var d = FamilyDistribution;
            yield return new InspectableNode("family distribution (all two-term pairs)",
                summary: $"Truly(P1)={Count(d, QFamily.P1)}, Uniform={Count(d, QFamily.Uniform)}, " +
                         $"Alternating={Count(d, QFamily.Alternating)}, Continuous={Count(d, QFamily.Continuous)}, " +
                         $"None={Count(d, QFamily.None_)} over {PairCount} pairs");

            yield return new InspectableNode("authority agreement",
                summary: $"{MatchCount}/{PairCount} pairs: router fate == PauliPairTrichotomy.Classify " +
                         $"({(AllPairsMatchAuthority ? "all match" : "MISMATCH")})");

            var soft = SoftExample;
            yield return new InspectableNode("soft example (XX+XZ)",
                summary: $"routed {soft.RoutedFate}/{soft.Family}, authority {soft.AuthorityFate}: {soft.Reason}");

            var hard = HardExample;
            yield return new InspectableNode("hard example (XY+XZ)",
                summary: $"routed {hard.RoutedFate}/{hard.Family}, authority {hard.AuthorityFate}; Q7 reason: {hard.Reason}");
        }
    }

    private static PauliTerm T(string label) => new(PauliLabel.Parse(label), Complex.One);

    /// <summary>Build the all-pairs routing-vs-authority comparison over the upper triangle of the
    /// nine bilinears (including self-pairs). The authority runs on an internal N=4 chain (the router
    /// is N-independent); the registry chain is untouched here to keep registry build cheap.</summary>
    private static IReadOnlyList<PairComparison> BuildComparisons()
    {
        var authorityChain = new ChainSystem(AuthorityN, 1.0, 0.05);
        var rows = new List<PairComparison>(Bilinears.Length * (Bilinears.Length + 1) / 2);
        for (int i = 0; i < Bilinears.Length; i++)
        {
            for (int j = i; j < Bilinears.Length; j++)
            {
                string a = Bilinears[i], b = Bilinears[j];
                var routed = TwoTermPalindromeRouting.Classify(a, b);
                var authority = PauliPairTrichotomy.Classify(authorityChain, new[] { T(a), T(b) });
                rows.Add(new PairComparison(
                    Term1: a, Term2: b,
                    RoutedFate: routed.Fate, Family: routed.Family, Reason: routed.Reason,
                    AuthorityFate: authority));
            }
        }
        return rows;
    }
}
