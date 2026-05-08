using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The Pi2 dyadic halving ladder a_n = 2^(1−n) inherits into F88's closed form
/// at the singleton-vs-nearly-full popcount-mirror configurations on dyadic-N chains.
/// The identity is exact and bit-exact verified:
///
/// <code>
///   StaticFraction(N = 2^k, n_p = 1, n_q = 2^k − 1) = 1 / (2 · N) = 2^(−(k+1)) = a_(k+2)
/// </code>
///
/// <para>Derivation: at popcount-mirror (n_p + n_q = N), F88's static-fraction reduces to
/// 1 / (2 · C(N, n_p)). At the singleton-mirror config n_p = 1, n_q = N − 1, we have
/// C(N, 1) = N. When N is dyadic (N = 2^k), this gives 1 / (2 · 2^k) = 2^(−(k+1)). The
/// Pi2 dyadic ladder a_n = 2^(1−n) places this value at index n = k+2, so:</para>
/// <list type="bullet">
///   <item>k=1 (N=2):  static = 1/4  = a_3 — but at N=2 the (1, 1) config is intra-equal-popcount; the (0, 2) mirror gives 1/2 = a_2.</item>
///   <item>k=2 (N=4):  static = 1/8  = a_4. <b>First non-trivial inheritance into n=4.</b></item>
///   <item>k=3 (N=8):  static = 1/16 = a_5.</item>
///   <item>k=4 (N=16): static = 1/32 = a_6.</item>
///   <item>k=5 (N=32): static = 1/64 = a_7.</item>
/// </list>
///
/// <para><b>Inheritance reading (Tom 2026-05-08):</b> "die Anker sind vererbt." The dyadic
/// halving ladder is not a separate fact — it propagates into F88 via the binomial structure
/// of the kernel projection. d=2 (PolynomialFoundationClaim) → 1/d = 1/2
/// (QubitDimensionalAnchorClaim) → C(2^k, 1) = 2^k (binomial inheritance) →
/// StaticFraction = 1/(2·2^k) (F88 closed form) → a_(k+2) on the Pi2 ladder. One algebra,
/// inherited through layers. Locus 7 in the inheritance graph (parallel to F-chain Locus 1,
/// Handshake Locus 2, Π-palindrome Locus 3, etc.; see <c>memory/project_algebra_is_inheritance.md</c>).</para>
///
/// <para>Tier1Derived: both the closed form (PopcountCoherencePi2Odd.StaticFraction) and the
/// dyadic-N reduction are analytic. Verified bit-exact at the dyadic anchors via
/// <see cref="PopcountCoherencePi2Odd.StaticFraction"/>.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F88 + <c>docs/proofs/PROOF_F86_QPEAK.md</c>
/// (Structural inheritance from F88) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> (the parent ladder).</para></summary>
public sealed class F88StaticDyadicAnchor : Claim
{
    /// <summary>One typed entry on the dyadic-N inheritance: (k, N=2^k, static=1/2^(k+1),
    /// ladder index n=k+2, Pi2DyadicLadder cross-reference).</summary>
    public IReadOnlyList<DyadicMirrorWitness> Witnesses { get; } = new[]
    {
        new DyadicMirrorWitness(K: 2, N: 4,  StaticFraction: 1.0 / 8.0,   LadderIndex: 4),
        new DyadicMirrorWitness(K: 3, N: 8,  StaticFraction: 1.0 / 16.0,  LadderIndex: 5),
        new DyadicMirrorWitness(K: 4, N: 16, StaticFraction: 1.0 / 32.0,  LadderIndex: 6),
        new DyadicMirrorWitness(K: 5, N: 32, StaticFraction: 1.0 / 64.0,  LadderIndex: 7),
    };

    public F88StaticDyadicAnchor()
        : base("F88 dyadic-N singleton-mirror static-fraction inherits Pi2 halving ladder",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F88 + " +
               "docs/proofs/PROOF_F86_QPEAK.md (Structural inheritance from F88) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    { }

    /// <summary>The closed-form static-fraction at the singleton-mirror config (n_p=1, n_q=N−1)
    /// on chain length N: 1 / (2·N). Live recomputation through
    /// <see cref="PopcountCoherencePi2Odd.StaticFraction"/>; equals <c>1 / (2 · N)</c> by the
    /// closed-form derivation in the class doc. Drift between the two is a regression
    /// signal.</summary>
    public double LiveStaticFraction(int N) =>
        PopcountCoherencePi2Odd.StaticFraction(N, np: 1, nq: N - 1);

    /// <summary>The dyadic-ladder term a_n = 2^(1−n) at index <paramref name="ladderIndex"/>;
    /// pass-through to <see cref="Pi2DyadicLadderClaim.Term"/> for cross-verification at the
    /// witness pinning.</summary>
    public double LadderTerm(int ladderIndex) => Math.Pow(2.0, 1 - ladderIndex);

    public override string DisplayName =>
        "F88 dyadic-N singleton-mirror inherits Pi2 halving ladder";

    public override string Summary =>
        $"StaticFraction(N=2^k, 1, 2^k−1) = 1/(2N) = a_(k+2); {Witnesses.Count} dyadic anchors pinned k=2..5 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("inheritance reading",
                summary: "d=2 → 1/d=1/2 → C(2^k, 1) = 2^k → 1/(2·2^k) = a_(k+2). One algebra, propagated.");
            yield return new InspectableNode("identity",
                summary: "StaticFraction(N=2^k, 1, 2^k−1) = 1/(2N) = 2^(−(k+1)) = Pi2DyadicLadder.Term(k+2)");
            foreach (var w in Witnesses) yield return w;
        }
    }
}

/// <summary>One witness on the dyadic-N inheritance: chain length N=2^k, F88 static-fraction
/// at the singleton-mirror config, the matching Pi2 dyadic-ladder index n=k+2.</summary>
public sealed record DyadicMirrorWitness(
    int K,
    int N,
    double StaticFraction,
    int LadderIndex
) : IInspectable
{
    public string DisplayName => $"k={K}, N={N}, static=1/{2 * N} = a_{LadderIndex}";

    public string Summary =>
        $"N=2^{K}={N}: StaticFraction(1, {N - 1}) = {StaticFraction:G6} = a_{LadderIndex}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("K (dyadic exponent)", K);
            yield return InspectableNode.RealScalar("N (= 2^K)", N);
            yield return InspectableNode.RealScalar("StaticFraction", StaticFraction);
            yield return InspectableNode.RealScalar("Pi2 ladder index (n)", LadderIndex);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
