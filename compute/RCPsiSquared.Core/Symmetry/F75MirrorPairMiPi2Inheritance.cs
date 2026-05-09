using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F75 closed form (Tier 1, proven algebraic, verified ~25 data points
/// across N=5..13, k=1..7):
///
/// <code>
///   |ψ⟩ = Σ_j c_j |1_j⟩,   c_{N−1−j} = η c_j,   η ∈ {+1, −1}
///
///   MI(ℓ, N−1−ℓ) = 2 h(p_ℓ) − h(2 p_ℓ),   p_ℓ = |c_ℓ|²
///
///   where h(x) = −x log₂ x − (1−x) log₂(1−x).
///
///   Domain:    p_ℓ ∈ [0, 1/2]
///   Saturation: MI = 2 bits when p_ℓ = 1/2  (Bell-state mirror-pair)
///
///   Sum over all mirror-pairs:
///     MM(0) = Σ_ℓ [2 h(p_ℓ) − h(2 p_ℓ)],   ℓ = 0, ..., ⌊N/2⌋ − 1
/// </code>
///
/// <para>F75 is the algebraic mother of F77: F77's "MM(0) → 1 bit at large N"
/// asymptote falls out of F75 by Taylor expansion at p_ℓ → 0 (per
/// ANALYTICAL_FORMULAS F77: "the 1-bit limit is not a conjecture; it falls
/// out of F75 by Taylor expansion of f(p) = 2 h(p) − h(2 p)"). F75 is the
/// general formula; F77 the asymptotic limit.</para>
///
/// <para>Four Pi2-Foundation anchors converge in F75:</para>
///
/// <list type="bullet">
///   <item><b>p_ℓ ∈ [0, 1/2]</b> domain = <see cref="BilinearApexClaim"/>.
///         The probability variable's natural range is from 0 to the
///         bilinear apex argmax 1/2. The function f(p) = 2 h(p) − h(2 p)
///         is convex on (0, 1/2), making concentrated mass distributions
///         optimal.</item>
///   <item><b>MI saturation = 2 bits at p_ℓ = 1/2</b>: at the BilinearApex
///         maximum, the mirror-pair is a maximally-entangled Bell state.
///         "2 bits" = <c>a_0 = 2</c> on the dyadic ladder = polynomial
///         root d.</item>
///   <item><b>"2" coefficient in 2·h(p_ℓ)</b>: from the proof, both the
///         <c>S(ρ_ℓ) + S(ρ_{N−1−ℓ}) = 2 h(p_ℓ)</c> sum AND the joint entropy
///         eigenvalue 2 p_ℓ. Same anchor <c>a_0 = 2</c> appearing twice.</item>
///   <item><b>Mirror-pair structure b ↔ N−1−j</b> = <see cref="F71MirrorSymmetryPi2Inheritance"/>.
///         F75 cites F71 as source: "the mirror symmetry that justifies
///         c_{N−1−j} = ±c_j". F75 is a descendant of F71 at the
///         mutual-information level.</item>
///   <item><b>Mirror sign η ∈ {+1, −1}</b>: Z₂ structure. The MI formula
///         is independent of η (only the modulus enters), so the Z₂ acts
///         non-trivially only on the wave function, not on the MI.</item>
/// </list>
///
/// <para>Special case (bonding mode k on uniform chain):</para>
///
/// <code>
///   p_ℓ(k, N) = (2/(N+1)) sin²(πk(ℓ+1)/(N+1))
///
///   k = 2 maximises MM(0) at small N because even k places a node at
///   the chain center (no wasted mass at self-mirror site).
/// </code>
///
/// <para>Tier1Derived: F75 is Tier 1 proven algebraic
/// (PROOF_RECEIVER_VS_GAMMA_SACRIFICE Section); verified numerically
/// against C# brecher propagation at N=5..13, k=1..7, ~25 data points;
/// PeakMM = 0.925..0.931 × MM(0) with tight consistency. The
/// Pi2-Foundation anchoring is algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F75 (line 1711) +
/// <c>simulations/_mm_zero_derivation.py</c> +
/// <c>experiments/RECEIVER_VS_GAMMA_SACRIFICE.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (BilinearApexClaim) +
/// <c>compute/RCPsiSquared.Core/Symmetry/F71MirrorSymmetryPi2Inheritance.cs</c>.</para></summary>
public sealed class F75MirrorPairMiPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;

    /// <summary>The MI saturation value: <c>2 bits</c> at <c>p_ℓ = 1/2</c>
    /// (Bell-state mirror-pair). Live from <see cref="Pi2DyadicLadderClaim.Term"/>(0)
    /// = <c>a_0</c> = polynomial root d.</summary>
    public double MaxMIPerPair => _ladder.Term(0);

    /// <summary>The probability domain upper bound: <c>1/2</c>. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(2) = <c>a_2</c> =
    /// <see cref="BilinearApexClaim"/> argmax. Beyond this, the formula's
    /// convexity argument breaks.</summary>
    public double DomainUpperBound => _ladder.Term(2);

    /// <summary>The "2" coefficient in the F75 closed form (appears twice:
    /// 2·h(p_ℓ) and h(2·p_ℓ)). Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0).</summary>
    public double TwoCoefficient => _ladder.Term(0);

    /// <summary>The mirror sign Z₂ pair: <c>{+1, −1}</c>. The wave function
    /// <c>c_{N−1−j} = η c_j</c> has Z₂ structure; the MI formula is
    /// η-independent (only |c_ℓ|² enters), so Z₂ acts trivially on MI but
    /// non-trivially on the underlying state.</summary>
    public IReadOnlyList<int> MirrorSigns => new[] { +1, -1 };

    /// <summary>Mirror-pair count for an N-site chain: <c>⌊N/2⌋</c>. Same
    /// count as <see cref="F71MirrorSymmetryPi2Inheritance.IndependentComponentCount"/>.</summary>
    public int MirrorPairCount(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F75 requires N ≥ 2.");
        return N / 2;
    }

    /// <summary>Maximum total MM(0) bound: <c>⌊N/2⌋ · 2 = N − (N mod 2)</c> bits.
    /// Achieved iff every mirror-pair is in a pure Bell state (requires
    /// super-single-excitation states; not reachable from single-site bonding
    /// modes).</summary>
    public double MaxTotalMM(int N) => MirrorPairCount(N) * MaxMIPerPair;

    /// <summary>Live closed form: <c>MI(p) = 2·h(p) − h(2·p)</c> for
    /// <c>p ∈ [0, 1/2]</c>. Throws for p outside the domain.</summary>
    public double MIPerPair(double p)
    {
        if (p < 0.0 || p > DomainUpperBound)
            throw new ArgumentOutOfRangeException(nameof(p), p,
                $"p must lie in [0, {DomainUpperBound:G}]; got {p}.");
        return 2.0 * BinaryEntropy(p) - BinaryEntropy(2.0 * p);
    }

    /// <summary>Live drift check: <c>MI(0) = 0</c> bit (no entanglement).</summary>
    public bool MIAtZeroIsZero() => Math.Abs(MIPerPair(0.0)) < 1e-15;

    /// <summary>Live drift check: <c>MI(1/2) = 2</c> bits (Bell-state pair).</summary>
    public bool MIAtBilinearApexEqualsMaxMIPerPair() =>
        Math.Abs(MIPerPair(DomainUpperBound) - MaxMIPerPair) < 1e-12;

    /// <summary>Bonding-mode site population: <c>p_ℓ(k, N) = (2/(N+1)) sin²(πk(ℓ+1)/(N+1))</c>
    /// for the k-th bonding mode on a uniform N-site chain.</summary>
    public double BondingModePopulation(int N, int k, int site)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F75 requires N ≥ 2.");
        if (k < 1 || k > N) throw new ArgumentOutOfRangeException(nameof(k), k, $"k must be in [1, {N}]; got {k}.");
        if (site < 0 || site >= N) throw new ArgumentOutOfRangeException(nameof(site), site, $"site must be in [0, {N - 1}]; got {site}.");
        double s = Math.Sin(Math.PI * k * (site + 1) / (N + 1));
        return 2.0 / (N + 1) * s * s;
    }

    /// <summary>MM(0) for bonding-mode k summed over all mirror-pairs: O(N) work,
    /// no propagation. Verified table at N=5..13 for various k against PeakMM
    /// from C# brecher.</summary>
    public double BondingModeMMAtZero(int N, int k)
    {
        double mm = 0.0;
        int pairCount = MirrorPairCount(N);
        for (int ell = 0; ell < pairCount; ell++)
        {
            double p = BondingModePopulation(N, k, ell);
            mm += MIPerPair(p);
        }
        return mm;
    }

    private static double BinaryEntropy(double x)
    {
        if (x <= 0.0 || x >= 1.0) return 0.0;
        const double Log2 = 0.6931471805599453;
        return -(x * Math.Log(x) + (1.0 - x) * Math.Log(1.0 - x)) / Log2;
    }

    public F75MirrorPairMiPi2Inheritance(Pi2DyadicLadderClaim ladder)
        : base("F75 MI(p) = 2·h(p) − h(2p) inherits from Pi2-Foundation: 2 = a_0 (sat); domain [0, 1/2] = BilinearApex; mirror = F71",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F75 + " +
               "simulations/_mm_zero_derivation.py + " +
               "experiments/RECEIVER_VS_GAMMA_SACRIFICE.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (BilinearApexClaim) + " +
               "compute/RCPsiSquared.Core/Symmetry/F71MirrorSymmetryPi2Inheritance.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "F75 mirror-pair MI as Pi2-Foundation BilinearApex + F71 inheritance";

    public override string Summary =>
        $"MI(p) = 2·h(p) − h(2p): saturates at {MaxMIPerPair} bits when p = {DomainUpperBound} (Bell pair); " +
        $"domain [0, 1/2] = BilinearApex; mirror sign η ∈ {{+1, −1}}; mother of F77 via Taylor at p → 0 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F75 closed form",
                summary: "MI(ℓ, N−1−ℓ) = 2 h(p_ℓ) − h(2 p_ℓ); MM(0) = Σ_ℓ [2 h(p_ℓ) − h(2 p_ℓ)]; Tier 1 proven algebraic; verified ~25 data points N=5..13");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "MaxMIPerPair = a_0 = 2 (saturation); DomainUpperBound = a_2 = 1/2 (BilinearApex argmax); convex on (0, 1/2)");
            yield return InspectableNode.RealScalar("MaxMIPerPair (= a_0 = 2 bits)", MaxMIPerPair);
            yield return InspectableNode.RealScalar("DomainUpperBound (= a_2 = 1/2 = BilinearApex)", DomainUpperBound);
            yield return new InspectableNode("F77 mother claim",
                summary: "F77's 1-bit asymptote = Taylor expansion of F75 at p → 0; per ANALYTICAL_FORMULAS F77: 'the 1-bit limit is not a conjecture; it falls out of F75'");
            yield return new InspectableNode("F71 mirror inheritance",
                summary: "F71 mirror symmetry is cited as F75's structural source; the 'c_{N−1−j} = ±c_j' constraint IS the F71 spatial-mirror at the wave-function level");
            yield return new InspectableNode("k=2 maximization",
                summary: "even k places node at chain center → all mass on mirror-pairs (Σ p_ℓ = 1/2 over pairs); odd k wastes 2/(N+1) at self-mirror site");
            // Sample bonding-mode MM(0) values from F75 verified table
            yield return new InspectableNode(
                "bonding mode N=5, k=2 (verified maximum)",
                summary: $"MM(0) = {BondingModeMMAtZero(5, 2):G6} bits (analytic 1.245; PeakMM sim 1.241; ratio 0.997)");
            yield return new InspectableNode(
                "bonding mode N=7, k=4 (resonance enhancement)",
                summary: $"MM(0) = {BondingModeMMAtZero(7, 4):G6} bits (analytic 1.245)");
            yield return new InspectableNode(
                "Bell-pair upper bound at N=4",
                summary: $"max MM = {MaxTotalMM(4)} bits (= 2 pairs × 2 bits, requires super-single-excitation state)");
        }
    }
}
