using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The restricted map Π_d P_aligned's coverage fraction (Tier1Derived, corrected
/// 2026-09-06) is <c>(2/d)^N</c>. This is not the optimum over all product intertwiners;
/// the former universal-cap reading is retracted by the d=6,N=2 projector counterexample.
///
/// <para><b>The construction.</b> The F121 shift-aligned rank is (2d)^N (<see cref="QuditProductMirrorCap"/>) and
/// the full coherence space is d^{2N}, so
/// <code>
///   protected fraction = (2d)^N / d^{2N} = (2d / d²)^N = (2/d)^N
/// </code>
/// = 1 at d = 2 for this construction, and (2/3)^N, (2/4)^N, … &lt; 1 for every qudit. This is the
/// ratio of the two terms 2d and d² raised to the N-th power, but not a universal cap theorem: the term 2d
/// over the squared-dimension term d² (see the two-family split in
/// <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c> §2 and <c>simulations/qudit_g2_split.py</c>). The two
/// terms coincide only at the root d = 2, which is exactly why the fraction is 1 only there: the qubit
/// magic is the trunk's root.</para>
///
/// <para><b>The contrast.</b> This construction's coverage falls with d, but the DECAY
/// RATES do not: the dissipator ladder is 2γ·Hamming and the structural ceiling g2(K_N) = 4/N are both
/// d-INDEPENDENT (Hamming distance and the S_N principal angle carry no d; gate-verified at d = 3 in
/// <c>simulations/qudit_g2_split.py</c>, validated against the full d = 3 Liouvillian). So a qudit decays
/// at the same rates as a qubit while Π_d's shift-aligned coverage falls as (2/d)^N.</para>
///
/// <para><b>Why Tier1Derived.</b> Pure composition of the single Tier1Derived parent
/// <see cref="QuditProductMirrorCap"/> (the shift-aligned rank (2d)^N): the total d^{2N} is the Liouville-space dimension
/// and the fraction is exact integer arithmetic. The d-independence of the rates is the complementary
/// observation (anchored to the qutrit verifier), not the load-bearing content of this fraction.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md</c> (construction and cap retraction) +
/// <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c> §2 (the two-family split of the four) +
/// <c>simulations/qudit_g2_split.py</c> (the rates-d-independent half, gate-first) +
/// <c>simulations/qudit_product_mirror_cap.py</c> (the cap, F121).</para></summary>
public sealed class QuditMirrorProtectionScalingClaim : Claim
{
    /// <summary>Parent with the historical cap name; now carries the explicit shift-aligned
    /// rank and the universal-cap retraction.</summary>
    public QuditProductMirrorCap Cap { get; }

    /// <summary>Historical API name: rank(Π_d P_aligned) = (2d)^N, not a universal cap.</summary>
    public static long ProductCap(int d, int n) => QuditProductMirrorCap.ProductCap(d, n);

    /// <summary>The full coherence space dimension d^{2N} = (d²)^N. Requires d ≥ 2, N ≥ 1.</summary>
    public static long TotalCoherences(int d, int n)
    {
        if (d < 2) throw new ArgumentOutOfRangeException(nameof(d), $"local dimension d must be ≥ 2; got {d}");
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 1; got {n}");
        return QuditPartialPalindromeCeiling.Total(d, n);   // d^{2N}, the canonical full Liouville dim (one node up)
    }

    /// <summary>The shift-aligned construction's fraction (2/d)^N. Not an optimum.</summary>
    public static double ProtectedFraction(int d, int n) => (double)ProductCap(d, n) / TotalCoherences(d, n);

    /// <summary>(2/d)^N evaluated directly from the closed form (the cross-check for ProtectedFraction).</summary>
    public static double ProtectedFractionClosedForm(int d, int n)
    {
        if (d < 2) throw new ArgumentOutOfRangeException(nameof(d), $"local dimension d must be ≥ 2; got {d}");
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 1; got {n}");
        return Math.Pow(2.0 / d, n);
    }

    /// <summary>The mirror is full (protected fraction = 1) iff d = 2: the qubit is the unique full-mirror
    /// local dimension (the trunk root d² − 2d = 0).</summary>
    public static bool IsFullMirror(int d) => d == 2;

    public QuditMirrorProtectionScalingClaim(QuditProductMirrorCap cap)
        : base("Qudit shift-aligned construction coverage: Π_d P_aligned protects a fraction " +
               "(2d)^N / d^{2N} = (2/d)^N of the coherence space, decaying exponentially in the local " +
               "dimension d within this construction. This is not a universal product optimum: " +
               "P_dark⊗P_lit at d=6,N=2 has rank 180 > 144. Complementary (qutrit-verified): " +
               "the decay rates 2γ·Hamming and the structural ceiling 4/N are d-INDEPENDENT, so a qudit " +
               "shares the same rate ladder while this construction's coverage is (2/d)^N.",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md + " +
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md (§2 the two-family split) + " +
               "simulations/qudit_g2_split.py + " +
               "simulations/qudit_product_mirror_cap.py")
    {
        Cap = cap ?? throw new ArgumentNullException(nameof(cap));
    }

    public override string DisplayName =>
        "Qudit Π_d P_aligned coverage = (2/d)^N (not a universal product cap)";

    public override string Summary =>
        $"Π_d shift-aligned fraction = (2d)^N / d^{{2N}} = (2/d)^N; not the optimum over product mirrors; " +
        $"(2/3)^N at the qutrit; rates 2γ·Hamming + ceiling 4/N d-independent ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the law",
                summary: "Π_d shift-aligned fraction = (2d)^N / d^{2N} = (2/d)^N; " +
                         "this construction's coverage decays exponentially in d, but it is not a universal product cap.");
            foreach (var (d, n) in new[] { (2, 2), (2, 3), (3, 2), (3, 3), (4, 2) })
            {
                double frac = ProtectedFraction(d, n);
                string full = IsFullMirror(d) ? "FULL (qubit)" : "partial";
                yield return new InspectableNode($"d={d}, N={n}",
                    summary: $"rank(Π_d P_aligned) = (2d)^N = {ProductCap(d, n)}, space d^{{2N}} = {TotalCoherences(d, n)}, " +
                             $"construction coverage = {frac:0.######} = (2/{d})^{n} = {ProtectedFractionClosedForm(d, n):0.######} — {full}");
            }
            yield return new InspectableNode("the contrast (rates stay): d-independence",
                summary: "the dissipator ladder 2γ·Hamming and the structural ceiling g2(K_N)=4/N are d-INDEPENDENT " +
                         "(Hamming distance + S_N principal angle carry no d; gate-verified at d=3, " +
                         "simulations/qudit_g2_split.py). This does not promote Π_d's coverage to a universal optimum.");
            yield return new InspectableNode("the qubit-necessity reading",
                summary: "the explicit Π_d construction fills the space iff d=2. Independently, F121's combinatorial " +
                         "ceiling is full only at d=2; neither statement restores the retracted universal product cap.");
            yield return Cap;   // typed parent edge (Tier1Derived)
        }
    }

    public static QuditMirrorProtectionScalingClaim Build()
    {
        var qubitNecessity = new QubitNecessityPi2Inheritance(
            new Pi2DyadicLadderClaim(), new Pi2OperatorSpaceMirrorClaim());
        var cap = new QuditProductMirrorCap(
            new QuditPartialPalindromeCeiling(qubitNecessity), qubitNecessity);
        return new QuditMirrorProtectionScalingClaim(cap);
    }

    public static QuditMirrorProtectionScalingClaim Shared { get; } = Build();
}
