using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The qudit mirror-protection scaling law (Tier1Derived, 2026-06-17): the fraction of the
/// d^{2N} coherence space that a per-site product mirror can palindrome-protect is
/// P(d, N) / d^{2N}, which is exactly <c>(2/d)^N</c> for every d ≤ 5, so the protection falls with
/// the local dimension d as the N-th power of 2/d, exponentially in the system size N, and the qubit
/// (d = 2) is the UNIQUE dimension with full (= 1) open-system mirror symmetry.
///
/// <para><b>The law.</b> The F121 product-mirror cap is P(d, N) = max_m (2d)^(N−2m)·(d³ − d²)^m
/// (<see cref="QuditProductMirrorCap"/>), which is (2d)^N for d ≤ 5, and the full coherence space is
/// d^{2N}, so
/// <code>
///   protected fraction = (2d)^N / d^{2N} = (2d / d²)^N = (2/d)^N     (d ≤ 5)
/// </code>
/// = 1 ⟺ d = 2 (the cap saturates the space), and (2/3)^N, (2/4)^N, (2/5)^N &lt; 1 for the qudits
/// below d = 6. This is the ratio of the two TERMS of the trunk polynomial d² − 2d = 0 raised to the
/// N-th power: the swap term 2d over the squared-dimension term d² (see the two-family split in
/// <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c> §2 and <c>simulations/qudit_g2_split.py</c>). The two
/// terms coincide only at the root d = 2, which is exactly why the fraction is 1 only there: the qubit
/// magic is the trunk's root. From d = 6 on, pairing a dark-only site with a lit-only site beats the
/// swap, and the fraction becomes ((d − 1)/d²)^(N/2) at even N and (2/d)·((d − 1)/d²)^((N−1)/2) at odd
/// N: still exponential in N, still below 1, and still falling with d.</para>
///
/// <para><b>The contrast (what makes the law physical).</b> The protection falls with d, but the DECAY
/// RATES do not: the dissipator ladder is 2γ·Hamming and the structural ceiling g2(K_N) = 4/N are both
/// d-INDEPENDENT (Hamming distance and the S_N principal angle carry no d; gate-verified at d = 3 in
/// <c>simulations/qudit_g2_split.py</c>, validated against the full d = 3 Liouvillian). So a qudit decays
/// at the same rates as a qubit but loses product-mirror protection as (2/d)^N (d ≤ 5; faster from
/// d = 6). The qubit is the
/// only carrier whose open-system mirror is complete, for product mirrors and for any mirror of the
/// full-Cartan dephasing dissipator (F121's ceiling is full only at d = 2 too).</para>
///
/// <para><b>Why Tier1Derived.</b> Pure composition of the single Tier1Derived parent
/// <see cref="QuditProductMirrorCap"/> (the cap P(d, N)): the total d^{2N} is the Liouville-space
/// dimension and the fraction is exact integer arithmetic. The d-independence of the rates is the
/// complementary observation (anchored to the qutrit verifier), not the load-bearing content of this
/// fraction.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md</c> §6 (the cap) +
/// <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c> §2 (the two-family split of the four) +
/// <c>simulations/qudit_g2_split.py</c> (the rates-d-independent half, gate-first) +
/// <c>simulations/qudit_product_mirror_cap.py</c> (the cap, F121).</para></summary>
public sealed class QuditMirrorProtectionScalingClaim : Claim
{
    /// <summary>Parent: the F121 product-mirror cap P(d, N). The protected fraction is this cap divided
    /// by the full coherence space d^{2N}. Cited, not re-derived.</summary>
    public QuditProductMirrorCap Cap { get; }

    /// <summary>The product-mirror cap P(d, N) (= (2d)^N for d ≤ 5) for local dimension d, system
    /// size N. Requires d ≥ 2, N ≥ 1.</summary>
    public static long ProductCap(int d, int n) => QuditProductMirrorCap.ProductCap(d, n);

    /// <summary>The full coherence space dimension d^{2N} = (d²)^N. Requires d ≥ 2, N ≥ 1.</summary>
    public static long TotalCoherences(int d, int n)
    {
        if (d < 2) throw new ArgumentOutOfRangeException(nameof(d), $"local dimension d must be ≥ 2; got {d}");
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 1; got {n}");
        return QuditPartialPalindromeCeiling.Total(d, n);   // d^{2N}, the canonical full Liouville dim (one node up)
    }

    /// <summary>The product-mirror-protected fraction P(d, N) / d^{2N}; (2/d)^N for d ≤ 5.
    /// Requires d ≥ 2, N ≥ 1.</summary>
    public static double ProtectedFraction(int d, int n) => (double)ProductCap(d, n) / TotalCoherences(d, n);

    /// <summary>(2/d)^N evaluated directly from the closed form (the d ≤ 5 face of ProtectedFraction).</summary>
    public static double ProtectedFractionClosedForm(int d, int n)
    {
        if (d < 2) throw new ArgumentOutOfRangeException(nameof(d), $"local dimension d must be ≥ 2; got {d}");
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 1; got {n}");
        return Math.Pow(2.0 / d, n);
    }

    /// <summary>Exact integer check of the d ≤ 5 law: P(d, N)·d^N == 2^N·d^{2N}, i.e.
    /// P / d^{2N} = (2/d)^N with no floating division.</summary>
    public static bool FractionIsTwoOverDToTheN(int d, int n) =>
        checked(ProductCap(d, n) * QuditPartialPalindromeCeiling.IntPow(d, n))
        == checked(QuditPartialPalindromeCeiling.IntPow(2, n) * TotalCoherences(d, n));

    /// <summary>The mirror is full (protected fraction = 1) iff d = 2: the qubit is the unique full-mirror
    /// local dimension (the trunk root d² − 2d = 0).</summary>
    public static bool IsFullMirror(int d) => d == 2;

    public QuditMirrorProtectionScalingClaim(QuditProductMirrorCap cap)
        : base("Qudit mirror-protection scaling: the per-site product mirror protects a fraction " +
               "P(d, N) / d^{2N} of the coherence space, = (2d)^N / d^{2N} = (2/d)^N for d ≤ 5, falling " +
               "as the N-th power of 2/d, exponentially in N; = 1 ⟺ d = 2, so the qubit is the UNIQUE dimension " +
               "with full open-system mirror symmetry. The swap term 2d and the space term d² are the two " +
               "terms of the trunk d² − 2d = 0 (raised to N), equal only at the root d = 2; from d = 6 on the " +
               "cap is larger than (2d)^N and the fraction is ((d − 1)/d²)^(N/2) at even N, still exponential. " +
               "Complementary (qutrit-verified): the decay rates 2γ·Hamming and the structural ceiling 4/N " +
               "are d-INDEPENDENT, so a qudit decays like a qubit but loses palindrome protection as (2/d)^N.",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md + " +
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md (§2 the two-family split) + " +
               "simulations/qudit_g2_split.py + " +
               "simulations/qudit_product_mirror_cap.py")
    {
        Cap = cap ?? throw new ArgumentNullException(nameof(cap));
    }

    public override string DisplayName =>
        "Qudit mirror-protection scaling: protected fraction = (2/d)^N for d ≤ 5, full iff d=2 (the qubit alone)";

    public override string Summary =>
        $"protected fraction = P(d, N) / d^{{2N}} = (2/d)^N for d ≤ 5; = 1 ⟺ d=2 (qubit = unique full mirror), " +
        $"(2/3)^N at the qutrit; rates 2γ·Hamming + ceiling 4/N d-independent ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the law",
                summary: "protected fraction = cap P(d, N) / space d^{2N} = (2/d)^N for d ≤ 5 = (swap term 2d / " +
                         "squared-dim term d²)^N; = 1 ⟺ d=2 (the trunk root d²−2d=0). Mirror protection falls " +
                         "as the N-th power of 2/d, exponentially in N; from d = 6 the cap exceeds (2d)^N and the fraction is " +
                         "((d−1)/d²)^(N/2) at even N.");
            foreach (var (d, n) in new[] { (2, 2), (2, 3), (3, 2), (3, 3), (4, 2), (6, 2) })
            {
                double frac = ProtectedFraction(d, n);
                string full = IsFullMirror(d) ? "FULL (qubit)" : "partial";
                string law = d <= 5
                    ? $"= (2/{d})^{n} = {ProtectedFractionClosedForm(d, n):0.######}"
                    : $"> (2/{d})^{n} = {ProtectedFractionClosedForm(d, n):0.######} (the (0, 2) pairing beats the swap)";
                yield return new InspectableNode($"d={d}, N={n}",
                    summary: $"cap P = {ProductCap(d, n)}, space d^{{2N}} = {TotalCoherences(d, n)}, " +
                             $"protected = {frac:0.######} {law}; {full}");
            }
            yield return new InspectableNode("the contrast (rates stay): d-independence",
                summary: "the dissipator ladder 2γ·Hamming and the structural ceiling g2(K_N)=4/N are d-INDEPENDENT " +
                         "(Hamming distance + S_N principal angle carry no d; gate-verified at d=3, " +
                         "simulations/qudit_g2_split.py). A qudit decays like a qubit but loses mirror protection as (2/d)^N.");
            yield return new InspectableNode("the qubit-necessity reading",
                summary: "full mirror ⟺ P(d, N) = d^{2N} ⟺ d² − 2d = 0 ⟺ d = 2: the open-system palindrome is " +
                         "complete ONLY for qubits, for product mirrors and (F121's ceiling) for any mirror of the " +
                         "full-Cartan dephasing dissipator. " +
                         "Two-level carriers are privileged (the trunk's root is the qubit).");
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
