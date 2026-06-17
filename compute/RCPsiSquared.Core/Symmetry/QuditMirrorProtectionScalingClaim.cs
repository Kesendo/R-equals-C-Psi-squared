using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The qudit mirror-protection scaling law (Tier1Derived, 2026-06-17): the fraction of the
/// d^{2N} coherence space that a per-site product mirror can palindrome-protect is exactly
/// <c>(2/d)^N</c>, so the protection decays EXPONENTIALLY in the local dimension d and the qubit (d = 2)
/// is the UNIQUE dimension with full (= 1) open-system mirror symmetry.
///
/// <para><b>The law.</b> The F121 product-mirror cap is (2d)^N (<see cref="QuditProductMirrorCap"/>) and
/// the full coherence space is d^{2N}, so
/// <code>
///   protected fraction = (2d)^N / d^{2N} = (2d / d²)^N = (2/d)^N
/// </code>
/// = 1 ⟺ d = 2 (the cap saturates the space), and (2/3)^N, (2/4)^N, … &lt; 1 for every qudit. This is the
/// ratio of the two TERMS of the trunk polynomial d² − 2d = 0 raised to the N-th power: the cap term 2d
/// over the squared-dimension term d² (see the two-family split in
/// <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c> §2 and <c>simulations/qudit_g2_split.py</c>). The two
/// terms coincide only at the root d = 2, which is exactly why the fraction is 1 only there: the qubit
/// magic is the trunk's root.</para>
///
/// <para><b>The contrast (what makes the law physical).</b> The protection falls with d, but the DECAY
/// RATES do not: the dissipator ladder is 2γ·Hamming and the structural ceiling g2(K_N) = 4/N are both
/// d-INDEPENDENT (Hamming distance and the S_N principal angle carry no d; gate-verified at d = 3 in
/// <c>simulations/qudit_g2_split.py</c>, validated against the full d = 3 Liouvillian). So a qudit decays
/// at the same rates as a qubit but loses palindrome protection exponentially as (2/d)^N. The qubit is the
/// only carrier whose open-system mirror is complete.</para>
///
/// <para><b>Why Tier1Derived.</b> Pure composition of the single Tier1Derived parent
/// <see cref="QuditProductMirrorCap"/> (the cap (2d)^N): the total d^{2N} is the Liouville-space dimension
/// and the fraction is exact integer arithmetic. The d-independence of the rates is the complementary
/// observation (anchored to the qutrit verifier), not the load-bearing content of this fraction.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md</c> (the cap (2d)^N) +
/// <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c> §2 (the two-family split of the four) +
/// <c>simulations/qudit_g2_split.py</c> (the rates-d-independent half, gate-first) +
/// <c>simulations/qudit_product_mirror_cap.py</c> (the cap, F121).</para></summary>
public sealed class QuditMirrorProtectionScalingClaim : Claim
{
    /// <summary>Parent: the F121 product-mirror cap (2d)^N. The protected fraction is this cap divided by
    /// the full coherence space d^{2N}. Cited, not re-derived.</summary>
    public QuditProductMirrorCap Cap { get; }

    /// <summary>The product-mirror cap (2d)^N for local dimension d, system size N. Requires d ≥ 2, N ≥ 1.</summary>
    public static long ProductCap(int d, int n) => QuditProductMirrorCap.ProductCap(d, n);

    /// <summary>The full coherence space dimension d^{2N} = (d²)^N. Requires d ≥ 2, N ≥ 1.</summary>
    public static long TotalCoherences(int d, int n)
    {
        if (d < 2) throw new ArgumentOutOfRangeException(nameof(d), $"local dimension d must be ≥ 2; got {d}");
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 1; got {n}");
        return QuditPartialPalindromeCeiling.IntPow((long)d * d, n);
    }

    /// <summary>The product-mirror-protected fraction (2/d)^N = (2d)^N / d^{2N}. Requires d ≥ 2, N ≥ 1.</summary>
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
        : base("Qudit mirror-protection scaling: the per-site product mirror protects a fraction " +
               "(2d)^N / d^{2N} = (2/d)^N of the coherence space, decaying exponentially in the local " +
               "dimension d; = 1 ⟺ d = 2, so the qubit is the UNIQUE dimension with full open-system " +
               "mirror symmetry. The cap (2d)^N and the space d^{2N} are the two terms of the trunk " +
               "d² − 2d = 0 (raised to N), equal only at the root d = 2. Complementary (qutrit-verified): " +
               "the decay rates 2γ·Hamming and the structural ceiling 4/N are d-INDEPENDENT, so a qudit " +
               "decays like a qubit but loses palindrome protection as (2/d)^N.",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md + " +
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md (§2 the two-family split) + " +
               "simulations/qudit_g2_split.py + " +
               "simulations/qudit_product_mirror_cap.py")
    {
        Cap = cap ?? throw new ArgumentNullException(nameof(cap));
    }

    public override string DisplayName =>
        "Qudit mirror-protection scaling: protected fraction = (2/d)^N, full iff d=2 (the qubit alone)";

    public override string Summary =>
        $"protected fraction = (2d)^N / d^{{2N}} = (2/d)^N; = 1 ⟺ d=2 (qubit = unique full mirror), " +
        $"(2/3)^N at the qutrit; rates 2γ·Hamming + ceiling 4/N d-independent ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the law",
                summary: "protected fraction = (2d)^N / d^{2N} = (2/d)^N = (cap term 2d / squared-dim term d²)^N; " +
                         "= 1 ⟺ d=2 (the trunk root d²−2d=0). Mirror protection decays exponentially in d.");
            foreach (var (d, n) in new[] { (2, 2), (2, 3), (3, 2), (3, 3), (4, 2) })
            {
                double frac = ProtectedFraction(d, n);
                string full = IsFullMirror(d) ? "FULL (qubit)" : "partial";
                yield return new InspectableNode($"d={d}, N={n}",
                    summary: $"cap (2d)^N = {ProductCap(d, n)}, space d^{{2N}} = {TotalCoherences(d, n)}, " +
                             $"protected = {frac:0.######} = (2/{d})^{n} = {ProtectedFractionClosedForm(d, n):0.######} — {full}");
            }
            yield return new InspectableNode("the contrast (rates stay): d-independence",
                summary: "the dissipator ladder 2γ·Hamming and the structural ceiling g2(K_N)=4/N are d-INDEPENDENT " +
                         "(Hamming distance + S_N principal angle carry no d; gate-verified at d=3, " +
                         "simulations/qudit_g2_split.py). A qudit decays like a qubit but loses mirror protection (2/d)^N.");
            yield return new InspectableNode("the qubit-necessity reading",
                summary: "full mirror ⟺ (2d)^N = d^{2N} ⟺ d² − 2d = 0 ⟺ d = 2: the open-system palindrome is " +
                         "complete ONLY for qubits. Two-level carriers are privileged (the trunk's root is the qubit).");
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
