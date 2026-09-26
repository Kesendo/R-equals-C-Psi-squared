using System.Collections.Generic;
using System.Globalization;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The qudit mirror-protection scaling law computed live (typed home:
/// <c>QuditMirrorProtectionScalingClaim</c>). The per-site product mirror protects a fraction
/// P(d, N) / d^{2N} of the coherence space, which is (2d)^N / d^{2N} = (2/d)^N for every d ≤ 5,
/// the N-th power of 2/d, exponential in N; = 1 ⟺ d = 2, so the qubit is the unique
/// dimension with full open-system mirror symmetry. This witness recomputes the cap P(d, N) live via
/// <see cref="QuditProductMirrorCap.ProductCap"/>, checks the integer identity P·d^N = 2^N·d^{2N}
/// (the fraction is (2/d)^N with no floating division) on the d ≤ 5 grid, shows the d = 6 row where
/// the (0, 2) pairing beats the swap, and confirms the full-iff-d=2 gate. The complementary half (the
/// decay rates 2γ·Hamming and the structural ceiling 4/N are d-INDEPENDENT) is gate-verified in
/// <c>simulations/qudit_g2_split.py</c>.</summary>
public sealed class QuditMirrorProtectionWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private static readonly (int D, int N)[] Grid =
        { (2, 2), (2, 3), (2, 4), (3, 2), (3, 3), (4, 2), (5, 2) };

    public string DisplayName => "QuditMirrorProtectionWitness (protected fraction = (2/d)^N for d ≤ 5, full iff d=2)";

    public string Summary =>
        "the qudit mirror-protection scaling law computed live (typed home: QuditMirrorProtectionScalingClaim): " +
        "the per-site product mirror protects P(d, N) / d^{2N} = (2/d)^N of the coherence space for d ≤ 5, " +
        "falling as the N-th power of 2/d; = 1 ⟺ d=2, so the qubit is the unique full-mirror dimension. From d = 6 " +
        "the cap exceeds (2d)^N (180 > 144 at N = 2) and the fraction stays exponential. Rates stay d-independent.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("the law",
                summary: "protected fraction = cap P(d, N) / space d^{2N} = (2/d)^N for d ≤ 5 = (swap term 2d / " +
                         "squared-dim term d²)^N; = 1 ⟺ d=2 (the trunk root d²−2d=0). Mirror protection falls as " +
                         "the N-th power of 2/d, exponentially in N.");

            bool allMatch = true;
            foreach (var (d, n) in Grid)
            {
                long cap = QuditMirrorProtectionScalingClaim.ProductCap(d, n);
                long total = QuditMirrorProtectionScalingClaim.TotalCoherences(d, n);
                double frac = QuditMirrorProtectionScalingClaim.ProtectedFraction(d, n);
                bool law = QuditMirrorProtectionScalingClaim.FractionIsTwoOverDToTheN(d, n);
                bool full = QuditMirrorProtectionScalingClaim.IsFullMirror(d);
                bool match = law && (cap == total) == full;
                allMatch &= match;
                yield return new InspectableNode($"d={d}, N={n}",
                    summary: $"cap P = (2d)^N = {cap.ToString(Inv)}, space d^{{2N}} = {total.ToString(Inv)}, " +
                             $"protected = {frac.ToString("0.######", Inv)} = (2/{d})^{n} " +
                             $"({(match ? "exact integer match" : "MISMATCH")}); {(full ? "FULL (qubit, the unique full mirror)" : "partial")}");
            }

            long cap6 = QuditMirrorProtectionScalingClaim.ProductCap(6, 2);
            long total6 = QuditMirrorProtectionScalingClaim.TotalCoherences(6, 2);
            bool beyond = cap6 == 180 && !QuditMirrorProtectionScalingClaim.FractionIsTwoOverDToTheN(6, 2);
            allMatch &= beyond;
            yield return new InspectableNode("d=6, N=2 (past the threshold)",
                summary: $"cap P = {cap6.ToString(Inv)} > (2d)^N = 144 (a dark-only site beside a lit-only site: " +
                         $"6·30), space {total6.ToString(Inv)}, protected = " +
                         $"{QuditMirrorProtectionScalingClaim.ProtectedFraction(6, 2).ToString("0.######", Inv)} = (d−1)/d² " +
                         $"> (2/6)^2 ({(beyond ? "as the lemma says" : "MISMATCH")}); partial");

            yield return new InspectableNode("gate: fraction == (2/d)^N for d ≤ 5, the d = 6 row above it, full ⟺ d=2",
                summary: allMatch
                    ? "PASS: every d ≤ 5 grid point satisfies P·d^N = 2^N·d^{2N} exactly, the fraction is 1 exactly " +
                      "at d=2 (and below 1 for every qudit), and the d = 6 row sits above (2/d)^N."
                    : "FAIL: a grid point disagrees (see the MISMATCH row above).");

            yield return new InspectableNode("the contrast (rates stay): d-independence",
                summary: "the dissipator ladder 2γ·Hamming and the structural ceiling g2(K_N)=4/N are d-INDEPENDENT " +
                         "(Hamming distance + S_N principal angle carry no d; gate-verified at d=3, " +
                         "simulations/qudit_g2_split.py). A qudit decays like a qubit but loses mirror protection as (2/d)^N.");

            yield return new InspectableNode("the qubit-necessity reading",
                summary: "full mirror ⟺ P(d, N) = d^{2N} ⟺ d² − 2d = 0 ⟺ d = 2: the open-system palindrome is " +
                         "complete ONLY for qubits, for product mirrors and (F121's ceiling) for any mirror of the " +
                         "full-Cartan dephasing dissipator; two-level carriers are privileged (the trunk's root is the qubit).");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
