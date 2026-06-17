using System.Collections.Generic;
using System.Globalization;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The qudit mirror-protection scaling law computed live (typed home:
/// <c>QuditMirrorProtectionScalingClaim</c>). The per-site product mirror protects a fraction
/// (2d)^N / d^{2N} = (2/d)^N of the coherence space, decaying exponentially in the local dimension d;
/// = 1 ⟺ d = 2, so the qubit is the unique dimension with full open-system mirror symmetry. This witness
/// recomputes the cap (2d)^N live via <see cref="QuditProductMirrorCap.ProductCap"/>, divides by the full
/// space d^{2N}, and confirms the closed form (2/d)^N and the full-iff-d=2 gate. The complementary half (the
/// decay rates 2γ·Hamming and the structural ceiling 4/N are d-INDEPENDENT) is gate-verified in
/// <c>simulations/qudit_g2_split.py</c>.</summary>
public sealed class QuditMirrorProtectionWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private static readonly (int D, int N)[] Grid =
        { (2, 2), (2, 3), (2, 4), (3, 2), (3, 3), (4, 2), (5, 2) };

    public string DisplayName => "QuditMirrorProtectionWitness (protected fraction = (2/d)^N, full iff d=2)";

    public string Summary =>
        "the qudit mirror-protection scaling law computed live (typed home: QuditMirrorProtectionScalingClaim): " +
        "the per-site product mirror protects (2d)^N / d^{2N} = (2/d)^N of the coherence space, decaying " +
        "exponentially in d; = 1 ⟺ d=2, so the qubit is the unique full-mirror dimension. Rates stay d-independent.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("the law",
                summary: "protected fraction = cap (2d)^N / space d^{2N} = (2/d)^N = (cap term 2d / squared-dim " +
                         "term d²)^N; = 1 ⟺ d=2 (the trunk root d²−2d=0). Exponential decay of mirror protection in d.");

            bool allMatch = true;
            foreach (var (d, n) in Grid)
            {
                long cap = QuditMirrorProtectionScalingClaim.ProductCap(d, n);
                long total = QuditMirrorProtectionScalingClaim.TotalCoherences(d, n);
                double frac = QuditMirrorProtectionScalingClaim.ProtectedFraction(d, n);
                double cf = QuditMirrorProtectionScalingClaim.ProtectedFractionClosedForm(d, n);
                bool full = QuditMirrorProtectionScalingClaim.IsFullMirror(d);
                bool match = System.Math.Abs(frac - cf) < 1e-12 && (System.Math.Abs(frac - 1.0) < 1e-12) == full;
                allMatch &= match;
                yield return new InspectableNode($"d={d}, N={n}",
                    summary: $"cap (2d)^N = {cap.ToString(Inv)}, space d^{{2N}} = {total.ToString(Inv)}, " +
                             $"protected = {frac.ToString("0.######", Inv)} = (2/{d})^{n} = {cf.ToString("0.######", Inv)} " +
                             $"({(match ? "match" : "MISMATCH")}); {(full ? "FULL (qubit, the unique full mirror)" : "partial")}");
            }

            yield return new InspectableNode("gate: fraction == (2/d)^N and full ⟺ d=2",
                summary: allMatch
                    ? "PASS — every grid point: (2d)^N / d^{2N} equals (2/d)^N to machine precision, and the " +
                      "fraction is 1 exactly at d=2 (and below 1 for every qudit)."
                    : "FAIL — a grid point disagrees (see the MISMATCH row above).");

            yield return new InspectableNode("the contrast (rates stay): d-independence",
                summary: "the dissipator ladder 2γ·Hamming and the structural ceiling g2(K_N)=4/N are d-INDEPENDENT " +
                         "(Hamming distance + S_N principal angle carry no d; gate-verified at d=3, " +
                         "simulations/qudit_g2_split.py). A qudit decays like a qubit but loses mirror protection (2/d)^N.");

            yield return new InspectableNode("the qubit-necessity reading",
                summary: "full mirror ⟺ (2d)^N = d^{2N} ⟺ d² − 2d = 0 ⟺ d = 2: the open-system palindrome is " +
                         "complete ONLY for qubits; two-level carriers are privileged (the trunk's root is the qubit).");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
