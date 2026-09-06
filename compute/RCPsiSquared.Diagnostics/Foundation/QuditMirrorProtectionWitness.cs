using System.Collections.Generic;
using System.Globalization;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The explicit Π_d P_aligned coverage computed live (typed home:
/// <c>QuditMirrorProtectionScalingClaim</c>). This witness recomputes its rank (2d)^N,
/// divides by d^{2N}, and confirms (2/d)^N. It is not the optimum over product intertwiners;
/// the d=6,N=2 rank-180 projector retracts that reading. The complementary half (the
/// decay rates 2γ·Hamming and the structural ceiling 4/N are d-INDEPENDENT) is gate-verified in
/// <c>simulations/qudit_g2_split.py</c>.</summary>
public sealed class QuditMirrorProtectionWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private static readonly (int D, int N)[] Grid =
        { (2, 2), (2, 3), (2, 4), (3, 2), (3, 3), (4, 2), (5, 2) };

    public string DisplayName => "QuditMirrorProtectionWitness (Π_d shift coverage = (2/d)^N; universal cap retracted)";

    public string Summary =>
        "The restriction Π_d P_aligned gives coverage (2d)^N / d^{2N} = (2/d)^N. This is not a universal " +
        "product optimum; P_dark⊗P_lit at d=6,N=2 has rank 180 > 144. Rates stay d-independent.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("the construction",
                summary: "Π_d P_aligned coverage = (2d)^N / d^{2N} = (2/d)^N; not a universal product cap.");

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
                    summary: $"rank(Π_d P_aligned) = (2d)^N = {cap.ToString(Inv)}, space d^{{2N}} = {total.ToString(Inv)}, " +
                             $"coverage = {frac.ToString("0.######", Inv)} = (2/{d})^{n} = {cf.ToString("0.######", Inv)} " +
                             $"({(match ? "match" : "MISMATCH")}); {(full ? "FULL (qubit, the unique full mirror)" : "partial")}");
            }

            yield return new InspectableNode("gate: Π_d construction coverage == (2/d)^N",
                summary: allMatch
                    ? "PASS — every grid point: (2d)^N / d^{2N} equals (2/d)^N to machine precision, and the " +
                      "construction fills the space at d=2 (and is partial for every tested qudit)."
                    : "FAIL — a grid point disagrees (see the MISMATCH row above).");

            yield return new InspectableNode("the contrast (rates stay): d-independence",
                summary: "the dissipator ladder 2γ·Hamming and the structural ceiling g2(K_N)=4/N are d-INDEPENDENT " +
                         "(Hamming distance + S_N principal angle carry no d; gate-verified at d=3, " +
                         "simulations/qudit_g2_split.py). The rate statement does not make the coverage an optimum.");

            yield return new InspectableNode("the qubit-necessity reading",
                summary: "Π_d fills the space at d=2; independently F121's combinatorial ceiling is full only at d=2. " +
                         "The former universal product-cap inference is retracted.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
