using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live witness for the UNIQUENESS_PROOF argument (2026-07-02): CΨ = ¼ is the UNIQUE
/// bifurcation boundary, and α = 2 is the UNIQUE Rényi order whose fold threshold is
/// state-independent. The value ¼ is already typed (QuarterAsBilinearMaxvalClaim,
/// PolynomialDiscriminantAnchorClaim); this witness recomputes the two ARGUMENTS that single it out
/// (typed as <see cref="RCPsiSquared.Core.Symmetry.QuarterBoundaryUniquenessClaim"/>).
///
/// <para>Argument 1 (the load-bearing forcing, Step 6). The generalized recursion R = C(Ψ + R)^α has
/// fold-bifurcation threshold CΨ*_α = (α−1)^(α−1) / (α^α · Ψ^(α−2)). The Ψ^(α−2) factor makes this
/// threshold state-DEPENDENT for every α except α = 2, where the exponent α−2 = 0 kills the Ψ
/// dependence and the value is exactly ¼. So α = 2 is the unique Rényi order with a state-independent
/// threshold, and ¼ is what it pins. The witness sweeps α and reports the threshold at two distinct
/// states Ψ; the spread |CΨ*_α(Ψ_lo) − CΨ*_α(Ψ_hi)| is 0 ONLY at α = 2.</para>
///
/// <para>Argument 2 (the discriminant boundary, Steps 1-3). Expanding R = C(Ψ + R)² gives
/// C·R² + (2CΨ − 1)·R + CΨ² = 0, whose discriminant is D = (2CΨ − 1)² − 4C²Ψ² = 1 − 4CΨ. So D = 0
/// (the fold tangency) iff CΨ = ¼; D &gt; 0 (two real fixed points) iff CΨ &lt; ¼; D &lt; 0 (none) iff
/// CΨ &gt; ¼. The witness reports D across CΨ, its single sign change at ¼ being the unique boundary.</para>
///
/// <para>Both cruxes are elementary arithmetic (no matrices), so the witness is exact and instant.
/// Anchors: <c>docs/proofs/UNIQUENESS_PROOF.md</c> + <c>simulations/review2_A3_renyi.py</c>.</para></summary>
public sealed class QuarterBoundaryUniquenessWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The unique bifurcation value CΨ = 1/4, pinned by α = 2.</summary>
    public const double QuarterBoundary = 0.25;

    public double PsiLow { get; }
    public double PsiHigh { get; }

    /// <summary>Per-α reading: the Rényi fold threshold at the two probe states and their spread
    /// (the state-dependence, 0 iff the threshold is state-independent).</summary>
    public sealed record AlphaRow(double Alpha, double ThresholdLow, double ThresholdHigh)
    {
        /// <summary>State-dependence: |CΨ*_α(Ψ_lo) − CΨ*_α(Ψ_hi)|. Exactly 0 iff α = 2.</summary>
        public double Spread => System.Math.Abs(ThresholdLow - ThresholdHigh);
    }

    public IReadOnlyList<AlphaRow> Rows { get; }

    /// <summary>The Rényi-α fold-bifurcation (tangency) threshold of R = C(Ψ + R)^α:
    /// CΨ*_α = (α−1)^(α−1) / (α^α · Ψ^(α−2)). State-independent iff α = 2 (where it equals 1/4).</summary>
    public static double CriticalThreshold(double alpha, double psi) =>
        System.Math.Pow(alpha - 1.0, alpha - 1.0) / (System.Math.Pow(alpha, alpha) * System.Math.Pow(psi, alpha - 2.0));

    /// <summary>The discriminant of the α = 2 fixed-point quadratic C·R² + (2CΨ−1)·R + CΨ² = 0:
    /// D = 1 − 4CΨ. D = 0 iff CΨ = 1/4 (the fold tangency); &gt; 0 below, &lt; 0 above.</summary>
    public static double Discriminant(double cpsi) => 1.0 - 4.0 * cpsi;

    private static readonly double[] AlphaGrid = { 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 4.0 };

    public QuarterBoundaryUniquenessWitness(double psiLow = 0.3, double psiHigh = 0.7)
    {
        if (psiLow is <= 0.0 or >= 1.0)
            throw new ArgumentOutOfRangeException(nameof(psiLow), psiLow, "probe state Ψ must lie in (0, 1).");
        if (psiHigh is <= 0.0 or >= 1.0)
            throw new ArgumentOutOfRangeException(nameof(psiHigh), psiHigh, "probe state Ψ must lie in (0, 1).");
        if (System.Math.Abs(psiLow - psiHigh) < 1e-9)
            throw new ArgumentException("psiLow and psiHigh must be distinct so the state-dependence spread is meaningful.");

        PsiLow = psiLow;
        PsiHigh = psiHigh;
        Rows = AlphaGrid
            .Select(a => new AlphaRow(a, CriticalThreshold(a, psiLow), CriticalThreshold(a, psiHigh)))
            .ToList();
    }

    /// <summary>The row at α = 2 (the unique state-independent order).</summary>
    public AlphaRow AtTwo => Rows.Single(r => r.Alpha == 2.0);

    private static string F(double v) => v.ToString("0.######", Inv);
    private static string E(double v) => v.ToString("E2", Inv);

    public string DisplayName =>
        $"QuarterBoundaryUniquenessWitness (α-sweep at Ψ={F(PsiLow)}, {F(PsiHigh)})";

    public string Summary
    {
        get
        {
            var two = AtTwo;
            int stateIndependent = Rows.Count(r => r.Spread < 1e-12);
            return $"the Rényi fold threshold CΨ*_α is state-independent for exactly {stateIndependent} of " +
                   $"{Rows.Count} swept orders — only α=2, where it equals {F(two.ThresholdLow)} = 1/4 " +
                   $"(spread {E(two.Spread)}); every other α carries a Ψ^(α−2) factor and depends on the state. " +
                   $"The α=2 fixed-point discriminant D = 1−4CΨ has its single zero at CΨ = 1/4 " +
                   $"(D(0.2)={F(Discriminant(0.2))}>0, D(0.25)={F(Discriminant(0.25))}, D(0.3)={F(Discriminant(0.3))}<0): " +
                   $"the unique bifurcation boundary.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // 1. The α-uniqueness: the state-dependence spread per α (0 only at α=2).
            yield return new InspectableNode(
                displayName: "α=2 is the unique state-independent Rényi order",
                summary: $"CΨ*_α = (α−1)^(α−1)/(α^α·Ψ^(α−2)); the Ψ^(α−2) factor vanishes ONLY at α=2, so the " +
                         $"threshold is state-independent (spread=0) there and there alone: " +
                         string.Join("; ", Rows.Select(r =>
                             $"α={F(r.Alpha)} spread={E(r.Spread)}{(r.Alpha == 2.0 ? " (=1/4, unique)" : "")}")) + ".",
                provenance: NodeProvenance.Live);

            // 2. One node per swept α.
            foreach (var r in Rows)
                yield return new InspectableNode(
                    displayName: $"α={F(r.Alpha)}: CΨ*_α = {F(r.ThresholdLow)} @Ψ={F(PsiLow)}, {F(r.ThresholdHigh)} @Ψ={F(PsiHigh)}",
                    summary: $"fold threshold at the two states; spread = {E(r.Spread)} " +
                             $"({(r.Spread < 1e-12 ? "state-INDEPENDENT: this is α=2, CΨ*=1/4" : "state-dependent (Ψ^(α−2)≠1)")}).",
                    provenance: NodeProvenance.Live);

            // 3. The discriminant boundary: D = 1−4CΨ, single sign change at 1/4.
            yield return new InspectableNode(
                displayName: "the discriminant boundary D = 1−4CΨ (unique zero at 1/4)",
                summary: $"R = C(Ψ+R)² expands to C·R² + (2CΨ−1)·R + CΨ² = 0, discriminant D = 1−4CΨ: " +
                         $"D(0.20)={F(Discriminant(0.20))} (>0, two real fixed points), " +
                         $"D(0.25)={F(Discriminant(0.25))} (=0, fold tangency), " +
                         $"D(0.30)={F(Discriminant(0.30))} (<0, none). The single sign change at CΨ = 1/4 is the " +
                         "unique bifurcation boundary (Steps 1-3 of the proof).",
                provenance: NodeProvenance.Live);

            // 4. Payload curve: the state-dependence spread vs α (a V touching zero at α=2).
            var xs = Rows.Select(r => r.Alpha).ToArray();
            var ys = Rows.Select(r => r.Spread).ToArray();
            yield return new InspectableNode(
                displayName: "state-dependence spread vs α (zero only at α=2)",
                summary: "the spread |CΨ*_α(Ψ_lo) − CΨ*_α(Ψ_hi)| across α: touches 0 uniquely at α=2, " +
                         "the state-independence that singles out ¼.",
                payload: new InspectablePayload.Curve("spread vs α", xs, ys, "α", "state-dependence spread"),
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
