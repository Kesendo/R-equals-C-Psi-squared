using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The typed ARGUMENT of UNIQUENESS_PROOF (2026-07-02). The VALUE ¼ is already typed
/// (<see cref="QuarterAsBilinearMaxvalClaim"/>, PolynomialDiscriminantAnchorClaim); this claim types
/// the two reasons that single it out.
///
/// <para>Argument 1 (the load-bearing forcing, Step 6): the generalized recursion R = C(Ψ + R)^α has
/// fold threshold CΨ*_α = (α−1)^(α−1) / (α^α · Ψ^(α−2)). The Ψ^(α−2) factor makes it state-DEPENDENT
/// for every α except α = 2, where α−2 = 0 kills the Ψ dependence and the value is exactly ¼. So
/// α = 2 is the UNIQUE Rényi order with a state-independent threshold, and ¼ is what it pins.</para>
///
/// <para>Argument 2 (the discriminant boundary, Steps 1-3): R = C(Ψ + R)² expands to
/// C·R² + (2CΨ − 1)·R + CΨ² = 0, discriminant D = 1 − 4CΨ; D = 0 (fold tangency) iff CΨ = ¼, with a
/// single sign change there — the UNIQUE bifurcation boundary.</para>
///
/// <para>Live witness (recomputes both, elementary arithmetic, no matrices):
/// <c>compute/RCPsiSquared.Diagnostics/Foundation/QuarterBoundaryUniquenessWitness.cs</c>,
/// <c>inspect --root quarter-uniqueness</c>.</para></summary>
public sealed class QuarterBoundaryUniquenessClaim : Claim
{
    /// <summary>The typed parent: ¼ as the bilinear maxval, the value this argument singles out.</summary>
    public QuarterAsBilinearMaxvalClaim Quarter { get; }

    public QuarterBoundaryUniquenessClaim(QuarterAsBilinearMaxvalClaim quarter)
        : base("Quarter-boundary uniqueness: CΨ = 1/4 is the unique bifurcation, α = 2 the unique Rényi order with a state-independent fold threshold",
               Tier.Tier1Derived,
               "docs/proofs/UNIQUENESS_PROOF.md")
    {
        Quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
    }

    public override string DisplayName =>
        "Quarter-boundary uniqueness: α = 2 is the unique state-independent Rényi order; CΨ = 1/4 the unique bifurcation";

    public override string Summary =>
        "CΨ = 1/4 is singled out two ways. (1) The Rényi fold threshold CΨ*_α = (α−1)^(α−1)/(α^α·Ψ^(α−2)) is " +
        "state-independent iff α = 2 (the Ψ^(α−2) factor vanishes), where it equals 1/4 — the unique order. " +
        "(2) The α=2 fixed-point quadratic C·R²+(2CΨ−1)·R+CΨ²=0 has discriminant D = 1−4CΨ, whose single zero " +
        "at CΨ = 1/4 is the unique bifurcation boundary. Value is the parent ¼ (QuarterAsBilinearMaxval); this " +
        "node types the uniqueness argument (live: QuarterBoundaryUniquenessWitness).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return Quarter;
            yield return new InspectableNode("the α=2 forcing (state-independence)",
                summary: "CΨ*_α = (α−1)^(α−1)/(α^α·Ψ^(α−2)) is the fold-bifurcation threshold of R = C(Ψ+R)^α; " +
                         "the Ψ^(α−2) factor makes it state-dependent for every α EXCEPT α = 2 (exponent α−2 = 0), " +
                         "where it is exactly 1/4 regardless of Ψ. α = 2 is thus the unique Rényi order whose " +
                         "critical boundary is state-independent — the reason 1/4 is singled out (Step 6).");
            yield return new InspectableNode("the discriminant boundary (uniqueness of 1/4)",
                summary: "R = C(Ψ+R)² ⟹ C·R² + (2CΨ−1)·R + CΨ² = 0, discriminant D = (2CΨ−1)² − 4C²Ψ² = 1 − 4CΨ. " +
                         "D = 0 iff CΨ = 1/4 (fold tangency); D > 0 (two real fixed points) below; D < 0 (none) above. " +
                         "The single sign change at 1/4 is the unique bifurcation boundary (Steps 1-3).");
            yield return new InspectableNode("live witness",
                summary: "QuarterBoundaryUniquenessWitness (inspect --root quarter-uniqueness) sweeps α and reports " +
                         "the fold threshold at two probe states Ψ — the spread is 0 only at α = 2 (=1/4) — and the " +
                         "discriminant D = 1−4CΨ across CΨ, its single zero at 1/4. Elementary arithmetic, exact.");
        }
    }
}
