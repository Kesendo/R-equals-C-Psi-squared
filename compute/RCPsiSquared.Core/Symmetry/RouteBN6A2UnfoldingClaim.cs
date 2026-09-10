using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F163: the complex-q Route-B A2 layer of the uniform N=6 open XY chain,
/// with uniform local Z dephasing, Delta=0 and gamma=1. The degree-133 irreducible
/// layer and Hermitian real-t embedding carry semisimplicity to all 266 loci.
/// The nonzero effective discriminants make the local end-bond response orders
/// one/equal/opposite = 1/1/2. This is an N=6 local theorem.</summary>
public sealed class RouteBN6A2UnfoldingClaim : Claim
{
    public const int Sites = 6;
    public const int FullDimension = 90;
    public const int SectorDimension = 45;
    public const int LociPerParity = 133;
    public const int TotalLoci = 266;
    public const int OneEndOrder = 1;
    public const int EqualEndOrder = 1;
    public const int OppositeEndOrder = 2;

    /// <summary>F131 supplies reflection order sorting, not nonzero coefficients.</summary>
    public MirrorOrderSortingClaim OrderSorting { get; }
    /// <summary>F89d supplies only the explicitly claimed partner-block transport.</summary>
    public F89CrossFoldSimilarityClaim PartnerTransport { get; }

    public RouteBN6A2UnfoldingClaim(MirrorOrderSortingClaim orderSorting,
        F89CrossFoldSimilarityClaim partnerTransport)
        : base("For the uniform N=6 open XY chain under local Z dephasing (Delta=0, gamma=1), " +
               "the 90-dimensional (1,2) coherence block splits as 45+45. Each reflection sector " +
               "has 133 Route-B A2 loci, all semisimple doubles, for 266 complex-q loci. " +
               "Locally in the end-bond epsilon, one/equal/opposite end profiles open two distinct local EP2 branches " +
               "at orders 1/1/2. F131 supplies the odd first-order selection rule; exact nonzero effective " +
               "discriminants supply the leading orders. F89d carries the independently certified " +
               "character to the (1,4) partner block by (q,lambda) -> (conj(q),-conj(lambda)-12). " +
               "Complex-q analytic continuation, local epsilon only, not an all-N or physical measurement claim.",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F163 + docs/proofs/PROOF_ROUTE_B_N6_UNFOLDING.md + " +
               "simulations/route_b_n6_exact_unfolding.py + simulations/results/route_b_n6_exact_unfolding.json + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/RouteBN6A2UnfoldingWitness.cs (inspect --root n6unfolding)")
    {
        OrderSorting = orderSorting ?? throw new ArgumentNullException(nameof(orderSorting));
        PartnerTransport = partnerTransport ?? throw new ArgumentNullException(nameof(partnerTransport));
    }

    public override string DisplayName => "F163: N=6 Route-B A2 semisimplicity and local unfolding";
    public override string Summary =>
        "N=6, 90=45+45, 133+133=266 complex-q semisimple A2 loci; local end-bond " +
        "epsilon orders one/equal/opposite = 1/1/2; not an all-N or physical measurement claim (Tier 1, derived)";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the certified layer",
                summary: "t=i*qCSharp, Lambda=2*lambda; each 45-dimensional sector has residual degree 32 " +
                         "and AT degree 13. The bounded-CRT layer factorization is consumed; the live witness " +
                         "rebuilds the exact pencil, reads the degree-133 cycle modulo 367 and the effective matrices modulo 101.");
            yield return new InspectableNode("F131 premise and nonzero leading orders",
                summary: "On the reflection-symmetric base, R V_- R=-V_- gives P V_- P=0 for the same-parity " +
                         "double plane. Parity alone does not prove a nonzero quadratic response. The local " +
                         "effective discriminants establish epsilon, epsilon, epsilon^2 for one/equal/opposite ends.");
            yield return new InspectableNode("the end-bond epsilon convention",
                summary: "The interior bonds remain one. The (left,right) end weights are " +
                         "(1+epsilon,1), (1+epsilon/2,1+epsilon/2), and (1+epsilon/2,1-epsilon/2). " +
                         "The orders describe EP locations in complex t=i*q, for sufficiently small nonzero epsilon.");
            yield return new InspectableNode("F89d partner transport",
                summary: "The (1,4) partner is an antiunitarily similar 90-dimensional block: " +
                         "(q,lambda) -> (conj(q),-conj(lambda)-12). It inherits the certified Jordan character, " +
                         "and is not an additional independently enumerated 266-locus census.");
            yield return OrderSorting;
            yield return PartnerTransport;
        }
    }
}
