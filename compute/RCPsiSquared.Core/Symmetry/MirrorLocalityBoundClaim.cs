using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The mirror's locality bound (experiments/THE_CUBE_SITE_BY_SITE.md, finding 1, the Z corner;
/// 2026-10-30). For a Lindbladian L = −i[H, ·] + Σ_l D_l with on every site a thermal bath (cooling σ⁻ at γ↓_l,
/// heating σ⁺ at γ↑_l) and an H commuting with Z⊗N, the string Z⊗N, the letter cube's vertex at (0, N, N) and the
/// deepest cell under the baths, is an eigenvector of L at −c, c = Σ_l (γ↓_l + γ↑_l). Any element E with
/// E·L·E⁻¹ = −L − c carries the −c mode onto the 0 mode, so E(Z⊗N) ∝ ρ_ss when the steady state is unique; and a
/// map of operator Schmidt rank r across a cut turns a product across that cut into an operator of rank at most r.
/// Hence, across every cut,
///
/// <para>Schmidt-rank(E) ≥ Schmidt-rank(ρ_ss).</para>
///
/// <para>At one temperature the steady state is the Gibbs product and the bound is 1, which the site-local M_β of
/// KMS Question 4 attains within its scope (XXZ bonds with axial DM, no fields, a finite nonzero ratio). For XXZ
/// bonds (two-body hopping) whose cooling-to-heating ratio differs across a hopping bond the steady state is no
/// product over the sites, exactly: a product steady state would be diagonal site by site by U(1) covariance, the
/// dissipator keeps it diagonal and the hopping takes it off the diagonal, so each factor is the site's own thermal
/// state, and the hopping commutes with that product only at equal ratios. So no invertible product map over the
/// sites, sandwich or not, is an element of the shape −L − c there: finding 1's site-local NONE on a second route.
/// Across a single cut the bound is the steady state's rank across that cut, witnessed full on the chains computed
/// (including a cut crossing only an equal-ratio bond, (3, 3, 2) at N = 3) and not derived. Existence of some element is free wherever a
/// reflection carries the pairing, as S_XY does for XXZ bonds at any rates, by <see cref="PalindromeTwoEndCountClaim"/>
/// §(f8)'s similarity argument (L real in the Pauli basis, every matrix similar to its transpose); where the
/// spectrum pairs and no reflection is known, existence stays open. What the bound measures is how far from a product
/// an element must sit; the least rank attained is not known (at least 4 at N = 2, the trivial ceiling 16).</para>
///
/// <para>Parents: <see cref="PalindromeTwoEndCountClaim"/> (the existence by similarity the bound presupposes)
/// and <see cref="F84ThermalAmplitudeDampingPi2Inheritance"/> (the thermal channel σ⁻/σ⁺ at unequal rates as a
/// typed object; the heights the Z corner adds up are read off that channel here, F84 types its I → Z arrow). Live:
/// the witness <c>MirrorLocalityBoundWitness</c> (Diagnostics), <c>inspect --root locality</c>, which builds the
/// integer generator in the Pauli-string basis, recovers the steady state exactly over ℚ and reads its Schmidt
/// rank as an exact rank, no eigensolver.</para></summary>
public sealed class MirrorLocalityBoundClaim : Claim
{
    public PalindromeTwoEndCountClaim TwoEndCount { get; }
    public F84ThermalAmplitudeDampingPi2Inheritance ThermalChannel { get; }

    public MirrorLocalityBoundClaim(PalindromeTwoEndCountClaim twoEndCount, F84ThermalAmplitudeDampingPi2Inheritance thermalChannel)
        : base("The mirror's locality bound: under thermal baths on every site and an H commuting with Z^N, Z^N is an " +
               "eigenvector of L at -c, every element E with E L E^-1 = -L - c maps it onto the unique steady state's " +
               "line, so the operator Schmidt rank of E across any cut is at least that of rho_ss; at one temperature " +
               "rho_ss is a product (bound 1, attained by M_beta within its scope), across a temperature gradient along a " +
               "hopping bond of XXZ bonds it is no product over the sites, so no invertible product map over the sites is an " +
               "element there",
               Tier.Tier1Derived,
               "experiments/THE_CUBE_SITE_BY_SITE.md")
    {
        TwoEndCount = twoEndCount ?? throw new ArgumentNullException(nameof(twoEndCount));
        ThermalChannel = thermalChannel ?? throw new ArgumentNullException(nameof(thermalChannel));
    }

    public override string DisplayName =>
        "The mirror's locality bound (Schmidt-rank(E) ≥ Schmidt-rank(ρ_ss), the Z corner Z⊗N)";

    public override string Summary =>
        "Z⊗N sits at −c = −Σ(γ↓ + γ↑) for every H commuting with Z⊗N; an element of the shape −L − c carries it onto " +
        "the steady state's line; so the element is at least as entangled as the bath's state across every cut. One " +
        "temperature: a product, bound 1. A gradient along a hopping bond (XXZ): no product over the sites, no site-local " +
        $"element; the least rank is open ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode(
                "the three steps",
                "1. L(Z⊗N) = −c·Z⊗N (heights add, [H, Z⊗N] = 0). 2. E L = (−L − c) E sends an eigenvector at λ to one at " +
                "−λ − c, so E(Z⊗N) ∈ ker L = span(ρ_ss). 3. E = Σ_{k≤r} A_k ⊗ B_k across a cut gives rank(E(P ⊗ Q)) ≤ r.");
            yield return new InspectableNode(
                "where the bound bites",
                "A product steady state is diagonal site by site (U(1) covariance), hence ⊗ of the one-site thermal states, " +
                "and the hopping commutes with it only at equal ratios: a gradient along a hop leaves no product over the " +
                "sites, so no invertible product map over the sites is an element. Exact ranks at the witness: N=2 (3,2 | 1,1) " +
                "rank 4 of 4; N=3 (3,2,1 | 1,1,1, Δ=1) rank 4 across both cuts, (3,3,2 | 1,1,1) rank 4 across the cut crossing " +
                "only the equal-ratio bond; N=4 (6,5,4,3 | 2,2,2,2) rank 16 of 16 at the middle; across single cuts witnessed, not derived.");
            yield return new InspectableNode(
                "what stays open",
                "The least Schmidt rank an element attains (at least 4 at N = 2, the trivial ceiling 16); whether an entangled " +
                "sandwich A·ρᵀ·B exists off the product line; whether an element exists where the spectrum pairs and no " +
                "reflection is known (where S_XY reflects, existence is free: parent F158 §(f8)).");
        }
    }
}
