using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>Tests for the F86 schema-level Tier1Derived claim that asserts
/// Q_peak ∈ {2 − 1/2, 2 + 1/2} = {1.5, 2.5} via composition of QEpLaw idealised
/// (Tier1Derived) + HalfAsStructuralFixedPointClaim (Tier1Derived). Bit-exact value
/// deviations are Tier2 (PolarityInheritanceLink); closed-form derivation of those
/// deviations is blocked by PROOF_F86B_OBSTRUCTION.</summary>
public class PolarityPairQPeakDecompositionClaimTests
{
    private readonly ITestOutputHelper _out;

    public PolarityPairQPeakDecompositionClaimTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Claim_IsTier1Derived_AtSchemaLevel()
    {
        var claim = PolarityPairQPeakDecompositionClaim.Build();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void QEpCentral_IsTwo_FromQEpLawAtGEffOne()
    {
        Assert.Equal(2.0, PolarityPairQPeakDecompositionClaim.QEpCentral);
    }

    [Fact]
    public void PolarityHalfMagnitude_IsHalf_FromHalfAsStructuralFixedPoint()
    {
        Assert.Equal(0.5, PolarityPairQPeakDecompositionClaim.PolarityHalfMagnitude);
    }

    [Fact]
    public void EndpointSchema_IsTwoPlusHalf_PositivePolarityPole()
    {
        Assert.Equal(2.5, PolarityPairQPeakDecompositionClaim.EndpointQPeakSchema);
    }

    [Fact]
    public void InteriorSchema_IsTwoMinusHalf_NegativePolarityPole()
    {
        Assert.Equal(1.5, PolarityPairQPeakDecompositionClaim.InteriorQPeakSchema);
    }

    [Fact]
    public void SchemaSum_IsTwiceCentral_PairAroundQEpIdealized()
    {
        // The pair {Interior=1.5, Endpoint=2.5} sums to 4 = 2·Q_EP_central = 2·2.
        // This is the polarity-pair invariant: r and −r cancel, leaving 2·central.
        var sum = PolarityPairQPeakDecompositionClaim.InteriorQPeakSchema
                + PolarityPairQPeakDecompositionClaim.EndpointQPeakSchema;
        Assert.Equal(2.0 * PolarityPairQPeakDecompositionClaim.QEpCentral, sum);
    }

    [Fact]
    public void ParentClaims_AreAllTier1Derived()
    {
        var claim = PolarityPairQPeakDecompositionClaim.Build();
        Assert.Equal(Tier.Tier1Derived, claim.Half.Tier);
        Assert.Equal(Tier.Tier1Derived, claim.PolarityOrigin.Tier);
        Assert.Equal(Tier.Tier1Derived, claim.DimensionalAnchor.Tier);
    }

    [Fact]
    public void Build_NeverThrows_AndPopulatesAllParents()
    {
        var claim = PolarityPairQPeakDecompositionClaim.Build();
        Assert.NotNull(claim.Half);
        Assert.NotNull(claim.PolarityOrigin);
        Assert.NotNull(claim.DimensionalAnchor);
    }

    [Fact]
    public void SchemaValues_MatchQAnchorMapEntries()
    {
        // The schema values 1.5 and 2.5 should equal QAnchorMap.CanonicalAnchors entries
        // for F86 Q_peak c=2 and Endpoint orbit candidate respectively.
        var map = new QAnchorMap();
        var interior = map.AnkerAt(PolarityPairQPeakDecompositionClaim.InteriorQPeakSchema);
        var endpoint = map.AnkerAt(PolarityPairQPeakDecompositionClaim.EndpointQPeakSchema);
        Assert.NotNull(interior);
        Assert.NotNull(endpoint);
        Assert.Equal(QBand.Peak, interior!.Band);
        Assert.Equal(QBand.EndpointOrbit, endpoint!.Band);
    }

    [Fact]
    public void EmpiricalDeviation_StaysWithinTenPercent_PolarityInheritanceLinkWitnesses()
    {
        // PolarityInheritanceLink's envelope table, c=2 N=5..8: Endpoint Q_peak 2.517..2.581
        // (at most 3.2% from the schema 2.5), Interior 1.483..1.611 (at most 7.4% from 1.5,
        // at N=8). The schema (Tier1Derived) is exact; the empirical values deviate. The 10%
        // is the claim's own stated finite-size bracket, checked here, not a law.
        const double relativeBudget = 0.10;
        foreach (var w in PolarityInheritanceLink.Build().Witnesses)
        {
            double endpointDev = Math.Abs(w.QPeakEndpoint / PolarityPairQPeakDecompositionClaim.EndpointQPeakSchema - 1.0);
            double interiorDev = Math.Abs(w.QPeakInterior / PolarityPairQPeakDecompositionClaim.InteriorQPeakSchema - 1.0);
            Assert.True(endpointDev <= relativeBudget,
                $"N={w.N} Endpoint {w.QPeakEndpoint} deviates from schema 2.5 by {endpointDev:P1}");
            Assert.True(interiorDev <= relativeBudget,
                $"N={w.N} Interior {w.QPeakInterior} deviates from schema 1.5 by {interiorDev:P1}");
        }
    }

    [Fact]
    public void Children_IncludeParentClaimsForGraphWalking()
    {
        IInspectable inspectable = PolarityPairQPeakDecompositionClaim.Build();
        var children = inspectable.Children.ToList();
        // Tier + Anchor + ExtraChildren (scalars + parent claims)
        Assert.True(children.Count >= 9);  // tier + anchor + 4 scalars + 3 notes + 3 parents = 12+
    }

    [Fact]
    public void Render_EmitsSchemaSummary()
    {
        var claim = PolarityPairQPeakDecompositionClaim.Build();
        var summary = claim.Summary;
        Assert.Contains("Tier1Derived", summary);
        Assert.Contains("schema", summary);
        Assert.Contains("2.5", summary);
        Assert.Contains("1.5", summary);
        _out.WriteLine($"DisplayName: {claim.DisplayName}");
        _out.WriteLine($"Summary: {claim.Summary}");
    }
}
